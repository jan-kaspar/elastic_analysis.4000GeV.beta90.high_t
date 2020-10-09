#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"

#include <string>
#include <map>
#include <cstdio>

#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"


using namespace std;

//----------------------------------------------------------------------------------------------------

void replace(string &str, const string &f, const string &r)
{
	size_t pos = str.find(f);
	if (pos == string::npos)
		return;

	str.replace(pos, f.length(), r);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

// global variables used in minimisation
TMatrixD S;				// smearing (response) matrix
TVectorD v_data;		// data vector (smeared)
TVectorD v_data_unc;	// vector of data uncertainties
TVectorD v_data_x;		// vector of data point abscissas
TVectorD v_true;		// true vector (unsmeared)
int bi_min, bi_max;		// bin range (inclusive) used in minimisation

//----------------------------------------------------------------------------------------------------

class S2_FCN : public ROOT::Minuit2::FCNBase
{
	public:
		S2_FCN() {}

  		double operator() (const std::vector<double> &) const;
  		double Up() const { return 1.; }
};

//----------------------------------------------------------------------------------------------------

int S2_debug = 0;

double S2_alpha = 0.;

double S2_FCN::operator() (const std::vector<double> &par) const
{
	if (S2_debug > 3)
		printf("\n\n--------------------\n");

	int dim = v_data.GetNrows();

	// ----- matrix-inversion component -----
	for (int i = 0; i < dim; i++)
		v_true(i) = par[i];

	TVectorD v_sm = S * v_true;

	// (inclusive indeces) starting with 0
	int i_min = bi_min - 1, i_max = bi_max - 1;
	int n_point = i_max - i_min + 1;

	double S2_mat_inv = 0.;
	for (int i = i_min; i <= i_max; i++)
	{
		double diff = (v_sm(i) - v_data(i)) / v_data_unc(i);
		S2_mat_inv += diff*diff;

		if (S2_debug > 3)
			printf("v_sm = %.2E, v_data = %.2E, unc = %.2E | rel diff = %.2f\n", v_sm(i), v_data(i), v_data_unc(i), diff);
	}

	// ----- smoothness component -----

	// check of positiveness
	int N_non_pos = 0;
	for (int i = i_min; i <= i_max; i++)
	{
		if (v_true(i) <= 0.)
			N_non_pos++;
	}

	double S2_smooth = 0.;
	double S2_smooth_I=0., S2_smooth_II=0., S2_smooth_III=0.;

	if (N_non_pos > 0)
	{
		S2_smooth = 1E6 * N_non_pos;
	} else {
		// integrate second derivatives of log(ds/dt)
		S2_smooth = 0.;
		for (int i = i_min + 1; i <= i_max - 1; i++)
		{
			double xm = v_data_x(i-1), ym = log(v_true(i-1));
			double x0 = v_data_x(i), y0 = log(v_true(i));
			double xp = v_data_x(i+1), yp = log(v_true(i+1));
	
			double fdm = (yp - y0) / (xp - x0);
			double fdp = (y0 - ym) / (x0 - xm);
	
			double w = (xp - xm) / 2.;
			double sd = (fdp - fdm) / w;
	
			double con = fabs(sd) * w;
	
			S2_smooth += con * con;

			if (x0 <= 0.4) S2_smooth_I += con * con;
			if (x0 > 0.4 && x0 <= 0.55) S2_smooth_II += con * con;
			if (x0 > 0.55) S2_smooth_III += con * con;
	
			if (S2_debug > 3)
				printf("%3i: w = %.2E, sd = %.2E, con = %.2E\n", i, w, sd, con);
		}

		// TODO: test only, remove; or keep ??
		S2_smooth = (S2_smooth_I + S2_smooth_II / 30. + S2_smooth_III) * 3.;
	}

	// ----- combine S2's -----

	double S2_comb = S2_mat_inv + S2_alpha * S2_smooth;

	// ----- debug print -----
	if (S2_debug > 0)
	{
		printf("mat_inv: S2=%.2E, n_pnt=%i, S2/n_pnt=%.2f | smooth: S2=%.2E, al*S2=%.2f | comb: S2=%.2f | ratio: sm/mi=%.2f\n",
			S2_mat_inv, n_point, S2_mat_inv/n_point,
			S2_smooth, S2_alpha*S2_smooth,
			S2_comb,
			(S2_alpha * S2_smooth) / (S2_mat_inv)
			);
	}

	if (S2_debug > 0)
		printf("smooth al*S2: I=%.3f, II=%.3f, III=%.3f\n", S2_alpha*S2_smooth_I, S2_alpha*S2_smooth_II, S2_alpha*S2_smooth_III);

	return S2_comb;
}

//----------------------------------------------------------------------------------------------------

void BuildInitialTrueVector(TH1D *h_data, const TVectorD &C)
{
	int dim = C.GetNrows();

	h_data->SetLineColor(1);
	h_data->SetName("h_data_input");
	
	// apply correction
	TH1D *h_data_corr = new TH1D(*h_data);
	for (int i = 0; i < dim; i++)
	{
		int bi = i + 1;
		h_data_corr->SetBinContent(bi, h_data->GetBinContent(bi) * C(i));
	}
	
	h_data_corr->SetLineColor(2);
	h_data_corr->SetName("h_data_corr");

	// extrapolate outside acceptance
	TH1D *h_data_corr_ext = new TH1D(*h_data_corr);
	h_data_corr_ext->SetLineColor(4);
	h_data_corr_ext->SetName("h_data_corr_ext");

	TF1 *ff = new TF1("ff", "[0]*exp([1]*x)");
	
	ff->SetParameters(530., -19.7);
	h_data_corr->Fit(ff, "Q+", "", anal.t_min, 0.1);
	for (int bi = 1; bi <= h_data_corr_ext->GetNbinsX(); bi++)
	{
		double c = h_data_corr_ext->GetBinCenter(bi);
		double v = h_data_corr_ext->GetBinContent(bi);

		if (v > 0.)
		{
			bi_min = bi;
			break;
		}

		double v_ext = ff->Eval(c);
		h_data_corr_ext->SetBinContent(bi, v_ext);
		h_data_corr_ext->SetBinError(bi, v_ext * 0.05);
	}

	ff->SetParameters(3., -5.);
	h_data_corr->Fit(ff, "Q+", "", 1.4, anal.t_max);
	for (int bi = 1; bi <= h_data_corr_ext->GetNbinsX(); bi++)
	{
		double c = h_data_corr_ext->GetBinCenter(bi);

		if (c < anal.t_max)
			bi_max = bi;

		if (c > anal.t_max && c < anal.t_max_full)
		{
			double v_ext = ff->Eval(c);
			h_data_corr_ext->SetBinContent(bi, v_ext);
			h_data_corr_ext->SetBinError(bi, v_ext * 0.7);
		}
	}

	delete ff;

	// regularise empty bins
	for (int bi = 1; bi <= h_data_corr_ext->GetNbinsX(); ++bi)
	{
		double v = h_data_corr_ext->GetBinContent(bi);
		double u = h_data_corr_ext->GetBinError(bi);

		if (v == 0. || u == 0.)
		{
			double v_left = 1., u_left = 1.;
			for (int bj = bi-1; bj >= 1; --bj)
			{
				v = h_data_corr_ext->GetBinContent(bj);
				u = h_data_corr_ext->GetBinError(bj);
				if (v > 0.)
				{
					v_left = v;
					u_left = u;
					break;
				}
			}
			
			double v_right = 1., u_right = 1.;
			for (int bj = bi+1; bj <= h_data_corr_ext->GetNbinsX(); ++bj)
			{
				v = h_data_corr_ext->GetBinContent(bj);
				u = h_data_corr_ext->GetBinError(bj);
				if (v > 0.)
				{
					v_right = v;
					u_right = u;
					break;
				}
			}

			v = (v_left + v_right) / 2.;
			u = (u_left + u_right) / 2.;
			h_data_corr_ext->SetBinContent(bi, v);
			h_data_corr_ext->SetBinError(bi, u);

			printf("bin %i regularised: value = %E, uncertainty = %E\n", bi, v, u);
		}
	}

	// comparison canvas
	TCanvas *c = new TCanvas("h_data cmp");
	c->SetLogy(1);
	h_data_corr_ext->Draw("");
	h_data_corr->Draw("same");
	h_data->Draw("same");
	c->Write();

	h_data_corr_ext->Write("h_unfold_init");

	// save to v_true and v_data_unc
	v_true.ResizeTo(dim);
	v_data_unc.ResizeTo(dim);
	for (int i = 0; i < dim; i++)
	{
		int bi = i + 1;
		v_true(i) = h_data_corr_ext->GetBinContent(bi);
		v_data_unc(i) = h_data_corr_ext->GetBinError(bi);
	}

	// print info
	printf("* bi_min = %i, bi_max = %i\n", bi_min, bi_max);
}

//----------------------------------------------------------------------------------------------------

void PrintMinuit2State(TFitterMinuit *minuit)
{
	for (int pi = 0; pi < minuit->GetNumberTotalParameters(); pi++)
	{
		char name[20];
		double value, verr, vlow, vhigh;
		minuit->GetParameter(pi, name, value, verr, vlow, vhigh);
		printf("\t%3i, %7s, value=%8.2E, unc=%8.2E, min=%+7.3f, max=%+7.3f", pi, minuit->GetParName(pi), value, verr, vlow, vhigh);

		if (vhigh > vlow)
			printf(", limited");

		if (minuit->IsFixed(pi))
			printf(", fixed");

		printf("\n");
	}
}

//----------------------------------------------------------------------------------------------------

void DoUnfolding(TH1D *h_data, const TMatrixD &S_in, const TVectorD &C_in, double alpha)
{
	// dimensionality check
	int dim_data = h_data->GetNbinsX();
	int dim_S_r = S_in.GetNrows(), dim_S_c = S_in.GetNcols();
	int dim_C_r = S_in.GetNrows();

	printf("h_data: %i\n", dim_data);
	printf("S: %i x %i\n", dim_S_r, dim_S_c);
	printf("C: %i\n", dim_C_r);

	if (dim_data != dim_S_r || dim_data != dim_S_c || dim_data != dim_C_r || dim_S_r != dim_S_c || dim_S_r != dim_C_r || dim_S_c != dim_C_r)
	{
		printf("ERROR: dimension mismatch, stop.\n");
		return;
	}

	int dim = dim_data;

	// copy smearing matrix
	S.ResizeTo(S_in);
	S = S_in;

	// build data vector
	v_data.ResizeTo(dim);
	v_data_x.ResizeTo(dim);
	for (int i = 0; i < dim; i++)
	{
		int bi = i + 1;		// histogram indeces start with 1
		v_data(i) = h_data->GetBinContent(bi);
		v_data_x(i) = h_data->GetBinCenter(bi);
	}

	// build initial true vector (output to v_true),
	// regularise uncertainties (extrapolation, empty bins) and
	// determine bin range
	BuildInitialTrueVector(h_data, C_in);

	// run minimisation
	S2_alpha = alpha;

	TFitterMinuit *minuit = new TFitterMinuit();
	S2_FCN fcn;
	minuit->SetMinuitFCN(&fcn);
	minuit->SetPrintLevel(0);
	minuit->SetStrategy(2);
	minuit->CreateMinimizer(TFitterMinuit::kMigrad);
	//minuit->CreateMinimizer(TFitterMinuit::kCombined);

	for (int i = 0; i < dim; i++)
	{
		char buf[30];
		sprintf(buf, "bin %i", i+1);
		minuit->SetParameter(i, buf, v_true(i), v_data_unc(i), 0., 0.);
	}

	int i_min = bi_min - 1, i_max = bi_max - 1;
	int margin = 1;

	for (int i = 0; i < i_min - margin; i++)
		minuit->FixParameter(i);

	for (int i = i_max + 1 + margin; i < dim; i++)
		minuit->FixParameter(i);

	printf("\n* state before fit\n");
	PrintMinuit2State(minuit);

	printf("\n* fit\n");
	double minuit_tolerance = 10.;	// without factor 1E-3
	S2_debug = 0;
	minuit->Minimize(0, minuit_tolerance);

	printf("\n* state after fit\n");
	PrintMinuit2State(minuit);

	// print details
	vector<double> param(dim);
	for (int i = 0; i < dim; i++)
		param[i] = minuit->GetParameter(i);

	printf("\n* S2 details at solution\n");
	S2_debug = 1;
	fcn(param);

	// get range of free parameters
	int pfi_min = 12345, pfi_max = 0;
	for (int i = 0; i < dim; i++)
	{
		if (! minuit->IsFixed(i))
		{
			pfi_min = min(pfi_min, i);
			pfi_max = max(pfi_max, i);
		}
	}

	printf("* range of free parameters: %i to %i\n", pfi_min, pfi_max);

	// correction
	for (int i = 0; i < dim; i++)
		v_true(i) = minuit->GetParameter(i);
	TVectorD v_sm = S * v_true;

	TMatrixD V_true(dim, dim);
	TVectorD corr(dim);
	TMatrixD M(dim, dim);
	for (int i = 0; i < dim; i++)
	{
		int bi = i + 1;
		double t = h_data->GetBinCenter(bi);
		corr(i) = (v_sm(i) > 0. && t > anal.t_min && t < anal.t_max) ? v_true(i) / v_sm(i) : 0.;

		for (int j = 0; j < dim; j++)
		{
			V_true(i, j) = 0.;
			if (pfi_min <= i && i <= pfi_max && pfi_min <= j && j <= pfi_max)
				V_true(i, j) = minuit->GetCovarianceMatrixElement(i - pfi_min, j - pfi_min); // this function only considers non-fixed parameters

			//if (i == j)
			//.	printf("i = %i, unc = %E, sqrt(V) = %E\n", i, minuit->GetParError(i), sqrt(V_true(i, j)));

			double R = v_sm(i);

			M(i, j) = 0.;
			if (R > 0.)
			{
				if (i == j) M(i, j) += 1./R;
				M(i, j) += - v_true(i) * S(i, j) / R/R;
			}

			//if (i == j)
			//.	printf("i = %3u | V(i, j) = %E | R = %E, T = %E, S(i, j) = %E | M(i, j) = %E\n", i, V_true(i, j), R, v_true(i), S(i, j), M(i, j));
		}
	}
	
	TMatrixD MT(TMatrixD::kTransposed, M);
	TMatrixD V_corr = M * V_true * MT;

	// visualise result
	TH1D *h_unfold = new TH1D(*h_data);
	h_unfold->SetName("h_unfold");
	h_unfold->SetLineColor(6);

	TH1D *h_unfold_resmear = new TH1D(*h_data);
	h_unfold_resmear->SetName("h_unfold_resmear");
	h_unfold_resmear->SetLineColor(8);

	TH1D *h_corr_raw = new TH1D(*h_data);
	h_corr_raw->SetName("h_corr_raw");
	h_corr_raw->SetLineColor(1);

	for (int i = 0; i < dim; i++)
	{
		int bi = i + 1;
		h_unfold->SetBinContent(bi, minuit->GetParameter(i));
		h_unfold->SetBinError(bi, minuit->GetParError(i));

		h_unfold_resmear->SetBinContent(bi, v_sm(i));
		h_unfold_resmear->SetBinError(bi, 0.);

		h_corr_raw->SetBinContent(bi, corr(i));
		h_corr_raw->SetBinError(bi, sqrt(V_corr(i, i)));
	}

	h_unfold->Write();
	h_unfold_resmear->Write();
	h_corr_raw->Write();

	// regularise high-|t| bins
	double av = 0.;
	TH1D *h_corr = new TH1D(*h_corr_raw);
	h_corr->SetName("h_corr");

	for (int bi = 1; bi <= h_corr->GetNbinsX(); ++bi)
	{
		if (h_corr->GetBinCenter(bi) < 1.7)
			av = h_corr->GetBinContent(bi);
		else
			h_corr->SetBinContent(bi, av);
	}

	h_corr->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	if (argc < 2)
	{
		printf("ERROR: Wrong number of parameters.\n");
		return 1;
	}
	
	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// binnings
	vector<string> binnings;
	vector< vector<double> > alphasPerBinning;
	binnings.push_back("ob-1-30-0.10"); alphasPerBinning.push_back({1E0, 5E-1, 1E-1, 5E-2, 1E-3});
	//binnings.push_back("ob-2-20-0.20"); alphasPerBinning.push_back({1E-3});1E0, 1E-1, 5E-2, 
	//binnings.push_back("ob-3-10-0.30"); alphasPerBinning.push_back({5E-1, 5E-2, 5E-3});
	binnings.push_back("bt1"); alphasPerBinning.push_back({1E0, 5E-1, 1E-1, 5E-2, 1E-2, 5E-3, 1E-3});
	binnings.push_back("bt2"); alphasPerBinning.push_back({1E0, 5E-1, 1E-1, 5E-2, 1E-2, 5E-3, 1E-3});

	// unsmearing models
	vector< vector<string> > models;
	models.push_back({"smearing_matrix_mc_<dgn>", "exp3-intf-exp1", "<binning>"});
	models.push_back({"smearing_matrix_mc_<dgn>", "p1*exp3+p2*exp2", "<binning>"});

	// open input file
	string fn_in = string("distributions_") + argv[1] + ".root";
	TFile *f_in = new TFile(fn_in.c_str());
	if (f_in->IsZombie())
	{
		printf("ERROR: can't open file `%s'.\n", fn_in.c_str());
		return 1;
	}
	
	// prepare output
	string fn_out = string("unfolding_gr_") + argv[1] + ".root";
	TFile *f_out = new TFile(fn_out.c_str(), "recreate");
	if (f_out->IsZombie())
	{
		printf("ERROR: can't open file `%s' for writing.\n", fn_out.c_str());
		return 3;
	}
	
	// run method for all combinations
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		TDirectory *binDir = f_out->mkdir(binnings[bi].c_str());

		for (unsigned int mi = 0; mi < models.size(); mi++)
		{
			// compile model file and tag strings
			string mf = models[mi][0] + ".root";
			replace(mf, "<dgn>", argv[1]);
			replace(mf, "<binning>", binnings[bi]);

			string mt = models[mi][1] + "/" + models[mi][2];
			replace(mt, "<dgn>", argv[1]);
			replace(mt, "<binning>", binnings[bi]);

			printf("\n\n");
			printf("* binning = %s\n", binnings[bi].c_str());
			printf("* model file = %s\n", mf.c_str());
			printf("* model tag = %s\n", mt.c_str());

			// directory name
			string dn = mf + "," + mt;
			replace(dn, "/", ",");
			TDirectory *modDir = binDir->mkdir(dn.c_str());

			// load objects
			TFile *f_mod_in = TFile::Open(mf.c_str());
			if (f_mod_in == NULL)
			{
				printf("ERROR: can't open file `%s', skipping.\n", mf.c_str());
				continue;
			}

			TMatrixD *S = (TMatrixD *) f_mod_in->Get((mt + "/S_norm_mat").c_str());
			TVectorD *C = (TVectorD *) f_mod_in->Get((mt + "/C_vec").c_str());

			TH1D *h_data = (TH1D *) f_in->Get(("normalization/" + binnings[bi] + "/h_t_normalized").c_str());

			if (!S || !C || !h_data)
			{
				printf("ERROR: can't load input data: S=%p, C=%p, h_data=%p, skipping\n", S, C, h_data);
				continue;
			}

			// try different values of alpha
			for (unsigned int ai = 0; ai < alphasPerBinning[bi].size(); ai++)
			{
				double alpha = alphasPerBinning[bi][ai];
				printf("* alpha = %.2E\n", alpha);

				char buf[50];
				sprintf(buf, "alpha=%.2E", alpha);
				gDirectory = modDir->mkdir(buf);
				
				// run unfolding
				DoUnfolding(h_data, *S, *C, alpha);
			}

			delete f_mod_in;
		}
	}

		
	delete f_out;

	return 0;
}
