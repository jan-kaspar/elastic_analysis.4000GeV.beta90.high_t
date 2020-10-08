#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

struct HistGroup
{
	TH1D *T;	///< histogram of t_true, with simulation event weight
	TH1D *R;	///< histogram of t_reco, with weight w_simu * w_reco

	TH2D *M;	///< histogram of t_reco vs. t_true, with weight w_simu * w_reco

	HistGroup() : T(NULL), R(NULL), M(NULL) {}

	void Init(const string &binning)
	{
		unsigned int N_bins;
		double *bin_edges;
		BuildBinning(anal, binning, bin_edges, N_bins, "./");

		T = new TH1D("", ";t true", N_bins, bin_edges);
		R = new TH1D("", ";t reco", N_bins, bin_edges);

		M = new TH2D("", ";t true;t reco", N_bins, bin_edges, N_bins, bin_edges);
	}

	void FillTrue(double t_true, double w_simu)
	{
		T->Fill(t_true, w_simu);
	}

	void FillReco(double t_true, double w_simu, double t_reco, double w_reco)
	{
		R->Fill(t_reco, w_simu * w_reco);

		M->Fill(t_true, t_reco, w_simu * w_reco);
	}

	void Write()
	{
		T->SetLineColor(1);
		T->Write("T");

		TH1D *T_norm = new TH1D(*T);
		T_norm->SetLineColor(2);
		T_norm->Scale(1., "width");
		T_norm->Write("T_norm");

		R->SetLineColor(4);
		R->Write("R");

		TH1D *R_norm = new TH1D(*R);
		R_norm->SetLineColor(6);
		R_norm->Scale(1., "width");
		R_norm->Write("R_norm");

		TH1D *C = new TH1D(*T);
		TVectorD C_vec(M->GetNbinsX());
		C->SetLineColor(1);
		for (int bi = 1; bi <= C->GetNbinsX(); bi++)
		{
			double t = T->GetBinCenter(bi);
			double v_T = T->GetBinContent(bi);
			double v_R = R->GetBinContent(bi);
			double corr = (v_R > 0. && t > anal.t_min && t < anal.t_max) ? v_T / v_R : 0.;

			C->SetBinContent(bi, corr);
			C_vec(bi-1) = corr;		// matrix indeces start with 0
		}
		C->Write("C");
		C_vec.Write("C_vec");

		M->Write("M");

		TH2D *S = new TH2D(*M);
		TH2D *S_norm = new TH2D(*M);
		TMatrixD S_mat(M->GetNbinsX(), M->GetNbinsX());
		TMatrixD S_norm_mat(M->GetNbinsX(), M->GetNbinsX());
		for (int bi = 1; bi <= M->GetNbinsX(); bi++)
		{
			for (int bj = 1; bj <= M->GetNbinsY(); bj++)
			{
				double v = M->GetBinContent(bj, bi);	// M_ij
				double n = T->GetBinContent(bj);		// T_j
				double e = (n > 0.) ? v / n : 0.;

				S->SetBinContent(bj, bi, e);			// S_ij = M_ij / T_j
				S_mat(bi-1, bj-1) = e;					// matrix indeces start with 0, note the i <--> j swap

				double w_j = T->GetBinWidth(bj);
				double w_i = T->GetBinWidth(bi);
				S_norm->SetBinContent(bj, bi, e * w_j / w_i);
				S_norm_mat(bi-1, bj-1) = e * w_j / w_i;
			}
		}
		S->Write("S");
		S_mat.Write("S_mat");
		S_norm->Write("S_norm");
		S_norm_mat.Write("S_norm_mat");

		// verification: calculate R_i = Sum_j S_ij T_j
		TH1D *R_test = new TH1D(*R);
		TH1D *R_norm_test = new TH1D(*R_norm);

		for (int bi = 1; bi <= R->GetNbinsX(); bi++)
		{
			double sum = 0.;		// R_i
			double sum_norm = 0.;	// R_norm_i 
			
			for (int bj = 1; bj <= T->GetNbinsX(); bj++)
			{
				// S_ij * T_j
				sum += S->GetBinContent(bj, bi) * T->GetBinContent(bj);

				sum_norm += S_norm->GetBinContent(bj, bi) * T_norm->GetBinContent(bj);
			}

			R_test->SetBinContent(bi, sum);

			R_norm_test->SetBinContent(bi, sum_norm);
		}

		R_test->SetLineColor(8);
		R_test->Write("R_test");

		R_norm_test->SetLineColor(8);
		R_norm_test->Write("R_norm_test");

		TCanvas *c = new TCanvas("R cmp");
		c->SetLogy(1);
		R->Draw("");
		R_test->Draw("same");
		c->Write();

		c = new TCanvas("R_norm cmp");
		c->SetLogy(1);
		R_norm->Draw("");
		R_norm_test->Draw("same");
		c->Write();

		// verification with matrices: calculate R_i = Sum_j S_ij T_j
		int dim = T->GetNbinsX();

		TVectorD T_vec(dim), T_norm_vec(dim);
		for (int bi = 1; bi <= dim; bi++)
		{
			int i = bi - 1;
			T_vec(i) = T->GetBinContent(bi);
			T_norm_vec(i) = T_norm->GetBinContent(bi);
		}

		TVectorD R_vec = S_mat * T_vec, R_norm_vec = S_norm_mat * T_norm_vec;

		TH1D *R_test_mat = new TH1D(*R); R_test_mat->SetLineColor(12);
		TH1D *R_norm_test_mat = new TH1D(*R); R_norm_test_mat->SetLineColor(12);
		
		TH1D *C_test_mat = new TH1D(*C); C_test_mat->SetLineColor(6);
		TH1D *C_norm_test_mat = new TH1D(*C); C_norm_test_mat->SetLineColor(8);

		for (int i = 0; i < dim; i++)
		{
			int bi = i + 1;

			R_test_mat->SetBinContent(bi, R_vec(i));
			R_norm_test_mat->SetBinContent(bi, R_norm_vec(i));

			double C_test = (R_vec(i) > 0.) ? T_vec(i) / R_vec(i) : 0.;
			double C_norm_test = (R_norm_vec(i) > 0.) ? T_norm_vec(i) / R_norm_vec(i) : 0.;

			C_test_mat->SetBinContent(bi, C_test);
			C_norm_test_mat->SetBinContent(bi, C_norm_test);
		}

		R_test_mat->Write("R_test_mat");
		R_norm_test_mat->Write("R_norm_test_mat");

		C_test_mat->Write("C_test_mat");
		C_norm_test_mat->Write("C_norm_test_mat");
	}
};

//----------------------------------------------------------------------------------------------------
	
TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// settings
	unsigned int N_ev = (unsigned int) 1E9;
	unsigned int seed = 1;

	printf("* settings:\n");
	printf("	N_ev = %.2E\n", (double) N_ev);
	printf("	seed = %u\n", seed);

	// binnings
	vector<string> binnings;
	//binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	//binnings.push_back("ob-2-20-0.20");
	//binnings.push_back("ob-3-10-0.30");
	binnings.push_back("bt1");
	binnings.push_back("bt2");

	// models
	vector<string> models;
	models.push_back("exp3-intf-exp1");
	models.push_back("p1*exp3+p2*exp2");

	// prepare output
	TFile *f_out = new TFile((string("smearing_matrix_mc_") + argv[1] + ".root").c_str(), "recreate");

	// build matrices
	for (unsigned int mi = 0; mi < models.size(); ++mi)
	{
		TDirectory *modDir = f_out->mkdir(models[mi].c_str());

		// load model
		TFile *f_in_model = new TFile(("../models/" + models[mi] + ".root").c_str());
		TGraph *g_dsdt = (TGraph *) f_in_model->Get("g_dsdt");
		TSpline *s_dsdt = BuildSpline(g_dsdt);

		// re-initialize seed
		gRandom->SetSeed(seed);
		
		// book histograms
		vector<HistGroup> hists;
		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			HistGroup hg;
			hg.Init(binnings[bi]);
			hists.push_back(hg);
		}
	
		// run simulation
		for (unsigned int evi = 0; evi < N_ev; evi++)
		{
			// generate true event
			Kinematics k_tr;
			k_tr.t = gRandom->Rndm() * (anal.t_max_full - anal.t_min_full) + anal.t_min_full;
			k_tr.phi = th_y_sign * gRandom->Rndm() * M_PI * th_y_sign;	// just one half
			double w = s_dsdt->Eval(k_tr.t);	// event weight
	
			k_tr.th = sqrt(k_tr.t) / env.p;
			k_tr.th_x = k_tr.th_x_L = k_tr.th_x_R = k_tr.th * cos(k_tr.phi);
			k_tr.th_y = k_tr.th_y_L = k_tr.th_y_R = k_tr.th * sin(k_tr.phi);
	
			for (unsigned int bi = 0; bi < binnings.size(); bi++)
				hists[bi].FillTrue(k_tr.t, w);
	
			// apply smearing
			// IMPOARTANT: make sure that this is compatible with the reconstruction formula in `common_algorithms.h'
			Kinematics k_sm = k_tr;
			k_sm.th_x_L += gRandom->Gaus() * anal.si_th_x_1arm;
			k_sm.th_x_R += gRandom->Gaus() * anal.si_th_x_1arm;
			k_sm.th_x = (k_sm.th_x_L + k_sm.th_x_R) / 2.;
	
			k_sm.th_y_L += gRandom->Gaus() * anal.si_th_y_1arm;
			k_sm.th_y_R += gRandom->Gaus() * anal.si_th_y_1arm;
			k_sm.th_y = (k_sm.th_y_L + k_sm.th_y_R) / 2.;
	
			k_sm.th = sqrt(k_sm.th_x*k_sm.th_x + k_sm.th_y*k_sm.th_y);
			k_sm.t = env.p*env.p * k_sm.th * k_sm.th;
	
			// acceptance cut and correction
			double phi_corr = 0., div_corr = 0.;
			bool skip = CalculateAcceptanceCorrections(+1., k_sm, anal, phi_corr, div_corr);
			if (!skip)
			{
				double corr = div_corr * phi_corr / 2.; // only one half in generation
	
				for (unsigned int bi = 0; bi < binnings.size(); bi++)
					hists[bi].FillReco(k_tr.t, w, k_sm.t, corr);
			}
		}

		// save histograms
		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			gDirectory = modDir->mkdir(binnings[bi].c_str());
			hists[bi].Write();	
		}
	}

	delete f_out;

	return 0;
}
