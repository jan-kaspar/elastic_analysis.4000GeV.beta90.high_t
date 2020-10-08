#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"

#include <string>
#include <map>
#include <cstdio>
#include <cassert>

#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

#include "NumericalIntegration.h"


using namespace std;

//----------------------------------------------------------------------------------------------------

unsigned long int_ws_phi_size, int_ws_MX_size, int_ws_MY_size, int_ws_DY_size;
gsl_integration_workspace *int_ws_phi, *int_ws_MX, *int_ws_MY, *int_ws_DY;

//----------------------------------------------------------------------------------------------------

double n_si;

double si_de_th_y_L, si_de_th_y_R;
double si_m_x, si_m_y, si_De_y; 

double th_x_lcut;
double th_x_hcut;
double th_y_lcut;
double th_y_hcut;

bool includeAcceptanceEffects = true;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

double dist_m_x(double m)
{
	double si = si_m_x;
	return exp( - (m * m) / 2. / si / si ) / sqrt(2. * M_PI ) / si;
}

//----------------------------------------------------------------------------------------------------

double dist_de_th_y_L(double de_th_y)
{
	double si = si_de_th_y_L;
	return exp( - (de_th_y * de_th_y) / 2. / si / si ) / sqrt(2. * M_PI ) / si;
}

//----------------------------------------------------------------------------------------------------

double dist_de_th_y_R(double de_th_y)
{
	double si = si_de_th_y_R;
	return exp( - (de_th_y * de_th_y) / 2. / si / si ) / sqrt(2. * M_PI ) / si;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TF1 *source_dist_t = NULL;

// true distribution of t
double dist_t_true(double t)
{
	return source_dist_t->Eval(t);
}

//----------------------------------------------------------------------------------------------------

// true distribution of th_x, th_y
double dist_th_x_th_y_true(double th_x, double th_y)
{
	double t = (th_x*th_x + th_y*th_y) * env.p*env.p;
	return dist_t_true(t);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

double IntegOverDY(double x, double *p, const void *)
{
	double m_y = p[2];
	double d_y = x;

	double de_th_y_L = m_y - d_y/2.;
	double de_th_y_R = m_y + d_y/2.;

	return dist_de_th_y_L(de_th_y_L) * dist_de_th_y_R(de_th_y_R);
}

//----------------------------------------------------------------------------------------------------

double IntegOverMX(double x, double *p, const void *)
{

	double th_x = p[0];
	double th_y = p[1];
	double m_x = x;
	double m_y = p[2];

	return dist_th_x_th_y_true(th_x - m_x, th_y - m_y) * dist_m_x(m_x);
}

//----------------------------------------------------------------------------------------------------

double IntegOverMY(double x, double *p, const void *)
{
	double param[] = { p[0], p[1], x };

	// ----- integral over m_x -----

	double range_x = n_si * si_m_x;
	double precision = 1E-2;
	double imx = RealIntegrate(IntegOverMX, param, NULL, -range_x, +range_x, 0., precision,
		int_ws_MX_size, int_ws_MX, "IntegOverMY-imx");

	// ----- integral over d_y -----

	double s_th_y_hcut_L, s_th_y_lcut_R, s_th_y_lcut_L, s_th_y_hcut_R;
	if (th_y_sign > 0)
	{
		s_th_y_lcut_L = anal.th_y_lcut_L;
		s_th_y_hcut_L = anal.th_y_hcut_L;
		s_th_y_lcut_R = anal.th_y_lcut_R;
		s_th_y_hcut_R = anal.th_y_hcut_R;
	} else {
		s_th_y_lcut_L = -anal.th_y_hcut_L;
		s_th_y_hcut_L = -anal.th_y_lcut_L;
		s_th_y_lcut_R = -anal.th_y_hcut_R;
		s_th_y_hcut_R = -anal.th_y_lcut_R;
	}

	// range for integration (d has sign)
	double th_y = p[1];
	double d_min = 2. * max(th_y - s_th_y_hcut_L, s_th_y_lcut_R - th_y);
	double d_max = 2. * min(th_y - s_th_y_lcut_L, s_th_y_hcut_R - th_y);

	double range = n_si * si_De_y;
	double UB = +range, LB = -range;

	if (includeAcceptanceEffects)
	{
		UB = min(+range, d_max);
		LB = max(-range, d_min);
	}

	double idy = 0.;
	if (UB > LB)
	{
		idy = RealIntegrate(IntegOverDY, param, NULL, LB, UB, 0., precision,
			int_ws_DY_size, int_ws_DY, "IntegOverMY-idy");
	}
	
	// ----- "divergence" acceptance correction -----

	double div_corr = 1.;

	if (includeAcceptanceEffects)
	{
		double LB_y = d_min / 2.;
		double UB_y = d_max / 2.;

		//printf("%E, %E\n", LB_y, UB_y);

		double F_y = (UB_y > LB_y) ? ( TMath::Erf(UB_y / anal.si_th_y_1arm) - TMath::Erf(LB_y / anal.si_th_y_1arm) ) / 2. : 0.;
		div_corr = 1. / F_y;
		//printf("F_y = %E, div_corr = %E\n", F_y, div_corr);
	}

	return idy * imx * div_corr;
}

//----------------------------------------------------------------------------------------------------

// distribution of th'_x and th'_y after smearing, including acceptance cuts and possibly "div" acceptance correction
double dist_th_x_th_y_sm_acc(double th_x, double th_y)
{
	double param[] = { th_x, th_y };
	double range = n_si * si_m_y;
	double precision = 1E-2;
	return RealIntegrate(IntegOverMY, param, NULL, -range, +range, 0., precision,
			int_ws_MY_size, int_ws_MY, "dist_th_x_th_y_sm");
}

//----------------------------------------------------------------------------------------------------

double IntegOverPhi(double x, double *p, const void *)
{
	double th = p[0];
	double phi = x;

	double th_x = th * cos(phi);
	double th_y = th * sin(phi);

	return dist_th_x_th_y_sm_acc(th_x, th_y);
}

//----------------------------------------------------------------------------------------------------

set<double> GetPhiLimits(double th)
{
	// get all intersections of const-th circle with acceptance boundaries
	set<double> phis;

	if (th > th_y_lcut)
	{
		double phi = asin(th_y_lcut / th);
		double ta_x = th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}
	
	if (th > th_y_hcut)
	{
		double phi = asin(th_y_hcut / th);
		double ta_x = th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}

	if (th > fabs(th_x_hcut))
	{
		double phi = acos(fabs(th_x_hcut) / th);
		double ta_y = th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(phi);
	}

	if (th > fabs(th_x_lcut))
	{
		double phi = acos(fabs(th_x_lcut) / th);
		double ta_y = th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(M_PI - phi);
	}

	return phis;
}

//----------------------------------------------------------------------------------------------------

// "reconstructed" distribution of t', after applying "phi" divergence correction
double dist_t_re(double t)
{
	double th = sqrt(t) / env.p;

	set<double> phis;
	if (includeAcceptanceEffects)
		phis = GetPhiLimits(th);
	else {
		phis.insert(0.);
		phis.insert(2.*M_PI);
	}

	// the number of intersections must be even
	if ((phis.size() % 2) == 1)
	{
		printf("ERROR: odd number of intersections in acceptance calculation\n");
	}

	// no intersection => no acceptances
	if (phis.size() == 0)
		return 0.;

	// calculate integrals over phi sections
	double precision = 1E-2;
	double param[] = { th };

	double phiSum = 0.;
	double integralSum = 0.;
	for (set<double>::iterator it = phis.begin(); it != phis.end(); ++it)
	{
		double phi_start = *it;
		++it;
		double phi_end = *it;

		phiSum += phi_end - phi_start;

		if (th_y_sign == +1)
			integralSum += RealIntegrate(IntegOverPhi, param, NULL, phi_start, phi_end, 0., precision, int_ws_phi_size, int_ws_phi, "dist_t_re");
		else
			integralSum += RealIntegrate(IntegOverPhi, param, NULL, -phi_end, -phi_start, 0., precision, int_ws_phi_size, int_ws_phi, "dist_t_re");
	}

	return integralSum / phiSum;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

vector<TF1 *> parameterizations;

TF1* InitFitFunctions(double norm = 5.)
{
	printf(">> InitFitFunctions(%E)\n", norm);

	TF1 *ff;

	// tried and non-working
	// 	exp3+exp2: too little flexibility at high-|t|
	// 	exp3+exp3: too little flexibility at high-|t|

	ff = new TF1("exp3+exp4", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x) + [4]*exp([5]*x + [6]*x*x + [7]*x*x*x + [8]*x*x*x*x)");
	ff->SetParameters(norm*1E2, -20., 0., 0., norm/1E1, -3., 0., 0., 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);

	// bad results
	/*
	ff = new TF1("exp5+erf*exp2", "[0] * exp(-[1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x) + TMath::Erf((x - [6])/[7]) * [8] * exp(-[9]*x + [10]*x*x)");
	ff->SetParameter(0, 539.);
	ff->SetParameter(1, 20.6);
	ff->SetParameter(2, 13.2);
	ff->SetParameter(3, -46.);
	ff->SetParameter(4, 63.6);
	ff->SetParameter(5, -26.2);

	ff->SetParameter(6, 0.55);
	ff->SetParameter(7, 0.3);

	ff->SetParameter(8, 1.4);
	ff->SetParameter(9, 4.8);
	ff->SetParameter(10, 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);
	*/
	
	ff = new TF1("p1*exp3+p1*exp1", "([0] + [1]*x) * exp([2]*x + [3]*x*x + [4]*x*x*x) + ([5] + [6]*x) * exp([7]*x)");
	ff->SetParameter(0, 544.);
	ff->SetParameter(1, 0.);
	ff->SetParameter(2, -20.1);
	ff->SetParameter(3, 8.37);
	ff->SetParameter(4, -11.1);
	ff->SetParameter(5, -6.99);
	ff->SetParameter(6, 12.85);
	ff->SetParameter(7, -6.17);
	parameterizations.push_back(ff);
	
	ff = new TF1("p1*exp3+p2*exp2", "([0] + [1]*x) * exp([2]*x + [3]*x*x + [4]*x*x*x) + ([5] + [6]*x + [7]*x*x) * exp([8]*x + [9]*x*x)");
	ff->SetParameter(0, 544.);
	ff->SetParameter(1, 0.);
	ff->SetParameter(2, -20.1);
	ff->SetParameter(3, 8.37);
	ff->SetParameter(4, -11.1);
	ff->SetParameter(5, -6.99);
	ff->SetParameter(6, 12.85);
	ff->SetParameter(7, 0.);
	ff->SetParameter(8, -6.17);
	ff->SetParameter(9, 0.);
	parameterizations.push_back(ff);

	ff = new TF1("exp3-intf-exp1", "[1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x) + [5]*[5]*exp(2*[6]*x)");
	ff->SetParameter(0, 2.77);
	ff->SetParameter(1, 24.33);
	ff->SetParameter(2, -9.71);
	ff->SetParameter(3, 4.52);
	ff->SetParameter(4, -3.32);
	ff->SetParameter(5, 1.33);
	ff->SetParameter(6, -2.45);
	parameterizations.push_back(ff);
	
	ff = new TF1("(exp3-intf-exp1)*expG", "( [1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x) + [5]*[5]*exp(2*[6]*x) ) *  exp([7] * exp(-(x-[9])^2/2/[8]^2) )");
	ff->SetParameter(0, 2.77);
	ff->SetParameter(1, 24.33);
	ff->SetParameter(2, -9.71);
	ff->SetParameter(3, 4.52);
	ff->SetParameter(4, -3.32);
	ff->SetParameter(5, 1.33);
	ff->SetParameter(6, -2.45);

	ff->SetParameter(7, -0.25);
	ff->SetParameter(8, 0.022); ff->SetParLimits(8, 0.020, 0.1);
	ff->SetParameter(9, 0.53);
	parameterizations.push_back(ff);

	return ff;
}

//----------------------------------------------------------------------------------------------------

void MakeFit(TH1D *input, TH1D *corr, TF1 *ff)
{
	// correct the input with current correction
	TH1D *input_corr = new TH1D(*input);
	input_corr->Multiply(corr);
	
	input_corr->Fit(ff, "Q",  "", anal.t_min, anal.t_max);
	input_corr->Fit(ff, "Q",  "", anal.t_min, anal.t_max);
	input_corr->Fit(ff, "IQ",  "", anal.t_min, anal.t_max);
	//input_corr->Fit(ff, "I",  "", anal.t_min, anal.t_max);

	printf("residual chi^2 = %.1E, number of points = %i, chi^2/npx = %.2f\n",
		ff->GetChisquare(), ff->GetNpx(), ff->GetChisquare() / ff->GetNpx());

	input_corr->SetName("input_corr");
	input_corr->Write();
	
	ff->Write();
}

//----------------------------------------------------------------------------------------------------

void CalculateCorrection(TF1 *ff, TH1D *corr, int color = 1)
{
	includeAcceptanceEffects = true;

	source_dist_t = ff;

	// sample distributions
	int N_per_bin = 5;

	TGraph *g_tr = new TGraph(); g_tr->SetName("g_tr"); g_tr->SetLineColor(1);
	TGraph *g_re = new TGraph(); g_re->SetName("g_re"); g_re->SetLineColor(2);

	TGraph *g_corr = new TGraph(); g_corr->SetName("g_corr"); g_corr->SetLineColor(color);

	for (int bi = 1; bi <= corr->GetNbinsX(); ++bi)
	{
		double le = corr->GetBinCenter(bi);
		double w = corr->GetBinWidth(bi);

		int limit = N_per_bin;
		if (bi == corr->GetNbinsX())
			limit++;

		for (int i = 0; i < limit; i++)
		{
			double t = le + double(i) / N_per_bin * w;

			double v_true = dist_t_true(t);
			double v_reco = dist_t_re(t);

			int idx = g_tr->GetN();
			g_tr->SetPoint(idx, t, v_true);
			g_re->SetPoint(idx, t, v_reco);

			g_corr->SetPoint(idx, t, (v_reco > 0.) ? v_true / v_reco : 0.);
		}
	}

	TCanvas *c = new TCanvas();
	c->SetName("cmp");
	c->SetLogy(1);
	g_tr->Draw("al");
	g_re->Draw("l");
	c->Write();

	g_corr->Write();

	// calculate correction
	for (int bin = 1; bin <= corr->GetNbinsX(); bin++)
	{
		double l = corr->GetBinLowEdge(bin);
		double w = corr->GetBinWidth(bin);
		//double r = l + w;

		//printf("bin = %3u, l = %.5f, r = %.5f\n", bin, l, r);

		// integrate distributions over the bin
		double i_tr = 0., i_re = 0.;
		unsigned int N_div = 10;
		for (unsigned int di = 0; di < N_div; di++)
		{
			double t = l + w/N_div * (di + 0.5);
			i_tr += g_tr->Eval(t);
			i_re += g_re->Eval(t);
		}

		double v_corr = (i_re > 0.) ? i_tr / i_re : 0.;

		//printf("\t1-corr = %.3E | %.3E\n", 1.-corr, 1.-g_unsm_corr->Eval((l+r)/2));

		corr->SetBinContent(bin, v_corr);
	}

	corr->SetLineColor(color);
	corr->Write();
}

//----------------------------------------------------------------------------------------------------

void SaveDetails(TF1 *ff)
{
	gDirectory = gDirectory->mkdir("details");
	
	source_dist_t = ff;


	TGraph *g_tr = new TGraph(); g_tr->SetName("g_tr"); g_tr->SetLineColor(1);
	TGraph *g_re = new TGraph(); g_re->SetName("g_re"); g_re->SetLineColor(2);
	TGraph *g_id = new TGraph(); g_id->SetName("g_id"); g_id->SetLineColor(4);
	
	TGraph *g_re_to_id = new TGraph(); g_re_to_id->SetName("g_re_to_id"); g_re_to_id->SetLineColor(6);

	for (double t = 0; t < 2.; t += 0.02)
	{
		double v_true = dist_t_true(t);
	
		includeAcceptanceEffects = true;
		double v_reco = dist_t_re(t);

		includeAcceptanceEffects = false;
		double v_reco_id = dist_t_re(t);

		int idx = g_tr->GetN();
		g_tr->SetPoint(idx, t, v_true);
		g_re->SetPoint(idx, t, v_reco);
		g_id->SetPoint(idx, t, v_reco_id);

		g_re_to_id->SetPoint(idx, t, (v_reco_id > 0.) ? v_reco / v_reco_id : 0.);
	}

	TCanvas *c = new TCanvas();
	c->SetName("cmp");
	c->SetLogy(1);
	g_tr->Draw("al");
	g_re->Draw("l");
	g_id->Draw("l");
	c->Write();

	g_re_to_id->Write();
}

//----------------------------------------------------------------------------------------------------

void DoUnfolding(TH1D *h_input, TF1 *parameterization)
{
	// initialize fit function
	TF1 *ff = new TF1(*parameterization);
	ff->SetName("ff");

	// make initial correction histogram
	TH1D *corr = new TH1D(*h_input);
	corr->SetName("corr");
	for (int i = 1; i <= corr->GetNbinsX(); i++)
	{
		corr->SetBinContent(i, 1.);
		corr->SetBinError(i, 0.);
	}

	// save current directory
	TDirectory *topDir = gDirectory;

	// run iterations
	unsigned int N_it = 3;
	int colors[] = { 1, 2, 4, 6, 8 };
	for (unsigned int it = 0; it < N_it; it++)
	{
		printf("* iteration: %i\n", it);
	
		char buf[100];
		sprintf(buf, "iteration %i", it);
		gDirectory = topDir->mkdir(buf);

		MakeFit(h_input, corr, ff);

		CalculateCorrection(ff, corr, colors[it]);
	}

	// save final result
	gDirectory = topDir;

	corr->SetName("corr_final");
	corr->Write();

	// save more information
	SaveDetails(ff);
}

//----------------------------------------------------------------------------------------------------

void RegulariseHighT(TH1D *h_dsdt, TH1D *h_dsdt_comb)
{
	double t_min = 1.;

	int bi_min = h_dsdt->FindBin(t_min);
	for (int bi = bi_min; bi <= h_dsdt->GetNbinsX(); ++bi)
	{
		h_dsdt->SetBinContent(bi, h_dsdt_comb->GetBinContent(bi));
		h_dsdt->SetBinError(bi, h_dsdt_comb->GetBinError(bi));
	}
}

//----------------------------------------------------------------------------------------------------

struct SigmaAdjustment
{
	double x, y;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	if (argc < 2)
	{
		printf("ERROR: Wrong number of parameters.\n");
		return 1;
	}
	
	// binnings
	vector<string> binnings;
	binnings.push_back("ub");
	binnings.push_back("NPB");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30");
	binnings.push_back("bt1");
	binnings.push_back("bt2");

	// fit parametrisations
	InitFitFunctions();

	// sigma adjustments
	vector<SigmaAdjustment> sigmaAdjustments;
	sigmaAdjustments.push_back({0., 0.});
	sigmaAdjustments.push_back({+1., 0.});
	sigmaAdjustments.push_back({0., +1.});
	//sigmaAdjustments.push_back({+1., +1.});	// linearity OK

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	th_x_lcut = anal.th_x_lcut;
	th_x_hcut = anal.th_x_hcut;
	th_y_lcut = anal.th_y_lcut;
	th_y_hcut = anal.th_y_hcut;

	// default smearing sigmas
	double si_th_x_1arm_def = anal.si_th_x_1arm;
	double si_th_y_1arm_def = anal.si_th_y_1arm;

	si_de_th_y_L = si_th_y_1arm_def;
	si_de_th_y_R = si_th_y_1arm_def;

	si_m_x = si_th_x_1arm_def / sqrt(2.);
	si_m_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R) / 2.;
	si_De_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R);

	n_si = 5.;

	printf("n_si = %.f\n", n_si);
	printf("si_m_x = %E\n", si_m_x);
	printf("si_m_y = %E\n", si_m_y);
	printf("si_De_y = %E\n", si_De_y);

	// print info
	printf("\n");
	printf("------------------------------ environment ------------------------------\n");
	env.Print();
	printf("\n");
	printf("------------------------------- analysis --------------------------------\n");
	anal.Print();
	printf("\n");

	// open input file
	string fn_in = string("distributions_") + argv[1] + ".root";
	TFile *f_in = new TFile(fn_in.c_str());
	if (f_in->IsZombie())
	{
		printf("ERROR: can't open file `%s'.\n", fn_in.c_str());
		return 1;
	}
	
	string fn_in_comb = "combine_distributions.root";
	TFile *f_in_comb = TFile::Open(fn_in_comb.c_str());
	assert(f_in_comb != NULL);
	
	// prepare output
	string fn_out = string("unfolding_cf_") + argv[1] + ".root";
	TFile *f_out = new TFile(fn_out.c_str(), "recreate");
	if (f_out->IsZombie())
	{
		printf("ERROR: can't open file `%s' for writing.\n", fn_out.c_str());
		return 3;
	}
	
	// initialise integration workspaces
	int_ws_phi_size = 100;
	int_ws_phi = gsl_integration_workspace_alloc(int_ws_phi_size);

	int_ws_MX_size = 100;
	int_ws_MX = gsl_integration_workspace_alloc(int_ws_MX_size);

	int_ws_MY_size = 100;
	int_ws_MY = gsl_integration_workspace_alloc(int_ws_MY_size);

	int_ws_DY_size = 100;
	int_ws_DY = gsl_integration_workspace_alloc(int_ws_DY_size);

	// run unfolding for all binnings
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		TDirectory *binningDir = f_out->mkdir(binnings[bi].c_str());

		string path = "normalization/" + binnings[bi] + "/h_t_normalized";
		TH1D *h_dsdt = (TH1D *) f_in->Get(path.c_str());
		
		TH1D *h_dsdt_comb = (TH1D *) f_in_comb->Get(path.c_str());

		// regularisation of h_dsdt at high |t|
		RegulariseHighT(h_dsdt, h_dsdt_comb);

		// try different fit parametrisations 
		for (unsigned int pi = 0; pi < parameterizations.size(); pi++)
		{
			TDirectory *paramDir = binningDir->mkdir(parameterizations[pi]->GetName());

			for (vector<SigmaAdjustment>::iterator si = sigmaAdjustments.begin(); si != sigmaAdjustments.end(); ++si)
			{
				char buf[50];
				sprintf(buf, "%+.0f,%+.0f", si->x, si->y);
				TDirectory *sigmaDir = paramDir->mkdir(buf);
				gDirectory = sigmaDir;
				
				printf("\n\n>> binning %s, parametrisation %s, sigma adjustments %+.0f %+.0f\n",
					binnings[bi].c_str(), parameterizations[pi]->GetName(), si->x, si->y);

				// adjust sigmas
				si_m_x = (si_th_x_1arm_def + si->x * anal.si_th_x_1arm_unc) / sqrt(2.);
				si_de_th_y_L = si_de_th_y_R = anal.si_th_y_1arm = si_th_y_1arm_def + si->y * anal.si_th_y_1arm_unc;

				printf("\tanal.si_th_y_1arm = %.3E\n", anal.si_th_y_1arm);
				printf("\tsi_m_x = %.3E\n", si_m_x);

				DoUnfolding(h_dsdt, parameterizations[pi]);

				// performance optimization
				if (binnings[bi] != "bt1")
					break;
			}
		}
	}
		
	delete f_out;

	return 0;
}
