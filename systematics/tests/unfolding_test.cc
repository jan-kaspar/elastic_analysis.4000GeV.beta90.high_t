#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSpline.h"

#include <string>
#include <map>
#include <cstdio>

#include "../../common_definitions.h"
#include "../../common_algorithms.h"
#include "../../DS4/parameters.h"
#include "../../common.h"

#include "../../NumericalIntegration.h"


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

TSpline *source_dist_t = NULL;

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

// distribution of t' after smearing
double dist_t_sm(double t)
{
	// TODO: phi integration of dist_th_x_th_y_sm
	return 0.;
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

void CalculateCorrection()
{
	includeAcceptanceEffects = true;

	// sample distributions
	TGraph *g_tr = new TGraph(); g_tr->SetName("g_tr"); g_tr->SetLineColor(1);
	TGraph *g_re = new TGraph(); g_re->SetName("g_re"); g_re->SetLineColor(2);

	TGraph *g_corr = new TGraph(); g_corr->SetName("g_corr"); g_corr->SetLineColor(4);

	for (double t = 0; t <= 2.1; t += 0.01)
	{
		double v_true = dist_t_true(t);
		double v_reco = dist_t_re(t);

		int idx = g_tr->GetN();
		g_tr->SetPoint(idx, t, v_true);
		g_re->SetPoint(idx, t, v_reco);

		g_corr->SetPoint(idx, t, (v_reco > 0.) ? v_true / v_reco : 0.);
		
	}

	TCanvas *c = new TCanvas();
	c->SetName("cmp");
	c->SetLogy(1);
	g_tr->Draw("al");
	g_re->Draw("l");
	c->Write();

	g_corr->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// init diagonal settings
	Init("45b_56t");
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// make sure that this is compatible with common_algorithms.h !!!
	th_x_lcut = -1E3;
	th_x_hcut = +1E3;
	th_y_lcut = max(anal.th_y_lcut_L, anal.th_y_lcut_R) + 1.5E-6;
	th_y_hcut = min(anal.th_y_hcut_L, anal.th_y_hcut_R) - 2.0E-6;

	// adjust smearing sigmas
	si_de_th_y_L = anal.si_th_y_1arm;
	si_de_th_y_R = anal.si_th_y_1arm;

	si_m_x = anal.si_th_x_2arm;
	si_m_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R) / 2.;
	si_De_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R);

	n_si = 5.;

	printf("n_si = %.f\n", n_si);
	printf("si_m_x = %E\n", si_m_x);
	printf("si_m_y = %E\n", si_m_y);
	printf("si_De_y = %E\n", si_De_y);

	/*
	anal.si_th_x_2arm += 0. * anal.si_th_x_2arm_unc;
	anal.si_th_y_1arm += 0. * anal.si_th_y_1arm_unc;
	*/

	// print info
	/*
	printf("\n");
	printf("------------------------------ environment ------------------------------\n");
	env.Print();
	printf("\n");
	printf("------------------------------- analysis --------------------------------\n");
	anal.Print();
	printf("\n");
	*/

	// open input file
	string fn_in = "test.root";
	TFile *f_in = new TFile(fn_in.c_str());
	if (f_in->IsZombie())
	{
		printf("ERROR: can't open file `%s'.\n", fn_in.c_str());
		return 1;
	}

	// get input graphs
	map<string, TGraph *> graphs;
	graphs["none"] = (TGraph *) f_in->Get("none/data fit 1/g_h_p");
	//graphs["beam-mom1"] = (TGraph *) f_in->Get("beam-mom1/data fit 1/g_h_p");
	graphs["beam-mom3"] = (TGraph *) f_in->Get("beam-mom3/data fit 1/g_h_p");
	graphs["opt-m1*3"] = (TGraph *) f_in->Get("opt-m1*3/data fit 1/g_h_p");
	graphs["opt-m2*3"] = (TGraph *) f_in->Get("opt-m2*3/data fit 1/g_h_p");
	
	// prepare output
	string fn_out = "unfolding_test.root";
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
	for (map<string, TGraph *>::iterator it = graphs.begin(); it != graphs.end(); ++it)
	{
		printf(">> %s : %p\n", it->first.c_str(), it->second);

		TDirectory *scDir = f_out->mkdir(it->first.c_str());
		gDirectory = scDir;

		TGraph *g = it->second;
		source_dist_t = new TSpline3("", g->GetX(), g->GetY(), g->GetN());

		CalculateCorrection();	
	}


		
	delete f_out;

	return 0;
}
