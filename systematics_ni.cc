#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSpline.h"

#include <string>
#include <map>
#include <cstdio>

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
// systematic biases

/**
 * For the calculation of acceptance effects, this program assumes that systematic biases in left
 * and right arms are identical.
 **/

// reconstruction shifts (rad)
double sh_th_x = 0., sh_th_y = 0.;

// reconstruction tilt (rad)
double tilt_th_x_th_y;

// reconstruction scales, th_y_reco = sc_th_y * th_y_true
double sc_th_x = 1., sc_th_y = 1.;

// relative error of beam momentum
double de_rel_p;
		
// uncertainty of the 3-out-of-4 efficiency (plateau intercept and slope)
double eff_intercept_error;
double eff_slope_error;
double eff_slope_fix_point;

// normalisation error
double de_norm;

bool useNonGaussianBeamDivergence;

bool useBiases = false;

//----------------------------------------------------------------------------------------------------

double n_si;	// for range of smearing integrations

double si_de_th_y_L, si_de_th_y_R;
double si_m_x, si_m_y, si_De_y; 

// local copy of cuts (can be modified by systematic scenarios)
double th_x_lcut, th_x_hcut, th_y_lcut, th_y_hcut;

bool useAcceptanceCuts = true;
bool useAcceptanceCorrections = true;

// non-gaussian distribution of De^{R-L} th_y
TF1 *dist_diffRL_th_y = NULL;

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
	if (useBiases && useNonGaussianBeamDivergence)
	{
		return sqrt(2.) * dist_diffRL_th_y->Eval(de_th_y * sqrt(2.));
	} else {
		double si = (useBiases) ? si_de_th_y_L : anal.si_th_y_1arm;
		return exp( - (de_th_y * de_th_y) / 2. / si / si ) / sqrt(2. * M_PI ) / si;
	}
}

//----------------------------------------------------------------------------------------------------

double dist_de_th_y_R(double de_th_y)
{
	if (useBiases && useNonGaussianBeamDivergence)
	{
		return sqrt(2.) * dist_diffRL_th_y->Eval(de_th_y * sqrt(2.));
	} else {
		double si = (useBiases) ? si_de_th_y_R : anal.si_th_y_1arm;
		return exp( - (de_th_y * de_th_y) / 2. / si / si ) / sqrt(2. * M_PI ) / si;
	}
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
	// observed angles
	double th_x_p = p[0];
	double th_y_p = p[1];

	// smearing parameters
	double m_x = x;
	double m_y = p[2];

	// true angles
	double th_x, th_y;

	double corr = 1.;

	if (useBiases)
	{
		double eff_corr = eff_intercept_error  +  eff_slope_error * (fabs(th_y_p) - eff_slope_fix_point);
	
		double de_p_norm_corr = 1. + de_rel_p;
	
		// transformation from reconstructed (primed) to original/true (non-primed) angles
		// Th' = M Th + De Th  ==>  Th = M^-1 (Th' - De The)
		// Mi = M^-1
		double th_x_ps = th_x_p - sh_th_x;
		double th_y_ps = th_y_p - sh_th_y;
	
		double Mi_xx = 1. / sc_th_x, Mi_xy = - tilt_th_x_th_y;
		double Mi_yx = 0., Mi_yy = 1. / sc_th_y;
	
		double D = Mi_xx * Mi_yy - Mi_xy * Mi_yx;
	
		th_x = Mi_xx * th_x_ps + Mi_xy * th_y_ps - m_x;
		th_y = Mi_yx * th_x_ps + Mi_yy * th_y_ps - m_y;

		corr = D * (1. + eff_corr) * (1. + de_norm) / de_p_norm_corr / de_p_norm_corr;
	} else {
		th_x = th_x_p - m_x;
		th_y = th_y_p - m_y;
	}

	return dist_th_x_th_y_true(th_x, th_y) * dist_m_x(m_x) * corr;
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

	if (useAcceptanceCuts)
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

	if (useAcceptanceCorrections)
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

/// returns phis always in the upper half-plane (between 0 and pi)
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
	double p_p = (useBiases) ? env.p*(1. + de_rel_p) : env.p;
	double th = sqrt(t) / p_p;

	set<double> phis;
	if (useAcceptanceCuts)
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

		if (th_y_sign == +1.)
			integralSum += RealIntegrate(IntegOverPhi, param, NULL, phi_start, phi_end, 0., precision, int_ws_phi_size, int_ws_phi, "dist_t_re");
		else
			integralSum += RealIntegrate(IntegOverPhi, param, NULL, -phi_end, -phi_start, 0., precision, int_ws_phi_size, int_ws_phi, "dist_t_re");
	}

	// "phi" acceptance correction, multiplied by Jacobian 1/(2 pi)
	double corr = (useAcceptanceCorrections) ? 1./phiSum : 1./(2.*M_PI);

	return integralSum * corr;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void CalculateOne()
{
	TGraph *g_h = new TGraph(); g_h->SetName("g_h"); g_h->SetLineColor(1);
	TGraph *g_h_p = new TGraph(); g_h_p->SetName("g_h_p"); g_h_p->SetLineColor(2);
	TGraph *g_h_obs = new TGraph(); g_h_obs->SetName("g_h_obs"); g_h_obs->SetLineColor(6);	// what is observed with no (acceptance) correction
	TGraph *g_r = new TGraph(); g_r->SetName("g_r"); g_r->SetLineColor(8);
	
	useAcceptanceCuts = true;

	double de_t = 1.;
	for (double t = 0.02; t <= 2.1; t += de_t)
	{
		useAcceptanceCorrections = true;
		useBiases = false;
		double h = dist_t_re(t);

		useAcceptanceCorrections = true;
		useBiases = true;
		double h_p = dist_t_re(t);

		useAcceptanceCorrections = false;
		useBiases = false;
		double h_obs = dist_t_re(t);

		double r = (h != 0.) ? h_p / h : 0.;

		int idx = g_h->GetN();

		g_h->SetPoint(idx, t, h);
		g_h_p->SetPoint(idx, t, h_p);
		g_h_obs->SetPoint(idx, t, h_obs);
		g_r->SetPoint(idx, t, r);

		de_t = 0.002;
		if (t > 0.2)
			de_t = 0.01;
	}

	g_h->Write();
	g_h_p->Write();
	g_h_obs->Write();
	g_r->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int SetScenario(const string &scenario)
{
	// reset to default cut positions (as in main analysis)
	th_x_lcut = anal.th_x_lcut;
	th_x_hcut = anal.th_x_hcut;
	th_y_lcut = anal.th_y_lcut;
	th_y_hcut = anal.th_y_hcut;

	// default smearing sigmas
	si_de_th_y_L = anal.si_th_y_1arm;
	si_de_th_y_R = anal.si_th_y_1arm;

	// defaults (no systematics)
	sh_th_x = sh_th_y = 0E-6;
	sc_th_x = sc_th_y = 1.;

	tilt_th_x_th_y = 0.;

	de_rel_p = 0.;

	eff_intercept_error = 0.;
	eff_slope_error = 0.;	
	eff_slope_fix_point = 40E-6;

	de_norm = 0.;

	useNonGaussianBeamDivergence = false;

	// set scenarios
	if (scenario.compare("none") == 0)
	{
		return 0;
	}

	if (scenario.compare("alig-sh-thx") == 0)
	{
		// paper: 0.8 urad
		sh_th_x = 0.76E-6;	// rad, 0.377m^-1 * si(de_x_1RP), si(de_x_1RP) = 2 um
		// TODO: adjust cuts ??
		return 0;
	}

	if (scenario.compare("alig-sh-thy") == 0)
	{
		// paper: 0.2 urad
		sh_th_y = 0.2E-6;	// rad
		// TODO: adjust cuts ??
		return 0;
	}

	if (scenario.compare("tilt-thx-thy") == 0)
	{
		// paper: th_x --> th_x + c*th_y, si[c] = 0.02
		tilt_th_x_th_y = 0.02;
		// TODO: adjust cuts ??
		return 0;
	}

	if (scenario.compare("opt-m1") == 0)
	{
		// paper: mode1, x = -0.182 %, y = +0.235 %
		sc_th_x = 1. +1.822E-03;
		sc_th_y = 1. -2.353E-03;
		// TODO: adjust cuts ??
		return 0;
	}

	if (scenario.compare("opt-m2") == 0)
	{
		// paper: mode2, x = -0.096 %, y = -0.074 %
		sc_th_x = 1. +9.579E-04;
		sc_th_y = 1. +7.418E-04;
		// TODO: adjust cuts ??
		return 0;
	}

	if (scenario.compare("eff-slp") == 0)
	{
		// no mentioned in the paper, value from "simu_direct.cc": 16 rad^-1
		//    the factor of sqrt(2) not used here, this will be accounted for in diagonal combination
		// summing in quadrature the uncertainties from "eff3outof4_details_fits.pdf" gives
		//   17.1 and 16.3 rad^-1 for the two diagonals
		eff_slope_error = 16.5;	// rad^-1
		return 0;
	}

	if (scenario.compare("eff-int") == 0)
	{
		// 3/4 efficiencies: 4x 0.1%
		// 2/4 efficiencies: 2x 0.7%
		// summed in quadrature: 1%
		eff_intercept_error = 0.01;
		return 0;
	}

	if (scenario.compare("acc-corr-sigma-unc") == 0)
	{
		// not in the paper, value taken from "simu_direct.cc"
		// compatible with "resolutions_vs_time_per_bunch.pdf"
		double de_si_th_y_1arm = 0.05E-6 / sqrt(2.);

		si_de_th_y_R -= de_si_th_y_1arm;
		si_de_th_y_L -= de_si_th_y_1arm;
		return 0;
	}

	if (scenario.compare("acc-corr-sigma-asym") == 0)
	{
		// paper: 15%, "simu_direct.cc": 10%
		double asymmetry = 0.15;	// 15% difference between L and R resolutions

		si_de_th_y_R = si_de_th_y_R / sqrt(1. + asymmetry);
		si_de_th_y_L = si_de_th_y_L * sqrt(1. + asymmetry);
		return 0;
	}

	if (scenario.compare("acc-corr-non-gauss") == 0)
	{
		useNonGaussianBeamDivergence = true;
		return 0;
	}

	if (scenario.compare("beam-mom") == 0)
	{
		// paper: 0.1 %
		de_rel_p = -0.001;
		return 0;
	}

	if (scenario.compare("norm") == 0)
	{
		// paper: 0.042
		de_norm = 0.042;
		return 0;
	}

	printf("ERROR: unknown scenario `%s'.\n", scenario.c_str());
	return 1;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------

TObject* GetObject(const string &file, const string &path)
{
	TFile *f = TFile::Open(file.c_str());
	if (!f)
	{
		printf("ERROR: can't load file `%s'.\n", file.c_str());
		return NULL;
	}

	TObject *o = f->Get(path.c_str());
	if (!o)
	{
		printf("ERROR: can't load object `%s' from file `%s'.\n", path.c_str(), file.c_str());
		return NULL;
	}

	return o;
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
	if (diagonal == dCombined || diagonal == ad45b_56b || diagonal == ad45t_56t)
		return rcIncompatibleDiagonal;

	// default smearing sigmas
	si_de_th_y_L = anal.si_th_y_1arm;
	si_de_th_y_R = anal.si_th_y_1arm;

	// ranges for integrals
	si_m_x = anal.si_th_x_2arm;
	si_m_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R) / 2.;
	si_De_y = sqrt(si_de_th_y_L*si_de_th_y_L + si_de_th_y_R*si_de_th_y_R);

	n_si = 5.;

	printf("n_si = %.f\n", n_si);
	printf("si_m_x = %E\n", si_m_x);
	printf("si_m_y = %E\n", si_m_y);
	printf("si_De_y = %E\n", si_De_y);

	// distribution for non-gaussian beam divergence
	string expr_diffRL_th_y = "[0] * exp( -x^2 / 2. / [1]^2) + [2] * exp( -x^2 / 2. / [3]^2)";

	TF1 *mom2_diffRL_th_y = new TF1("dist_diffRL_th_y", ("("+expr_diffRL_th_y+") * x*x").c_str());
	mom2_diffRL_th_y->SetParameters(9.15E9/7.853E+04, 2.28E-6, 9.15E9*0.26/7.853E+04, 4.4E-6);
	
	double de_th_y_range = 50E-6;
	double var_diffRL_th_y = mom2_diffRL_th_y->Integral(-de_th_y_range, +de_th_y_range);
	printf("* non-gaussian distribution of De^{R-L} th_y\n");
	printf("\tbefore width adjustment\n");
	printf("\t\tsi_diffRL_th_y = %.3E\n", sqrt(var_diffRL_th_y));
	printf("\t\tanal.si_th_y_1arm * sqrt(2) = %E\n", anal.si_th_y_1arm * sqrt(2.));
	double de_th_y_stretch = anal.si_th_y_1arm * sqrt(2.) / sqrt(var_diffRL_th_y);
	printf("\t\tde_th_y_stretch = %f\n", de_th_y_stretch);

	dist_diffRL_th_y = new TF1("dist_diffRL_th_y", expr_diffRL_th_y.c_str());
	dist_diffRL_th_y->SetParameters(9.15E9/7.853E+04/de_th_y_stretch, 2.28E-6*de_th_y_stretch,
		9.15E9*0.26/7.853E+04/de_th_y_stretch, 4.4E-6*de_th_y_stretch);

	mom2_diffRL_th_y->SetParameters(9.15E9/7.853E+04/de_th_y_stretch, 2.28E-6*de_th_y_stretch,
		9.15E9*0.26/7.853E+04/de_th_y_stretch, 4.4E-6*de_th_y_stretch);

	printf("\tafter width adjustment\n");
	var_diffRL_th_y = mom2_diffRL_th_y->Integral(-de_th_y_range, +de_th_y_range);
	printf("\t\tsi_diffRL_th_y = %.3E\n", sqrt(var_diffRL_th_y));

	double norm_diffRL_th_y = dist_diffRL_th_y->Integral(-de_th_y_range, +de_th_y_range);
	printf("\t\tnorm_diffRL_th_y = %.3E\n", norm_diffRL_th_y);

	// print info
	printf("\n");
	printf("------------------------------ environment ------------------------------\n");
	env.Print();
	printf("\n");
	printf("------------------------------- analysis --------------------------------\n");
	anal.Print();
	printf("\n");

	// scenarios
	vector<string> scenarios;
	scenarios.push_back("none");
	scenarios.push_back("alig-sh-thx");
	scenarios.push_back("alig-sh-thy");
	scenarios.push_back("tilt-thx-thy");
	scenarios.push_back("opt-m1");
	scenarios.push_back("opt-m2");
	scenarios.push_back("eff-slp");
	//scenarios.push_back("eff-int");
	scenarios.push_back("acc-corr-sigma-unc");
	scenarios.push_back("acc-corr-sigma-asym");
	scenarios.push_back("acc-corr-non-gauss");
	scenarios.push_back("beam-mom");
	scenarios.push_back("norm");

	// load input dsigma/dt distribution
	map<string, TSpline*> models;
	string mod_dir = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta90,high_t/models/";
	models["exp3-intf-exp1"] = BuildSpline((TGraph *) GetObject(mod_dir + "exp3-intf-exp1.root", "g_dsdt"));
	models["p1*exp3+p2*exp2"] = BuildSpline((TGraph *) GetObject(mod_dir + "p1*exp3+p2*exp2.root", "g_dsdt"));
	
	// prepare output
	string fn_out = string("systematics_ni_") + argv[1] + ".root";
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

	// evaluate effects
	for (unsigned int sci = 0; sci < scenarios.size(); sci++)
	{
		const string &scenario = scenarios[sci];

		if (SetScenario(scenario) != 0)
			continue;
	
		printf("\tscenario: %s\n", scenario.c_str());
		TDirectory *sc_dir = f_out->mkdir(scenario.c_str());

		for (map<string, TSpline*>::iterator mit = models.begin(); mit != models.end(); ++mit)
		{
			printf("\t\tmodel: %s\n", mit->first.c_str());
			TDirectory *mod_dir = sc_dir->mkdir(mit->first.c_str());
			gDirectory = mod_dir;

			// sample distributions
			source_dist_t = mit->second;
			CalculateOne();
		}
	}
		
	delete f_out;

	return 0;
}
