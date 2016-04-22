#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

#include <string>
#include <vector>
#include <set>

#include "numerical_integration_common.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

double p = 4E3;	// GeV

double th_y_sign;

// acceptance cuts in th_x and th_y (rad)
double th_x_lcut = 0.;
double th_x_hcut = 0.;
double th_y_lcut = 0.;
double th_y_hcut = 0.;

// reconstruction shifts (rad)
double sh_th_x = 0., sh_th_y = 0.;

// reconstruction tilt (rad)
double thx_thy_tilt;

// reconstruction scales, th_y_reco = sc_th_y * th_y_true
double sc_th_x = 1., sc_th_y = 1.;

// relative error of beam momentum
double de_rel_p;
		
// uncertainty of the 3-out-of-4 efficiency plateau slope
double eff_slope_error;
double eff_slope_fix_point;

double de_norm;

//----------------------------------------------------------------------------------------------------

double dist_true_th_x_th_y(double th_x, double th_y)
{
	double t = (th_x*th_x + th_y*th_y) * p*p;
	return dist_true_t(t);
}

//----------------------------------------------------------------------------------------------------

bool use_reco_dist_th_x_th_y;

double dist_reco_th_x_th_y(double th_x_p, double th_y_p)
{
	double eff_slp_corr = eff_slope_error * (fabs(th_y_p) - eff_slope_fix_point);

	double de_p_norm_corr = (use_reco_dist_th_x_th_y) ? (1. + de_rel_p) : 1.;

	// transformation from reconstructed (primed) to original/true (non-primed) angles
	// Th' = M Th + De Th  ==>  Th = M^-1 (Th' - De The)
	// Mi = M^-1
	double th_x_ps = th_x_p - sh_th_x;
	double th_y_ps = th_y_p - sh_th_y;

	double Mi_xx = 1. / sc_th_x, Mi_xy = -thx_thy_tilt;
	double Mi_yx = 0., Mi_yy = 1. / sc_th_y;

	double D = Mi_xx * Mi_yy - Mi_xy * Mi_yx;

	double th_x = Mi_xx * th_x_ps + Mi_xy * th_y_ps;
	double th_y = Mi_yx * th_x_ps + Mi_yy * th_y_ps;

	return dist_true_th_x_th_y(th_x, th_y)
		* D
		* (1. + eff_slp_corr)
		* (1. + de_norm)
		/ de_p_norm_corr / de_p_norm_corr;
}

//----------------------------------------------------------------------------------------------------

double IntegOverPhi(double x[], double par[])
{
	double phi = x[0];
	double th_p = par[0];

	double th_x_p = th_p * cos(phi);
	double th_y_p = th_p * sin(phi);

	double result = (use_reco_dist_th_x_th_y) ? dist_reco_th_x_th_y(th_x_p, th_y_p) : dist_true_th_x_th_y(th_x_p, th_y_p);
	//printf("phi = %.3f => %.5E\n", phi, result);

	return result;
}

//----------------------------------------------------------------------------------------------------

bool applyAcceptanceCorrection = false;

double dist_reco_t(double t_p)
{
	double p_p = (use_reco_dist_th_x_th_y) ? p*(1. + de_rel_p) : p;
	double th_p = sqrt(t_p) / p_p;
	
	double param[1] = { th_p };

	// get all intersections of const-th circle with acceptance boundaries
	set<double> phis;

	if (th_p > th_y_lcut)
	{
		double phi = asin(th_y_lcut / th_p);
		double ta_x = th_p * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}
	
	if (th_p > th_y_hcut)
	{
		double phi = asin(th_y_hcut / th_p);
		double ta_x = th_p * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}

	if (th_p > fabs(th_x_hcut))
	{
		double phi = acos(fabs(th_x_hcut) / th_p);
		double ta_y = th_p * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(phi);
	}

	if (th_p > fabs(th_x_lcut))
	{
		double phi = acos(fabs(th_x_lcut) / th_p);
		double ta_y = th_p * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(M_PI - phi);
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
	double precision = 1E-5;
	double phiSum = 0.;
	double integralSum = 0.;
	for (set<double>::iterator it = phis.begin(); it != phis.end(); ++it)
	{
		double phi_start = *it;
		++it;
		double phi_end = *it;

		phiSum += phi_end - phi_start;

		if (th_y_sign == +1)
			integralSum += DoubleInt(IntegOverPhi, phi_start, phi_end, param, precision);
		else
			integralSum += DoubleInt(IntegOverPhi, -phi_end, -phi_start, param, precision);
	}

	if (applyAcceptanceCorrection)
		return integralSum / phiSum;
	else
		return integralSum;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

enum DSType { dsUnknown, ds2a, ds2b } ds = dsUnknown; 
enum DGNType { dgnUnknown, dgn_45b, dgn_45t } dgn = dgnUnknown;

//----------------------------------------------------------------------------------------------------

int SetScenario(const string &scenario)
{
	// defaults (no systematics)
	sh_th_x = sh_th_y = 0E-6;
	sc_th_x = sc_th_y = 1.;

	thx_thy_tilt = 0.;

	de_rel_p = 0.;

	eff_slope_error = 0.;	
	eff_slope_fix_point = 0E-6;

	de_norm = 0.;

	if (scenario.compare("none") == 0)
	{
		return 0;
	}

	if (scenario.compare("alig-sh-thx") == 0)
	{
		sh_th_x = 0.28E-6;
		return 0;
	}

	if (scenario.compare("alig-sh-thy:D+0,R+1") == 0)
	{
		sh_th_y = 0.18E-6;
		return 0;
	}

	if (scenario.compare("alig-sh-thy:D+1,R+0") == 0)
	{
		sh_th_y = -0.05E-6 / 2.;
		return 0;
	}
	
	if (scenario.compare("thx-thy-tilt") == 0)
	{
		thx_thy_tilt = 0.005;
		return 0;
	}

	if (scenario.compare("opt-m1") == 0)
	{
		sc_th_x = 1. - -1.822E-03;
		sc_th_y = 1. - +2.353E-03;

		return 0;
	}

	if (scenario.compare("opt-m1*3") == 0)
	{
		sc_th_x = 1. - -1.822E-03*3.;
		sc_th_y = 1. - +2.353E-03*3.;

		return 0;
	}

	if (scenario.compare("opt-m2") == 0)
	{
		sc_th_x = 1. - -9.579E-04;
		sc_th_y = 1. - -7.418E-04;

		return 0;
	}

	if (scenario.compare("opt-m2*3") == 0)
	{
		sc_th_x = 1. - -9.579E-04*3.;
		sc_th_y = 1. - -7.418E-04*3.;

		return 0;
	}

	if (scenario.compare("eff-slp") == 0)
	{
		eff_slope_fix_point = 40E-6;		// |th_y|, rad
		if (ds == ds2a)
			eff_slope_error = 90.;			// 1/rad
		if (ds == ds2b)
			eff_slope_error = 60.;			// 1/rad

		return 0;
	}

	if (scenario.compare("beam-mom-1") == 0)
	{
		de_rel_p = -0.001 * -1.;
		return 0;
	}

	if (scenario.compare("beam-mom1") == 0)
	{
		de_rel_p = -0.001;
		return 0;
	}

	if (scenario.compare("beam-mom2") == 0)
	{
		de_rel_p = -0.001 * 2.;
		return 0;
	}

	if (scenario.compare("beam-mom3") == 0)
	{
		de_rel_p = -0.001 * 3.;
		return 0;
	}

	if (scenario.compare("norm") == 0)
	{
		de_norm = 0.042;

		return 0;
	}

	printf("ERROR: unknown scenario `%s'.\n", scenario.c_str());
	return 1;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int Setup(const string &s_dataset, const string &s_diagonal)
{
	// parse dataset and diagonal
	if (s_dataset.compare("DS2a") == 0) ds = ds2a;
	if (s_dataset.compare("DS2b") == 0) ds = ds2b;

	if (ds == dsUnknown)
	{
		printf("ERROR: dataset `%s' not recognized.\n", s_dataset.c_str());
		return 1;
	}

	if (s_diagonal.compare("45b_56t") == 0) dgn = dgn_45b;
	if (s_diagonal.compare("45t_56b") == 0) dgn = dgn_45t;

	if (dgn == dgnUnknown)
	{
		printf("ERROR: diagonal `%s' not recognized.\n", s_diagonal.c_str());
		return 2;
	}

	//--------------------------------------------------

	// sign of th_y
	if (dgn == dgn_45b)
		th_y_sign = +1.;
	else
		th_y_sign = -1.;
	
	th_x_lcut = -10000E-6;
	th_x_hcut = +10000E-6;

	if (dgn == dgn_45b)
	{
		th_y_lcut = 13.8E-6 + 0.2E-6;
		th_y_hcut = 100E-6 - 1.0E-6;
	}

	if (dgn == dgn_45t)
	{
		th_y_lcut = 5.8E-6 + 0.2E-6;
		th_y_hcut = 100E-6 - 1.0E-6;
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: numerical_integration <dataset> <diagonal> <output file name>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	if (argc != 4)
	{
		printf("ERROR: wrong number of parameters\n");
		PrintUsage();
		return 1;
	}

	// get command-line arguments
	string dataset = argv[1];
	string diagonal = argv[2];
	string outFileName = argv[3];

	// set-up for chosen dataset and diagonal
	if (Setup(dataset, diagonal) != 0)
	{
		printf("ERROR: Setup failed.\n");
		return 2;
	}

	// scenarios
	vector<string> scenarios;
	scenarios.push_back("none");
	scenarios.push_back("alig-sh-thx");
	scenarios.push_back("alig-sh-thy:D+0,R+1");
	scenarios.push_back("alig-sh-thy:D+1,R+0");
	scenarios.push_back("thx-thy-tilt");
	scenarios.push_back("opt-m1");
	scenarios.push_back("opt-m1*3");
	scenarios.push_back("opt-m2");
	scenarios.push_back("opt-m2*3");
	scenarios.push_back("eff-slp");
	scenarios.push_back("beam-mom-1");
	scenarios.push_back("beam-mom1");
	scenarios.push_back("beam-mom2");
	scenarios.push_back("beam-mom3");
	scenarios.push_back("norm");

	// load input dsigma/dt distribution
	if (LoadTDistributions() != 0)
		return 1;

	// prepare output
	TFile *f_out = new TFile(outFileName.c_str(), "recreate");

	WriteTDistributions();

	// evaluate effects
	for (unsigned int sci = 0; sci < scenarios.size(); sci++)
	{
		const string &scenario = scenarios[sci];

		if (SetScenario(scenario) != 0)
			continue;
	
		printf("\tmode: %s\n", scenario.c_str());
		TDirectory *mode_dir = f_out->mkdir(scenario.c_str());

		for (unsigned int dti = 2; dti < 3; dti++)
		{
			t_dist_type = dti;
			string dist_label = TDistTypeName(t_dist_type);
		
			printf("\t\tdist: %s\n", dist_label.c_str());
			TDirectory *dist_dir = mode_dir->mkdir(dist_label.c_str());
			gDirectory = dist_dir;

			// sample distributions
			TGraph *g_h = new TGraph(); g_h->SetName("g_h"); g_h->SetLineColor(1);
			TGraph *g_h_p = new TGraph(); g_h_p->SetName("g_h_p"); g_h_p->SetLineColor(2);
			TGraph *g_h_obs = new TGraph(); g_h_obs->SetName("g_h_obs"); g_h_obs->SetLineColor(6);
			TGraph *g_r = new TGraph(); g_r->SetName("g_r"); g_r->SetLineColor(8);
			for (double t = 0.02; t <= 2.1; t += 0.001)
			{
				applyAcceptanceCorrection = true;
				use_reco_dist_th_x_th_y = false;
				double h = dist_reco_t(t);

				applyAcceptanceCorrection = true;
				use_reco_dist_th_x_th_y = true;
				double h_p = dist_reco_t(t);

				applyAcceptanceCorrection = false;
				use_reco_dist_th_x_th_y = false;
				double h_obs = dist_reco_t(t);

				double r = (h != 0.) ? h_p / h : 0.;

				int idx = g_h->GetN();

				g_h->SetPoint(idx, t, h);
				g_h_p->SetPoint(idx, t, h_p);
				g_h_obs->SetPoint(idx, t, h_obs);
				g_r->SetPoint(idx, t, r);
			}

			g_h->Write();
			g_h_p->Write();
			g_h_obs->Write();
			g_r->Write();
		}
	}


	delete f_out;

	return 0;
}
