#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

//----------------------------------------------------------------------------------------------------

double B = 20.;
double p = 4E3;

double th_y_lcut = 33.5E-6;
double th_y_hcut = 103E-6;

//----------------------------------------------------------------------------------------------------

double cs_obs(double *x, double *par)
{
	double t = x[0];

	double th = sqrt(t) / p;

	if (th < th_y_lcut)
		return 0.;

	double phi_max = (th > th_y_hcut) ? asin(th_y_hcut / th) : M_PI/2.;
	double phi_min = asin(th_y_lcut / th);
	double A = (phi_max - phi_min) * 2. / M_PI;

	double cs = exp(-B*t);

	return A * cs;
}

//----------------------------------------------------------------------------------------------------

void const_N_binning()
{
	// --------------------

	TF1 *f_cs_obs = new TF1("f_cs_obs", cs_obs, 0, 1, 0);

	f_cs_obs->SetNpx(1000);
	f_cs_obs->Draw();
	
	// --------------------

	double t_min = pow(th_y_lcut * p, 2);	
	double t_max = 0.2;

	TGraph *g_cs_obs_int = new TGraph();
	double S = 0.;
	double dt = 0.0001;
	for (double t = t_min; t <= 1.; t += dt)
	{
		int idx = g_cs_obs_int->GetN();
		g_cs_obs_int->SetPoint(idx, t, S);

		double I = f_cs_obs->Integral(t, t+dt);
		S += I;
	}

	new TCanvas();
	g_cs_obs_int->Draw("alp");
	
	// --------------------

	double cs_int_min = g_cs_obs_int->Eval(t_min);
	double cs_int_max = g_cs_obs_int->Eval(t_max);

	double p_max;

	TGraph *g_int_inv = new TGraph();
	for (int i = 0; i < g_cs_obs_int->GetN(); i++)
	{
		double t, I;
		g_cs_obs_int->GetPoint(i, t, I);
		double Ir = (I - cs_int_min) / (cs_int_max - cs_int_min);

		g_int_inv->SetPoint(i, Ir, t);
		p_max = Ir;
	}

	new TCanvas();
	g_int_inv->Draw("alp");
	
	// --------------------

	for (double t = 0; t < t_min; t += t_min/10)
	{
		printf("%.5f\n", t);
	}

	printf("\n");

	double de_p = 0.001;	// 1/N_bins
	for (double p = 0.; p <= p_max; p += de_p)
	{
		double t = g_int_inv->Eval(p);

		//printf("p = %.3f, t = %.4f\n", p, t);
		printf("%.5f\n", t);
	}
}
