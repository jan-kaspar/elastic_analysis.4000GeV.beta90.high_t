#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSpline.h"

#include <cmath>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------

TSpline *dsdt_spline = NULL;

double dsdt_model(double t)
{
	if (dsdt_spline)
		return dsdt_spline->Eval(t);
	{
		double B = 20.;
		return exp(-B * t);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	double th_y_sign = -1.;

	double gamma = 0.02;
	
	double p = 4E3;

	unsigned int N_ev = (unsigned int) 1E8;

	//string model = "exp";
	string model = "mod";

	string label = "";

	// get input
	if (model.compare("exp") == 0)
		dsdt_spline = NULL;
	else {
		string mod_dir = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta90,high_t/models/";
		TFile *f_in = TFile::Open((mod_dir + "exp3-intf-exp1.root").c_str());
		dsdt_spline = BuildSpline((TGraph *) f_in->Get("g_dsdt"));
	}

	printf("dsdt_spline = %p\n", dsdt_spline);

	// prepare output
	char buf[100];
	sprintf(buf, "simu,%s,N_ev=%.1E,ga=%+.3f,thysign=%+.0f%s.root", model.c_str(), (double) N_ev, gamma, th_y_sign, label.c_str());
	TFile *f_out = new TFile(buf, "recreate");

	TH1D *h_t_true = new TH1D("h_t_true", "", 100, 0., 2.);
	TH1D *h_t_reco = new TH1D("h_t_reco", "", 100, 0., 2.);

	TProfile *p_de_t_vs_t = new TProfile("p_de_t_vs_t", ";|t|;t' - t", 100, 0., 2.);
	TProfile2D *p_de_t_vs_th_x_th_y = new TProfile2D("p_de_t_vs_th_x_th_y", ";#theta_{x};#theta_{y};t' - t", 100, -200E-6, +200E-6, 100, -200E-6, +200E-6);
	TProfile2D *p_coef_vs_th_x_th_y = new TProfile2D("p_coef_vs_th_x_th_y", ";#theta_{x};#theta_{y};coef", 100, -200E-6, +200E-6, 100, -200E-6, +200E-6);

	// event loop
	for (unsigned int ev = 0; ev < N_ev; ev++)
	{
		// generate true event
		double t = gRandom->Rndm() * 2.;
		double w = dsdt_model(t);
		double phi = gRandom->Rndm() * M_PI * th_y_sign;

		double th = sqrt(t) / p;
		double th_x = th * cos(phi);
		double th_y = th * sin(phi);

		// simulate systematic effect
		double th_x_p = th_x * 1. + th_y * gamma;
		double th_y_p = th_x * 0. + th_y * 1.;

		double t_p = p*p * (th_x_p*th_x_p + th_y_p*th_y_p);

		// acceptance
		double th_y_min = 30E-6;
		double th_y_max = 100E-6;
		double th_x_min = -10000E-6;
		double th_x_max = +10000E-6;

		bool acc_skip = (fabs(th_y) < th_y_min || fabs(th_y) > th_y_max || th_x < th_x_min || th_x > th_x_max);
		
		bool acc_skip_p = (fabs(th_y_p) < th_y_min || fabs(th_y_p) > th_y_max || th_x_p < th_x_min || th_x_p > th_x_max);

		// fill histograms
		if (!acc_skip)
		{
			h_t_true->Fill(t, w);

			p_de_t_vs_t->Fill(t, t_p - t);
			p_de_t_vs_th_x_th_y->Fill(th_x, th_y, t_p - t);
			
			if (fabs(th_x) > 5E-6 && fabs(th_y) > 5E-6)
				p_coef_vs_th_x_th_y->Fill(th_x, th_y, (t_p - t) / th_x / th_y);
		}

		if (!acc_skip_p)
			h_t_reco->Fill(t_p, w);	
	} 

	TH1D *h_t_ratio = new TH1D(*h_t_true);
	h_t_ratio->SetName("h_t_ratio");
	for (int bi = 1; bi <= h_t_ratio->GetNbinsX(); bi++)
	{
		double n = h_t_reco->GetBinContent(bi);
		double n_u = h_t_reco->GetBinError(bi);

		double d = h_t_true->GetBinContent(bi);

		double r = (d > 0.) ? n / d : 0.;
		double r_u = (d > 0.) ? n_u / d : 0.;

		h_t_ratio->SetBinContent(bi, r);
		h_t_ratio->SetBinError(bi, r_u);
	}

	h_t_true->SetLineColor(2);
	h_t_reco->SetLineColor(4);

	h_t_true->Sumw2();
	h_t_reco->Sumw2();

	h_t_true->Write();
	h_t_reco->Write();

	h_t_ratio->Write();

	p_de_t_vs_t->Write();
	p_de_t_vs_th_x_th_y->Write();
	p_coef_vs_th_x_th_y->Write();

	delete f_out;
	return 0;
}
