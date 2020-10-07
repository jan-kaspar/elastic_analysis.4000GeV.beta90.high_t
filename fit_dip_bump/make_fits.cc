#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

#include <TDirectory.h>
#include <memory>
#include <vector>
#include <map>

using namespace std;

//----------------------------------------------------------------------------------------------------

std::unique_ptr<TH1D> ApplySystematicMode(const TH1D *h_input, const TGraph *g_mode)
{
	std::unique_ptr<TH1D> h_output(new TH1D(*h_input));

	for (int bi = 1; bi < h_input->GetNbinsX(); ++bi)
	{
		const double t = h_input->GetBinCenter(bi);
		const double v = h_input->GetBinContent(bi);

		const double f = 1. + g_mode->Eval(t);

		h_output->SetBinContent(bi, f * v);
	}

	return h_output;
}

//----------------------------------------------------------------------------------------------------

struct FitResults
{
	double t_dip=0., t_dip_unc=0.;
	double dsdt_dip=0., dsdt_dip_unc=0.;
	double t_bmp=0., t_bmp_unc=0.;
	double dsdt_bmp=0., dsdt_bmp_unc=0.;
	double R=0., R_unc=0.;

	void Print() const
	{
		printf("  dip: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_dip, t_dip_unc, dsdt_dip, dsdt_dip_unc);
		printf("  bump: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_bmp, t_bmp_unc, dsdt_bmp, dsdt_bmp_unc);
		printf("  R = dsdt_bump / dsdt_dip = %.3f +- %.3f\n", R, R_unc);
	}
};

//----------------------------------------------------------------------------------------------------

FitResults MakeFit(TH1D *h_dsdt, const std::string &model, bool save=true)
{
	FitResults r;

	if (model == "local")
	{
		unique_ptr<TF1> ff2(new TF1("ff2", "[0] + [1] * pow(x - [2], 2)"));
		ff2->SetParameters(0.0147, 1., 0.523);
		h_dsdt->Fit(ff2.get(), "Q", "", 0.47, 0.57);
		if (save)
			ff2->Write("f_dip");

		r.t_dip = ff2->GetParameter(2);
		r.t_dip_unc = ff2->GetParError(2);
		r.dsdt_dip = ff2->GetParameter(0);
		r.dsdt_dip_unc = ff2->GetParError(0);

		unique_ptr<TF1> ff3(new TF1("ff3", "[0] + [1] * pow(x - [2], 2)"));
		ff3->SetParameters(0.0295, -1., 0.69);
		h_dsdt->Fit(ff3.get(), "Q", "", 0.55, 0.85);
		if (save)
			ff3->Write("f_bump");

		r.t_bmp = ff3->GetParameter(2);
		r.t_bmp_unc = ff3->GetParError(2);
		r.dsdt_bmp = ff3->GetParameter(0);
		r.dsdt_bmp_unc = ff3->GetParError(0);

		r.R = r.dsdt_bmp / r.dsdt_dip;
		r.R_unc = r.R * sqrt(pow(r.dsdt_dip_unc / r.dsdt_dip, 2.) + pow(r.dsdt_bmp_unc / r.dsdt_bmp, 2.));

		return r;
	}

	printf("ERROR in MakeFit: model `%s` not known.\n", model.c_str());

	return r;
}

//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> models = {
		"local"
	};

	// get input
	TFile *f_in = TFile::Open("../DS-merged/merged.root");
	TH1D *h_dsdt = (TH1D *) f_in->Get("bt1/DS4-sc/combined/h_dsdt");

	// load systematics
	map<string, TGraph*> syst_modes;
	unique_ptr<TFile> f_in_syst(TFile::Open("../DS4/systematics_matrix.root"));
	for (const string &m : {"alig-sh-thx", "alig-sh-thy", "tilt-thx-thy", "opt-m1",
		"opt-m2", "acc-corr-sigma-unc", "acc-corr-sigma-asym", "acc-corr-non-gauss", "eff-slp", "beam-mom", "norm"})
	{
		TGraph *c1 = (TGraph *) f_in_syst->Get(("contributions/" + m + "/g_eff_comb1").c_str());
		if (c1)
			syst_modes[m + ":1"] = c1;

		TGraph *c2 = (TGraph *) f_in_syst->Get(("contributions/" + m + "/g_eff_comb2").c_str());
		if (c2)
			syst_modes[m + ":2"] = c2;
	}

	// prepare output
	TFile *f_out = TFile::Open("make_fits.root", "recreate");

	h_dsdt->Write("h_dsdt");

	// make fits for all models
	for (const auto &model : models)
	{
		printf("------ %s ------\n", model.c_str());

		TDirectory *d_model = f_out->mkdir(model.c_str());

		// central fit
		gDirectory = d_model->mkdir("central");
		FitResults r_central = MakeFit(h_dsdt, model);
		printf("* central fit:\n");
		r_central.Print();

		// systematics
		vector<double> vsum(5, 0.);
		for (const auto &syst_mode : syst_modes)
		{
			gDirectory = d_model->mkdir(syst_mode.first.c_str());
			auto h_dsdt_mod = ApplySystematicMode(h_dsdt, syst_mode.second);
			FitResults r_syst = MakeFit(h_dsdt_mod.get(), model);

			vsum[0] += pow(r_syst.t_dip - r_central.t_dip, 2);
			vsum[1] += pow(r_syst.dsdt_dip - r_central.dsdt_dip, 2);
			vsum[2] += pow(r_syst.t_bmp - r_central.t_bmp, 2);
			vsum[3] += pow(r_syst.dsdt_bmp - r_central.dsdt_bmp, 2);
			vsum[4] += pow(r_syst.R - r_central.R, 2);
		}

		FitResults r_syst_unc;
		r_syst_unc.t_dip_unc = sqrt(vsum[0]);
		r_syst_unc.dsdt_dip_unc = sqrt(vsum[1]);
		r_syst_unc.t_bmp_unc = sqrt(vsum[2]);
		r_syst_unc.dsdt_bmp_unc = sqrt(vsum[3]);
		r_syst_unc.R_unc = sqrt(vsum[4]);

		printf("* systematics:\n");
		r_syst_unc.Print();
	}

	// clean up
	delete f_out;
	delete f_in;

	return 0;
}