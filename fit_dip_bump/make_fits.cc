#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <TDirectory.h>
#include <TFitResultPtr.h>
#include <TRandom3.h>
#include "TMatrixDSymEigen.h"
#include "TFitResult.h"

#include <cstdio>
#include <memory>
#include <vector>
#include <map>

#include "../stat.h"

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

	bool valid = true;

	void Print() const
	{
		printf("  dip: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_dip, t_dip_unc, dsdt_dip, dsdt_dip_unc);
		printf("  bump: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_bmp, t_bmp_unc, dsdt_bmp, dsdt_bmp_unc);
		printf("  R = dsdt_bump / dsdt_dip = %.3f +- %.3f\n", R, R_unc);
	}
};

//----------------------------------------------------------------------------------------------------

void FindMinimum(TF1 *f, double x_min, double x_max, double &x_found, double &y_found)
{
	y_found = 1E100;

	for (double x = x_min; x <= x_max; x += 0.0001)
	{
		const double y = f->Eval(x);
		if (y < y_found)
		{
			x_found = x;
			y_found = y;
		}
	}
}

//----------------------------------------------------------------------------------------------------

void FindMaximum(TF1 *f, double x_min, double x_max, double &x_found, double &y_found)
{
	y_found = -1E100;

	for (double x = x_min; x <= x_max; x += 0.0001)
	{
		const double y = f->Eval(x);
		if (y > y_found)
		{
			x_found = x;
			y_found = y;
		}
	}
}

//----------------------------------------------------------------------------------------------------

struct RangeData
{
	string label;
	vector<double> data;
	bool enabled;
};

//----------------------------------------------------------------------------------------------------

vector<RangeData> GetRanges(const std::string &model)
{
	// NB: first range is the central one

	const vector<string> tags = { "C", "L", "H" };

	vector<RangeData> list;

	if (model == "local")
	{
		const vector<double> dip_min = { 0.470, 0.455, 0.485 };
		const vector<double> dip_max = { 0.555, 0.530, 0.580 };

		const vector<double> bmp_min = { 0.570, 0.540, 0.600 };
		const vector<double> bmp_max = { 0.850, 0.780, 0.900 };

		for (unsigned int dmi = 0; dmi < 3; ++dmi)
			for (unsigned int dma = 0; dma < 3; ++dma)
				for (unsigned int bmi = 0; bmi < 3; ++bmi)
					for (unsigned int bma = 0; bma < 3; ++bma)
					{
						RangeData rd;
						string l_dip = tags[dmi] + tags[dma];
						string l_bmp = tags[bmi] + tags[bma];
						rd.label = l_dip + "," + l_bmp;
						rd.data = { dip_min[dmi], dip_max[dma], bmp_min[bmi], bmp_max[bma] };
						rd.enabled = true;

						if (l_bmp == "HL")
							rd.enabled = false;

						list.push_back(rd);
					}
	}

	if (model == "exp2+exp3")
	{
		//const vector<double> min = { 0.35, 0.33, 0.36 };
		const vector<double> min = { 0.42, 0.41, 0.43 };
		const vector<double> max = { 1.00, 0.90, 1.20 };

		for (unsigned int mi = 0; mi < 3; ++mi)
			for (unsigned int ma = 0; ma < 3; ++ma)
			{
				RangeData rd;
				rd.label = tags[mi] + tags[ma];
				rd.data = { min[mi], max[ma] };
				rd.enabled = true;
				list.push_back(rd);
			}
	}

	if (list.empty())
		printf("ERROR in GetRanges: model `%s` not known.\n", model.c_str());

	return list;
}

//----------------------------------------------------------------------------------------------------

struct FitData
{
	shared_ptr<TF1> ff;
	TFitResultPtr fr;
};

//----------------------------------------------------------------------------------------------------

vector<FitData> MakeFit(TH1D *h_dsdt, const std::string &model, const vector<double> &rangeData)
{
	const bool save = true;

	vector<FitData> fits;

	if (model == "local")
	{
		{
			const double x_min = rangeData[0], x_max = rangeData[1];

			const int idx_min = h_dsdt->FindBin(x_min);
			const int idx_max = h_dsdt->FindBin(x_max);

			const double x_min_edge = h_dsdt->GetBinLowEdge(idx_min);
			const double x_max_edge = h_dsdt->GetBinLowEdge(idx_max) + h_dsdt->GetBinWidth(idx_max);

			shared_ptr<TF1> ff2(new TF1("ff2", "[0] + [1] * pow(x - [2], 2)"));
			ff2->SetParameters(0.0147, 1., 0.523);
			ff2->SetRange(x_min_edge, x_max_edge);

			TFitResultPtr result = h_dsdt->Fit(ff2.get(), "RQS", "");

			// TODO: remove
			h_dsdt->Write("hist_with_fit_dip");

			if (save)
				ff2->Write("f_dip");

			fits.push_back({ff2, result});
		}

		{
			const double x_min = rangeData[2], x_max = rangeData[3];

			const int idx_min = h_dsdt->FindBin(x_min);
			const int idx_max = h_dsdt->FindBin(x_max);

			const double x_min_edge = h_dsdt->GetBinLowEdge(idx_min);
			const double x_max_edge = h_dsdt->GetBinLowEdge(idx_max) + h_dsdt->GetBinWidth(idx_max);

			shared_ptr<TF1> ff3(new TF1("ff3", "[0] + [1] * pow(x - [2], 2)"));
			ff3->SetParameters(0.0295, -1., 0.69);
			ff3->SetRange(x_min_edge, x_max_edge);

			TFitResultPtr result = h_dsdt->Fit(ff3.get(), "RQS", "");

			// TODO: remove
			h_dsdt->Write("hist_with_fit_bmp");

			if (save)
				ff3->Write("f_bump");

			fits.push_back({ff3, result});
		}

		return fits;
	}

	if (model == "exp2+exp3")
	{
		double x_min = rangeData[0], x_max = rangeData[1];

		int idx_min = h_dsdt->FindBin(x_min);
		int idx_max = h_dsdt->FindBin(x_max);

		double x_min_edge = h_dsdt->GetBinLowEdge(idx_min);
		double x_max_edge = h_dsdt->GetBinLowEdge(idx_max) + h_dsdt->GetBinWidth(idx_max);

		shared_ptr<TF1> ff(new TF1("ff", "[0] * exp([1]*(x-0.4) + [2]*(x-0.4)*(x-0.4)) + [3]*exp([4]*(x-0.7) + [5]*(x-0.7)*(x-0.7) + [6]*(x-0.7)*(x-0.7)*(x-0.7))"));
		ff->SetParameters(0.144, -26.2, -31.2, 0.03, 0., -19.5, 39.7);
		ff->SetRange(x_min_edge, x_max_edge);

		TFitResultPtr result = h_dsdt->Fit(ff.get(), "RQS", "");

		// TODO: remove
		h_dsdt->Write("hist_with_fit");

		if (save)
			ff->Write("f_global");

		fits.push_back({ff, result});

		return fits;
	}

	printf("ERROR in MakeFit: model `%s` not known.\n", model.c_str());

	return fits;
}

//----------------------------------------------------------------------------------------------------

FitResults AnalyzeFit(const vector<FitData> &fits, const std::string &model)
{
	//printf("* AnalyzeFit\n");

	FitResults r;

	if (model == "local")
	{
		const auto f_dip = fits[0].ff;
		const auto f_bmp = fits[1].ff;

		r.t_dip = f_dip->GetParameter(2);
		r.t_dip_unc = f_dip->GetParError(2);
		r.dsdt_dip = f_dip->GetParameter(0);
		r.dsdt_dip_unc = f_dip->GetParError(0);

		r.t_bmp = f_bmp->GetParameter(2);
		r.t_bmp_unc = f_bmp->GetParError(2);
		r.dsdt_bmp = f_bmp->GetParameter(0);
		r.dsdt_bmp_unc = f_bmp->GetParError(0);

		r.R = r.dsdt_bmp / r.dsdt_dip;
		r.R_unc = r.R * sqrt(pow(r.dsdt_dip_unc / r.dsdt_dip, 2.) + pow(r.dsdt_bmp_unc / r.dsdt_bmp, 2.));

		return r;
	}

	if (model == "exp2+exp3")
	{
		const auto ff = fits[0].ff;

		// TODO:
		//for (int pi = 0; pi < ff->GetNpar(); ++pi)
		//	printf(" - pi = %u: par = %.2E\n", pi, ff->GetParameter(pi));

		FindMinimum(ff.get(), 0.49, 0.56, r.t_dip, r.dsdt_dip);
		r.t_dip_unc = 0.;
		r.dsdt_dip_unc = 0.;

		FindMaximum(ff.get(), 0.60, 0.80, r.t_bmp, r.dsdt_bmp);
		r.t_bmp_unc = 0.;
		r.dsdt_bmp_unc = 0.;

		r.R = r.dsdt_bmp / r.dsdt_dip;
		r.R_unc = 0.;

		if (fabs(r.dsdt_bmp - 0.03) > 0.02)
			r.valid = false;

		return r;
	}

	printf("ERROR in AnalyzeFit: model `%s` not known.\n", model.c_str());

	return r;
}

//----------------------------------------------------------------------------------------------------

void GeneratePerturbation(const vector<FitData> &input, const vector<shared_ptr<TMatrixD>> &genMat, vector<FitData> &output,
	vector<shared_ptr<Stat>> &stats)
{
	//printf("* GeneratePerturbation\n");

	for (unsigned int fi = 0; fi < input.size(); ++fi)
	{
		//printf("  - fi = %u\n", fi);

		const TF1 *f_input = input[fi].ff.get();
		const TMatrixD &m_gen = * genMat[fi].get();
		TF1 *f_output = output[fi].ff.get();
		const auto st = stats[fi];

		TVectorD delta(m_gen.GetNrows()), rdm(m_gen.GetNrows());
		for (int i = 0; i < m_gen.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TVectorD params_new(m_gen.GetNrows());

		for (int pi = 0; pi < m_gen.GetNrows(); ++pi)
		{
			f_output->SetParameter(pi, f_input->GetParameter(pi) + delta(pi));

			params_new(pi) = f_output->GetParameter(pi);

			//printf("    - pi = %u: input = %.2E, delta = %.2E, output = %.2E\n", pi, f_input->GetParameter(pi), delta(pi), f_output->GetParameter(pi));
		}

		st->Fill(params_new);
	}
}

//----------------------------------------------------------------------------------------------------

FitResults AnalyzeFitUncertainty(const vector<FitData> &fits, const std::string &model)
{
	//printf("* AnalyzeFitUncertainty\n");

	// prepare perturbation generation
	vector<FitData> fits_perturbed;
	vector<shared_ptr<TMatrixD>> m_gens;
	vector<shared_ptr<Stat>> stats;

	for (const auto fd : fits)
	{
		const TF1 *ff = fd.ff.get(); 
		shared_ptr<TF1> nf(new TF1(* ff));
		TFitResultPtr fake_ptr;
		fits_perturbed.push_back({nf, fake_ptr});

		const auto &V = fd.fr->GetCovarianceMatrix();

		TMatrixDSymEigen eig_decomp(V);
		TVectorD eig_values(eig_decomp.GetEigenValues());
		TMatrixDSym S(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;

		m_gens.push_back(make_shared<TMatrixD>(eig_decomp.GetEigenVectors() * S));

		stats.push_back(make_shared<Stat>(V.GetNrows()));
	}

	// evaluate random perturbations
	Stat st(5);

	for (unsigned int ev = 0; ev < 10000; ++ev)
	{
		GeneratePerturbation(fits, m_gens, fits_perturbed, stats);

		const FitResults &r_per = AnalyzeFit(fits_perturbed, model);

		if (!r_per.valid)
			continue;

		st.Fill(
			r_per.t_dip,
			r_per.dsdt_dip,
			r_per.t_bmp,
			r_per.dsdt_bmp,
			r_per.R);
	}

	// print stats
	/*
	for (unsigned int fi = 0; fi < fits.size(); ++fi)
	{
		printf(" - fi = %u\n", fi);

		const auto *ff = fits[fi].ff.get();

		for (int pi = 0; pi < ff->GetNpar(); ++pi)
		{
			printf("    - pi = %u: mean = %.3E, unc = %.3E\n", pi, ff->GetParameter(pi), ff->GetParError(pi));
		}

		stats[fi]->PrintMeanAndStdDev();
	}
	*/

	// build output
	FitResults r;
	r.t_dip_unc = st.GetStdDev(0);
	r.dsdt_dip_unc = st.GetStdDev(1);
	r.t_bmp_unc = st.GetStdDev(2);
	r.dsdt_bmp_unc = st.GetStdDev(3);
	r.R_unc = st.GetStdDev(4);

	return r;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> models = {
		"local",
		"exp2+exp3"
	};

	// get input
	TFile *f_in = TFile::Open("../DS-merged/merged.root");
	TH1D *h_dsdt = (TH1D *) f_in->Get("bt1/DS4-sc/combined/h_dsdt");

	// load systematics
	const vector<string> syst_mode_names = {
		"alig-sh-thx", "alig-sh-thy", "tilt-thx-thy",
		"opt-m1", "opt-m2",
		"acc-corr-sigma-unc", "acc-corr-sigma-asym", "acc-corr-non-gauss",
		"eff-slp", "beam-mom",
		"unsm-sigma-x", "unsm-sigma-y", "unsm-model",
		"norm"
	};

	map<string, TGraph*> syst_modes;
	unique_ptr<TFile> f_in_syst(TFile::Open("../DS4/systematics_matrix.root"));
	for (const string &m : syst_mode_names)
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
		printf("\n------ %s ------\n", model.c_str());

		TDirectory *d_model = f_out->mkdir(model.c_str());

		// get list of ranges
		const auto &ranges = GetRanges(model);
		const auto range_central = ranges[0].data;

		// central fit
		gDirectory = d_model->mkdir("central");
		const auto &fits_central = MakeFit(h_dsdt, model, range_central);
		FitResults r_central = AnalyzeFit(fits_central, model);
		printf("* central fit:\n");
		r_central.Print();

		// central fit uncertainty
		FitResults r_central_unc = AnalyzeFitUncertainty(fits_central, model);
		printf("* central fit uncertainty:\n");
		r_central_unc.Print();

		// systematics
		TDirectory *d_systematics = d_model->mkdir("systematics");

		vector<double> vsum(5, 0.);
		for (const auto &syst_mode : syst_modes)
		{
			gDirectory = d_systematics->mkdir(syst_mode.first.c_str());

			auto h_dsdt_mod = ApplySystematicMode(h_dsdt, syst_mode.second);
			h_dsdt_mod->Write("h_dsdt_mod");

			FitResults r_syst = AnalyzeFit(MakeFit(h_dsdt_mod.get(), model, range_central), model);

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

		// fit range dependence
		TDirectory *d_fit_range = d_model->mkdir("fit range");

		Stat st(5);

		for (const auto &range : ranges)
		{
			gDirectory = d_fit_range->mkdir(range.label.c_str());

			FitResults r_range = AnalyzeFit(MakeFit(h_dsdt, model, range.data), model);

			if (!range.enabled)
			{
				TGraph *dummy = new TGraph();
				dummy->Write("disabled");
				continue;
			}

			st.Fill(
				r_range.t_dip - r_central.t_dip,
				r_range.dsdt_dip - r_central.dsdt_dip,
				r_range.t_bmp - r_central.t_bmp,
				r_range.dsdt_bmp - r_central.dsdt_bmp,
				r_range.R - r_central.R);
		}

		FitResults r_fit_range_unc;
		r_fit_range_unc.t_dip_unc = st.GetStdDev(0);
		r_fit_range_unc.dsdt_dip_unc = st.GetStdDev(1);
		r_fit_range_unc.t_bmp_unc = st.GetStdDev(2);
		r_fit_range_unc.dsdt_bmp_unc = st.GetStdDev(3);
		r_fit_range_unc.R_unc = st.GetStdDev(4);

		printf("* fit range uncertainty:\n");
		r_fit_range_unc.Print();

		// print summary
		printf("* summary:\n");
		printf("  t_dip    = %.4f +- %.4f (stat) +- %.4f (syst) +- %.4f (fit range)\n",
			r_central.t_dip, r_central_unc.t_dip_unc, r_syst_unc.t_dip_unc, r_fit_range_unc.t_dip_unc);
		printf("  dsdt_dip = %.4f +- %.4f (stat) +- %.4f (syst) +- %.4f (fit range)\n",
			r_central.dsdt_dip, r_central_unc.dsdt_dip_unc, r_syst_unc.dsdt_dip_unc, r_fit_range_unc.dsdt_dip_unc);
		printf("  t_bmp    = %.4f +- %.4f (stat) +- %.4f (syst) +- %.4f (fit range)\n",
			r_central.t_bmp, r_central_unc.t_bmp_unc, r_syst_unc.t_bmp_unc, r_fit_range_unc.t_bmp_unc);
		printf("  dsdt_bmp = %.4f +- %.4f (stat) +- %.4f (syst) +- %.4f (fit range)\n",
			r_central.dsdt_bmp, r_central_unc.dsdt_bmp_unc, r_syst_unc.dsdt_bmp_unc, r_fit_range_unc.dsdt_bmp_unc);
		printf("  R        = %.4f +- %.4f (stat) +- %.4f (syst) +- %.4f (fit range)\n",
			r_central.R, r_central_unc.R_unc, r_syst_unc.R_unc, r_fit_range_unc.R_unc);
	}

	// clean up
	delete f_out;
	delete f_in;

	return 0;
}