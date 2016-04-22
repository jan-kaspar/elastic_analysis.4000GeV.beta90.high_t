#include "TDirectory.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"

#include <vector>
#include <string>
#include <cmath>

#include "../../stat.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void AnalyzeOne(const vector<string> &files, const string &obj_path, TDirectory *d_out)
{
	TH1D *h_ref = NULL;

	int Nb = 0;
	Stat st_v, st_u;

	// loop over datasets
	for (unsigned int fi = 0; fi < files.size(); fi++)
	{
		string file = files[fi];
		TFile *f_in = new TFile(file.c_str());

		TH1D *h = (TH1D *) f_in->Get(obj_path.c_str());

		if (!h)
		{
			printf("ERROR: can't load object `%s' from file `%s'. Skipping.\n", obj_path.c_str(), file.c_str());
			continue;
		}

		if (h_ref == NULL)
		{
			gDirectory = d_out;
			h_ref = new TH1D(*h);
			Nb = h_ref->GetNbinsX();
			st_v.Init(Nb);
			st_u.Init(Nb);
		}

		// loop over bins
		vector<double> vec_v(Nb), vec_u(Nb);
		for (int bi = 1; bi <= Nb; bi++)
		{
			int i = bi - 1;

			double v = h->GetBinContent(bi);
			double u = h->GetBinError(bi);

			vec_v[i] = v;
			vec_u[i] = u;
		}

		st_v.Fill(vec_v);
		st_u.Fill(vec_u);
	
		delete f_in;
	}

	if (h_ref == NULL)
	{
		printf("ERROR: no data to process. Skipping.\n");
		return;
	}
	
	// make final histograms
	gDirectory = d_out;

	TH1D *h_v_mean = new TH1D(*h_ref); h_v_mean->SetName("h_v_mean"); h_v_mean->SetLineColor(1);
	TH1D *h_v_mean_sqrt = new TH1D(*h_ref); h_v_mean_sqrt->SetName("h_v_mean_sqrt"); h_v_mean_sqrt->SetLineColor(1);
	TH1D *h_v_rms = new TH1D(*h_ref); h_v_rms->SetName("h_v_rms"); h_v_rms->SetLineColor(2);
	TH1D *h_u_mean = new TH1D(*h_ref); h_u_mean->SetName("h_u_mean"); h_u_mean->SetLineColor(4);

	for (int i = 0; i < Nb; i++)
	{
		int bi = i + 1;

		double v_m = st_v.GetMean(i);
		double v_m_u = st_v.GetMeanUnc(i);
		double v_s = st_v.GetStdDev(i);
		double v_s_u = st_v.GetStdDevUnc(i);

		// debug
		/*
		double v_s_ug = st_v.GetStdDevUncGauss(i);
		double v_s_u_2d = (v_s > 0.) ? st_v.GetCovarianceUnc(i, i) / 2. / v_s : 0.;
		printf("idx=%3i | S1=%5.0f | %E, %E, %E | %E, %E\n", i, st_v.S1,
			v_s_u, v_s_ug, v_s_u_2d,
			v_s, sqrt(st_v.GetCovariance(i, i)));
		*/

		double u_m = st_u.GetMean(i);
		double u_m_u = st_u.GetMeanUnc(i);

		double v_m_sqrt = sqrt(v_m);
		double v_m_sqrt_u = (v_m > 0.) ? v_m_u / 2. / sqrt(v_m) : 0.;
		//printf("%i, %E, %E\n", bi, v_m_sqrt, v_m_sqrt_u);

		h_v_mean->SetBinContent(bi, v_m); h_v_mean->SetBinError(bi, v_m_u);
		h_v_mean_sqrt->SetBinContent(bi, v_m_sqrt); h_v_mean_sqrt->SetBinError(bi, v_m_sqrt_u);
		h_v_rms->SetBinContent(bi, v_s); h_v_rms->SetBinError(bi, v_s_u);
		h_u_mean->SetBinContent(bi, u_m); h_u_mean->SetBinError(bi, u_m_u);
	}

	h_v_mean->Write();
	h_v_mean_sqrt->Write();
	h_v_rms->Write();
	h_u_mean->Write();

	// make ratio
	TH1D *h_v_rms_over_u_mean = new TH1D(*h_v_rms);
	h_v_rms_over_u_mean->SetName("h_v_rms_over_u_mean");
	h_v_rms_over_u_mean->Divide(h_u_mean);
	h_v_rms_over_u_mean->Fit("pol1", "Q", "", 0.02, 0.3);
	h_v_rms_over_u_mean->Write();

	h_v_rms_over_u_mean->Fit("pol1", "WQ", "", 0.02, 0.3);
	h_v_rms_over_u_mean->Write("h_v_rms_over_u_mean, pol1, ew");
	
	TF1 *f1 = new TF1("f1", "[0]-[0]+1");
	h_v_rms_over_u_mean->Fit(f1, "Q", "", 0.02, 0.3);
	h_v_rms_over_u_mean->Write("h_v_rms_over_u_mean, f1");

	/*
	TF1 *ff = new TF1("ff", "[0]/abs(x-[1])^[2] + [3] + [4]*x + [5]*x^2 + [6]*x^3");
	ff->SetParameters(1.39438e-04, 4.13146e-03, 1.97862e+00, 3.40574e-01, 1.24682e+01, -8.42749e+01, 1.77299e+02);
	h_v_rms_over_u_mean->Fit(ff, "Q", "", 0.01, 0.2);
	h_v_rms_over_u_mean->Fit(ff, "Q", "", 0.01, 0.2);
	h_v_rms_over_u_mean->Write("h_v_rms_over_u_mean, ff");
	*/

	// correlation plots
	const double *binEdges = h_ref->GetXaxis()->GetXbins()->GetArray();

	TH2D *h_v_cov = new TH2D("h_v_cov", "", Nb, binEdges, Nb, binEdges);
	TH2D *h_v_cor = new TH2D("h_v_cor", "", Nb, binEdges, Nb, binEdges);
	
	TGraphErrors *g_v_cov = new TGraphErrors(); g_v_cov->SetName("g_v_cov");
	TGraphErrors *g_v_cor = new TGraphErrors(); g_v_cor->SetName("g_v_cor");

	for (int i = 0; i < Nb; i++)
	{
		for (int j = 0; j < Nb; j++)
		{
			int bi = i + 1;
			int bj = j + 1;

			double ti = h_v_cov->GetXaxis()->GetBinCenter(bi);
			double tj = h_v_cov->GetXaxis()->GetBinCenter(bj);

			double t_max = 0.2;
			if (ti > t_max && tj > t_max)
				continue;

			double c = st_v.GetCovariance(i, j);
			double c_u = st_v.GetCovarianceUnc(i, j);

			double corr = st_v.GetCorrelation(i, j);
			double corr_u = st_v.GetCorrelationUnc(i, j);

			if (i == j)
				corr = 0;

			//if (j > i)
			//{
				int idx = g_v_cov->GetN();
				g_v_cov->SetPoint(idx, idx, c);
				g_v_cov->SetPointError(idx, 0., c_u);

				g_v_cor->SetPoint(idx, idx, corr);
				g_v_cor->SetPointError(idx, 0., corr_u);
			//}

			//printf("%i, %E, %E\n", idx, corr, corr_u);

			h_v_cov->SetBinContent(bi, bj, c);
			h_v_cor->SetBinContent(bi, bj, corr);
		}
	}

	h_v_cov->Write();
	h_v_cor->Write();
	g_v_cov->Write();
	g_v_cor->Write();
}

