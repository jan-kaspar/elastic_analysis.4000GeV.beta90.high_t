#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF1.h"

#include <cmath>

//----------------------------------------------------------------------------------------------------

TGraph* MakeFixStatUncBinSizeGraph(TProfile *p_acc_corr, TF1 *ff_cs, double lumi, double u, const char *name)
{
	TGraph *g_bs = new TGraph();
	g_bs->SetName(name);

	for (int i = 0; i < p_acc_corr->GetNbinsX(); i++)
	{
		double t = p_acc_corr->GetBinCenter(i);
		double acc_corr = p_acc_corr->GetBinContent(i);

		if (acc_corr == 0)
			continue;

		double A = 1./acc_corr;

		double dsdt = ff_cs->Eval(t);
		double bin_size = 1. / (u*u * lumi * dsdt * A);

		g_bs->SetPoint(g_bs->GetN(), t, bin_size);
	}

	return g_bs;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_acc_corr = new TFile("acceptance_correction_DS4.root");
	TProfile *p_acc_corr = (TProfile *) f_acc_corr->Get("p_t_eb_full_corr");

	TFile *f_cs = new TFile("dsigma_dt_fit_DS4.root");
	TF1 *ff_cs = (TF1 *) f_cs->Get("ff");
	
	TFile *f_simu = new TFile("simu.root");
	TGraph *g_rms_t = (TGraph *) f_simu->Get("g_RMS_de_t_vs_t");

	// prepare output
	TFile *f_out = new TFile("generators.root", "recreate");

	g_rms_t->SetName("g_rms_t");
	g_rms_t->Write();

	MakeFixStatUncBinSizeGraph(p_acc_corr, ff_cs, 740.7E3, 0.10, "g_bs_stat_unc_10")->Write();
	MakeFixStatUncBinSizeGraph(p_acc_corr, ff_cs, 740.7E3, 0.20, "g_bs_stat_unc_20")->Write();
	MakeFixStatUncBinSizeGraph(p_acc_corr, ff_cs, 740.7E3, 0.30, "g_bs_stat_unc_30")->Write();

	delete f_out;

	return 0;
}
