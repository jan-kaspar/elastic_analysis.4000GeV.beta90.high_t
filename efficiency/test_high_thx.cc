#include <cstdio>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include "../stat.h"

using namespace std;

int main()
{
	vector<string> diagonals = {
		"45b_56t",
		"45t_56b",
	};

	vector<string> RPs = {
		"L_F", 
		"L_N", 
		"R_N", 
		"R_F", 
	};

	TF1 *ff = new TF1("ff", "[0]");
	
	TF1 *ffs = new TF1("ff", "[0] - [1]/2 * (TMath::Erf((abs(x) - [2])/sqrt(2)/[3]) + 1)");

	Stat st(1);

	TFile *f_out = TFile::Open("test_high_thx.root", "recreate");
	assert(f_out != NULL);

	for (auto diagonal : diagonals)
	{
		printf("\n\n* %s\n", diagonal.c_str());

		double S_diff = 0.;

		TFile *f_in = TFile::Open(("../DS4/eff3outof4_more_cuts_"+diagonal+".root").c_str());
		assert(f_in != NULL);

		TDirectory *dgnDir = f_out->mkdir(diagonal.c_str());

		for (auto RP : RPs)
		{
			TDirectory *rpDir = dgnDir->mkdir(RP.c_str());
			gDirectory = rpDir;

			TH1D *h = (TH1D *) f_in->Get(("excluded RPs "+RP+"/n_si 3.0/th_x dependence/h_refined_ratio").c_str());
			assert(h != NULL);

			h->Fit(ff, "Q", "", -80E-6, +80E-6);
			double eff_c = ff->GetParameter(0) * 100.;

			//h->Fit(ff, "Q", "", 170E-6, +300E-6);
			h->Fit(ff, "Q", "", -300E-6, -170E-6);
			double eff_r = ff->GetParameter(0) * 100.;
			double eff_r_u = ff->GetParError(0) * 100.;

			ffs->SetParameters(eff_c/100., (eff_c - eff_r)/100., 150E-6, 20E-6);
			ffs->SetParLimits(2, 130E-6, 160E-6);
			ffs->SetParLimits(3, 15E-6, 30E-6);
			h->Fit(ffs, "Q", "", 0E-6, 300E-6);
			h->Write();

			double diff = eff_r - eff_c;

			S_diff += diff;
			st.Fill(diff);

			printf("%.1f, %.1f => %+.1f +- %.1f (%.1f si)\n", eff_c, eff_r, eff_r - eff_c, eff_r_u, fabs(eff_r - eff_c) / eff_r_u);

			printf("\tc=%.1f, de=%.1f, m=%.1f, si=%.1f\n", ffs->GetParameter(0)*100., ffs->GetParameter(1)*100.,
				ffs->GetParameter(2)*1E6, ffs->GetParameter(3)*1E6);
		}

		printf("S_diff = %.1f\n", S_diff);

		delete f_in;
	}

	st.PrintMeanAndStdDev();

	delete f_out;
	return 0;
}
