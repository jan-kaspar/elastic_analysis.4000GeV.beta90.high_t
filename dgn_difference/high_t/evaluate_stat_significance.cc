#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOne(TH1D *h_45b, TH1D *h_45t)
{
	double S2 = 0.;
	int N = 0;
	for (int bi = 1; bi <= h_45b->GetNbinsX(); bi++)
	{
		double c = h_45b->GetBinCenter(bi);

		if (c < 1.1 || c > 1.9)
			continue;

		N++;

		double v_45b = h_45b->GetBinContent(bi), u_45b = h_45b->GetBinError(bi);
		double v_45t = h_45t->GetBinContent(bi), u_45t = h_45t->GetBinError(bi);
		double u = sqrt(u_45b*u_45b + u_45t*u_45t);
		double rel = (v_45b - v_45t) / u;
		S2 += rel * rel;
	}

	double p_value = TMath::Prob(S2, N);
	double sig = sqrt(2.) * TMath::ErfcInverse(p_value);

	//printf("\tS2 / N = %.1f / %i = %.2f, p_value = %.2f, significance = %.1f\n", S2, N, S2/N, p_value, sig);
	printf("%.1f / %i = %.2f & %.2f & %.1f\n", S2, N, S2/N, p_value, sig);
}

//----------------------------------------------------------------------------------------------------

int main()
{
	string dataset = "DS4";

	vector<string> binnings = {
		"ob-1-30-0.10",
		"ob-2-20-0.20",
		"ob-3-10-0.30",
	};

	//--------------------

	// get files
	TFile *f_in_45b = TFile::Open(("../"+dataset+"/distributions_45b_56t.root").c_str());
	TFile *f_in_45t = TFile::Open(("../"+dataset+"/distributions_45t_56b.root").c_str());

	assert(f_in_45b != NULL);
	assert(f_in_45t != NULL);

	for (auto binning : binnings)
	{
		printf("* %s\n", binning.c_str());

		// get histograms
		TH1D *h_45b = (TH1D *) f_in_45b->Get(("normalization/"+binning+"/h_t_normalized").c_str());
		TH1D *h_45t = (TH1D *) f_in_45t->Get(("normalization/"+binning+"/h_t_normalized").c_str());

		assert(h_45b != NULL);
		assert(h_45t != NULL);

		printf("\\hbox{%s} & ", binning.c_str());
		ProcessOne(h_45b, h_45t);
	}

	return 0;
}
