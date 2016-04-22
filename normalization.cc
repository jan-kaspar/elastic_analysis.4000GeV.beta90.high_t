#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TH1D.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void MakeSum(TH1D *h, int bin_min, int bin_max, double &S, double &uS)
{
	S = 0.;
	double Suu = 0.;
	for (int bi = bin_min; bi <= bin_max; bi++)
	{
		double c = h->GetBinContent(bi);
		double w = h->GetBinWidth(bi);
		double u = h->GetBinError(bi);
		
		S += c*w;
		Suu += u*u*w*w;
	}
	uS = sqrt(Suu);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	string dgn_str = argv[1];

	Init(dgn_str);
	if (diagonal == ad45b_56b || diagonal == ad45t_56t)
		return rcIncompatibleDiagonal;

	// get input
	string fn_in = (diagonal == dCombined) ? "combine_distributions.root" : "distributions_" + dgn_str + ".root";
	TFile *f_in = new TFile(fn_in.c_str());
	TH1D *h_in = (TH1D *) f_in->Get("normalization+unfolding/ub/h_t_normalized_unsmeared");
	
	TFile *f_ref = new TFile("../../4000GeV,beta90/DS-merged/merged.root");
	TH1D *h_ref = (TH1D *) f_ref->Get(("ub/DS2/" + dgn_str + "/h_dsdt").c_str());

	// parameters
	int bi_min = 5;
	int bi_max = 12;

	// integration range
	printf("\n");
	printf("*  left edge of bin %2u: %.4f (this), %.4f (ref)\n", bi_min, h_in->GetBinLowEdge(bi_min), h_ref->GetBinLowEdge(bi_min));
	printf("* right edge of bin %2u: %.4f (this), %.4f (ref)\n", bi_max,
		h_in->GetBinLowEdge(bi_max) + h_in->GetBinWidth(bi_max),
		h_ref->GetBinLowEdge(bi_max) + h_ref->GetBinWidth(bi_max));

	// make sums
	printf("\n");

	double S_this = 0., uS_this = 0.;
	MakeSum(h_in, bi_min, bi_max, S_this, uS_this);
	printf("* this: S = %.3f +- %.3f\n", S_this, uS_this);

	double S_ref = 0., uS_ref = 0.;
	MakeSum(h_ref, bi_min, bi_max, S_ref, uS_ref);
	printf("* ref : S = %.3f +- %.3f\n", S_ref, uS_ref);

	// ratios
	printf("\n");

	double c = S_ref / S_this;
	double uc = sqrt(  uS_ref*uS_ref/S_this/S_this + S_ref*S_ref/S_this/S_this/S_this/S_this * uS_this*uS_this  );

	double d = S_this / S_ref;
	double ud = sqrt(  uS_this*uS_this/S_ref/S_ref + S_this*S_this/S_ref/S_ref/S_ref/S_ref * uS_ref*uS_ref  );
	
	printf("* ref / this = %.4f +- %.4f\n", c, uc);
	
	printf("* this / ref = %.4f +- %.4f\n", d, ud);

	delete f_ref;
	delete f_in;
	return 0;
}
