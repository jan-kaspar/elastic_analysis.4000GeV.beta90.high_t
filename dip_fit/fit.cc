#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
//#include "TGraph.h"

#include <string>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

TF1 *ff = NULL;

void Fit(TH1D *h_dsdt)
{
	ff->SetParameters(200., 0.52, -4.);

	h_dsdt->Fit(ff, "", "", 0.48, 0.57);

	printf("dip = %.3f +- %.3f\n", ff->GetParameter(1), ff->GetParError(1));

	h_dsdt->Write();
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// prepare analysis
	ff = new TF1("ff", "exp([0]*(x-[1])^2 + [2])");

	// prepare input
	TFile *f_in = new TFile("../DS-merged/merged.root");

	// prepare output
	TFile *f_out = new TFile("fit.root", "recreate");

	{
		TH1D *h_dsdt = (TH1D *) f_in->Get("ob-2-20-0.20/DS4-sc/combined/h_dsdt");

		//gDirectory = f_out->mkdir(diagonals[dgni].c_str());
		
		Fit(h_dsdt);
	}

	delete f_out;
	return 0;
}
