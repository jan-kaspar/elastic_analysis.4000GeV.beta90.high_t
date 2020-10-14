#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include <memory>
#include <vector>
#include <map>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_in = TFile::Open("../DS-merged/merged.root");
	TH1D *h_dsdt = (TH1D *) f_in->Get("bt1/DS4-sc/combined/h_dsdt");

	// prepare output
	TFile *f_out = TFile::Open("make_fits.root", "recreate");

	h_dsdt->Write("h_dsdt");

	// make fits
	TF1 *ff;

	ff = new TF1("p1*exp3+p2*exp2", "([0] + [1]*x) * exp([2]*x + [3]*x*x + [4]*x*x*x) + ([5] + [6]*x + [7]*x*x) * exp([8]*x + [9]*x*x)");
	ff->SetParameter(0, 544.);
	ff->SetParameter(1, 0.);
	ff->SetParameter(2, -20.1);
	ff->SetParameter(3, 8.37);
	ff->SetParameter(4, -11.1);
	ff->SetParameter(5, -6.99);
	ff->SetParameter(6, 12.85);
	ff->SetParameter(7, 0.);
	ff->SetParameter(8, -6.17);
	ff->SetParameter(9, 0.);
	ff->SetRange(0, 2.);
	h_dsdt->Fit(ff, "", "", 0., 1.85);
	ff->Write();

	ff = new TF1("exp3-intf-exp1", "[1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x) + [5]*[5]*exp(2*[6]*x)");
	ff->SetParameter(0, 2.77);
	ff->SetParameter(1, 24.33);
	ff->SetParameter(2, -9.71);
	ff->SetParameter(3, 4.52);
	ff->SetParameter(4, -3.32);
	ff->SetParameter(5, 1.33);
	ff->SetParameter(6, -2.45);
	ff->SetRange(0, 2.);
	h_dsdt->Fit(ff, "", "", 0., 1.85);
	ff->Write();

	// clean up
	delete f_out;
	delete f_in;

	return 0;
}