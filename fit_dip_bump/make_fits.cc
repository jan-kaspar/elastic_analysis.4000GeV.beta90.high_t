#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

int main()
{
	TFile *f_in = TFile::Open("../DS-merged/merged.root");
	TH1D *h_dsdt = (TH1D *) f_in->Get("bt1/DS4-sc/combined/h_dsdt");
	//TH1D *h_dsdt = (TH1D *) f_in->Get("ob-2-20-0.20/DS4-sc/combined/h_dsdt");

	TFile *f_out = TFile::Open("make_fits.root", "recreate");

	TF1 *ff0 = new TF1("ff0", "[1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x) + [5]*[5]*exp(2*[6]*x)");
	ff0->SetParameter(0, 2.77);
	ff0->SetParameter(1, 24.33);
	ff0->SetParameter(2, -9.71);
	ff0->SetParameter(3, 4.52);
	ff0->SetParameter(4, -3.32);
	ff0->SetParameter(5, 1.33);
	ff0->SetParameter(6, -2.45);

	h_dsdt->Fit(ff0, "+", "", 0.0, 1.9);

	ff0->Write("ff0");

	TF1 *ff1 = new TF1("ff1", "[1]*[1]*exp(2*[2]*x + 2*[3]*x*x + 2*[4]*x*x*x) + 2 * cos([0]) * [1]*exp([2]*x + [3]*x*x + [4]*x*x*x) * [5]*exp([6]*x + [7]*x*x) + [5]*[5]*exp(2*[6]*x + 2*[7]*x*x)");
	ff1->SetParameter(0, 2.77);
	ff1->SetParameter(1, 24.33);
	ff1->SetParameter(2, -9.71);
	ff1->SetParameter(3, 4.52);
	ff1->SetParameter(4, -3.32);
	ff1->SetParameter(5, 1.33);
	ff1->SetParameter(6, -2.45);
	ff1->SetParameter(7, 0.);
	
	h_dsdt->Fit(ff1, "+", "", 0.0, 1.9);

	ff1->Write("ff1");
	
	TF1 *ff2 = new TF1("ff2", "[0] + [1] * pow(x - [2], 2)");
	ff2->SetParameters(0.0147, 1., 0.523);
	h_dsdt->Fit(ff2, "+", "", 0.47, 0.57);

	const double t_dip = ff2->GetParameter(2), t_dip_unc = ff2->GetParError(2);
	const double dsdt_dip = ff2->GetParameter(0), dsdt_dip_unc = ff2->GetParError(0);
	printf("* dip: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_dip, t_dip_unc, dsdt_dip, dsdt_dip_unc);

	TF1 *ff3 = new TF1("ff3", "[0] + [1] * pow(x - [2], 2)");
	ff3->SetParameters(0.0295, -1., 0.69);
	h_dsdt->Fit(ff3, "+", "", 0.55, 0.85);

	const double t_bump = ff3->GetParameter(2), t_bump_unc = ff3->GetParError(2);
	const double dsdt_bump = ff3->GetParameter(0), dsdt_bump_unc = ff3->GetParError(0);
	printf("* bump: t = %.4f +- %.4f, dsdt = %.4f +- %.4f\n", t_bump, t_bump_unc, dsdt_bump, dsdt_bump_unc);

	const double R = dsdt_bump / dsdt_dip, R_unc = R * sqrt(pow(dsdt_dip_unc / dsdt_dip, 2.) + pow(dsdt_bump_unc / dsdt_bump, 2.));

	printf("\n");
	printf("* R = dsdt_bump / dsdt_dip = %.3f +- %.3f\n", R, R_unc);

	h_dsdt->Write("h_dsdt");

	delete f_out;
	delete f_in;

	return 0;
}
