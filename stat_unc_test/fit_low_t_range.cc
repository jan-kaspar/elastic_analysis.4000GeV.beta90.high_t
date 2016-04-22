#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

TGraphErrors* BuildGraph(TH1D *h, TF1 *f_fit, const string &name)
{
	TGraphErrors *g = new TGraphErrors();
	g->SetName(name.c_str());

	for (int bi = 1; bi <= h->GetNbinsX(); bi++)
	{
		double l = h->GetBinLowEdge(bi);
		double w = h->GetBinWidth(bi);
		double r = l + w;

		double V = f_fit->Integral(l, r) / w;

		double xr;
		while (r - l > w/10000)
		{
			xr = (r+l)/2.;
			if (f_fit->Eval(xr) < V)
			{
				r -= (r-l)/2.;
			} else {
				l += (r-l)/2.;
			}
		}
		xr = (r+l)/2.;

		//printf("\t\ti=%2u: c = %.5f, xr = %.5f\n", bi, h->GetBinCenter(bi), xr);

		int idx = g->GetN();
		g->SetPoint(idx, xr, h->GetBinContent(bi));
		g->SetPointError(idx, 0., h->GetBinError(bi));
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// configuration
	string dataset = "DS4";

	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	vector<string> binnings;
	//binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30");

	// fit range
	double t_min = 0.;
	double t_max = 0.2;

	// fitting function
	TF1 *ff = new TF1("ff", "[0] * exp([1]*x + [2]*x*x + [3]*x*x*x)");

	// prepare output
	TFile *f_out = TFile::Open("fit_low_t_range.root", "recreate");

	// make fits
	for (unsigned int dgni = 0; dgni < diagonals.size(); dgni++)
	{
		printf("* %s\n", diagonals[dgni].c_str());
		
		TDirectory *dgnDir = f_out->mkdir(diagonals[dgni].c_str());

		TFile *f_in = TFile::Open(("../"+dataset+"/distributions_" + diagonals[dgni] + ".root").c_str());

		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			printf("\t%s\n", binnings[bi].c_str());

			TDirectory *binDir = dgnDir->mkdir(binnings[bi].c_str());
			gDirectory = binDir;

			TH1D *h_in = (TH1D *) f_in->Get(("normalization/" + binnings[bi] + "/h_t_normalized").c_str());

			// initial parameters
			ff->SetParameters(535., -20.3, 7., -15.);

			// run fit
			h_in->Fit(ff, "Q", "", t_min, t_max);
			h_in->Fit(ff, "IQ", "", t_min, t_max);

			double chisq = ff->GetChisquare();
			double ndf = ff->GetNDF();

			printf("\t\thist : chi^2 / ndf = %.3f / %.0f = %.3f\n", chisq, ndf, chisq / ndf);

			h_in->Write("h_dsdt");

			char buf[100];
			sprintf(buf, "$\\ch^2$ / ndf = %.3f / %.0f = %.3f\n", chisq, ndf, chisq / ndf);
			ff->SetTitle(buf);
			ff->Write();

			// build graph with representative points
			TGraphErrors *g_dsdt = BuildGraph(h_in, ff, "g_dsdt");
			g_dsdt->Write();

			// fit graph
			g_dsdt->Fit(ff, "Q", "", t_min, t_max);
			
			chisq = ff->GetChisquare();
			ndf = ff->GetNDF();

			printf("\t\tgraph: chi^2 / ndf = %.3f / %.0f = %.3f\n", chisq, ndf, chisq / ndf);
		}

		delete f_in;
	}
	
	delete f_out;
	return 0;
}
