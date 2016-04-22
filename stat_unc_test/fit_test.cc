#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <set>

using namespace std;

//----------------------------------------------------------------------------------------------------

double ChiSq(double x, int ki)
{
	double k = ki;
	return pow(x, k/2. - 1.) * exp(-x/2.) / (pow(2., k/2.) * TMath::Gamma(k/2.));
}

//----------------------------------------------------------------------------------------------------

TGraph* SampleChiSqDist(int k, const char *name)
{
	TGraph *g = new TGraph();
	g->SetName(name);

	for (double x = 0; x <= 10.; x += 0.01)
	{
		double y = ChiSq(x, k);

		int idx = g->GetN();
		g->SetPoint(idx, x, y);
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

TGraph* SampleChiSqCumul(int ki, const char *name)
{
	TGraph *g = new TGraph();
	g->SetName(name);

	double k = ki;

	for (double x = 0; x <= 10.; x += 0.01)
	{
		double y = TMath::Gamma(k/2., x/2.);

		int idx = g->GetN();
		g->SetPoint(idx, x, y);
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

TH1D* ChiSqHist(int k, const char *name)
{
	TH1D *h = new TH1D(name, ";chi^{2}", 10, 0., 10.);

	for (int bi = 0; bi <= h->GetNbinsX(); bi++)
	{
		int N_div = 10;
		double l = h->GetBinLowEdge(bi);
		double w = h->GetBinWidth(bi);

		double S = 0.;
		for (int i = 0; i < N_div; i++)
		{
			double x = l + w / N_div * (0.5 + i);
			S += ChiSq(x, k);
		}
		double v = S / N_div;

		h->SetBinContent(bi, v);
	}

	return h;
}

//----------------------------------------------------------------------------------------------------

TF1 *f_exp1;
TF1 *f_exp2;

//----------------------------------------------------------------------------------------------------

TGraph* BuildKolmogorovGraph(set<double> vals, const string &name)
{
	TGraph *g = new TGraph();
	g->SetName(name.c_str());

	double S = 0.;
	for (set<double>::iterator it = vals.begin(); it != vals.end(); ++it)
	{
		int idx = g->GetN();
		g->SetPoint(idx, *it, S);
		S += 1. / vals.size();
		g->SetPoint(idx+1, *it, S);
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

void OneFitTest(TH1D *h, int bi_min, int bi_max, int fit_bins, TF1 *ff, const string &name)
{
	TGraphErrors *g_chisq = new TGraphErrors();
	g_chisq->SetName(("g_"+name).c_str());
	
	TH1D *h_chisq = new TH1D("", ";chi^{2}", 10, 0., 10.);
	h_chisq->SetName(("h_"+name).c_str());

	set<double> chisqs;

	for (int bi = bi_min; bi <= bi_max;)
	{
		int bi_beg = bi;
		int bi_end = bi + fit_bins - 1;

		if (bi_end > bi_max)
			break;

		//printf("\t%u --> %u\n", bi_beg, bi_end);

		double x_beg = h->GetBinLowEdge(bi_beg) + 0.1 * h->GetBinWidth(bi_beg);
		double x_end = h->GetBinLowEdge(bi_end) + 0.9 * h->GetBinWidth(bi_end);

		h->Fit(ff, "Q", "", x_beg, x_end);
		h->Fit(ff, "IQ", "", x_beg, x_end);

		double chisq = ff->GetChisquare();

		chisqs.insert(chisq);

		h_chisq->Fill(chisq);

		double x_mean = (h->GetBinLowEdge(bi_end) + h->GetBinWidth(bi_end) + h->GetBinLowEdge(bi_beg)) / 2.;
		double x_hw = (h->GetBinLowEdge(bi_end) + h->GetBinWidth(bi_end) - h->GetBinLowEdge(bi_beg)) / 2.;

		int idx = g_chisq->GetN();
		g_chisq->SetPoint(idx, x_mean, chisq);
		g_chisq->SetPointError(idx, x_hw, 0.);

		bi = bi_end + 1;
	}

	h_chisq->Scale(1./h_chisq->GetEntries());

	g_chisq->Write();
	h_chisq->Write();

	BuildKolmogorovGraph(chisqs, "g_kol_"+name)->Write();
}

//----------------------------------------------------------------------------------------------------

void AnalyzeOne(TH1D *h)
{
	printf(">> AnalyzeOne %p\n", h);

	// determine bin range
	int bi_min = -1, bi_max = -1;
	for (int bi = 1; bi <= h->GetNbinsX(); bi++)
	{
		double v = h->GetBinContent(bi);

		if (bi_min == -1 && v > 0)
			bi_min = bi;

		if (bi_min > -1 && v <= 0)
		{
			bi_max = bi - 1;
			break;
		}
	}

	if (bi_max == -1)
		bi_max = h->GetNbinsX();

	// run tests
	OneFitTest(h, bi_min, bi_max, 3, f_exp1, "3,exp1");
	OneFitTest(h, bi_min, bi_max, 4, f_exp1, "4,exp1");
	OneFitTest(h, bi_min, bi_max, 5, f_exp1, "5,exp1");

	OneFitTest(h, bi_min, bi_max, 4, f_exp1, "4,exp2");
	OneFitTest(h, bi_min, bi_max, 5, f_exp1, "5,exp2");
	OneFitTest(h, bi_min, bi_max, 6, f_exp1, "6,exp2");
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// init fit functions
	f_exp1 = new TF1("f_exp1", "[0] * exp([1]*x)");
	f_exp2 = new TF1("f_exp2", "[0] * exp([1]*x + [2]*x*x)");

	// settings
	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	vector<string> binnings;
	//binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30");

	// prepare output
	TFile *f_out = TFile::Open("fit_test.root", "recreate");

	// process data
	for (unsigned int dgni = 0; dgni < diagonals.size(); dgni++)
	{
		TDirectory *dgnDir = f_out->mkdir(diagonals[dgni].c_str());

		TFile *f_in = TFile::Open(("../DS4/distributions_" + diagonals[dgni] + ".root").c_str());

		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			TDirectory *binDir = dgnDir->mkdir(binnings[bi].c_str());
			gDirectory = binDir;

			TH1D *h_in = (TH1D *) f_in->Get(("normalization/" + binnings[bi] + "/h_t_normalized").c_str());

			AnalyzeOne(h_in);
		}

		delete f_in;
	}

	gDirectory = f_out;

	SampleChiSqDist(1, "g_chisq_ndf1")->Write();
	SampleChiSqDist(2, "g_chisq_ndf2")->Write();
	SampleChiSqDist(3, "g_chisq_ndf3")->Write();

	SampleChiSqCumul(1, "g_cum_chisq_ndf1")->Write();
	SampleChiSqCumul(2, "g_cum_chisq_ndf2")->Write();
	SampleChiSqCumul(3, "g_cum_chisq_ndf3")->Write();
	
	ChiSqHist(1, "h_chisq_ndf1")->Write();
	ChiSqHist(2, "h_chisq_ndf2")->Write();
	ChiSqHist(3, "h_chisq_ndf3")->Write();

	delete f_out;
	return 0;
}
