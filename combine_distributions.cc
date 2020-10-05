#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

/*
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"

*/

#include <cassert>

using namespace std;

//----------------------------------------------------------------------------------------------------

void MkdirP(TFile *f_out, const string &path, string &objName)
{
	TDirectory *d = f_out;

	size_t ps = 0;
	size_t pe = path.find("/", ps);
	while (pe != string::npos)
	{
		string d_name = path.substr(ps, pe-ps);

		TDirectory *d_next = d->GetDirectory(d_name.c_str());
		if (!d_next)
			d_next = d->mkdir(d_name.c_str());

		d = d_next;

		ps = pe + 1;
		pe = path.find("/", ps);
	}

	gDirectory = d;

	objName = path.substr(ps, string::npos);
}

//----------------------------------------------------------------------------------------------------

TH1D* Merge(const vector<TH1D *> &hists, TH1D *m)
{
	for (int bi = 1; bi <= m->GetNbinsX(); bi++)
	{
		double Svw = 0., Sw = 0.;
		for (unsigned int hi = 0; hi < hists.size(); hi++)
		{
			double v = hists[hi]->GetBinContent(bi);
			double e = hists[hi]->GetBinError(bi);
			double w = (e > 0.) ? 1./e/e : 0.;

			Sw += w;
			Svw += v * w;
		}

		double v = (Sw > 0.) ? Svw / Sw : 0.;
		double e = (Sw > 0.) ? 1. / sqrt(Sw) : 0.;

		m->SetBinContent(bi, v);
		m->SetBinError(bi, e);
	}

	return m;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	// binnings
	vector<string> binnings;
	binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30"); 
	binnings.push_back("bt1");
	binnings.push_back("bt2");

	// histograms to merge
	vector<string> sources;
	for (auto binning : binnings)
	{
		sources.push_back("normalization/"+binning+"/h_t_normalized");
		sources.push_back("normalization+unfolding/"+binning+"/h_t_normalized_unsmeared");
	}

	// input files
	TFile *f_in_45b = TFile::Open("distributions_45b_56t.root");
	assert(f_in_45b != NULL);

	TFile *f_in_45t = TFile::Open("distributions_45t_56b.root");
	assert(f_in_45t != NULL);

	// output file
	TFile *f_out = TFile::Open("combine_distributions.root", "recreate");
	assert(f_out != NULL);

	// process
	for (auto src : sources)
	{
		printf("* %s\n", src.c_str());

		TH1D *h_45b = (TH1D *) f_in_45b->Get(src.c_str());
		TH1D *h_45t = (TH1D *) f_in_45t->Get(src.c_str());
		assert(h_45b != NULL);
		assert(h_45t != NULL);

		string objName;
		MkdirP(f_out, src, objName);

		TH1D *h = new TH1D(*h_45b);
		h->SetName(objName.c_str());

		Merge({h_45b, h_45t}, h);

		h->Write();
	}

	delete f_out;
	delete f_in_45b;
	delete f_in_45t;
	return 0;
}
