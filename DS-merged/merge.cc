#include "TFile.h"
#include "TH1D.h"

#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct shist
{
	TH1D *hist;
	double scale;
	shist(TH1D *_h, double _s) : hist(_h), scale(_s) {}
};

//----------------------------------------------------------------------------------------------------

TH1D* Merge(const vector<shist> &hists)
{
	if (hists.empty())
		return NULL;

	TH1D *m = new TH1D(*hists[0].hist);
	m->SetName("h_dsdt");

	for (int bi = 1; bi <= m->GetNbinsX(); bi++)
	{
		double Svw = 0., Sw = 0.;
		for (unsigned int hi = 0; hi < hists.size(); hi++)
		{
			double v = hists[hi].hist->GetBinContent(bi);
			double e = hists[hi].hist->GetBinError(bi) * hists[hi].scale;
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

struct Entry
{
	string input;
	double stat_unc_scale;
	string output;
	bool merge;

	Entry(const string &_i, double _sus, const string &_o, bool _m) : input(_i), stat_unc_scale(_sus), output(_o), merge(_m) {}
};

//----------------------------------------------------------------------------------------------------

int main()
{
	vector<Entry> entries;
	entries.push_back(Entry("DS4", 1., "DS4", false));
	
	entries.push_back(Entry("DS4", 1./1.17596, "DS4-sc", false));
	
	/*
	entries.push_back(Entry("DS4-b26", 1., "DS4-b26", true));
	entries.push_back(Entry("DS4-b2990", 1., "DS4-b2990", true));
	entries.push_back(Entry("DS4-b648", 1., "DS4-b648", true));
	*/
	
	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	vector<string> binnings;
	binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30");

	// prepare output
	TFile *f_out = new TFile("merged.root", "recreate");

	// loop over binnings
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		printf("\t%s\n", binnings[bi].c_str());

		TDirectory *binningDir = f_out->mkdir(binnings[bi].c_str());

		vector<shist> full_list;				// list of histograms for final merge
		map<string, vector<shist> > full_map;	// map: diagonal --> list of inputs

		for (unsigned int ei = 0; ei < entries.size(); ei++)
		{
			printf("\t\t%s\n", entries[ei].output.c_str());

			TDirectory *datasetDir = binningDir->mkdir(entries[ei].output.c_str());

			vector<shist> ds_list;
			for (unsigned int dgni = 0; dgni < diagonals.size(); dgni++)
			{
				printf("\t\t\t%s\n", diagonals[dgni].c_str());

				string fn = "../" + entries[ei].input + "/distributions_"+diagonals[dgni]+".root";
				TFile *f_in = TFile::Open(fn.c_str());
				if (!f_in)
				{
					printf("ERROR: File `%s' cannot be opened.\n", fn.c_str());
					return 1;
				}

				string on = "normalization+unfolding/"+binnings[bi]+"/h_t_normalized_unsmeared";
				//string on = "normalization/"+binnings[bi]+"/h_t_normalized";
				TH1D *h_dsdt = (TH1D *) f_in->Get(on.c_str());
				if (!h_dsdt)
				{
					printf("ERROR: Object `%s' cannot be loaded from file `%s'.\n", on.c_str(), fn.c_str());
					return 2;
				}
				
				h_dsdt->SetName("h_dsdt");

				shist sc_hist(h_dsdt, entries[ei].stat_unc_scale);

				ds_list.push_back(sc_hist);

				if (entries[ei].merge)
				{
					full_list.push_back(sc_hist);
					full_map[diagonals[dgni]].push_back(sc_hist);
				}

				gDirectory = datasetDir->mkdir(diagonals[dgni].c_str());
				h_dsdt->Write();
			}

			// save histogram merged from all inputs of a dataset
			gDirectory = datasetDir->mkdir("combined");
			TH1D *h_m = Merge(ds_list);
			if (h_m)
				h_m->Write();
		}

		// save merged histograms
		TDirectory *mergedDir = binningDir->mkdir("merged");

		for (map<string, vector<shist> >::iterator it = full_map.begin(); it != full_map.end(); ++it)
		{
			gDirectory = mergedDir->mkdir(it->first.c_str());
			TH1D *h_m = Merge(it->second);
			if (h_m)
				h_m->Write();
		}
		
		gDirectory = mergedDir->mkdir("combined");
		TH1D *h_m = Merge(full_list);
		if (h_m)
			h_m->Write();
	}

	delete f_out;

	return 0;
}
