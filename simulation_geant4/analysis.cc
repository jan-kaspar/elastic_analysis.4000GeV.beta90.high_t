#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <cmath>

#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

// map: element set, condition -> count 
typedef map<string, map<string, unsigned long > > CounterMap;

//----------------------------------------------------------------------------------------------------

#if 0
/**
NON-PARALLEL pattern recognition run with these parameters:
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99
**/

#endif


//----------------------------------------------------------------------------------------------------

struct RPStruct
{
	string name;

	RPRootDumpDigiInfo *digi;
	RPRootDumpPatternInfo *pat;
	RPRootDumpTrackInfo *tr;
	
	RPRootDumpTrackInfo *tr_id;

	bool prot;
	bool pl_suff;
	bool pl_too_full_u, pl_too_full_v;
	bool pat_suff;
	bool pat_more;
	//bool tr_any;
	bool tr_val;

	double mco;

	void AssignBranches(const string &_n, TChain *ch_full, TChain *ch_ideal, unsigned int id)
	{
		name = _n;

		char buf[100];

		digi = new RPRootDumpDigiInfo;
		sprintf(buf, "digi_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch_full->SetBranchAddress(buf, &digi);

		pat = new RPRootDumpPatternInfo();
		sprintf(buf, "nonpar_patterns_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "nonpar_patterns_rp_%u.", id); ch_full->SetBranchAddress(buf, &pat);

		tr = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_full->SetBranchAddress(buf, &tr);

		tr_id = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch_ideal->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_ideal->SetBranchAddress(buf, &tr_id);
	}
	
	bool RPTooFullU()
	{
		// count planes that could have (in principle) contributed to pattern search
		unsigned N = 0;
		for (unsigned int i = 1; i < 10; i += 2)
			if (digi->numberOfClusters[i] <= 5)
				N++;
	
		return (N < 3);
	}
	
	bool RPTooFullV()
	{
		unsigned N = 0;
		for (unsigned int i = 0; i < 10; i += 2)
			if (digi->numberOfClusters[i] <= 5)
				N++;
	
		return (N < 3);
	}

	void Analyze(CounterMap &c);
};

//----------------------------------------------------------------------------------------------------

void RPStruct::Analyze(CounterMap &c)
{
	prot = tr_id->valid;

	pl_suff = (digi->uPlanesOn >= 3 && digi->vPlanesOn >= 3);
	pl_too_full_u = RPTooFullU();
	pl_too_full_v = RPTooFullV();
	pat_suff = (pat->u.size() > 0 || pat->v.size() > 0) || (pl_too_full_u || pl_too_full_v);
	pat_more = (pat->u.size() > 1 && pat->v.size() > 1) || (pl_too_full_u || pl_too_full_v);
	//tr_any = (tr->valid || mtr->size() > 0);

	tr_val = tr->valid;

	if (prot) c[name]["prot"]++;

	if (prot & pat_suff) c[name]["prot & pat_suff"]++;
	if (prot & !pat_suff) c[name]["prot & !pat_suff"]++;
	if (!prot & pat_suff) c[name]["!prot & pat_suff"]++;
	if (!prot & !pat_suff) c[name]["!prot & !pat_suff"]++;

	if (prot & tr_val) c[name]["prot & tr_val"]++;
	if (prot & !tr_val) c[name]["prot & !tr_val"]++;
	if (!prot & tr_val) c[name]["!prot & tr_val"]++;
	if (!prot & !tr_val) c[name]["!prot & !tr_val"]++;

	// determine mean cluster occupancy
	mco = 0.;
	for (unsigned int i = 0; i < 10; i++)
		mco += digi->numberOfClusters[i];
	mco /= 10;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct ActivityCorrelationStruct
{
	string name;

	TH1D *h_N, *h_F, *h_NO, *h_FO;
	TH2D *h2_N_F, *h2_N_FO, *h2_F_FO;

	void Init(const string &_n)
	{
		name = _n;
		h_N = new TH1D("", ";mco_N", 30, 0., 30.);
		h_F = new TH1D("", ";mco_F", 30, 0., 30.);
		h_NO = new TH1D("", ";mco_NO", 30, 0., 30.);
		h_FO = new TH1D("", ";mco_FO", 30, 0., 30.);

		h2_N_F = new TH2D("", ";mso_F;mco_N", 30, 0., 30., 30, 0., 30.);
		h2_N_FO = new TH2D("", ";mco_FO;mco_N", 30, 0., 30., 30, 0., 30.);
		h2_F_FO = new TH2D("", ";mco_FO;mco_F", 30, 0., 30., 30, 0., 30.);
	}

	void Fill(double mco_N, double mco_F, double mco_NO, double mco_FO)
	{
		h_N->Fill(mco_N);
		h_F->Fill(mco_F);
		h_NO->Fill(mco_NO);
		h_FO->Fill(mco_FO);

		h2_N_F->Fill(mco_F, mco_N);
		h2_N_FO->Fill(mco_FO, mco_N);
		h2_F_FO->Fill(mco_FO, mco_F);
	}

	void Write()
	{
		TDirectory *topDir = gDirectory;
		gDirectory = topDir->mkdir(name.c_str());

		h_N->Write("h_N");
		h_F->Write("h_F");
		h_NO->Write("h_NO");
		h_FO->Write("h_FO");

		h2_N_F->Write("h2_N_F");
		h2_N_FO->Write("h2_N_FO");
		h2_F_FO->Write("h2_F_FO");

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

TH1D* MakeOneRatio(TH1D *h_n, TH1D *h_d)
{
	TH1D *h = new TH1D(*h_n);
	h->Sumw2();

	for (int bi = 1; bi <= h->GetNbinsX(); bi++)
	{
		double n = h->GetBinContent(bi);
		//double n_u = h->GetBinError(bi);
		double d = h_d->GetBinContent(bi);

		double r = (d > 0) ? n / d : 0.;
		double r_u = (d > 0) ? sqrt(n * (1. - n/d)) / d : 0.;

		h->SetBinContent(bi, r);
		h->SetBinError(bi, r_u);
	}

	return h;
}

//----------------------------------------------------------------------------------------------------

TH2D* MakeOneRatio2D(TH2D *h_n, TH2D *h_d)
{
	TH2D *h = new TH2D(*h_n);
	h->Sumw2();

	for (int bi = 1; bi <= h->GetNbinsX(); bi++)
	{
		for (int bj = 1; bj <= h->GetNbinsY(); bj++)
		{
			double n = h->GetBinContent(bi, bj);
			double d = h_d->GetBinContent(bi, bj);
	
			double r = (d > 0) ? n / d : 0.;
			double r_u = (d > 0) ? sqrt(n * (1. - n/d)) / d : 0.;
	
			h->SetBinContent(bi, bj, r);
			h->SetBinError(bi, bj, r_u);
		}
	}

	return h;
}

//----------------------------------------------------------------------------------------------------

struct Hists
{
	TH1D *h_th_x = NULL, *h_th_y = NULL;
	TH2D *h_th_x_th_y = NULL;

	void Init()
	{
		h_th_x = new TH1D("", ";th_x", 100, -350E-6, +350E-6);
		h_th_y = new TH1D("", ";th_y", 120, -120E-6, +120E-6);

		h_th_x_th_y = new TH2D("", ";th_x;th_y", 35, -350E-6, +350E-6, 24, -120E-6, +120E-6);
	}

	void Fill(double th_x, double th_y)
	{
		if (fabs(th_y) > 35E-6 && fabs(th_y) < 100E-6)
			h_th_x->Fill(th_x);

		h_th_y->Fill(th_y);

		h_th_x_th_y->Fill(th_x, th_y);
	}

	void Write()
	{
		h_th_x->Write("h_th_x");
		h_th_y->Write("h_th_y");

		h_th_x_th_y->Write("h_th_x_th_y");
	}

	static Hists MakeRatio(Hists &num, Hists &den)
	{
		Hists rat;
		rat.h_th_x = MakeOneRatio(num.h_th_x, den.h_th_x);
		rat.h_th_y = MakeOneRatio(num.h_th_y, den.h_th_y);
		rat.h_th_x_th_y = MakeOneRatio2D(num.h_th_x_th_y, den.h_th_x_th_y);

		// determine fit range
		double m = rat.h_th_y->GetBinCenter(rat.h_th_y->GetMaximumBin());
		if (m > 0)
			rat.h_th_y->Fit("pol1", "Q", "", +30E-6, +100E-6);
		else
			rat.h_th_y->Fit("pol1", "Q", "", -100E-6, -30E-6);

		return rat;
	}
};

//----------------------------------------------------------------------------------------------------

struct StationStruct
{
	string name;

	RPStruct N_T, N_B, F_T, F_B;

	ActivityCorrelationStruct acs_T, acs_B;

	unordered_map<string, Hists> hists;

	void AssignBranches(const string &_n, TChain *ch_full, TChain *ch_ideal,
			unsigned int nt, unsigned int nb, unsigned int ft, unsigned int fb)
	{
		name = _n;

		N_T.AssignBranches("N_T", ch_full, ch_ideal, nt);
		N_B.AssignBranches("N_B", ch_full, ch_ideal, nb);
		F_T.AssignBranches("F_T", ch_full, ch_ideal, ft);
		F_B.AssignBranches("F_B", ch_full, ch_ideal, fb);

		acs_T.Init("acs_T");
		acs_B.Init("acs_B");
	}

	void FillThHist(const string &desc, double th_x, double th_y)
	{
		auto it = hists.find(desc);
		if (it == hists.end())
		{
			it = hists.insert({desc, Hists()}).first;
			it->second.Init();
		}

		it->second.Fill(th_x, th_y);
	}

	void MakeHistRatio(const string &desc_num, const string &desc_den)
	{
		if (hists.find(desc_num) == hists.end())
			return;

		if (hists.find(desc_den) == hists.end())
			return;

		const Hists &h = Hists::MakeRatio(hists[desc_num], hists[desc_den]);

		hists[desc_num + " OVER " + desc_den] = h;
	}

	void Analyze(unsigned long ev, double th_x, double th_y, CounterMap &c);

	void MakeRatios();

	void Write()
	{
		TDirectory *topDir = gDirectory;
		TDirectory *stDir = topDir->mkdir(name.c_str());

		for (auto p : hists)
		{
			gDirectory = stDir->mkdir(p.first.c_str());	
			p.second.Write();
		}

		gDirectory = stDir;
		acs_T.Write();
		acs_B.Write();

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

void StationStruct::Analyze(unsigned long /*ev*/, double th_x, double th_y, CounterMap &c)
{
	c["total"]["total"]++;

	N_T.Analyze(c);
	N_B.Analyze(c);
	F_T.Analyze(c);
	F_B.Analyze(c);

	if (N_T.prot && F_T.prot) c["N_T,F_T"]["prot & prot"]++;
	if (N_T.prot && F_T.prot && N_T.tr_val && F_T.tr_val) c["N_T,F_T"]["prot & prot & tr_val & tr_val"]++;
	if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val) c["N_T,F_T"]["prot & prot & tr_val & !tr_val"]++;
	if (N_T.prot && F_T.prot && !N_T.tr_val && F_T.tr_val) c["N_T,F_T"]["prot & prot & !tr_val & tr_val"]++;
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val) c["N_T,F_T"]["prot & prot & !tr_val & !tr_val"]++;

	if (N_B.prot && F_B.prot) c["N_B,F_B"]["prot & prot"]++;
	if (N_B.prot && F_B.prot && N_B.tr_val && F_B.tr_val) c["N_B,F_B"]["prot & prot & tr_val & tr_val"]++;
	if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val) c["N_B,F_B"]["prot & prot & tr_val & !tr_val"]++;
	if (N_B.prot && F_B.prot && !N_B.tr_val && F_B.tr_val) c["N_B,F_B"]["prot & prot & !tr_val & tr_val"]++;
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val) c["N_B,F_B"]["prot & prot & !tr_val & !tr_val"]++;
	
	if (N_T.prot && F_T.prot && N_B.tr_val) c["N_T,F_T,N_B"]["prot & prot & opp-tr_val"]++;
	if (N_T.prot && F_T.prot && F_B.tr_val) c["N_T,F_T,F_B"]["prot & prot & opp-tr_val"]++;
	if (N_T.prot && F_T.prot && N_B.tr_val && F_B.tr_val) c["N_T,F_T,N_B,F_B"]["prot & prot & opp-tr_val & opp-tr_val"]++;
	
	if (N_T.prot && F_T.prot && N_B.pl_suff) c["N_T,F_T,N_B"]["prot & prot & opp-pl_suff"]++;
	if (N_T.prot && F_T.prot && F_B.pl_suff) c["N_T,F_T,F_B"]["prot & prot & opp-pl_suff"]++;
	if (N_T.prot && F_T.prot && N_B.pl_suff && F_B.pl_suff) c["N_T,F_T,N_B,F_B"]["prot & prot & opp-pl_suff & opp-pl_suff"]++;
	
	if (N_B.prot && F_B.prot && N_T.tr_val) c["N_B,F_B,N_T"]["prot & prot & opp-tr_val"]++;
	if (N_B.prot && F_B.prot && F_T.tr_val) c["N_B,F_B,F_T"]["prot & prot & opp-tr_val"]++;
	if (N_B.prot && F_B.prot && N_T.tr_val && F_T.tr_val) c["N_B,F_B,N_T,F_T"]["prot & prot & opp-tr_val & opp-tr_val"]++;
	
	if (N_B.prot && F_B.prot && N_T.pl_suff) c["N_B,F_B,N_T"]["prot & prot & opp-pl_suff"]++;
	if (N_B.prot && F_B.prot && F_T.pl_suff) c["N_B,F_B,F_T"]["prot & prot & opp-pl_suff"]++;
	if (N_B.prot && F_B.prot && N_T.pl_suff && F_T.pl_suff) c["N_B,F_B,N_T,F_T"]["prot & prot & opp-pl_suff & opp-pl_suff"]++;

	//--------------------
	
	if (N_T.prot && F_T.prot) FillThHist("N_T.prot & N_F.prot", th_x, th_y);
	if (N_T.prot && F_T.prot && N_T.tr_val && F_T.tr_val) FillThHist("N_T.prot & F_T.prot & N_T.tr_val & F_T.tr_val", th_x, th_y);
	if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val) FillThHist("N_T.prot & F_T.prot & N_T.tr_val & !F_T.tr_val", th_x, th_y);
	if (N_T.prot && F_T.prot && !N_T.tr_val && F_T.tr_val) FillThHist("N_T.prot & F_T.prot & !N_T.tr_val & F_T.tr_val", th_x, th_y);
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val) FillThHist("N_T.prot & F_T.prot & !N_T.tr_val & !F_T.tr_val", th_x, th_y);

	if (N_B.prot && F_B.prot) FillThHist("N_B.prot & N_F.prot", th_x, th_y);
	if (N_B.prot && F_B.prot && N_B.tr_val && F_B.tr_val) FillThHist("N_B.prot & F_B.prot & N_B.tr_val & F_B.tr_val", th_x, th_y);
	if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val) FillThHist("N_B.prot & F_B.prot & N_B.tr_val & !F_B.tr_val", th_x, th_y);
	if (N_B.prot && F_B.prot && !N_B.tr_val && F_B.tr_val) FillThHist("N_B.prot & F_B.prot & !N_B.tr_val & F_B.tr_val", th_x, th_y);
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val) FillThHist("N_B.prot & F_B.prot & !N_B.tr_val & !F_B.tr_val", th_x, th_y);

	// activity correlations
	if (fabs(th_y) > 35E-6 && fabs(th_y) < 90E-6)
	{
		if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val)
			acs_T.Fill(N_T.mco, F_T.mco, N_B.mco, F_B.mco);

		if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val)
			acs_B.Fill(N_B.mco, F_B.mco, N_T.mco, F_T.mco);
	}

	// event print
	/*
	if (!N_T.prot && N_T.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, N_T\n", ev, name.c_str());

	if (!F_T.prot && F_T.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, F_T\n", ev, name.c_str());
	*/
	
	/*
	//if (name.compare("R") == 0 && fabs(th_y) > 30E-6 && fabs(th_y) < 35E-6)
	if (name.compare("R") == 0 && fabs(th_y) > 95E-6 && fabs(th_y) < 100E-6)
	{
		if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val)
			printf("\tevent_selection.insert({%lu, \"%s, N_T & !F_T\"});\n", ev, name.c_str());
		//if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val)
		//	printf("\tevent_selection.insert({%lu, \"%s, N_B & F_B\"});\n", ev, name.c_str());
	}
	*/

	//if (name.compare("R") == 0 && fabs(th_y) > 30E-6 && fabs(th_y) < 35E-6)
	//if (name.compare("R") == 0 && fabs(th_y) > 95E-6 && fabs(th_y) < 100E-6)
	//{
	//if (name.compare("L") == 0 && fabs(th_y) > 35E-6 && fabs(th_y) < 90E-6)
	//	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val)
	//		printf("\tevent_selection.insert({%lu, \"%s, !N_B & !F_B\"});\n", ev, name.c_str());
	//}

	/*
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_T & F_T\"});\n", ev, name.c_str());
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_B & F_B\"});\n", ev, name.c_str());
	*/
}

//----------------------------------------------------------------------------------------------------

void StationStruct::MakeRatios()
{
	MakeHistRatio("N_T.prot & F_T.prot & N_T.tr_val & F_T.tr_val", "N_T.prot & N_F.prot");
	MakeHistRatio("N_T.prot & F_T.prot & N_T.tr_val & !F_T.tr_val", "N_T.prot & N_F.prot");
	MakeHistRatio("N_T.prot & F_T.prot & !N_T.tr_val & F_T.tr_val", "N_T.prot & N_F.prot");
	MakeHistRatio("N_T.prot & F_T.prot & !N_T.tr_val & !F_T.tr_val", "N_T.prot & N_F.prot");

	MakeHistRatio("N_B.prot & F_B.prot & N_B.tr_val & F_B.tr_val", "N_B.prot & N_F.prot");
	MakeHistRatio("N_B.prot & F_B.prot & N_B.tr_val & !F_B.tr_val", "N_B.prot & N_F.prot");
	MakeHistRatio("N_B.prot & F_B.prot & !N_B.tr_val & F_B.tr_val", "N_B.prot & N_F.prot");
	MakeHistRatio("N_B.prot & F_B.prot & !N_B.tr_val & !F_B.tr_val", "N_B.prot & N_F.prot");
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

string replace(const string &src, const string &fnd, const string &rep)
{
	string n = src;

	size_t pos = n.find(fnd);
	if (pos != string::npos)
		n.replace(pos, fnd.size(), rep);

	return n;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// select input
	vector<string> tags;
	/*
	for (unsigned int i = 1; i <= 7; i++)
	{
		char buf[100];
		sprintf(buf, "simulations/elastic,no_details/<TYPE>_2E5_seed%i.root", i);
		tags.push_back(buf);
	}
	*/

	for (unsigned int i = 1; i <= 10; i++)
	{
		char buf[100];
		sprintf(buf, "simulations/flat,no_details/2E5_%i_<TYPE>.root", i);
		//sprintf(buf, "simulations/flat,no_details,short/1E4_%i_<TYPE>.root", i);
		tags.push_back(buf);
	}

	//tags.push_back("th_y_low_3E4_seed2");
	//tags.push_back("th_y_high_2E4_seed2");
	//tags.push_back("flat_1E5_seed1");

	// get input
	TChain *ch_full = new TChain("TotemNtuple");
	TChain *ch_ideal = new TChain("TotemNtuple");
	for (unsigned int i = 0; i < tags.size(); i++)
	{
		ch_full->Add(replace(tags[i], "<TYPE>", "ntuple_full").c_str());
		ch_ideal->Add(replace(tags[i], "<TYPE>", "ntuple_ideal").c_str());
	}

	// sanity check
	printf("full ntuple: %llu\n", ch_full->GetEntries());
	printf("ideal ntuple: %llu\n", ch_ideal->GetEntries());
	if (ch_full->GetEntries() != ch_ideal->GetEntries())
	{
		printf("ERROR: full and ideal ntuples are not compatible.\n");
		return 1;
	}
	
	// prepare output
	TFile *f_out = new TFile("analysis.root", "recreate");

	// select and link input branches
	ch_full->SetBranchStatus("*", 0);
	ch_ideal->SetBranchStatus("*", 0);

	EventMetaData *metaData = new EventMetaData();
	ch_full->SetBranchStatus("event_info.*", 1);
	ch_full->SetBranchAddress("event_info.", &metaData);

	//TriggerData *triggerData = new TriggerData();
	//ch_full->SetBranchStatus("trigger_data.*", 1);
	//ch_full->SetBranchAddress("trigger_data.", &triggerData);

	RPRootDumpReconstructedProton *simProton_R = new RPRootDumpReconstructedProton();
	ch_full->SetBranchStatus("sim_prot_right.*", 1);
	ch_full->SetBranchAddress("sim_prot_right.", &simProton_R);

	StationStruct st_L, st_R;
	st_L.AssignBranches("L", ch_full, ch_ideal, 20, 21, 24, 25);
	st_R.AssignBranches("R", ch_full, ch_ideal, 120, 121, 124, 125);

	// prepare counters and histograms
	map<string, CounterMap> counters;	// map: station label -> CounterMap
	
	// loop over events
	for (unsigned int en = 0; en < ch_full->GetEntries(); en++)
	{
		// load all data for enent en
		ch_full->GetEvent(en);
		ch_ideal->GetEvent(en);

		// get event number
		unsigned long ev = metaData->event_no;

		double th_x = simProton_R->thx;
		double th_y = simProton_R->thy;
		
		// run analysis
		st_L.Analyze(ev, th_x, th_y, counters["L"]);
		st_R.Analyze(ev, th_x, th_y, counters["R"]);
	}

	// print results
	for (map<string, CounterMap>::iterator stit = counters.begin(); stit != counters.end(); ++stit)
	{
		printf("\n\n==================== %s ====================\n", stit->first.c_str());

		for (map<string, map<string, unsigned long > >::iterator elit = stit->second.begin(); elit != stit->second.end(); ++elit)
		{
			printf("* %s\n", elit->first.c_str());

			for (map<string, unsigned long >::iterator cit = elit->second.begin(); cit != elit->second.end(); ++cit)
			{
				printf("\t%s: %lu\n", cit->first.c_str(), cit->second);
			}
		}

		CounterMap &c = stit->second;

		printf("\n");
		printf("N_T: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["N_T"]["prot & tr_val"] / c["N_T"]["prot"], 100. * c["N_T"]["prot & !tr_val"] / c["N_T"]["prot"]);
		printf("N_B: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["N_B"]["prot & tr_val"] / c["N_B"]["prot"], 100. * c["N_B"]["prot & !tr_val"] / c["N_B"]["prot"]);
		printf("F_T: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["F_T"]["prot & tr_val"] / c["F_T"]["prot"], 100. * c["F_T"]["prot & !tr_val"] / c["F_T"]["prot"]);
		printf("F_B: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["F_B"]["prot & tr_val"] / c["F_B"]["prot"], 100. * c["F_B"]["prot & !tr_val"] / c["F_B"]["prot"]);

		printf("\n");
		printf("N_T, F_T: N(prot & prot & tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & tr_val & tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & tr_val & !tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & !tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & !tr_val & tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & !tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & !tr_val & !tr_val"] / c["N_T,F_T"]["prot & prot"]);

		printf("\n");
		printf("N_B, F_B: N(prot & prot & tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & tr_val & tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & tr_val & !tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & !tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & !tr_val & tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & !tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & !tr_val & !tr_val"] / c["N_B,F_B"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_T.prot & F_T.prot & N_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B"]["prot & prot & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & F_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,F_B"]["prot & prot & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & N_B.tr_val & F_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B,F_B"]["prot & prot & opp-tr_val & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_T.prot & F_T.prot & N_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B"]["prot & prot & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & F_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,F_B"]["prot & prot & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & N_B.pl_suff & F_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B,F_B"]["prot & prot & opp-pl_suff & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_B.prot & F_B.prot & N_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T"]["prot & prot & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & F_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,F_T"]["prot & prot & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & N_T.tr_val & F_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T,F_T"]["prot & prot & opp-tr_val & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_B.prot & F_B.prot & N_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T"]["prot & prot & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & F_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,F_T"]["prot & prot & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & N_T.pl_suff & F_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T,F_T"]["prot & prot & opp-pl_suff & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
	}

	// save results
	st_L.MakeRatios();
	st_R.MakeRatios();

	st_L.Write();
	st_R.Write();

	delete f_out;
	return 0;
}
