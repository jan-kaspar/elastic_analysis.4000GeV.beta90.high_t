#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"

#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct CutPlots
{
	TH1D *h_dth_y, *h_y_dy, *h_vtx_x;

	void Init()
	{
		h_dth_y = new TH1D("", "", 100, -5.*5.6E-6, +5.*5.6E-6);
		h_y_dy = new TH1D("", "", 100, -5.*0.022, +5.*0.022);
		h_vtx_x = new TH1D("", "", 100, -5.*0.274, +5.*0.274);
	}

	void Fill(double dth_y, double y_dy, double vtx_x)
	{
		h_dth_y->Fill(dth_y);
		h_y_dy->Fill(y_dy);
		h_vtx_x->Fill(vtx_x);
	}

	void Write(int color=1)
	{
		h_dth_y->SetName("h_dth_y"); h_dth_y->SetLineColor(color); h_dth_y->Write();
		h_y_dy->SetName("h_y_dy"); h_y_dy->SetLineColor(color); h_y_dy->Write();
		h_vtx_x->SetName("h_vtx_x"); h_vtx_x->SetLineColor(color); h_vtx_x->Write();
	}
};

//----------------------------------------------------------------------------------------------------

struct DistPlots
{
	TH1D *h_th_x, *h_th_y;
	TH2D *h_th_x_th_y, *hc_th_x_th_y;

	void Init()
	{
		if (th_x_binning_edges_1d == NULL)
			BuildThBinning();

		h_th_x = new TH1D("", ";#theta_{x}", th_x_binning_n_1d, th_x_binning_edges_1d);
		h_th_y = new TH1D("", ";#theta_{y}", th_y_binning_n_1d, th_y_binning_edges_1d);

		h_th_x_th_y = new TH2D("", ";#theta_{x};#theta_{y}", th_x_binning_n_2d, th_x_binning_edges_2d,
			th_y_binning_n_2d, th_y_binning_edges_2d);
		hc_th_x_th_y = new TH2D("", ";#theta_{x};#theta_{y}", th_x_binning_n_2d_coarse,
			th_x_binning_edges_2d_coarse, th_y_binning_n_2d_coarse, th_y_binning_edges_2d_coarse);
	}

	void Fill(double th_x, double th_y)
	{
		if (fabs(th_y) > anal.th_y_lcut && fabs(th_y) < anal.th_y_hcut)
			h_th_x->Fill(th_x);

		h_th_y->Fill(th_y);

		h_th_x_th_y->Fill(th_x, th_y);
		hc_th_x_th_y->Fill(th_x, th_y);
	}

	void Write()
	{
		h_th_x->SetName("h_th_x"); h_th_x->Write();
		h_th_y->SetName("h_th_y"); h_th_y->Write();
		h_th_x_th_y->SetName("h_th_x_th_y"); h_th_x_th_y->Write();
		hc_th_x_th_y->SetName("hc_th_x_th_y"); hc_th_x_th_y->Write();
	}
};

//----------------------------------------------------------------------------------------------------

struct Plots
{
	bool initialized;
	CutPlots cuts_full, cuts_sel;
	DistPlots dist_3t, dist_4t;

	Plots() : initialized(false) {}

	void Init()
	{
		cuts_full.Init();
		cuts_sel.Init();
		dist_3t.Init();
		dist_4t.Init();
		initialized = true;
	}

	void FillCuts(double dth_y, double y_dy, double vtx_x, bool selected)
	{
		if (!initialized)
			Init();

		cuts_full.Fill(dth_y, y_dy, vtx_x);
		if (selected)
			cuts_sel.Fill(dth_y, y_dy, vtx_x);
	}

	void FillDists(double th_x, double th_y, bool four_tracks)
	{
		if (!initialized)
			Init();

		dist_3t.Fill(th_x, fabs(th_y));
		if (four_tracks)
			dist_4t.Fill(th_x, fabs(th_y));
	}

	void Write()
	{
		TDirectory *topDir = gDirectory;

		gDirectory = topDir->mkdir("cuts_full");
		cuts_full.Write(1);

		gDirectory = topDir->mkdir("cuts_sel");
		cuts_sel.Write(2);

		gDirectory = topDir->mkdir("dist_3t");
		dist_3t.Write();

		gDirectory = topDir->mkdir("dist_4t");
		dist_4t.Write();
	
		// make ratios
		TF1 *ff = new TF1("ff", "[0]");
		
		gDirectory = topDir->mkdir("th_x dependence");
		printf("\tth_x\n");
		TH1D *h_simple_ratio_vs_th_x = MakeSimpleRatio(dist_4t.h_th_x, dist_3t.h_th_x, ff, -100E-6, 100E-6, true);
		h_simple_ratio_vs_th_x->SetName("h_simple_ratio");
		h_simple_ratio_vs_th_x->Write();

		TH1D *h_refined_ratio_vs_th_x = MakeRefinedRatio(dist_4t.h_th_x, dist_3t.h_th_x, ff, -100E-6, 100E-6, true);
		h_refined_ratio_vs_th_x->SetName("h_refined_ratio");
		h_refined_ratio_vs_th_x->Write();

		gDirectory = topDir->mkdir("th_y dependence");
		printf("\tth_y\n");
		TH1D *h_simple_ratio_vs_th_y = MakeSimpleRatio(dist_4t.h_th_y, dist_3t.h_th_y, ff, anal.eff_th_y_min, 100E-6, true);
		h_simple_ratio_vs_th_y->SetName("h_simple_ratio");
		h_simple_ratio_vs_th_y->Write();

		TH1D *h_refined_ratio_vs_th_y = MakeRefinedRatio(dist_4t.h_th_y, dist_3t.h_th_y, ff, anal.eff_th_y_min, 100E-6, true);
		h_refined_ratio_vs_th_y->SetName("h_refined_ratio");
		h_refined_ratio_vs_th_y->Write();

		gDirectory = topDir->mkdir("th_x, th_y dependence");
		TH2D *h_simple_ratio_vs_th_x_th_y = MakeSimpleRatio(dist_4t.h_th_x_th_y, dist_3t.h_th_x_th_y);
		h_simple_ratio_vs_th_x_th_y->SetName("h_simple_ratio");
		h_simple_ratio_vs_th_x_th_y->Write();

		gDirectory = topDir->mkdir("th_x, th_y dependence, coarse");
		TH2D *hc_simple_ratio_vs_th_x_th_y = MakeSimpleRatio(dist_4t.hc_th_x_th_y, dist_3t.hc_th_x_th_y);
		hc_simple_ratio_vs_th_x_th_y->SetName("hc_simple_ratio");
		hc_simple_ratio_vs_th_x_th_y->Write();

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

void AnalyzeOne(const EventRed &ev, const HitData &h_al, unsigned int rpi, const string &rps,
		unsigned int nsii, double n_si, vector< vector<Plots> > &plots)
{
	// which RPs selected?
	bool L_F = (rps.find("L_F") == string::npos);
	bool L_N = (rps.find("L_N") == string::npos);
	bool R_N = (rps.find("R_N") == string::npos);
	bool R_F = (rps.find("R_F") == string::npos);

	// determine reference arm (opposite to the tested RP)
	enum { aUndefined, aL, aR } ref_arm = aUndefined;
	if (L_F && L_N) ref_arm = aL;
	if (R_F && R_N) ref_arm = aR;

	if (ref_arm == aUndefined)
		return;

	// tracks in all selected RPs?
	if ( (L_F && !ev.h.v_L_F) || (L_N && !ev.h.v_L_N) || (R_N && !ev.h.v_R_N) || (R_F && !ev.h.v_R_F) ) 
		return;

	// tracks in all 4 diagonal RPs?
	bool allTracks = ev.h.v_L_F && ev.h.v_L_N && ev.h.v_R_N && ev.h.v_R_F;
	
	// reconstruction for the reference arm
	Kinematics k_std = DoReconstruction(h_al, env);
	Kinematics k_ref;
	if (ref_arm == aL) { k_ref.th_x = k_std.th_x_L; k_ref.th_y = k_std.th_y_L;	}
	if (ref_arm == aR) { k_ref.th_x = k_std.th_x_R; k_ref.th_y = k_std.th_y_R;	}

	// kinematics for cuts
	Kinematics k_cuts = k_std;

	// reconstruction for the tested arm
	// only th_y can be reconstructed from single RP
	Kinematics k_test;
	k_test.th_y = 0.;
	if (ref_arm == aR)
	{
		if (L_N) k_test.th_y = -h_al.y_L_N / env.L_y_L_N;
		if (L_F) k_test.th_y = -h_al.y_L_F / env.L_y_L_F;
		k_cuts.th_y_L = k_test.th_y;
	}
	if (ref_arm == aL)
	{
		if (R_N) k_test.th_y = +h_al.y_R_N / env.L_y_R_N;
		if (R_F) k_test.th_y = +h_al.y_R_F / env.L_y_R_F;
		k_cuts.th_y_R = k_test.th_y;
	}

	//printf("k_ref.th_y = %E, k_test.th_y = %E | th_y_cuts_L= %E, th_y_cuts_R = %E\n",
	//.	k_ref.th_y, k_test.th_y, k_cuts.th_y_L, k_cuts.th_y_R);

	// test cuts, th_y L-R comparison (2) and in reference arm: x* reasonability (3/4) and y-th_y correlation (5/6)
	CutData cd;
	anal.n_si = n_si;
	anal.EvaluateCuts(h_al, k_cuts, cd);

	bool select = cd.ct[2];
	if (ref_arm == aL)
		select &= cd.ct[4] & cd.ct[6];
	if (ref_arm == aR)
		select &= cd.ct[3] & cd.ct[5];

	// fill cut plots
	if (ref_arm == aL)
		plots[rpi][nsii].FillCuts(cd.cv[2], cd.cv[6], cd.cv[4], select);
	
	if (ref_arm == aR)
		plots[rpi][nsii].FillCuts(cd.cv[2], cd.cv[5], cd.cv[3], select);

	if (!select)
		return;
	
	// fill distributions
	plots[rpi][nsii].FillDists(k_ref.th_x, k_ref.th_y, allTracks);
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// alignment init
	for (unsigned int i = 0; i < alignmentSources.size(); ++i)
	{
		printf("\n---------- alignment source %u ----------\n", i);
		alignmentSources[i].Init();
	}
	printf("\n\n");

	// init cut parameterss
	anal.BuildCuts();

	// print info
	printf("\n");
	printf("------------------------------ environment ------------------------------\n");
	env.Print();
	printf("\n");
	printf("------------------------------- analysis --------------------------------\n");
	anal.Print();
	printf("\n");

	// init files
	TFile *inF;
	inF = new TFile((string("distill_") + argv[1] + ".root").c_str());
	TFile *outF = new TFile((string("eff3outof4_more_cuts_") + argv[1] + ".root").c_str(), "recreate");

	// get input data
	TTree *inT = (TTree *) inF->Get("distilled");
	EventRed ev;
	inT->SetBranchAddress("timestamp", &ev.timestamp);
	inT->SetBranchAddress("run_num", &ev.run_num);
	inT->SetBranchAddress("bunch_num", &ev.bunch_num);
	inT->SetBranchAddress("event_num", &ev.event_num);
	inT->SetBranchAddress("trigger_bits", &ev.trigger_bits);

	inT->SetBranchAddress("v_L_F", &ev.h.v_L_F); inT->SetBranchAddress("x_L_F", &ev.h.x_L_F); inT->SetBranchAddress("y_L_F", &ev.h.y_L_F);
	inT->SetBranchAddress("v_L_N", &ev.h.v_L_N); inT->SetBranchAddress("x_L_N", &ev.h.x_L_N); inT->SetBranchAddress("y_L_N", &ev.h.y_L_N);
	inT->SetBranchAddress("v_R_N", &ev.h.v_R_N); inT->SetBranchAddress("x_R_N", &ev.h.x_R_N); inT->SetBranchAddress("y_R_N", &ev.h.y_R_N);
	inT->SetBranchAddress("v_R_F", &ev.h.v_R_F); inT->SetBranchAddress("x_R_F", &ev.h.x_R_F); inT->SetBranchAddress("y_R_F", &ev.h.y_R_F);

	// vector of cut sigma multiples
	vector<double> n_si;
	n_si.push_back(1.);
	n_si.push_back(3.);
	n_si.push_back(5.);

	// studied RP combinations
	vector<string> rps;
	rps.push_back("L_F");
	rps.push_back("L_N");
	rps.push_back("R_N");
	rps.push_back("R_F");

	// book histograms
	vector< vector<Plots> > plots(rps.size(), vector<Plots>(n_si.size()));

	// build histograms
	for (int ev_idx = 0; ev_idx < inT->GetEntries(); ++ev_idx)
	{
		inT->GetEntry(ev_idx);
		
		// choose the desired trigger
		if ((ev.trigger_bits & 3) == 0)	// RP trigger only
			continue;

		// remove troublesome runs
		unsigned int run = ev.run_num / 10000;
		unsigned int file = ev.run_num % 10000;
		if (SkipRun(run, file, true))
			continue;

		// select the elastic-trigger bunch(es) only
		if (SkipBunch(run, ev.bunch_num))
			continue;

		// apply fine alignment
		HitData h_al = ev.h;
		for (unsigned int i = 0; i < alignmentSources.size(); ++i)
		{
			AlignmentData alData = alignmentSources[i].Eval(ev.timestamp);
			h_al = h_al.ApplyAlignment(alData);
		}

		// loop over tested RP combinations
		for (unsigned int rpi = 0; rpi < rps.size(); rpi++)
		{
			for (unsigned int nsii = 0; nsii < n_si.size(); nsii++)
				AnalyzeOne(ev, h_al, rpi, rps[rpi], nsii, n_si[nsii], plots);
		}
	}
	
	// process and save histograms
	for (unsigned int rpi = 0; rpi < rps.size(); rpi++)
	{
		printf("\n\n------------------------- excluded RPs: %s -------------------------\n", rps[rpi].c_str());

		char buf[100];
		sprintf(buf, "excluded RPs %s", rps[rpi].c_str());
		TDirectory *rpDir = outF->mkdir(buf);

		for (unsigned int nsii = 0; nsii < n_si.size(); nsii++)
		{
			printf("* n_si: %.1f\n", n_si[nsii]);

			sprintf(buf, "n_si %.1f", n_si[nsii]);
			TDirectory *nsiDir = rpDir->mkdir(buf);
			gDirectory = nsiDir;

			plots[rpi][nsii].Write();
		}
	}
	
	return 0;
}
