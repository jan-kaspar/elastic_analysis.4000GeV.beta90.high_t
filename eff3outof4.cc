#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"

#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct HistGroup
{
	TH1D *de_th_x, *de_th_y;
	TH1D *th_x;
	TH1D *th_y;
	TH2D *th_x_th_y;

	void Init()
	{
		de_th_x = new TH1D("", ";#Delta #theta_{x}", 1000, 0., 0.);
		de_th_y = new TH1D("", ";#Delta #theta_{y}", 1000, 0., 0.);

		if (th_x_binning_edges_1d == NULL)
			BuildThBinning();

		th_x = new TH1D("", ";#theta_{x}", th_x_binning_n_1d, th_x_binning_edges_1d);
		th_y = new TH1D("", ";#theta_{y}", th_y_binning_n_1d, th_y_binning_edges_1d);

		th_x_th_y = new TH2D("", ";#theta_{x};#theta_{y}", th_x_binning_n_2d, th_x_binning_edges_2d,
			th_y_binning_n_2d, th_y_binning_edges_2d);
	}

	void Fill(double thx, double thy)
	{
		if (!th_x)
			Init();
					
		if (fabs(thy) > anal.th_y_lcut && fabs(thy) < anal.th_y_hcut)
			th_x->Fill(thx);

		th_y->Fill(fabs(thy));
		
		th_x_th_y->Fill(thx, fabs(thy));
	}

	void FillDelta(double dthx, double dthy)
	{
		if (!de_th_x)
			Init();

		de_th_x->Fill(dthx);
		de_th_y->Fill(dthy);
	}
};

//----------------------------------------------------------------------------------------------------
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

	// init files
	TFile *inF;
	inF = new TFile((string("distill_") + argv[1] + ".root").c_str());
	TFile *outF = new TFile((string("eff3outof4_") + argv[1] + ".root").c_str(), "recreate");

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

	// tolerances (= 1 sigma of left-right difference)
	//double si_de_th_x = 2000E-6;
	double si_de_th_y = 3.3E-6;

	// vector of cut sigma multiples
	vector<double> n_si = { 1., 3., 5. };

	// studied RP combinations
	vector<string> rps;
	rps.push_back("L_F");
	rps.push_back("L_N");
	rps.push_back("R_N");
	rps.push_back("R_F");

	rps.push_back("L_F,R_N");
	rps.push_back("L_F,R_F");
	rps.push_back("L_N,R_N");
	rps.push_back("L_N,R_F");

	// book histograms
	vector< vector<HistGroup> > h_sel(rps.size(), vector<HistGroup>(n_si.size())), h_full(rps.size(), vector<HistGroup>(n_si.size()));

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
			bool L_F = (rps[rpi].find("L_F") == string::npos);
			bool L_N = (rps[rpi].find("L_N") == string::npos);
			bool R_N = (rps[rpi].find("R_N") == string::npos);
			bool R_F = (rps[rpi].find("R_F") == string::npos);

			// determine reference arm (opposite to the tested RP)
			enum { aUndefined, aL, aR } ref_arm = aUndefined;
			if (L_F && L_N) ref_arm = aL;
			if (R_F && R_N) ref_arm = aR;

			// tracks in all selected RPs?
			if ( (L_F && !ev.h.v_L_F) || (L_N && !ev.h.v_L_N) || (R_N && !ev.h.v_R_N) || (R_F && !ev.h.v_R_F) ) 
				continue;
			
			// tracks in all 4/4 RPs?
			bool allTracks = ev.h.v_L_F && ev.h.v_L_N && ev.h.v_R_N && ev.h.v_R_F;

			// calculate elastic kinematics for the selected tracks

			double th_y_L_sel = 0., th_x_L_sel = 0., norm_L = 0.;
			if (L_F) { norm_L += 1.; th_x_L_sel += -h_al.x_L_F / env.L_x_L_F; th_y_L_sel += -h_al.y_L_F / env.L_y_L_F; }
			if (L_N) { norm_L += 1.; th_x_L_sel += -h_al.x_L_N / env.L_x_L_N; th_y_L_sel += -h_al.y_L_N / env.L_y_L_N; }
			th_x_L_sel /= norm_L; th_y_L_sel /= norm_L;
			
			double th_y_R_sel = 0., th_x_R_sel = 0., norm_R = 0.;
			if (R_F) { norm_R += 1.; th_x_R_sel += +h_al.x_R_F / env.L_x_R_F; th_y_R_sel += +h_al.y_R_F / env.L_y_R_F; }
			if (R_N) { norm_R += 1.; th_x_R_sel += +h_al.x_R_N / env.L_x_R_N; th_y_R_sel += +h_al.y_R_N / env.L_y_R_N; }
			th_x_R_sel /= norm_R; th_y_R_sel /= norm_R;

			//double th_x_sel = (th_x_R_sel + th_x_L_sel) / 2.;
			double th_y_sel = (th_y_R_sel + th_y_L_sel) / 2.;

			double de_th_x_sel = th_x_R_sel - th_x_L_sel;
			double de_th_y_sel = th_y_R_sel - th_y_L_sel;

			// reference quantities
			double th_x_ref = 0.;
		
			if (ref_arm == aL)
			{
				double D_x_L = - env.L_x_L_N * env.v_x_L_F + env.L_x_L_F * env.v_x_L_N;
				th_x_ref = (env.v_x_L_F * h_al.x_L_N - env.v_x_L_N * h_al.x_L_F) / D_x_L;
			}
		
			if (ref_arm == aR)
			{
				double D_x_R = env.L_x_R_N * env.v_x_R_F - env.L_x_R_F * env.v_x_R_N;
				th_x_ref = (env.v_x_R_F * h_al.x_R_N - env.v_x_R_N * h_al.x_R_F) / D_x_R;
			}

			double th_y_ref = th_y_sel;
		
			// cuts on elastic kinematics at various sigma-levels
			for (unsigned int nsi = 0; nsi < n_si.size(); nsi++)
			{
				h_sel[rpi][nsi].FillDelta(de_th_x_sel, de_th_y_sel);

				//bool cut_th_x = (fabs(de_th_x_sel) < n_si[nsi] * si_de_th_x);
				bool cut_th_y = (fabs(de_th_y_sel) < n_si[nsi] * si_de_th_y);

				if (cut_th_y)
				{
					h_sel[rpi][nsi].Fill(th_x_ref, th_y_ref);
			
					if (allTracks)
						h_full[rpi][nsi].Fill(th_x_ref, th_y_ref);
				}
			}
		}
	}

	// process and save histograms
	TF1 *ff = new TF1("ff", "[0]");

	for (unsigned int rpi = 0; rpi < rps.size(); rpi++)
	{
		printf("\n\n------------------------- excluded RPs: %s -------------------------\n", rps[rpi].c_str());

		char buf[100];
		sprintf(buf, "excluded RPs %s", rps[rpi].c_str());
		TDirectory *rpDir = outF->mkdir(buf);

		for (unsigned int nsi = 0; nsi < n_si.size(); nsi++)
		{
			printf("* n_si: %.1f\n", n_si[nsi]);

			sprintf(buf, "n_si %.1f", n_si[nsi]);
			TDirectory *nsiDir = rpDir->mkdir(buf);

			gDirectory = nsiDir;
			h_sel[rpi][nsi].de_th_x->Write("de_th_x");
			h_sel[rpi][nsi].de_th_y->Write("de_th_y");

			TCanvas *c;

			//--------------------

			gDirectory = nsiDir->mkdir("th_x dependence");
			printf("\tth_x dependence\n");
			c = new TCanvas();
			h_sel[rpi][nsi].th_x->SetName("h_sel.th_x");
			h_sel[rpi][nsi].th_x->Sumw2();
			h_sel[rpi][nsi].th_x->SetLineColor(2);
			h_sel[rpi][nsi].th_x->Draw();
			h_full[rpi][nsi].th_x->SetName("h_full.th_x");
			h_full[rpi][nsi].th_x->Sumw2();
			h_full[rpi][nsi].th_x->SetLineColor(4);
			h_full[rpi][nsi].th_x->Draw("same");
			c->Write("th_x comparison");

			TH1D *h_simple_ratio_vs_th_x = MakeSimpleRatio(h_full[rpi][nsi].th_x, h_sel[rpi][nsi].th_x, ff, -20E-6, 20E-6, false);
			h_simple_ratio_vs_th_x->SetName("h_simple_ratio.th_x");
			h_simple_ratio_vs_th_x->Write();

			TH1D *h_refined_ratio_vs_th_x = MakeRefinedRatio(h_full[rpi][nsi].th_x, h_sel[rpi][nsi].th_x, ff, -20E-6, 20E-6, false);
			h_refined_ratio_vs_th_x->SetName("h_refined_ratio.th_x");
			h_refined_ratio_vs_th_x->Write();

			//--------------------

			gDirectory = nsiDir->mkdir("th_y dependence");
			printf("\tth_y dependence\n");
			c = new TCanvas();
			h_sel[rpi][nsi].th_y->SetName("h_sel.th_y");
			h_sel[rpi][nsi].th_y->Sumw2();
			h_sel[rpi][nsi].th_y->SetLineColor(2);
			h_sel[rpi][nsi].th_y->Draw();
			h_full[rpi][nsi].th_y->SetName("h_full.th_y");
			h_full[rpi][nsi].th_y->Sumw2();
			h_full[rpi][nsi].th_y->SetLineColor(4);
			h_full[rpi][nsi].th_y->Draw("same");
			c->Write("th_y comparison");

			TH1D *h_simple_ratio_vs_th_y = MakeSimpleRatio(h_full[rpi][nsi].th_y, h_sel[rpi][nsi].th_y, ff, anal.eff_th_y_min, 100E-6, true);
			h_simple_ratio_vs_th_y->SetName("h_simple_ratio.th_y");
			h_simple_ratio_vs_th_y->Write();

			TH1D *h_refined_ratio_vs_th_y = MakeRefinedRatio(h_full[rpi][nsi].th_y, h_sel[rpi][nsi].th_y, ff, anal.eff_th_y_min, 100E-6, true);
			h_refined_ratio_vs_th_y->SetName("h_refined_ratio.th_y");
			h_refined_ratio_vs_th_y->Write();

			//--------------------

			gDirectory = nsiDir->mkdir("th_x, th_y dependence");
			h_sel[rpi][nsi].th_x_th_y->SetName("h_sel.th_x_th_y");
			h_sel[rpi][nsi].th_x_th_y->Write();
			h_full[rpi][nsi].th_x_th_y->SetName("h_full.th_x_th_y");
			h_full[rpi][nsi].th_x_th_y->Write();

			TH2D *h_simple_ratio_vs_th_x_th_y = MakeSimpleRatio(h_full[rpi][nsi].th_x_th_y, h_sel[rpi][nsi].th_x_th_y);
			h_simple_ratio_vs_th_x_th_y->SetName("h_simple_ratio_vs_th_x_th_y");
			h_simple_ratio_vs_th_x_th_y->Write();
		}
	}
	
	return 0;
}
