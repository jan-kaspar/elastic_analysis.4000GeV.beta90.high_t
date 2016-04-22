#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1D.h"

#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	InitInputFiles();
	TChain *ch = new TChain(input_ntuple_name.c_str());
	for (unsigned int i = 0; i < input_files.size(); i++)
	{
		ch->Add(input_files[i].c_str());
		printf("+ %s\n", input_files[i].c_str());
	}
	printf(">> chain entries: %llu\n", ch->GetEntries());

	// select and link input branches
	ch->SetBranchStatus("*", 0);

	EventMetaData *metaData = new EventMetaData();
	ch->SetBranchStatus("event_info.*", 1);
	ch->SetBranchAddress("event_info.", &metaData);

	TriggerData *triggerData = new TriggerData();
	ch->SetBranchStatus("trigger_data.*", 1);
	ch->SetBranchAddress("trigger_data.", &triggerData);

	RPRootDumpTrackInfo *rp_L_F = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_L_N = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_R_N = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_R_F = new RPRootDumpTrackInfo();
	
	RPRootDumpTrackInfo *rp_L_FH = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_L_NH = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_R_NH = new RPRootDumpTrackInfo();
	RPRootDumpTrackInfo *rp_R_FH = new RPRootDumpTrackInfo();

	// verticals in 45
	if (diagonal == d45b_56t || diagonal == ad45b_56b)
	{
		ch->SetBranchStatus("track_rp_25.*", 1);
		ch->SetBranchAddress("track_rp_25.", &rp_L_F);

		ch->SetBranchStatus("track_rp_21.*", 1);
		ch->SetBranchAddress("track_rp_21.", &rp_L_N);
	}
	
	if (diagonal == d45t_56b || diagonal == ad45t_56t)
	{
		ch->SetBranchStatus("track_rp_24.*", 1);
		ch->SetBranchAddress("track_rp_24.", &rp_L_F);

		ch->SetBranchStatus("track_rp_20.*", 1);
		ch->SetBranchAddress("track_rp_20.", &rp_L_N);
	}
	
	// verticals in 56
	if (diagonal == d45t_56b || diagonal == ad45b_56b)
	{
		ch->SetBranchStatus("track_rp_121.*", 1);
		ch->SetBranchAddress("track_rp_121.", &rp_R_N);

		ch->SetBranchStatus("track_rp_125.*", 1);
		ch->SetBranchAddress("track_rp_125.", &rp_R_F);
	}

	if (diagonal == d45b_56t || diagonal == ad45t_56t)
	{
		ch->SetBranchStatus("track_rp_120.*", 1);
		ch->SetBranchAddress("track_rp_120.", &rp_R_N);

		ch->SetBranchStatus("track_rp_124.*", 1);
		ch->SetBranchAddress("track_rp_124.", &rp_R_F);
	}

	// horizontals
	ch->SetBranchStatus("track_rp_23.*", 1);
	ch->SetBranchAddress("track_rp_23.", &rp_L_FH);

	ch->SetBranchStatus("track_rp_22.*", 1);
	ch->SetBranchAddress("track_rp_22.", &rp_L_NH);

	ch->SetBranchStatus("track_rp_122.*", 1);
	ch->SetBranchAddress("track_rp_122.", &rp_R_NH);

	ch->SetBranchStatus("track_rp_123.*", 1);
	ch->SetBranchAddress("track_rp_123.", &rp_R_FH);

	// ouput file
	TFile *outF = new TFile((string("distill_") + argv[1] + "_new.root").c_str(), "recreate");

	// set up output tree
	EventRed ev;
	TTree *outT = new TTree("distilled", "bla");

	outT->Branch("v_L_F", &ev.h.v_L_F); outT->Branch("x_L_F", &ev.h.x_L_F); outT->Branch("y_L_F", &ev.h.y_L_F);
	outT->Branch("v_L_N", &ev.h.v_L_N); outT->Branch("x_L_N", &ev.h.x_L_N); outT->Branch("y_L_N", &ev.h.y_L_N);
	outT->Branch("v_R_N", &ev.h.v_R_N); outT->Branch("x_R_N", &ev.h.x_R_N); outT->Branch("y_R_N", &ev.h.y_R_N);
	outT->Branch("v_R_F", &ev.h.v_R_F); outT->Branch("x_R_F", &ev.h.x_R_F); outT->Branch("y_R_F", &ev.h.y_R_F);
	
	outT->Branch("v_L_FH", &ev.hH.v_L_F); outT->Branch("x_L_FH", &ev.hH.x_L_F); outT->Branch("y_L_FH", &ev.hH.y_L_F);
	outT->Branch("v_L_NH", &ev.hH.v_L_N); outT->Branch("x_L_NH", &ev.hH.x_L_N); outT->Branch("y_L_NH", &ev.hH.y_L_N);
	outT->Branch("v_R_NH", &ev.hH.v_R_N); outT->Branch("x_R_NH", &ev.hH.x_R_N); outT->Branch("y_R_NH", &ev.hH.y_R_N);
	outT->Branch("v_R_FH", &ev.hH.v_R_F); outT->Branch("x_R_FH", &ev.hH.x_R_F); outT->Branch("y_R_FH", &ev.hH.y_R_F);

	outT->Branch("timestamp", &ev.timestamp);
	outT->Branch("run_num", &ev.run_num);
	outT->Branch("bunch_num", &ev.bunch_num);
	outT->Branch("event_num", &ev.event_num);
	outT->Branch("trigger_num", &ev.trigger_num);
	outT->Branch("trigger_bits", &ev.trigger_bits);

	// loop over the chain entries
	long int evi = 0;
	for (; evi < ch->GetEntries(); evi++)
	{
		ch->GetEvent(evi);

		// copy RP data
		ev.h.v_L_F = rp_L_F->valid; ev.h.x_L_F = rp_L_F->x; ev.h.y_L_F = rp_L_F->y;
		ev.h.v_L_N = rp_L_N->valid; ev.h.x_L_N = rp_L_N->x; ev.h.y_L_N = rp_L_N->y;
		ev.h.v_R_N = rp_R_N->valid; ev.h.x_R_N = rp_R_N->x; ev.h.y_R_N = rp_R_N->y;
		ev.h.v_R_F = rp_R_F->valid; ev.h.x_R_F = rp_R_F->x; ev.h.y_R_F = rp_R_F->y;

		ev.hH.v_L_F = rp_L_FH->valid; ev.hH.x_L_F = rp_L_FH->x; ev.hH.y_L_F = rp_L_FH->y;
		ev.hH.v_L_N = rp_L_NH->valid; ev.hH.x_L_N = rp_L_NH->x; ev.hH.y_L_N = rp_L_NH->y;
		ev.hH.v_R_N = rp_R_NH->valid; ev.hH.x_R_N = rp_R_NH->x; ev.hH.y_R_N = rp_R_NH->y;
		ev.hH.v_R_F = rp_R_FH->valid; ev.hH.x_R_F = rp_R_FH->x; ev.hH.y_R_F = rp_R_FH->y;

		// event interesting?
		bool save = (ev.h.v_L_F || ev.h.v_L_N) && (ev.h.v_R_N || ev.h.v_R_F);
		//bool save = (ev.h.v_L_F && ev.h.v_L_N) && (ev.h.v_R_N && ev.h.v_R_F);
		if (!save)
			continue;

		// copy metadata
		ev.timestamp = metaData->timestamp - timestamp0;
		ev.run_num = metaData->run_no;
		ev.bunch_num = triggerData->bunch_num;
		ev.event_num = triggerData->event_num;
		ev.trigger_num = triggerData->trigger_num;
		ev.trigger_bits = triggerData->input_status_bits;

		// skip if bunch not selected
		unsigned int run = ev.run_num / 10000;
		if (SkipBunch(run, ev.bunch_num))
			continue;

		outT->Fill();
	}
	printf(">> last event number: %li\n", evi);

	// save output tree
	outT->Write();

	delete outF;
	return 0;
}
