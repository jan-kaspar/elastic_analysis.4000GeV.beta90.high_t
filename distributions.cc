#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

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
#include "TKey.h"

#include <cmath>
#include <cassert>

using namespace std;

//----------------------------------------------------------------------------------------------------

string replace(const string &orig, const string &w, const string &r)
{
	string text = orig;

	size_t pos = 0;
	pos = text.find(w, pos);
	while (pos != string::npos)
	{
		text.replace(pos, w.size(), r);

		pos = text.find(w, pos+1);
	}

	return text;
}

//----------------------------------------------------------------------------------------------------

void ProfileToRMSGraph(TProfile *p, TGraphErrors *g)
{
	for (int bi = 1; bi <= p->GetNbinsX(); ++bi)
	{
		double c = p->GetBinCenter(bi);

		double N = p->GetBinEntries(bi);
		double Sy = p->GetBinContent(bi) * N;
		double Syy = p->GetSumw2()->At(bi);
		
		double si_sq = Syy/N - Sy*Sy/N/N;
		double si = (si_sq >= 0.) ? sqrt(si_sq) : 0.;
		double si_unc_sq = si_sq / 2. / N;	// Gaussian approximation
		double si_unc = (si_unc_sq >= 0.) ? sqrt(si_unc_sq) : 0.;

		int idx = g->GetN();
		g->SetPoint(idx, c, si);
		g->SetPointError(idx, 0., si_unc);
	}
}

//----------------------------------------------------------------------------------------------------

TGraphErrors* MakeRelDiff(TH1D *input)
{
	TGraphErrors *g_rel_diff = new TGraphErrors();

	for (int i = 1; i <= input->GetNbinsX(); ++i)
	{
		double t = input->GetBinCenter(i);
		double v = input->GetBinContent(i);
		double v_u = input->GetBinError(i);

		if (t > 0.4)
			continue;

		double ref = 530. * exp(-19.6*t);

		double rv = v / ref - 1.;
		double rv_u = v_u / ref;

		int idx = g_rel_diff->GetN();
		g_rel_diff->SetPoint(idx, t, rv);
		g_rel_diff->SetPointError(idx, 0., rv_u);
	}

	return g_rel_diff;
}

//----------------------------------------------------------------------------------------------------

unsigned int SuppressLowStatisticsBins(TProfile *p, int threshold)
{
	unsigned int reasonableBins = 0;
	for (int bi = 1; bi <= p->GetNbinsX(); ++bi)
	{
		if (p->GetBinEntries(bi) < threshold)
		{
			p->SetBinContent(bi, 0.);
			p->SetBinError(bi, 0.);
		} else
			reasonableBins++;
	}

	return reasonableBins;
}

//----------------------------------------------------------------------------------------------------

void HideLowTBins(TH1D *h, double threshold)
{
	for (int bi = 1; bi <= h->GetNbinsX(); ++bi)
	{
		if (h->GetBinCenter(bi) < threshold)
		{
			h->SetBinContent(bi, 0.);
			h->SetBinError(bi, 0.);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void DecoupleResolution(TGraphErrors *gRMS_vtx_x_vs_time, TGraphErrors *gRMS_diffLR_th_x_vs_time,
		TGraphErrors *g_beam_div_x_vs_time, TGraphErrors *g_sensor_res_x_vs_time)
{
	g_beam_div_x_vs_time->SetName("g_beam_div_x_vs_time"); g_beam_div_x_vs_time->SetTitle(";timestamp;beam divergence in x");
	g_sensor_res_x_vs_time->SetName("g_sensor_res_x_vs_time"); g_sensor_res_x_vs_time->SetTitle(";timestamp;sensor resolution in x");

	for (int i = 0; i <= gRMS_vtx_x_vs_time->GetN(); ++i)
	{
		double time=0., si_diff=0., si_vtx=0.;
		gRMS_vtx_x_vs_time->GetPoint(i, time, si_vtx);
		double si_vtx_unc = gRMS_vtx_x_vs_time->GetErrorY(i);

		gRMS_diffLR_th_x_vs_time->GetPoint(i, time, si_diff);
		double si_diff_unc = gRMS_diffLR_th_x_vs_time->GetErrorY(i);

		double si_bdx = si_vtx * sqrt(2.) / 90. * 1E-3;	// in rad
		double si_bdx_unc = si_vtx_unc * sqrt(2.) / 90. * 1E-3;

		double si_srx = sqrt(si_diff*si_diff/2. - si_bdx*si_bdx);
		double si_srx_unc = sqrt(si_diff*si_diff/4.*si_diff_unc*si_diff_unc + si_bdx*si_bdx*si_bdx_unc*si_bdx_unc) / si_srx;

		g_beam_div_x_vs_time->SetPoint(i, time, si_bdx);
		g_beam_div_x_vs_time->SetPointError(i, 0., si_bdx_unc);

		g_sensor_res_x_vs_time->SetPoint(i, time, si_srx);
		g_sensor_res_x_vs_time->SetPointError(i, 0., si_srx_unc);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc < 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined || diagonal == ad45b_56b || diagonal == ad45t_56t)
		return rcIncompatibleDiagonal;
	
	// default parameters
	unsigned int detailsLevel = 10; 	// 0: no details, 1: some details, >= 2 all details
	bool overrideCutSelection = false;	// whether the default cut selection should be overriden by the command-line selection
	string cutSelectionString;
	string outputDir = ".";
	string inputDir = ".";
	double input_n_si = 4.0;
	int time_group_divisor = 0;
	int time_group_remainder = 0;
	int event_group_divisor = 0;
	int event_group_index = 0;

	// parse command line arguments, starting from index 2
	for (int i = 2; i < argc; i++)
	{
		//printf("%u => %s\n", i, argv[i]);

		if (strcmp(argv[i], "-no-details") == 0)
		{
			detailsLevel = 0;
			continue;
		}

		if (strcmp(argv[i], "-details") == 0)
		{
			if (argc-1 > i)
				detailsLevel = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-cuts") == 0)
		{
			if (argc-1 > i)
			{
				cutSelectionString = argv[++i];
				overrideCutSelection = true;
			}
			continue;
		}
		
		if (strcmp(argv[i], "-output-dir") == 0)
		{
			if (argc-1 > i)
				outputDir = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "-input-dir") == 0)
		{
			if (argc-1 > i)
				inputDir = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "-n-si") == 0)
		{
			if (argc-1 > i)
				input_n_si = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-tg-divisor") == 0)
		{
			if (argc-1 > i)
				time_group_divisor = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-tg-remainder") == 0)
		{
			if (argc-1 > i)
				time_group_remainder = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-eg-divisor") == 0)
		{
			if (argc-1 > i)
				event_group_divisor = (int) atof(argv[++i]);
			continue;
		}


		if (strcmp(argv[i], "-eg-index") == 0)
		{
			if (argc-1 > i)
				event_group_index = (int) atof(argv[++i]);
			continue;
		}


		printf("ERROR: unknown parameter `%s'.\n", argv[i]);
		return 3;
	}
	
	printf("* detailsLevel = %u\n", detailsLevel);
	printf("* outputDir = %s\n", outputDir.c_str());
	printf("* inputDir = %s\n", inputDir.c_str());
	printf("* input n_si = %.3f\n", input_n_si);
	printf("* time_group_divisor = %i\n", time_group_divisor);
	printf("* time_group_remainder = %i\n", time_group_remainder);
	printf("* event_group_divisor = %i\n", event_group_divisor);
	printf("* event_group_index = %i\n", event_group_index);

	// select cuts
	anal.BuildCuts(); 
	anal.n_si = input_n_si;

	if (overrideCutSelection)
	{
		anal.cuts.clear();
		char buf[100];
		strcpy(buf, cutSelectionString.c_str());

		// workaround for a Bash bug? --> keeping quotes around the parameter
		for (unsigned int i = 0; i < strlen(buf); i++)
			if (buf[i] == '"')
				buf[i] = ' ';

		printf("* setting cuts from selection string [%s]\n", buf);

		if (strcmp(buf, "  ") == 0 || strlen(buf) == 0)
			printf("* cut selection list empty\n");
		else {
			char *pch = strtok(buf, ",");
			while (pch != NULL)
			{
				unsigned int cut = atoi(pch);
				if (cut < 1 || cut > 7)
				{
					printf("ERROR: invalid cut number %u.\n", cut);
					return 2;
				}
				anal.cuts.push_back(cut);
				pch = strtok (NULL, ",");
			}
		}
	}

	// choice of reconstruction formulae	
	printf(">> useOldThetaXReconstruction = %i\n", useOldThetaXReconstruction);
	
	// print info
	printf("\n");
	printf("------------------------------ environment ------------------------------\n");
	env.Print();
	printf("\n");
	printf("------------------------------- analysis --------------------------------\n");
	anal.Print();
	printf("\n");

	// alignment init
	for (unsigned int i = 0; i < alignmentSources.size(); ++i)
	{
		printf("\n---------- alignment source %u ----------\n", i);
		int r = alignmentSources[i].Init(inputDir + "/");
		if (r != 0)
			return 1;
	}
	printf("\n\n");
	
	// binnings
	vector<string> binnings;
	binnings.push_back("ub");
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("ob-2-20-0.20");
	binnings.push_back("ob-3-10-0.30");

	// init files
	TFile *inF = TFile::Open((inputDir + "/distill_" + argv[1] + ".root").c_str());
	assert(inF != NULL);

	TFile *outF = TFile::Open((outputDir+"/distributions_" + argv[1] + ".root").c_str(), "recreate");
	assert(outF != NULL);

	// get input data
	TTree *inT = (TTree *) inF->Get("distilled");
	EventRed ev;
	inT->SetBranchAddress("v_L_F", &ev.h.v_L_F); inT->SetBranchAddress("x_L_F", &ev.h.x_L_F); inT->SetBranchAddress("y_L_F", &ev.h.y_L_F);
	inT->SetBranchAddress("v_L_N", &ev.h.v_L_N); inT->SetBranchAddress("x_L_N", &ev.h.x_L_N); inT->SetBranchAddress("y_L_N", &ev.h.y_L_N);
	inT->SetBranchAddress("v_R_N", &ev.h.v_R_N); inT->SetBranchAddress("x_R_N", &ev.h.x_R_N); inT->SetBranchAddress("y_R_N", &ev.h.y_R_N);
	inT->SetBranchAddress("v_R_F", &ev.h.v_R_F); inT->SetBranchAddress("x_R_F", &ev.h.x_R_F); inT->SetBranchAddress("y_R_F", &ev.h.y_R_F);
	
	inT->SetBranchAddress("v_L_FH", &ev.hH.v_L_F); inT->SetBranchAddress("x_L_FH", &ev.hH.x_L_F); inT->SetBranchAddress("y_L_FH", &ev.hH.y_L_F);
	inT->SetBranchAddress("v_L_NH", &ev.hH.v_L_N); inT->SetBranchAddress("x_L_NH", &ev.hH.x_L_N); inT->SetBranchAddress("y_L_NH", &ev.hH.y_L_N);
	inT->SetBranchAddress("v_R_NH", &ev.hH.v_R_N); inT->SetBranchAddress("x_R_NH", &ev.hH.x_R_N); inT->SetBranchAddress("y_R_NH", &ev.hH.y_R_N);
	inT->SetBranchAddress("v_R_FH", &ev.hH.v_R_F); inT->SetBranchAddress("x_R_FH", &ev.hH.x_R_F); inT->SetBranchAddress("y_R_FH", &ev.hH.y_R_F);

	inT->SetBranchAddress("timestamp", &ev.timestamp);
	inT->SetBranchAddress("run_num", &ev.run_num);
	inT->SetBranchAddress("bunch_num", &ev.bunch_num);
	inT->SetBranchAddress("event_num", &ev.event_num);
	inT->SetBranchAddress("trigger_num", &ev.trigger_num);
	inT->SetBranchAddress("trigger_bits", &ev.trigger_bits);

	// read the selection of events
	if (cmsZeroBiasSelection)
	{
		//LoadZeroBiasEventSelection((inputDir + "/CMSZeroBiasEvents").c_str());
		
		LoadZeroBiasEventSelection((inputDir + "/CMSZeroBiasEvents_8369").c_str());
		LoadZeroBiasEventSelection((inputDir + "/CMSZeroBiasEvents_8371").c_str());
		LoadZeroBiasEventSelection((inputDir + "/CMSZeroBiasEvents_8372").c_str());

		PrepareZeroBiasEventSelection();
	}

	// load pile-up corrections
	TGraph *corrg_pileup = NULL;
	if (anal.use_pileup_efficiency_fits)
	{
		string path = inputDir + "/pileup_fit_combined.root";
		TFile *puF = TFile::Open(path.c_str());
		if (!puF)
		{
			printf("ERROR: file `%s' cannot be opened.\n", path.c_str());
			return 1;
		}

		if (diagonal == d45b_56t)
			corrg_pileup = (TGraph *) puF->Get("45b_56t/dgn");
		if (diagonal == d45t_56b)
			corrg_pileup = (TGraph *) puF->Get("45t_56b/dgn");
	}

	// load resolution fits
	map<unsigned int, TF1 *> bunch_g_th_x_diffRL_RMS, bunch_g_th_y_diffRL_RMS;

	if (anal.use_time_dependent_resolutions)
	{
		string path = inputDir + "/resolution_fit_" + argv[1] + ".root";
		TFile *resFile = TFile::Open(path.c_str());
		if (!resFile)
		{
			printf("ERROR: file `%s' with resolution fits cannot be opened.\n", path.c_str());
			return 1;
		}

		TIter next(resFile->GetListOfKeys());
		TObject *o;
		while ((o = next()))
		{
			TKey *k = (TKey *) o;
			if (k->IsFolder())
			{
				TDirectory *bunchDir = (TDirectory *) k->ReadObj();
	
				string name = k->GetName();
				unsigned int bunch = atoi(name.substr(6).c_str());

				bunch_g_th_x_diffRL_RMS[bunch] = (TF1 *) bunchDir->Get("fit_x");
				bunch_g_th_y_diffRL_RMS[bunch] = (TF1 *) bunchDir->Get("fit_y");
			}
		}

		printf("\n>> using time-dependent resolutions:\n");
		for (map<unsigned int, TF1 *>::iterator it = bunch_g_th_x_diffRL_RMS.begin(); it != bunch_g_th_x_diffRL_RMS.end(); ++it)
		{
			unsigned int b = it->first;
			printf("\tbunch %u: %p, %p\n", b, bunch_g_th_x_diffRL_RMS[b], bunch_g_th_y_diffRL_RMS[b]);
		}
	}

	// get th_y* dependent efficiency correction
	TF1 *f_3outof4_efficiency_L_F = NULL;
	TF1 *f_3outof4_efficiency_L_N = NULL;
	TF1 *f_3outof4_efficiency_R_N = NULL;
	TF1 *f_3outof4_efficiency_R_F = NULL;
	if (anal.use_3outof4_efficiency_fits)
	{
		string path = inputDir + "/eff3outof4_details_fit_old.root";
		TFile *effFile = TFile::Open(path.c_str());
		if (!effFile)
		{
			printf("ERROR: file `%s' with efficiency fits cannot be opened.\n", path.c_str());
			return 1;
		}
		
		string diagonal = argv[1];
		f_3outof4_efficiency_L_F = (TF1 *) effFile->Get( (diagonal + "/L_F/fit").c_str() );
		f_3outof4_efficiency_L_N = (TF1 *) effFile->Get( (diagonal + "/L_N/fit").c_str() );
		f_3outof4_efficiency_R_N = (TF1 *) effFile->Get( (diagonal + "/R_N/fit").c_str() );
		f_3outof4_efficiency_R_F = (TF1 *) effFile->Get( (diagonal + "/R_F/fit").c_str() );

		printf("\n>> using 3-out-of-4 fits: %p, %p, %p, %p\n",
			f_3outof4_efficiency_L_F, f_3outof4_efficiency_L_N, f_3outof4_efficiency_R_N, f_3outof4_efficiency_R_F);
	}

	// get unsmearing corrections
	printf("\n>> unsmearing_file = %s\n", unsmearing_file.c_str());
	printf(">> unsmearing_object = %s\n", unsmearing_object.c_str());

	map<unsigned int, TH1D *> map_unsmearing_correction;

	TFile *unsmearing_correction_file = TFile::Open((inputDir + "/" + unsmearing_file).c_str());
	if (!unsmearing_correction_file)
	{
		printf("ERROR: file `%s' can not be opened.\n", unsmearing_file.c_str());
	} else {
		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			string path = replace(unsmearing_object, "<binning>", binnings[bi]);
			TH1D *obj = (TH1D *) unsmearing_correction_file->Get(path.c_str());

			if (!obj)
				printf("ERROR: unsmearing correction object `%s' cannot be loaded.\n", path.c_str());

			map_unsmearing_correction[bi] = obj;
		}

		printf(">> loaded unsmearing corrections:\n");
		for (map<unsigned int, TH1D *>::iterator it = map_unsmearing_correction.begin(); it != map_unsmearing_correction.end(); ++it)
			printf("\tbinning %s: %p\n", binnings[it->first].c_str(), it->second);
	}

	// book metadata histograms
	TH1D *h_timestamp_dgn = new TH1D("h_timestamp_dgn", ";timestamp;rate   (Hz)", 97001, 16E3-0.5, 113E3+0.5);
	TH1D *h_timestamp_B0 = new TH1D("h_timestamp_B0", ";timestamp;rate   (Hz)", 97001, 16E3-0.5, 113E3+0.5);
	TH1D *h_timestamp_sel = new TH1D("h_timestamp_sel", ";timestamp;rate   (Hz)", 97001, 16E3-0.5, 113E3+0.5);

	TGraph *g_run_vs_timestamp = new TGraph(); g_run_vs_timestamp->SetName("g_run_vs_timestamp"); g_run_vs_timestamp->SetTitle(";timestamp;run");
	TGraph *g_ev_num_vs_timestamp = new TGraph(); g_ev_num_vs_timestamp->SetName("g_ev_num_vs_timestamp"); g_ev_num_vs_timestamp->SetTitle(";timestamp;ev_num");
	TGraph *g_tr_num_vs_timestamp = new TGraph(); g_tr_num_vs_timestamp->SetName("g_tr_num_vs_timestamp"); g_tr_num_vs_timestamp->SetTitle(";timestamp;tr_num");
	TGraph *g_bunch_num_vs_timestamp = new TGraph(); g_bunch_num_vs_timestamp->SetName("g_bunch_num_vs_timestamp"); g_bunch_num_vs_timestamp->SetTitle(";timestamp;bunch");
	
	TGraph *g_timestamp_vs_ev_idx_dgn = new TGraph(); g_timestamp_vs_ev_idx_dgn->SetName("g_timestamp_vs_ev_idx_dgn"); g_timestamp_vs_ev_idx_dgn->SetTitle(";event index in distilled TTree;timestamp");
	TGraph *g_timestamp_vs_ev_idx_sel = new TGraph(); g_timestamp_vs_ev_idx_sel->SetName("g_timestamp_vs_ev_idx_sel"); g_timestamp_vs_ev_idx_sel->SetTitle(";event index in distilled TTree;timestamp");

	// book hit-distribution histograms
	TH2D *h_y_L_F_vs_x_L_F_al_nosel = new TH2D("h_y_L_F_vs_x_L_F_al_nosel", ";x^{L,F};y^{L,F}", 100, -20., 20., 400, -40., +40.);
	TH2D *h_y_L_N_vs_x_L_N_al_nosel = new TH2D("h_y_L_N_vs_x_L_N_al_nosel", ";x^{L,N};y^{L,N}", 100, -20., 20., 400, -40., +40.);
	TH2D *h_y_R_N_vs_x_R_N_al_nosel = new TH2D("h_y_R_N_vs_x_R_N_al_nosel", ";x^{R,N};y^{R,N}", 100, -20., 20., 400, -40., +40.);
	TH2D *h_y_R_F_vs_x_R_F_al_nosel = new TH2D("h_y_R_F_vs_x_R_F_al_nosel", ";x^{R,F};y^{R,F}", 100, -20., 20., 400, -40., +40.);
	
	TH2D *h_y_L_F_vs_x_L_F_noal_sel = new TH2D("h_y_L_F_vs_x_L_F_noal_sel", ";x^{L,F};y^{L,F}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_L_N_vs_x_L_N_noal_sel = new TH2D("h_y_L_N_vs_x_L_N_noal_sel", ";x^{L,N};y^{L,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_N_vs_x_R_N_noal_sel = new TH2D("h_y_R_N_vs_x_R_N_noal_sel", ";x^{R,N};y^{R,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_F_vs_x_R_F_noal_sel = new TH2D("h_y_R_F_vs_x_R_F_noal_sel", ";x^{R,F};y^{R,F}", 100, -10., +10., 300, -30., +30.);
	
	TH2D *h_y_L_F_vs_x_L_F_al_sel = new TH2D("h_y_L_F_vs_x_L_F_al_sel", ";x^{L,F};y^{L,F}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_L_N_vs_x_L_N_al_sel = new TH2D("h_y_L_N_vs_x_L_N_al_sel", ";x^{L,N};y^{L,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_N_vs_x_R_N_al_sel = new TH2D("h_y_R_N_vs_x_R_N_al_sel", ";x^{R,N};y^{R,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_F_vs_x_R_F_al_sel = new TH2D("h_y_R_F_vs_x_R_F_al_sel", ";x^{R,F};y^{R,F}", 100, -10., +10., 300, -30., +30.);

	TGraph *g_y_L_F_vs_x_L_F_al_sel = new TGraph(); g_y_L_F_vs_x_L_F_al_sel->SetName("g_y_L_F_vs_x_L_F_al_sel");
	TGraph *g_y_L_N_vs_x_L_N_al_sel = new TGraph(); g_y_L_N_vs_x_L_N_al_sel->SetName("g_y_L_N_vs_x_L_N_al_sel");
	TGraph *g_y_R_N_vs_x_R_N_al_sel = new TGraph(); g_y_R_N_vs_x_R_N_al_sel->SetName("g_y_R_N_vs_x_R_N_al_sel");
	TGraph *g_y_R_F_vs_x_R_F_al_sel = new TGraph(); g_y_R_F_vs_x_R_F_al_sel->SetName("g_y_R_F_vs_x_R_F_al_sel");
	
	TH2D *h_y_L_FH_vs_x_L_FH_al_nosel = new TH2D("h_y_L_FH_vs_x_L_FH_al_nosel", ";x^{L,F};y^{L,F}", 100, 0., +20., 300, -30., +30.);
	TH2D *h_y_L_NH_vs_x_L_NH_al_nosel = new TH2D("h_y_L_NH_vs_x_L_NH_al_nosel", ";x^{L,N};y^{L,N}", 100, 0., +20., 300, -30., +30.);
	TH2D *h_y_R_NH_vs_x_R_NH_al_nosel = new TH2D("h_y_R_NH_vs_x_R_NH_al_nosel", ";x^{R,N};y^{R,N}", 100, 0., +20., 300, -30., +30.);
	TH2D *h_y_R_FH_vs_x_R_FH_al_nosel = new TH2D("h_y_R_FH_vs_x_R_FH_al_nosel", ";x^{R,F};y^{R,F}", 100, 0., +20., 300, -30., +30.);
	
	TH2D *h_y_L_FH_vs_x_L_FH_al_sel = new TH2D("h_y_L_FH_vs_x_L_FH_al_sel", ";x^{L,F};y^{L,F}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_L_NH_vs_x_L_NH_al_sel = new TH2D("h_y_L_NH_vs_x_L_NH_al_sel", ";x^{L,N};y^{L,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_NH_vs_x_R_NH_al_sel = new TH2D("h_y_R_NH_vs_x_R_NH_al_sel", ";x^{R,N};y^{R,N}", 100, -10., +10., 300, -30., +30.);
	TH2D *h_y_R_FH_vs_x_R_FH_al_sel = new TH2D("h_y_R_FH_vs_x_R_FH_al_sel", ";x^{R,F};y^{R,F}", 100, -10., +10., 300, -30., +30.);
	

	// book alignment histograms
	map<signed int, TGraph *> g_y_L_N_vs_x_L_N_sel, g_y_L_F_vs_x_L_F_sel, g_y_R_N_vs_x_R_N_sel, g_y_R_F_vs_x_R_F_sel, g_w_vs_timestamp_sel;
	map<signed int, TH1D *> tm_h_th_x_L, tm_h_th_x_R;
	map<signed int, TProfile *> tm_p_diffLR_th_x, tm_p_diffLR_th_y;
	map<signed int, TProfile *> tm_p_x_L_F_vs_th_x_L, tm_p_x_R_F_vs_th_x_R;

	// book cut histograms
	map<unsigned int, TH1D *> h_cq;
	map<unsigned int, TH2D *> h2_cq, h2_cq_full;
	map<unsigned int, TGraph *> g_cq;
	map<unsigned int, TProfile *> p_cq, p_cq_time;
	for (unsigned int i = 1; i <= anal.N_cuts; ++i)
	{
		char name[100], title[100];

		int N2 = 0;
		double x_min=0., x_max=0., y_min=0., y_max = 0., q_min = 0., q_max = 0.;

		if (i == 1) { N2 = 100; x_min = -400E-6; x_max = +400E-6; y_min = -400E-6; y_max = 400E-6; q_min = -200E-6; q_max = +200E-6; }
		if (i == 2) { N2 = 200; x_min = -150E-6; x_max = +150E-6; y_min = -150E-6; y_max = 150E-6; q_min = -50E-6; q_max = +50E-6;}
		if (i == 3) { N2 = 100; x_min = -400E-6; x_max = +400E-6; y_min = -1.; y_max = +1.; }
		if (i == 4) { N2 = 100; x_min = -400E-6; x_max = +400E-6; y_min = -1.; y_max = +1.; }
		if (i == 5) { N2 = 200; x_min = -30.; x_max = +30.; y_min = -3.; y_max = +3.; q_min = -0.3; q_max = +0.3; }
		if (i == 6) { N2 = 200; x_min = -30.; x_max = +30.; y_min = -3.; y_max = +3.; q_min = -0.3; q_max = +0.3; }
		if (i == 7) { N2 = 100; x_min = -200E-6; x_max = +200E-6; y_min = -0.2; y_max = +0.2; q_min = -200E-3; q_max = +200E-3; }
		if (i == 9) { N2 = 100; x_min = -1.; x_max = +1.; y_min = -1.; y_max = +1.; q_min = -0.04; q_max = +0.04; }
		
		sprintf(name, "h_cq%i", i); sprintf(title, ";cq%i", i); h_cq[i] = new TH1D(name, title, 200, q_min, q_max);

		sprintf(name, "h2_cq%i", i); sprintf(title, ";%s;%s", anal.cqaN[i].c_str(), anal.cqbN[i].c_str()); h2_cq[i] = new TH2D(name, title, N2, x_min, x_max, N2, y_min, y_max);
		sprintf(name, "h2_cq_full%i", i); sprintf(title, ";%s;%s", anal.cqaN[i].c_str(), anal.cqbN[i].c_str()); h2_cq_full[i] = new TH2D(name, title, N2, x_min, x_max, N2, y_min, y_max);

		sprintf(name, "g_cq%i", i); sprintf(title, ";%s;%s", anal.cqaN[i].c_str(), anal.cqbN[i].c_str()); g_cq[i] = new TGraph(); g_cq[i]->SetName(name); g_cq[i]->SetTitle(title);
		sprintf(name, "p_cq%i", i); sprintf(title, ";%s;%s", anal.cqaN[i].c_str(), anal.cqbN[i].c_str()); p_cq[i] = new TProfile(name, title, 300, 0., 0.);
		sprintf(name, "p_cq_time%i", i); sprintf(title, ";time   (s);mean of cq%i", i); p_cq_time[i] = new TProfile(name, title, 808, 16E3, 113E3);	// 2-minutes blocks
	}
	
	// book histograms for selected hits
	TProfile *p_x_vs_y_L_F = new TProfile("p_x_vs_y_L_F", ";y^{L,F};x^{L,F};", 50, 0., 0.);
	TProfile *p_x_vs_y_L_N = new TProfile("p_x_vs_y_L_N", ";y^{L,N};x^{L,N};", 50, 0., 0.);
	TProfile *p_x_vs_y_R_N = new TProfile("p_x_vs_y_R_N", ";y^{R,N};x^{R,N};", 50, 0., 0.);
	TProfile *p_x_vs_y_R_F = new TProfile("p_x_vs_y_R_F", ";y^{R,F};x^{R,F};", 50, 0., 0.);

	TProfile *p_x_vs_y_L_F_noal = new TProfile("p_x_vs_y_L_F_noal", ";y^{L,F};x^{L,F};", 50, 0., 0.);
	TProfile *p_x_vs_y_L_N_noal = new TProfile("p_x_vs_y_L_N_noal", ";y^{L,N};x^{L,N};", 50, 0., 0.);
	TProfile *p_x_vs_y_R_N_noal = new TProfile("p_x_vs_y_R_N_noal", ";y^{R,N};x^{R,N};", 50, 0., 0.);
	TProfile *p_x_vs_y_R_F_noal = new TProfile("p_x_vs_y_R_F_noal", ";y^{R,F};x^{R,F};", 50, 0., 0.);

	TH2D *h_y_L_diffFN_vs_y_L_N = new TH2D("h_y_L_diffFN_vs_y_L_N", ";y^{LN};y^{LF} - y^{LN}", 300, -30., +30., 300, -3., +3.);
	TH2D *h_y_R_diffFN_vs_y_R_N = new TH2D("h_y_R_diffFN_vs_y_R_N", ";y^{RN};y^{RF} - y^{RN}", 300, -30., +30., 300, -3., +3.);

	TH2D *h_y_L_ratioFN_vs_y_L_N = new TH2D("h_y_L_ratioFN_vs_y_L_N", ";y^{LN};y^{LF} / y^{LN}", 300, -30., +30., 300, 1.08, 1.14);
	TH2D *h_y_R_ratioFN_vs_y_R_N = new TH2D("h_y_R_ratioFN_vs_y_R_N", ";y^{RN};y^{RF} / y^{RN}", 300, -30., +30., 300, 1.08, 1.14);

	// book angluar histograms
	TH1D *th_x_diffLR = new TH1D("th_x_diffLR", ";#theta_{x}^{R} - #theta_{x}^{L}", 1000, -500E-6, +500E-6); th_x_diffLR->Sumw2();
	TH1D *th_y_diffLR = new TH1D("th_y_diffLR", ";#theta_{y}^{R} - #theta_{y}^{L}", 1000, -50E-6, +50E-6); th_y_diffLR->Sumw2();

	TH1D *th_x_diffLF = new TH1D("th_x_diffLF", ";#theta_{x}^{L} - #theta_{x}", 400, -200E-6, +200E-6); th_x_diffLF->Sumw2();
	TH1D *th_x_diffRF = new TH1D("th_x_diffRF", ";#theta_{x}^{R} - #theta_{x}", 400, -200E-6, +200E-6); th_x_diffRF->Sumw2();
	
	TH2D *h_th_x_diffLR_vs_th_x = new TH2D("h_th_x_diffLR_vs_th_x", ";#theta_{x};#theta_{x}^{R} - #theta_{x}^{L}", 100, -200E-6, +200E-6, 100, -50E-6, +50E-6);
	TH2D *h_th_y_diffLR_vs_th_y = new TH2D("h_th_y_diffLR_vs_th_y", ";#theta_{y};#theta_{y}^{R} - #theta_{y}^{L}", 120, -120E-6, +120E-6, 150, -15E-6, +15E-6);
	
	TProfile *p_th_x_diffLR_vs_th_x = new TProfile("p_th_x_diffLR_vs_th_x", ";#theta_{x};#theta_{x}^{R} - #theta_{x}^{L}", 200, -200E-6, +200E-6);
	TProfile *p_th_y_diffLR_vs_th_y = new TProfile("p_th_y_diffLR_vs_th_y", ";#theta_{y};#theta_{y}^{R} - #theta_{y}^{L}", 120, -120E-6, +120E-6);
	TProfile *p_th_y_L_diffNF_vs_th_y_L = new TProfile("p_th_y_L_diffNF_vs_th_y_L", ";#theta_{y}^{L};#theta_{y}^{LF} - #theta_{y}^{LN}", 120, -120E-6, +120E-6);
	TProfile *p_th_y_R_diffNF_vs_th_y_R = new TProfile("p_th_y_R_diffNF_vs_th_y_R", ";#theta_{y}^{R};#theta_{y}^{RF} - #theta_{y}^{RN}", 120, -120E-6, +120E-6);
	
	TH1D *th_x_diffLR_safe = new TH1D("th_x_diffLR_safe", ";#theta_{x}^{R} - #theta_{x}^{L}", 1000, -500E-6, +500E-6); th_x_diffLR_safe->Sumw2();
	TH1D *th_y_diffLR_safe = new TH1D("th_y_diffLR_safe", ";#theta_{y}^{R} - #theta_{y}^{L}", 1000, -50E-6, +50E-6); th_y_diffLR_safe->Sumw2();

	TProfile *p_th_x_L_vs_th_y_L = new TProfile("p_th_x_L_vs_th_y_L", ";#theta_{y}^{L};#theta_{x}^{L}", 30, 0., 0.);
	TProfile *p_th_x_R_vs_th_y_R = new TProfile("p_th_x_R_vs_th_y_R", ";#theta_{y}^{R};#theta_{x}^{R}", 30, 0., 0.);

	TH2D *h_th_y_L_vs_th_x_L = new TH2D("h_th_y_L_vs_th_x_L", ";#theta_{x}^{L};#theta_{y}^{L}", 150, -150E-6, +150E-6, 750, -150E-6, +150E-6);
	TH2D *h_th_y_R_vs_th_x_R = new TH2D("h_th_y_R_vs_th_x_R", ";#theta_{x}^{R};#theta_{y}^{R}", 150, -150E-6, +150E-6, 750, -150E-6, +150E-6);
	TH2D *h_th_y_vs_th_x = new TH2D("h_th_y_vs_th_x", ";#theta_{x};#theta_{y}", 100, -150E-6, +150E-6, 100, -150E-6, +150E-6);
	
	TGraph *g_th_y_L_vs_th_x_L = new TGraph(); g_th_y_L_vs_th_x_L->SetName("g_th_y_L_vs_th_x_L"); g_th_y_L_vs_th_x_L->SetTitle(";#theta_{x}^{L};#theta_{y}^{L}");
	TGraph *g_th_y_R_vs_th_x_R = new TGraph(); g_th_y_R_vs_th_x_R->SetName("g_th_y_R_vs_th_x_R"); g_th_y_R_vs_th_x_R->SetTitle(";#theta_{x}^{R};#theta_{y}^{R}");
	TGraph *g_th_y_vs_th_x = new TGraph(); g_th_y_vs_th_x->SetName("g_th_y_vs_th_x"); g_th_y_vs_th_x->SetTitle(";#theta_{x}^{L};#theta_{y}^{L}");
	
	TH2D *h_th_y_L_vs_th_y_R = new TH2D("h_th_y_L_vs_th_y_R", ";#theta_{y}^{R};#theta_{y}^{L}", 300, -150E-6, +150E-6, 300, -150E-6, +150E-6);
	TGraph *g_th_y_L_vs_th_y_R = new TGraph(); g_th_y_L_vs_th_y_R->SetName("g_th_y_L_vs_th_y_R"); g_th_y_L_vs_th_y_R->SetTitle(";#theta_{y}^{R};#theta_{y}^{L}");
	
	TH1D *h_th_x_raw = new TH1D("h_th_x_raw", ";#theta_{x}", 200, -200E-6, +200E-6);
	TH1D *h_th_y_raw = new TH1D("h_th_y_raw", ";#theta_{y}", 960, -120E-6, +120E-6);
	TH1D *h_th_y_raw_flipped = new TH1D("h_th_y_raw_flipped", ";#theta_{y}", 960, -120E-6, +120E-6);
	
	TH1D *h_th_x_L = new TH1D("h_th_x_L", ";#theta_{x}^{L}", 200, -200E-6, +200E-6); h_th_x_L->SetLineColor(4);
	TH1D *h_th_x_R = new TH1D("h_th_x_R", ";#theta_{x}^{R}", 200, -200E-6, +200E-6); h_th_x_R->SetLineColor(2);
	
	TH1D *h_th_y_L = new TH1D("h_th_y_L", ";#theta_{y}^{L}", 960, -120E-6, +120E-6); h_th_y_L->SetLineColor(4);
	TH1D *h_th_y_R = new TH1D("h_th_y_R", ";#theta_{y}^{R}", 960, -120E-6, +120E-6); h_th_y_R->SetLineColor(2);

	// book alternative angular histograms
	TH1D *h_ta_th_x = new TH1D("h_ta_th_x", ";#tau_{x} - #theta_{x}", 100, 0., 0.);
	TProfile *p_ta_th_x_vs_th_x = new TProfile("p_ta_th_x_vs_th_x", ";#theta_{y};#tau_{y} - #theta_{y}", 100, 0., 0.);

	TH1D *h_ta_th_y_R = new TH1D("h_ta_th_y_R", ";#tau_{y}^{R} - #theta_{y}^{R}", 100, 0., 0.);
	TH1D *h_ta_th_y_L = new TH1D("h_ta_th_y_L", ";#tau_{y}^{L} - #theta_{y}^{L}", 100, 0., 0.);
	TH1D *h_ta_th_y = new TH1D("h_ta_th_y", ";#tau_{y} - #theta_{y}", 100, 0., 0.);
	TH1D *h_ta_y_diffLR = new TH1D("h_ta_y_diffLR", ";#Delta^{R-L} #tau_{y}", 100, 0., 0.);
	TProfile *p_ta_th_y_vs_th_y = new TProfile("p_ta_th_y_vs_th_y", ";#theta_{y};#tau_{y} - #theta_{y}", 100, 0., 0.);

	// vertex histograms
	TH1D *h_vtx_x_L = new TH1D("h_vtx_x_L", ";{x}^{*,L}", 100, 0., 0.);
	TH1D *h_vtx_x_R = new TH1D("h_vtx_x_R", ";{x}^{*,R}", 100, 0., 0.);
	TH1D *h_vtx_x = new TH1D("h_vtx_x", ";{x}^{*}", 100, 0., 0.);
	
	TProfile *p_x_L_F_vs_th_x = new TProfile("p_x_L_F_vs_th_x", ";#theta_{x};x^{L,F}", 100, 0., 0.);
	TProfile *p_x_L_N_vs_th_x = new TProfile("p_x_L_N_vs_th_x", ";#theta_{x};x^{L,N}", 100, 0., 0.);
	TProfile *p_x_R_F_vs_th_x = new TProfile("p_x_R_F_vs_th_x", ";#theta_{x};x^{R,F}", 100, 0., 0.);
	TProfile *p_x_R_N_vs_th_x = new TProfile("p_x_R_N_vs_th_x", ";#theta_{x};x^{R,N}", 100, 0., 0.);
	
	TProfile *p_x_L_F_vs_vtx_x = new TProfile("p_x_L_F_vs_vtx_x", ";vtx_{x};x^{L,F}", 100, 0., 0.);
	TProfile *p_x_L_N_vs_vtx_x = new TProfile("p_x_L_N_vs_vtx_x", ";vtx_{x};x^{L,N}", 100, 0., 0.);
	TProfile *p_x_R_F_vs_vtx_x = new TProfile("p_x_R_F_vs_vtx_x", ";vtx_{x};x^{R,F}", 100, 0., 0.);
	TProfile *p_x_R_N_vs_vtx_x = new TProfile("p_x_R_N_vs_vtx_x", ";vtx_{x};x^{R,N}", 100, 0., 0.);
	
	TProfile *p_vtx_x_L_vs_th_x = new TProfile("p_vtx_x_L_vs_th_x", ";#theta_{x};x^{*,L}", 100, 0., 0.);
	TProfile *p_vtx_x_L_vs_th_x_L = new TProfile("p_vtx_x_L_vs_th_x_L", ";#theta_{x}^{L};x^{*,L}", 100, 0., 0.);
	TProfile *p_vtx_x_R_vs_th_x = new TProfile("p_vtx_x_R_vs_th_x", ";#theta_{x};x^{*,R}", 100, 0., 0.);
	TProfile *p_vtx_x_R_vs_th_x_R = new TProfile("p_vtx_x_R_vs_th_x_R", ";#theta_{x}^{R};x^{*,R}", 100, 0., 0.);
	
	TProfile *p_vtx_x_L_vs_th_y_L = new TProfile("p_vtx_x_L_vs_th_y_L", ";#theta_{y}^{L};x^{*,L}", 100, 0., 0.);
	TProfile *p_vtx_x_R_vs_th_y_R = new TProfile("p_vtx_x_R_vs_th_y_R", ";#theta_{y}^{R};x^{*,R}", 100, 0., 0.);
	
	TH2D *h_vtx_x_L_vs_th_x_L = new TH2D("h_vtx_x_L_vs_th_x_L", ";#theta_{x}^{L};x^{*,L}", 100, 0., 0., 100, 0., 0.);
	TH2D *h_vtx_x_R_vs_th_x_R = new TH2D("h_vtx_x_R_vs_th_x_R", ";#theta_{x}^{R};x^{*,R}", 100, 0., 0., 100, 0., 0.);

	TH1D *h_vtx_x_diffLR = new TH1D("h_vtx_x_diffLR", ";x^{*,R} - x^{*,L}", 100, 0., 0.);
	TH2D *h_vtx_x_L_vs_vtx_x_R = new TH2D("h_vtx_x_L_vs_vtx_x_R", ";x^{*,R};x^{*,L}", 100, 0., 0., 100, 0., 0.);
	TProfile *p_vtx_x_diffLR_vs_th_y = new TProfile("p_vtx_x_diffLR_vs_th_y", ";#theta_{y};x^{*,R} - x^{*,L}", 100, 0., 0.);
	TProfile *p_vtx_x_diffLR_vs_th_x = new TProfile("p_vtx_x_diffLR_vs_th_x", ";#theta_{x};x^{*,R} - x^{*,L}", 100, 0., 0.);
	TH2D *h_vtx_x_diffLR_vs_th_y = new TH2D("h_vtx_x_diffLR_vs_th_y", ";#theta_{y};x^{*,R} - x^{*,L}", 100, 0., 0., 100, 0., 0.);
	TH2D *h_vtx_x_diffLR_vs_th_x = new TH2D("h_vtx_x_diffLR_vs_th_x", ";#theta_{x};x^{*,R} - x^{*,L}", 100, 0., 0., 100, 0., 0.);
	
	TH2D *h_vtx_y_L_vs_vtx_y_R = new TH2D("h_vtx_y_L_vs_vtx_y_R", ";y^{*,R};y^{*,L}", 100, 0., 0., 100, 0., 0.);
	
	// optics histograms
	TH2D *h_x_L_F_vs_th_x_L = new TH2D("h_x_L_F_vs_th_x_L", ";#theta_{x}^{*,L};x^{L,F}", 100, -200E-6, +200E-6, 100, -1.5, +1.5);
	TH2D *h_x_R_F_vs_th_x_R = new TH2D("h_x_R_F_vs_th_x_R", ";#theta_{x}^{*,R};x^{R,F}", 100, -200E-6, +200E-6, 100, -1.5, +1.5);

	TProfile *p_x_L_N_vs_th_x_L = new TProfile("p_x_L_N_vs_th_x_L", ";#theta_{x}^{*,L};x^{L,N}", 100, 0., 0.);
	TProfile *p_x_L_F_vs_th_x_L = new TProfile("p_x_L_F_vs_th_x_L", ";#theta_{x}^{*,L};x^{L,F}", 100, 0., 0.);
	TProfile *p_x_R_N_vs_th_x_R = new TProfile("p_x_R_N_vs_th_x_R", ";#theta_{x}^{*,R};x^{R,N}", 100, 0., 0.);
	TProfile *p_x_R_F_vs_th_x_R = new TProfile("p_x_R_F_vs_th_x_R", ";#theta_{x}^{*,R};x^{R,F}", 100, 0., 0.);

	TH1D *h_ratio_ta_x_to_th_x = new TH1D("h_ratio_ta_x_to_th_x", ";#tau_{x}/#theta_{x}", 4000, 0.8, 1.2);
	TH1D *h_rel_diff_ta_x_th_x = new TH1D("h_rel_diff_ta_x_th_x", ";(#tau_{x} - #theta_{x}) / (#tau_{x} + #theta_{x}) * 2", 4000, 0., 0.);

	// time-dependence histograms
	int N_time_bins = 72;	// 1 bin = 10 minutes
	double time_min = 70.2E3, time_max = 113.4E3;

	TProfile *p_diffLR_th_x_vs_time = new TProfile("p_diffLR_th_x_vs_time", ";timestamp;mean of #Delta^{R-L}#theta_{x}", N_time_bins, time_min, time_max);
	TGraphErrors *gRMS_diffLR_th_x_vs_time = new TGraphErrors; gRMS_diffLR_th_x_vs_time->SetName("gRMS_diffLR_th_x_vs_time"); gRMS_diffLR_th_x_vs_time->SetTitle(";timestamp;RMS of #Delta^{R-L}#theta_{x}");

	TProfile *p_diffLR_th_y_vs_time = new TProfile("p_diffLR_th_y_vs_time", ";timestamp;mean of #Delta^{R-L}#theta_{y}", N_time_bins, time_min, time_max);
	TGraphErrors *gRMS_diffLR_th_y_vs_time = new TGraphErrors; gRMS_diffLR_th_y_vs_time->SetName("gRMS_diffLR_th_y_vs_time"); gRMS_diffLR_th_y_vs_time->SetTitle(";timestamp;RMS of #Delta^{R-L}#theta_{y}");
	
	TProfile *p_diffNF_th_y_L_vs_time = new TProfile("p_diffNF_th_y_L_vs_time", ";timestamp;mean of #Delta^{F-N}#theta_{y}^{L}", N_time_bins, time_min, time_max);
	TGraphErrors *gRMS_diffNF_th_y_L_vs_time = new TGraphErrors; gRMS_diffNF_th_y_L_vs_time->SetName("gRMS_diffNF_th_y_L_vs_time"); gRMS_diffNF_th_y_L_vs_time->SetTitle(";timestamp;RMS of #Delta^{F-N}#theta_{y}^{L}");
	
	TProfile *p_diffNF_th_y_R_vs_time = new TProfile("p_diffNF_th_y_R_vs_time", ";timestamp;mean of #Delta^{F-N}#theta_{y}^{R}", N_time_bins, time_min, time_max);
	TGraphErrors *gRMS_diffNF_th_y_R_vs_time = new TGraphErrors; gRMS_diffNF_th_y_R_vs_time->SetName("gRMS_diffNF_th_y_R_vs_time"); gRMS_diffNF_th_y_R_vs_time->SetTitle(";timestamp;RMS of #Delta^{F-N}#theta_{y}^{R}");
	
	TProfile *p_vtx_x_vs_time = new TProfile("p_vtx_x_vs_time", ";timestamp;mean of x^{*}", N_time_bins, time_min, time_max);
	TGraphErrors *gRMS_vtx_x_vs_time = new TGraphErrors; gRMS_vtx_x_vs_time->SetName("gRMS_vtx_x_vs_time"); gRMS_vtx_x_vs_time->SetTitle(";timestamp;RMS of x^{*}");

	map<unsigned int, TProfile *> p_bunch_diffLR_th_x_vs_time, p_bunch_diffLR_th_y_vs_time, p_bunch_vtx_x_vs_time; // map: bunch -> profile

	TProfile *p_th_x_R_vs_time = new TProfile("p_th_x_R_vs_time", ";timestamp;#theta_{x}^{R}", N_time_bins, time_min, time_max);
	TProfile *p_th_y_R_vs_time = new TProfile("p_th_y_R_vs_time", ";timestamp;#theta_{y}^{R}", N_time_bins, time_min, time_max);
	TProfile *p_th_x_L_vs_time = new TProfile("p_th_x_L_vs_time", ";timestamp;#theta_{x}^{L}", N_time_bins, time_min, time_max);
	TProfile *p_th_y_L_vs_time = new TProfile("p_th_y_L_vs_time", ";timestamp;#theta_{y}^{L}", N_time_bins, time_min, time_max);
	
	TProfile *p_input_beam_div_x_vs_time = new TProfile("p_input_beam_div_x_vs_time", ";timestamp", N_time_bins, time_min, time_max);
	TProfile *p_input_beam_div_y_vs_time = new TProfile("p_input_beam_div_y_vs_time", ";timestamp", N_time_bins, time_min, time_max);
	
	// book acceptance-correction histograms
	TProfile *p_t_ub_div_corr = new TProfile("p_t_ub_div_corr", ";t_ub_{y}", 2000., 0., 0.2);

	map<unsigned int, TH1D*> bh_t_Nev_before, bh_t_Nev_after_no_corr;
	map<unsigned int, TH1D*> bh_t_before, bh_t_after_no_corr, bh_t_after;
	map<unsigned int, TProfile*> bp_t_phi_corr, bp_t_full_corr, bp_t_full_acc;

	for (unsigned int bi = 0; bi < binnings.size(); ++bi)
	{
		unsigned int N_bins;
		double *bin_edges;
		BuildBinning(anal, binnings[bi], bin_edges, N_bins, inputDir + "/");

		bh_t_Nev_before[bi] = new TH1D("h_t_Nev_before", ";|t|;events per bin", N_bins, bin_edges); bh_t_Nev_before[bi]->Sumw2();
		bh_t_Nev_after_no_corr[bi] = new TH1D("h_t_Nev_after_no_corr", ";|t|;events per bin", N_bins, bin_edges); bh_t_Nev_after_no_corr[bi]->Sumw2();
		bh_t_before[bi] = new TH1D("h_t_before", ";|t|", N_bins, bin_edges); bh_t_before[bi]->Sumw2();
		bh_t_after[bi] = new TH1D("h_t_after", ";|t|", N_bins, bin_edges); bh_t_after[bi]->Sumw2();
		bh_t_after_no_corr[bi] = new TH1D("h_t_after_no_corr", ";|t|", N_bins, bin_edges); bh_t_after_no_corr[bi]->Sumw2();
		bp_t_phi_corr[bi] = new TProfile("p_t_phi_corr", ";t", N_bins, bin_edges, "s");
		bp_t_full_corr[bi] = new TProfile("p_t_full_corr", ";t", N_bins, bin_edges, "s");
		bp_t_full_acc[bi] = new TProfile("p_t_full_acc", ";t", N_bins, bin_edges);
	}

	TH2D *h_th_y_vs_th_x_before = new TH2D("h_th_y_vs_th_x_before", ";#theta_{x};#theta_{y}", 150, -150E-6, +150E-6, 150, -150E-6, +150E-6); h_th_y_vs_th_x_before->Sumw2();
	TH2D *h_th_y_vs_th_x_after = new TH2D("h_th_y_vs_th_x_after", ";#theta_{x};#theta_{y}", 150, -150E-6, +150E-6, 150, -150E-6, +150E-6); h_th_y_vs_th_x_after->Sumw2();
	TH2D *h_th_vs_phi_after = new TH2D("h_th_vs_phi_after", ";#phi;#theta", 16*2*4, -M_PI, M_PI, 16*4, 0E-6, 320E-6); h_th_vs_phi_after->Sumw2();
	
	TGraph *g_weight_vs_th_y = new TGraph(); g_weight_vs_th_y->SetName("g_weight_vs_th_y"); g_weight_vs_th_y->SetTitle(";#theta_{y};weight");

	// book normalization histograms
	TProfile *p_norm_corr = new TProfile("p_norm_corr", ";timestamp", 1000, 16E3, 113E3);
	TProfile *p_3outof4_corr = new TProfile("p_3outof4_corr", ";#theta_{y}^{*}", 120, -120E-6, 120E-6);

	map<unsigned int, TH1D*> bh_t_normalized;
	for (unsigned int bi = 0; bi < binnings.size(); ++bi)
	{
		bh_t_normalized[bi] = new TH1D(* bh_t_after[bi]);
		bh_t_normalized[bi]->SetName("h_t_normalized");
	}
	
	TH2D *h_th_y_vs_th_x_normalized = new TH2D("h_th_y_vs_th_x_normalized", ";#theta_{x};#theta_{y}", 320, -320E-6, +320E-6, 120, -120E-6, +120E-6); h_th_y_vs_th_x_normalized->Sumw2();

	TGraph *g_norm_corr_vs_div_corr = new TGraph(); g_norm_corr_vs_div_corr->SetName("g_norm_corr_vs_div_corr"); g_norm_corr_vs_div_corr->SetTitle(";div_corr;norm_corr");

	// book background histograms
	map<unsigned int, TH1D *> hb_cq;
	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
	{
		char name[100], title[100];
		sprintf(name, "hb_cq%i", ci); sprintf(title, ";cq%i", ci); hb_cq[ci] = new TH1D(name, title, 6000, -1000, +1000); 
	}

	TH1D *hb_th_x_L = new TH1D("hb_th_x_L", ";-#theta_{x}^{L}", 100, 0, 0); hb_th_x_L->Sumw2();
	TH1D *hb_th_x_R = new TH1D("hb_th_x_R", ";#theta_{x}^{R}", 100, 0, 0); hb_th_x_R->Sumw2();
	TH1D *hb_th_y_L = new TH1D("hb_th_y_L", ";#theta_{y}^{L}", 100, 0, 0); hb_th_y_L->Sumw2();
	TH1D *hb_th_y_R = new TH1D("hb_th_y_R", ";#theta_{y}^{R}", 100, 0, 0); hb_th_y_R->Sumw2();

	TH2D *hb_th_y_L_vs_th_x_L = new TH2D("hb_th_y_L_vs_th_x_L", ";-#theta_{x}^{L};#theta_{y}^{L}", 100, -2E-3, +2E-3, 200, -120E-6, 120E-6);
	TH2D *hb_th_y_R_vs_th_x_R = new TH2D("hb_th_y_R_vs_th_x_R", ";#theta_{x}^{R};#theta_{y}^{R}", 100, -2E-3, +2E-3, 200, -120E-6, 120E-6);
	
	TH1D *hb_th_y_diffLR = new TH1D("hb_th_y_diffLR", ";#theta_{y}^{R} - #theta_{y}^{L}", 500, 0, 0);
	TH1D *hb_th_x_diffLR = new TH1D("hb_th_x_diffLR", ";#theta_{x}^{R} - #theta_{x}^{L}", 500, 0, 0);
	
	TH2D *hb_th_y_diffLR_vs_th_y = new TH2D("hb_th_y_diffLR_vs_th_y", ";#theta_{y};#theta_{y}^{R} - #theta_{y}^{L}", 100, 0, 0, 100, 0, 0);
	TH2D *hb_th_x_diffLR_vs_th_y = new TH2D("hb_th_x_diffLR_vs_th_y", ";#theta_{y};#theta_{x}^{R} - #theta_{x}^{L}", 100, 0, 0, 100, 0, 0);
	
	TH2D *hb_th_y_L_vs_th_y_R = new TH2D("hb_th_y_L_vs_th_y_R", ";-#theta_{y}^{R};#theta_{y}^{L}", 100, 0, 0, 100, 0, 0);
	TH2D *hb_th_x_L_vs_th_x_R = new TH2D("hb_th_x_L_vs_th_x_R", ";-#theta_{x}^{R};#theta_{x}^{L}", 100, 0, 0, 100, 0, 0);
	
	TH1D *hb_th_y_6cut = new TH1D("hb_th_y_6cut", ";|#theta_{y}|", 200, 0, 100E-6); hb_th_y_6cut->Sumw2();
	TH1D *hb_th_y_6cut_cut7fail = new TH1D("hb_th_y_6cut_cut7fail", ";|#theta_{y}|", 100, 0, 100E-6); hb_th_y_6cut_cut7fail->Sumw2();

	// zero counters
	unsigned long n_ev_full = 0;
	map<unsigned int, unsigned long> n_ev_cut;
	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
		n_ev_cut[ci] = 0;

	double th_min = 1E100;

	unsigned int N_anal=0, N_anal_zeroBias=0;
	unsigned int N_zeroBias_el=0, N_zeroBias_el_RP_trig=0;
	unsigned int N_4outof4=0, N_el=0;
	unsigned int N_el_T2trig=0, N_4outof4_T2trig=0;
	unsigned int N_el_raw=0;

	// build histograms
	for (int ev_idx = 0; ev_idx < inT->GetEntries(); ++ev_idx)
	{
		inT->GetEntry(ev_idx);

		// remove troublesome runs
		unsigned int run = ev.run_num / 10000;
		unsigned int file = ev.run_num % 10000;
		if (SkipRun(run, file, true))
			continue;

		// check time - selected?
		if (anal.SkipTime(ev.timestamp))
			continue;

		if (time_group_divisor != 0)
		{
			double time_group_interval = 1.;	// s
			int time_group = int(ev.timestamp / time_group_interval);
			if ( (time_group % time_group_divisor) != time_group_remainder)
				continue;
		}

		// diagonal cut
		bool allDiagonalRPs = (ev.h.v_L_F && ev.h.v_L_N && ev.h.v_R_N && ev.h.v_R_F);
		if (!allDiagonalRPs)
			continue;

		if ((ev.trigger_bits & 64) != 0)
			N_4outof4_T2trig++;
		
		h_timestamp_dgn->Fill(ev.timestamp);
		//g_timestamp_vs_ev_idx_dgn->SetPoint(g_timestamp_vs_ev_idx_dgn->GetN(), ev_idx, ev.timestamp);

		// select the elastic-trigger bunch(es) only
		if (SkipBunch(run, ev.bunch_num))
			continue;

		// zero bias event?
		bool zero_bias_event = IsZeroBias(ev.trigger_bits, ev.run_num, ev.event_num);

		N_anal++;
		if (zero_bias_event)
			N_anal_zeroBias++;

		// apply fine alignment
		HitData h_al = ev.h;
		HitData hH_al = ev.hH;
		for (unsigned int i = 0; i < alignmentSources.size(); ++i)
		{
			AlignmentData alData = alignmentSources[i].Eval(ev.timestamp);

			h_al = h_al.ApplyAlignment(alData);
			hH_al = hH_al.ApplyInterpolatedAlignment(alData, 215.078, 219.550);
		}

		// fill pre-selection histograms
		h_y_L_N_vs_x_L_N_al_nosel->Fill(h_al.x_L_N, h_al.y_L_N);
		h_y_L_F_vs_x_L_F_al_nosel->Fill(h_al.x_L_F, h_al.y_L_F);
		h_y_R_N_vs_x_R_N_al_nosel->Fill(h_al.x_R_N, h_al.y_R_N);
		h_y_R_F_vs_x_R_F_al_nosel->Fill(h_al.x_R_F, h_al.y_R_F);

		if (ev.hH.v_L_F) h_y_L_FH_vs_x_L_FH_al_nosel->Fill(hH_al.x_L_F, hH_al.y_L_F);
		if (ev.hH.v_L_N) h_y_L_NH_vs_x_L_NH_al_nosel->Fill(hH_al.x_L_N, hH_al.y_L_N);
		if (ev.hH.v_R_N) h_y_R_NH_vs_x_R_NH_al_nosel->Fill(hH_al.x_R_N, hH_al.y_R_N);
		if (ev.hH.v_R_F) h_y_R_FH_vs_x_R_FH_al_nosel->Fill(hH_al.x_R_F, hH_al.y_R_F);

		if (detailsLevel >= 2)
		{
			h_timestamp_B0->Fill(ev.timestamp);
			g_run_vs_timestamp->SetPoint(g_run_vs_timestamp->GetN(), ev.timestamp, ev.run_num);
			//g_ev_num_vs_timestamp->SetPoint(g_ev_num_vs_timestamp->GetN(), ev.timestamp, ev.event_num);
			//g_tr_num_vs_timestamp->SetPoint(g_tr_num_vs_timestamp->GetN(), ev.timestamp, ev.trigger_num);
		}

		Kinematics k = DoReconstruction(h_al, env);

		// alternative theta_x reconstruction
		// this is used in the old 90m analysis (for the PRL publication)
		double ta_x = (h_al.x_R_F - h_al.x_R_N - h_al.x_L_F + h_al.x_L_N) / (env.L_x_R_F - env.L_x_R_N) / 2.;

		// alternative theta_y reconstruction
		double ta_y_L = - (h_al.y_L_F - h_al.y_L_N) / (env.L_y_L_F - env.L_y_L_N);
		double ta_y_R = + (h_al.y_R_F - h_al.y_R_N) / (env.L_y_R_F - env.L_y_R_N);
		double ta_y = (ta_y_L + ta_y_R) / 2.;

		// cut evaluation
		CutData cd;
		bool select = anal.EvaluateCuts(h_al, k, cd);
		//bool bckg = !select;

		// increment counters
		n_ev_full++;
		for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
		{
			if (cd.ct[ci])
				n_ev_cut[ci]++;
		}
		
		// fill background distributions
		for (unsigned int qi = 1; qi <= anal.N_cuts; ++qi)
		{
			bool reduced_select = true;
			for (unsigned int i = 0; i < anal.cuts.size(); i++)
			{
				unsigned int ci = anal.cuts[i];
				if (ci != qi)
					reduced_select &= cd.ct[ci];
			}

			if (reduced_select) {
				hb_cq[qi]->Fill(cd.cv[qi] / anal.csi[qi]);
			}
		}

		/*
		if (bckg) {
			hb_th_y_L->Fill(th_y_L);
			hb_th_y_R->Fill(th_y_R);
			hb_th_x_L->Fill(-th_x_L);
			hb_th_x_R->Fill(th_x_R);
		
			hb_th_y_L_vs_th_x_L->Fill(-th_x_L, th_y_L);
			hb_th_y_R_vs_th_x_R->Fill(th_x_R, th_y_R);

			hb_th_y_diffLR->Fill(th_y_R - th_y_L);
			hb_th_x_diffLR->Fill(th_x_R - th_x_L);
	
			hb_th_y_diffLR_vs_th_y->Fill(th_y, th_y_R - th_y_L);
			hb_th_x_diffLR_vs_th_y->Fill(th_y, th_x_R - th_x_L);
	
			hb_th_y_L_vs_th_y_R->Fill(th_y_R, th_y_L);
			hb_th_x_L_vs_th_x_R->Fill(th_x_R, th_x_L);
		}
		*/
		
		// fill no-cut histograms
		for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
		{
			//h2_cq_full[ci]->Fill(ccb[ci]*cqa[ci] - cca[ci]*cqb[ci], cca[ci]*cqa[ci] + ccb[ci]*cqb[ci] + ccc[ci]);
			h2_cq_full[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
		}
		
		/*
		bool select_6cut = true;
		for (unsigned int i = 0; i < cuts.size(); i++)
			if (cuts[i] != 7)
				select_6cut &= ct[cuts[i]];

		if (select_6cut)
			hb_th_y_6cut->Fill(fabs(th_y));

		if (select_6cut && !ct[7])
			hb_th_y_6cut_cut7fail->Fill(fabs(th_y));
		*/

		N_4outof4++;

		// elastic cut
		if (!select)
			continue;

		N_el++;

		// fill zero-bias plots
		if (zero_bias_event)
		{
			N_zeroBias_el++;
			if ((ev.trigger_bits & 1) != 0)
				N_zeroBias_el_RP_trig++;
		}

		// elastic events with T2 trigger?
		if ((ev.trigger_bits & 64) != 0)
			N_el_T2trig++;
		
		// choose the desired trigger
		if ((ev.trigger_bits & 3) == 0)	// RP trigger only
			continue;
		
		// check event group
		N_el_raw++;

		if (event_group_divisor > 0)
		{
			int event_group = (N_el_raw-1) / event_group_divisor;
			if (event_group < event_group_index)
				continue;
			if (event_group > event_group_index)
				break;
		}

		// determine normalization factors (luminosity + corrections)
		double inefficiency_3outof4 = anal.inefficiency_3outof4;
		if (anal.use_3outof4_efficiency_fits)
		{
			inefficiency_3outof4 = 0.;
			inefficiency_3outof4 += 1. - f_3outof4_efficiency_L_F->Eval(k.th_y * 1E6);
			inefficiency_3outof4 += 1. - f_3outof4_efficiency_L_N->Eval(k.th_y * 1E6);
			inefficiency_3outof4 += 1. - f_3outof4_efficiency_R_N->Eval(k.th_y * 1E6);
			inefficiency_3outof4 += 1. - f_3outof4_efficiency_R_F->Eval(k.th_y * 1E6);
		}

		double inefficiency_shower_near = anal.inefficiency_shower_near;

		double inefficiency_pile_up = anal.inefficiency_pile_up;
		if (anal.use_pileup_efficiency_fits)
			inefficiency_pile_up = corrg_pileup->Eval(ev.timestamp);

		double norm_corr =
			1./(1. - (inefficiency_3outof4 + inefficiency_shower_near))
			* 1./(1. - inefficiency_pile_up)
			* 1./(1. - anal.inefficiency_trigger)
			* 1./(1. - anal.inefficiency_DAQ);

		p_norm_corr->Fill(ev.timestamp, norm_corr);
		p_3outof4_corr->Fill(k.th_y, inefficiency_3outof4);

		double normalization = anal.bckg_corr * norm_corr / anal.L_int;

		// data for alignment
		// (SHOULD use hit positions WITHOUT alignment corrections, i.e. ev.h)
		signed int period = int((ev.timestamp - anal.alignment_t0) / anal.alignment_ts);
		if (detailsLevel >= 2)
		{
			if (g_y_L_N_vs_x_L_N_sel.find(period) == g_y_L_N_vs_x_L_N_sel.end())
			{
				g_y_L_N_vs_x_L_N_sel[period] = new TGraph();
				g_y_L_F_vs_x_L_F_sel[period] = new TGraph();
				g_y_R_N_vs_x_R_N_sel[period] = new TGraph();
				g_y_R_F_vs_x_R_F_sel[period] = new TGraph();
				g_w_vs_timestamp_sel[period] = new TGraph();

				tm_h_th_x_L[period] = new TH1D("", ";#theta_{x}^{L}", 100, -150E-6, +150E-6);
				tm_h_th_x_R[period] = new TH1D("", ";#theta_{x}^{R}", 100, -150E-6, +150E-6);

				tm_p_diffLR_th_x[period] = new TProfile("", ";#theta_{y}   (#murad);#Delta^{R-L} #theta_{y}   (#murad)", 300, -300E-6, +300E-6);
				tm_p_diffLR_th_y[period] = new TProfile("", ";#theta_{y}   (#murad);#Delta^{R-L} #theta_{y}   (#murad)", 200, -100E-6, +100E-6);

				tm_p_x_L_F_vs_th_x_L[period] = new TProfile("", ";#theta_{x}^{L}   (#murad);x^{LF}   (mm)", 200, -200E-6, +200E-6);
				tm_p_x_R_F_vs_th_x_R[period] = new TProfile("", ";#theta_{x}^{R}   (#murad);x^{RF}   (mm)", 200, -200E-6, +200E-6);
			}
	
			g_y_L_N_vs_x_L_N_sel[period]->SetPoint(g_y_L_N_vs_x_L_N_sel[period]->GetN(), ev.h.x_L_N, ev.h.y_L_N);
			g_y_L_F_vs_x_L_F_sel[period]->SetPoint(g_y_L_F_vs_x_L_F_sel[period]->GetN(), ev.h.x_L_F, ev.h.y_L_F);
			g_y_R_N_vs_x_R_N_sel[period]->SetPoint(g_y_R_N_vs_x_R_N_sel[period]->GetN(), ev.h.x_R_N, ev.h.y_R_N);
			g_y_R_F_vs_x_R_F_sel[period]->SetPoint(g_y_R_F_vs_x_R_F_sel[period]->GetN(), ev.h.x_R_F, ev.h.y_R_F);
			g_w_vs_timestamp_sel[period]->SetPoint(g_w_vs_timestamp_sel[period]->GetN(), ev.timestamp, norm_corr);
		}

		// fill raw histograms
		h_timestamp_sel->Fill(ev.timestamp);
		if (detailsLevel >= 2)
		{
			g_timestamp_vs_ev_idx_sel->SetPoint(g_timestamp_vs_ev_idx_sel->GetN(), ev_idx, ev.timestamp);
			g_bunch_num_vs_timestamp->SetPoint(g_bunch_num_vs_timestamp->GetN(), ev.timestamp, ev.bunch_num);
		}
		
		for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
		{
			h_cq[ci]->Fill(cd.cv[ci]);
			h2_cq[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
			//h2_cq[ci]->Fill(ccb[ci]*cqa[ci] - cca[ci]*cqb[ci], cca[ci]*cqa[ci] + ccb[ci]*cqb[ci] + ccc[ci]);
			if (detailsLevel >= 2)
				g_cq[ci]->SetPoint(g_cq[ci]->GetN(), cd.cqa[ci], cd.cqb[ci]);
			p_cq[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
			p_cq_time[ci]->Fill(ev.timestamp, cd.cv[ci]);
		}

		h_y_L_N_vs_x_L_N_noal_sel->Fill(ev.h.x_L_N, ev.h.y_L_N);
		h_y_L_F_vs_x_L_F_noal_sel->Fill(ev.h.x_L_F, ev.h.y_L_F);
		h_y_R_N_vs_x_R_N_noal_sel->Fill(ev.h.x_R_N, ev.h.y_R_N);
		h_y_R_F_vs_x_R_F_noal_sel->Fill(ev.h.x_R_F, ev.h.y_R_F);

		h_y_L_N_vs_x_L_N_al_sel->Fill(h_al.x_L_N, h_al.y_L_N);
		h_y_L_F_vs_x_L_F_al_sel->Fill(h_al.x_L_F, h_al.y_L_F);
		h_y_R_N_vs_x_R_N_al_sel->Fill(h_al.x_R_N, h_al.y_R_N);
		h_y_R_F_vs_x_R_F_al_sel->Fill(h_al.x_R_F, h_al.y_R_F);

		{
			int idx = g_y_L_F_vs_x_L_F_al_sel->GetN();
			if (idx < 100000)
			{
				g_y_L_F_vs_x_L_F_al_sel->SetPoint(idx, h_al.x_L_F, h_al.y_L_F);
				g_y_L_N_vs_x_L_N_al_sel->SetPoint(idx, h_al.x_L_N, h_al.y_L_N);
				g_y_R_N_vs_x_R_N_al_sel->SetPoint(idx, h_al.x_R_N, h_al.y_R_N);
				g_y_R_F_vs_x_R_F_al_sel->SetPoint(idx, h_al.x_R_F, h_al.y_R_F);
			}
		}

		if (ev.hH.v_L_F) h_y_L_FH_vs_x_L_FH_al_sel->Fill(hH_al.x_L_F, hH_al.y_L_F);
		if (ev.hH.v_L_N) h_y_L_NH_vs_x_L_NH_al_sel->Fill(hH_al.x_L_N, hH_al.y_L_N);
		if (ev.hH.v_R_N) h_y_R_NH_vs_x_R_NH_al_sel->Fill(hH_al.x_R_N, hH_al.y_R_N);
		if (ev.hH.v_R_F) h_y_R_FH_vs_x_R_FH_al_sel->Fill(hH_al.x_R_F, hH_al.y_R_F);

		p_x_vs_y_L_F->Fill(h_al.y_L_F, h_al.x_L_F);
		p_x_vs_y_L_N->Fill(h_al.y_L_N, h_al.x_L_N);
		p_x_vs_y_R_N->Fill(h_al.y_R_N, h_al.x_R_N);
		p_x_vs_y_R_F->Fill(h_al.y_R_F, h_al.x_R_F);
		
		p_x_vs_y_L_F_noal->Fill(ev.h.y_L_F, ev.h.x_L_F);
		p_x_vs_y_L_N_noal->Fill(ev.h.y_L_N, ev.h.x_L_N);
		p_x_vs_y_R_N_noal->Fill(ev.h.y_R_N, ev.h.x_R_N);
		p_x_vs_y_R_F_noal->Fill(ev.h.y_R_F, ev.h.x_R_F);

		h_y_L_diffFN_vs_y_L_N->Fill(h_al.y_L_N, h_al.y_L_F - h_al.y_L_N);
		h_y_R_diffFN_vs_y_R_N->Fill(h_al.y_R_N, h_al.y_R_F - h_al.y_R_N);

		h_y_L_ratioFN_vs_y_L_N->Fill(h_al.y_L_N, h_al.y_L_F / h_al.y_L_N);
		h_y_R_ratioFN_vs_y_R_N->Fill(h_al.y_R_N, h_al.y_R_F / h_al.y_R_N);

		th_x_diffLR->Fill(k.th_x_R - k.th_x_L);
		th_y_diffLR->Fill(k.th_y_R - k.th_y_L);

		th_x_diffLF->Fill(k.th_x_L - k.th_x);
		th_x_diffRF->Fill(k.th_x_R - k.th_x);

		h_th_x_diffLR_vs_th_x->Fill(k.th_x, k.th_x_R - k.th_x_L);
		h_th_y_diffLR_vs_th_y->Fill(k.th_y, k.th_y_R - k.th_y_L);

		p_th_x_diffLR_vs_th_x->Fill(k.th_x, k.th_x_R - k.th_x_L);
		p_th_y_diffLR_vs_th_y->Fill(k.th_y, k.th_y_R - k.th_y_L);
		p_th_y_L_diffNF_vs_th_y_L->Fill(k.th_y_L, k.th_y_L_F - k.th_y_L_N);
		p_th_y_R_diffNF_vs_th_y_R->Fill(k.th_y_R, k.th_y_R_F - k.th_y_R_N);
		
		double safe_th_y_min = anal.th_y_lcut + 5E-6;
		double safe_th_y_max = anal.th_y_hcut - 5E-6;
		if (fabs(k.th_y) > safe_th_y_min && fabs(k.th_y) < safe_th_y_max)
		{
			th_y_diffLR_safe->Fill(k.th_y_R - k.th_y_L);
			th_x_diffLR_safe->Fill(k.th_x_R - k.th_x_L);
		}
	
		p_th_x_L_vs_th_y_L->Fill(k.th_y_L, k.th_x_L);
		p_th_x_R_vs_th_y_R->Fill(k.th_y_R, k.th_x_R);
	
		h_th_y_L_vs_th_x_L->Fill(k.th_x_L, k.th_y_L);
		h_th_y_R_vs_th_x_R->Fill(k.th_x_R, k.th_y_R);
		h_th_y_vs_th_x->Fill(k.th_x, k.th_y);
	
		if (detailsLevel >= 1)
		{
			g_th_y_L_vs_th_x_L->SetPoint(g_th_y_L_vs_th_x_L->GetN(), k.th_x_L, k.th_y_L);
			g_th_y_R_vs_th_x_R->SetPoint(g_th_y_R_vs_th_x_R->GetN(), k.th_x_R, k.th_y_R);
		}
	
		h_th_y_L_vs_th_y_R->Fill(k.th_y_R, k.th_y_L);
		if (detailsLevel >= 2)
			g_th_y_L_vs_th_y_R->SetPoint(g_th_y_L_vs_th_y_R->GetN(), k.th_y_R, k.th_y_L);

		h_th_x_raw->Fill(k.th_x);
		h_th_y_raw->Fill(k.th_y);
		h_th_y_raw_flipped->Fill(-k.th_y);
		
		h_th_x_L->Fill(k.th_x_L);
		h_th_x_R->Fill(k.th_x_R);
		
		h_th_y_L->Fill(k.th_y_L);
		h_th_y_R->Fill(k.th_y_R);

		h_vtx_x_L->Fill(k.vtx_x_L);
		h_vtx_x_R->Fill(k.vtx_x_R);
		h_vtx_x->Fill(k.vtx_x);

		p_x_L_F_vs_th_x->Fill(k.th_x, h_al.x_L_F);
		p_x_L_N_vs_th_x->Fill(k.th_x, h_al.x_L_N);
		p_x_R_F_vs_th_x->Fill(k.th_x, h_al.x_R_F);
		p_x_R_N_vs_th_x->Fill(k.th_x, h_al.x_R_N);
		
		p_x_L_F_vs_vtx_x->Fill(k.vtx_x, h_al.x_L_F);
		p_x_L_N_vs_vtx_x->Fill(k.vtx_x, h_al.x_L_N);
		p_x_R_F_vs_vtx_x->Fill(k.vtx_x, h_al.x_R_F);
		p_x_R_N_vs_vtx_x->Fill(k.vtx_x, h_al.x_R_N);
	
		p_vtx_x_L_vs_th_x->Fill(k.th_x, k.vtx_x_L);
		p_vtx_x_L_vs_th_x_L->Fill(k.th_x_L, k.vtx_x_L);
		p_vtx_x_R_vs_th_x->Fill(k.th_x, k.vtx_x_R);
		p_vtx_x_R_vs_th_x_R->Fill(k.th_x_R, k.vtx_x_R);
		
		p_vtx_x_L_vs_th_y_L->Fill(k.th_y_L, k.vtx_x_L);
		p_vtx_x_R_vs_th_y_R->Fill(k.th_y_R, k.vtx_x_R);
	
		h_vtx_x_L_vs_th_x_L->Fill(k.th_x_L, k.vtx_x_L);
		h_vtx_x_R_vs_th_x_R->Fill(k.th_x_R, k.vtx_x_R);

		h_vtx_x_diffLR->Fill(k.vtx_x_R - k.vtx_x_L);
		h_vtx_x_L_vs_vtx_x_R->Fill(k.vtx_x_R, k.vtx_x_L);
		p_vtx_x_diffLR_vs_th_y->Fill(k.th_y, k.vtx_x_R - k.vtx_x_L);
		p_vtx_x_diffLR_vs_th_x->Fill(k.th_x, k.vtx_x_R - k.vtx_x_L);
		h_vtx_x_diffLR_vs_th_y->Fill(k.th_y, k.vtx_x_R - k.vtx_x_L);
		h_vtx_x_diffLR_vs_th_x->Fill(k.th_x, k.vtx_x_R - k.vtx_x_L);
		
		h_vtx_y_L_vs_vtx_y_R->Fill(k.vtx_y_R, k.vtx_y_L);

		h_x_L_F_vs_th_x_L->Fill(k.th_x_L, h_al.x_L_F);
		h_x_R_F_vs_th_x_R->Fill(k.th_x_R, h_al.x_R_F);

		p_x_L_N_vs_th_x_L->Fill(k.th_x_L, h_al.x_L_N);
		p_x_L_F_vs_th_x_L->Fill(k.th_x_L, h_al.x_L_F);
		p_x_R_N_vs_th_x_R->Fill(k.th_x_R, h_al.x_R_N);
		p_x_R_F_vs_th_x_R->Fill(k.th_x_R, h_al.x_R_F);

		h_ratio_ta_x_to_th_x->Fill(ta_x/k.th_x);
		h_rel_diff_ta_x_th_x->Fill((ta_x - k.th_x) / (ta_x + k.th_x) * 2.);

		h_ta_th_x->Fill(ta_x - k.th_x);
		p_ta_th_x_vs_th_x->Fill(k.th_x, ta_x - k.th_x);

		h_ta_th_y_R->Fill(ta_y_R - k.th_y_R);
		h_ta_th_y_L->Fill(ta_y_L - k.th_y_L);
		h_ta_th_y->Fill(ta_y - k.th_y);
		h_ta_y_diffLR->Fill(ta_y_R - ta_y_L);
		p_ta_th_y_vs_th_y->Fill(k.th_y, ta_y - k.th_y);

		if (fabs(k.th_y) > safe_th_y_min && fabs(k.th_y) < safe_th_y_max)
		{
			p_diffLR_th_x_vs_time->Fill(ev.timestamp, k.th_x_R - k.th_x_L);
			p_diffLR_th_y_vs_time->Fill(ev.timestamp, k.th_y_R - k.th_y_L);

			p_diffNF_th_y_L_vs_time->Fill(ev.timestamp, k.th_y_L_F - k.th_y_L_N);
			p_diffNF_th_y_R_vs_time->Fill(ev.timestamp, k.th_y_R_F - k.th_y_R_N);

			p_vtx_x_vs_time->Fill(ev.timestamp, k.vtx_x);

			if (p_bunch_diffLR_th_x_vs_time.find(ev.bunch_num) == p_bunch_diffLR_th_x_vs_time.end())
			{
				p_bunch_diffLR_th_x_vs_time[ev.bunch_num] = new TProfile("p_diffLR_th_x_vs_time", ";timestamp;mean of #Delta^{R-L}#theta_{x}", N_time_bins, time_min, time_max);
				p_bunch_diffLR_th_y_vs_time[ev.bunch_num] = new TProfile("p_diffLR_th_y_vs_time", ";timestamp;mean of #Delta^{R-L}#theta_{y}", N_time_bins, time_min, time_max);
				p_bunch_vtx_x_vs_time[ev.bunch_num] = new TProfile("p_vtx_x_vs_time", ";timestamp;mean of x^{*}", N_time_bins, time_min, time_max);
			}
			
			p_bunch_diffLR_th_x_vs_time[ev.bunch_num]->Fill(ev.timestamp, k.th_x_R - k.th_x_L);
            p_bunch_diffLR_th_y_vs_time[ev.bunch_num]->Fill(ev.timestamp, k.th_y_R - k.th_y_L);
			p_bunch_vtx_x_vs_time[ev.bunch_num]->Fill(ev.timestamp, k.vtx_x);

			if (detailsLevel >= 2)
			{
				tm_h_th_x_L[period]->Fill(k.th_x_L);
				tm_h_th_x_R[period]->Fill(k.th_x_R);

				tm_p_diffLR_th_x[period]->Fill(k.th_x, k.th_x_R - k.th_x_L);
				tm_p_diffLR_th_y[period]->Fill(k.th_y, k.th_y_R - k.th_y_L);

				tm_p_x_L_F_vs_th_x_L[period]->Fill(k.th_x_L, h_al.x_L_F);
				tm_p_x_R_F_vs_th_x_R[period]->Fill(k.th_x_R, h_al.x_R_F);
			}
		}

		p_th_x_R_vs_time->Fill(ev.timestamp, k.th_x_R);
		p_th_y_R_vs_time->Fill(ev.timestamp, k.th_y_R);
		p_th_x_L_vs_time->Fill(ev.timestamp, k.th_x_L);
		p_th_y_L_vs_time->Fill(ev.timestamp, k.th_y_L);

		// set time-dependent resolutions
		if (anal.use_time_dependent_resolutions)
		{
			anal.si_th_x_1arm = bunch_g_th_x_diffRL_RMS[ev.bunch_num]->Eval(ev.timestamp) / sqrt(2.);
			anal.si_th_y_1arm = bunch_g_th_y_diffRL_RMS[ev.bunch_num]->Eval(ev.timestamp) / sqrt(2.);
		}

		p_input_beam_div_x_vs_time->Fill(ev.timestamp, anal.si_th_x_1arm);
		p_input_beam_div_y_vs_time->Fill(ev.timestamp, anal.si_th_y_1arm);

		// calculate acceptance divergence correction
		double phi_corr = 0., div_corr = 0.;
		
		bool skip = CalculateAcceptanceCorrections(th_y_sign, k, anal, phi_corr, div_corr);

		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			bh_t_Nev_before[bi]->Fill(k.t, 1.);
			bh_t_before[bi]->Fill(k.t, 1.);
		}

		h_th_y_vs_th_x_before->Fill(k.th_x, k.th_y, 1.);

		if (skip)
			continue;

		double corr = div_corr * phi_corr;

		th_min = min(th_min, k.th);

		// fill acceptance histograms
		p_t_ub_div_corr->Fill(k.t_y, div_corr);

		for (unsigned int bi = 0; bi < binnings.size(); bi++)
		{
			bh_t_Nev_after_no_corr[bi]->Fill(k.t, 1.);

			bh_t_after_no_corr[bi]->Fill(k.t, 1.);
			bh_t_after[bi]->Fill(k.t, corr);

			bp_t_phi_corr[bi]->Fill(k.t, phi_corr);
			bp_t_full_corr[bi]->Fill(k.t, corr);
			bp_t_full_acc[bi]->Fill(k.t, 1./corr);
		}
		
		h_th_y_vs_th_x_after->Fill(k.th_x, k.th_y, div_corr);
		h_th_vs_phi_after->Fill(k.phi, k.th, div_corr);

		if (detailsLevel >= 2)
			g_weight_vs_th_y->SetPoint(g_weight_vs_th_y->GetN(), k.th_y, div_corr);

		// apply normalization
		for (unsigned int bi = 0; bi < binnings.size(); bi++)
			bh_t_normalized[bi]->Fill(k.t, corr * normalization);
		
		h_th_y_vs_th_x_normalized->Fill(k.th_x, k.th_y, div_corr * normalization);

		if (detailsLevel >= 1)
		{
			g_th_y_vs_th_x->SetPoint(g_th_y_vs_th_x->GetN(), k.th_x, k.th_y);
			g_norm_corr_vs_div_corr->SetPoint(g_norm_corr_vs_div_corr->GetN(), div_corr, norm_corr);
		}
	}
	
	printf("---------------------------- after event loop ---------------------------\n");

	printf(">> th_min = %E\n", th_min);

	printf("\n");
	printf("N_anal = %u\n", N_anal);
	printf("N_anal_zeroBias = %u\n", N_anal_zeroBias);
	printf("N_zeroBias_el = %u\n", N_zeroBias_el);
	printf("N_zeroBias_el_RP_trig = %u\n", N_zeroBias_el_RP_trig);

	printf("N_4outof4 = %u\n", N_4outof4);
	printf("N_el = %u\n", N_el);
	printf("N_el_T2trig = %u\n", N_el_T2trig);
	printf("N_4outof4_T2trig = %u\n", N_4outof4_T2trig);

	// derived plots
	TGraphErrors *th_y_sigmaLR_vs_th_y = new TGraphErrors();
	th_y_sigmaLR_vs_th_y->SetName("th_y_sigmaLR_vs_th_y");
	th_y_sigmaLR_vs_th_y->SetTitle(";#theta_{y};RMS of #Delta^{R-L} #theta_{y}");
	ProfileToRMSGraph(p_th_y_diffLR_vs_th_y, th_y_sigmaLR_vs_th_y);

	TGraphErrors *th_x_sigmaLR_vs_th_x = new TGraphErrors();
	th_x_sigmaLR_vs_th_x->SetName("th_x_sigmaLR_vs_th_x");
	th_x_sigmaLR_vs_th_x->SetTitle(";#theta_{x};RMS of #Delta^{R-L} #theta_{x}");
	ProfileToRMSGraph(p_th_x_diffLR_vs_th_x, th_x_sigmaLR_vs_th_x);

	// normalize histograms
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		bh_t_before[bi]->Scale(1., "width");
		bh_t_after_no_corr[bi]->Scale(1., "width");
		bh_t_after[bi]->Scale(1., "width");
	
		bh_t_normalized[bi]->Scale(1., "width");
	}
	
	//h_th_y_vs_th_x_normalized->Scale(1., "width");
	
	hb_th_x_L->Scale(1., "width");
	hb_th_x_R->Scale(1., "width");
	hb_th_y_L->Scale(1., "width");
	hb_th_y_R->Scale(1., "width");
	
	th_y_diffLR->Scale(1., "width");
	th_x_diffLR->Scale(1., "width");
	th_y_diffLR_safe->Scale(1., "width");
	th_x_diffLR_safe->Scale(1., "width");

	// hide bins with high uncertainty
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
		HideLowTBins(bh_t_normalized[bi], anal.t_min_fit);

	// fit histograms
	double th_y_low_bound  = (diagonal == d45b_56t) ? anal.th_y_lcut + 5E-6 : -(anal.th_y_hcut - 5E-6);
	double th_y_high_bound = (diagonal == d45b_56t) ? anal.th_y_hcut - 5E-6 : -(anal.th_y_lcut + 5E-6);
	printf("\n* th_y fit bounds: from %E to %E\n", th_y_low_bound, th_y_high_bound);

	printf("\n* fitting p_th_x_diffLR_vs_th_x\n");
	p_th_x_diffLR_vs_th_x->Fit("pol1", "", "", -100E-6, +100E-6);
	printf("\n* fitting p_th_y_diffLR_vs_th_y\n");
	p_th_y_diffLR_vs_th_y->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	printf("\n* fitting p_th_y_L_diffNF_vs_th_y_L\n");
	p_th_y_L_diffNF_vs_th_y_L->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	printf("\n* fitting p_th_y_R_diffNF_vs_th_y_R\n");
	p_th_y_R_diffNF_vs_th_y_R->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	
	printf("\n* fitting p_x_L_F_vs_th_x\n");
	p_x_L_F_vs_th_x->Fit("pol1");
	printf("\n* fitting p_x_L_N_vs_th_x\n");
	p_x_L_N_vs_th_x->Fit("pol1");
	printf("\n* fitting p_x_R_F_vs_th_x\n");
	p_x_R_F_vs_th_x->Fit("pol1");
	printf("\n* fitting p_x_R_N_vs_th_x\n");
	p_x_R_N_vs_th_x->Fit("pol1");
	
	printf("\n* fitting p_x_L_F_vs_vtx_x\n");
	p_x_L_F_vs_vtx_x->Fit("pol1");
	printf("\n* fitting p_x_L_N_vs_vtx_x\n");
	p_x_L_N_vs_vtx_x->Fit("pol1");
	printf("\n* fitting p_x_R_F_vs_vtx_x\n");
	p_x_R_F_vs_vtx_x->Fit("pol1");
	printf("\n* fitting p_x_R_N_vs_vtx_x\n");
	p_x_R_N_vs_vtx_x->Fit("pol1");

	printf("\n* fitting p_vtx_x_L_vs_th_x_L\n");
	p_vtx_x_L_vs_th_x_L->Fit("pol1", "", "", -120E-6, +120E-6);
	printf("* fitting p_vtx_x_L_vs_th_x\n");
	p_vtx_x_L_vs_th_x->Fit("pol1", "", "", -120E-6, +120E-6);
	printf("* fitting p_vtx_x_R_vs_th_x_R\n");
	p_vtx_x_R_vs_th_x_R->Fit("pol1", "", "", -120E-6, +120E-6);
	printf("* fitting p_vtx_x_R_vs_th_x\n");
	p_vtx_x_R_vs_th_x->Fit("pol1", "", "", -120E-6, +120E-6);
	
	printf("* fitting p_vtx_x_diffLR_vs_th_x\n");
	p_vtx_x_diffLR_vs_th_x->Fit("pol1", "", "", -120E-6, +120E-6);

	printf("* fitting p_vtx_x_L_vs_th_y_L\n");
	p_vtx_x_L_vs_th_y_L->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	printf("* fitting p_vtx_x_R_vs_th_y_R\n");
	p_vtx_x_R_vs_th_y_R->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	
	printf("* fitting p_th_x_L_vs_th_y_L\n");
	p_th_x_L_vs_th_y_L->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);
	printf("* fitting p_th_x_R_vs_th_y_R\n");
	p_th_x_R_vs_th_y_R->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);

	printf("* fitting h_th_x_R\n");
	h_th_x_R->Fit("gaus", "", "", -50E-6, +50E-6);
	printf("* fitting h_th_x_L\n");
	h_th_x_L->Fit("gaus", "", "", -50E-6, +50E-6);

	printf("* fitting p_x_L_N_vs_th_x_L\n");
	p_x_L_N_vs_th_x_L->Fit("pol1");
	printf("* fitting p_x_L_F_vs_th_x_L\n");
	p_x_L_F_vs_th_x_L->Fit("pol1");
	printf("* fitting p_x_R_N_vs_th_x_R\n");
	p_x_R_N_vs_th_x_R->Fit("pol1");
	printf("* fitting p_x_R_F_vs_th_x_R\n");
	p_x_R_F_vs_th_x_R->Fit("pol1");

	printf("* fitting p_diffLR_th_x_vs_time\n");
	p_diffLR_th_x_vs_time->Fit("pol1");
	printf("* fitting p_diffLR_th_y_vs_time\n");
	p_diffLR_th_y_vs_time->Fit("pol1");

	printf("* fitting p_th_x_R_vs_time\n");
	p_th_x_R_vs_time->Fit("pol1");
	printf("* fitting p_th_y_R_vs_time\n");
	p_th_y_R_vs_time->Fit("pol1");
	printf("* fitting p_th_x_L_vs_time\n");
	p_th_x_L_vs_time->Fit("pol1");
	printf("* fitting p_th_y_L_vs_time\n");
	p_th_y_L_vs_time->Fit("pol1");

	// double-gauss fit for background histograms
	//TF1 *dg = new TF1("dg", "[0]*exp(-(x-[1])*(x-[1])/2/[2]/[2]) + [3]*exp(-(x-[4])*(x-[4])/2/[5]/[5])");
	for (map<unsigned int, TH1D *>::iterator it = hb_cq.begin(); it != hb_cq.end(); ++it) {
		/*
		printf("* Fitting background projection %u\n", it->first);
		dg->SetParameters(5E3, 0., 1., 1E2, 0., 3.);
		it->second->Fit(dg, "", "", -9., +9.);

		double A = dg->GetParameter(0), mu = dg->GetParameter(1), si = dg->GetParameter(2);
		double signal_int = A*si*sqrt(2.*M_PI)/2. * ( TMath::Erf((+3. - mu)/sqrt(2.)/si) - TMath::Erf((-3. - mu)/sqrt(2.)/si) );
		
		A = dg->GetParameter(3); mu = dg->GetParameter(4); si = dg->GetParameter(5);
		double background_int = A*si*sqrt(2.*M_PI)/2. * ( TMath::Erf((+3. - mu)/sqrt(2.)/si) - TMath::Erf((-3. - mu)/sqrt(2.)/si) );
		printf("\tsignal int = %E\n", signal_int);
		printf("\tbackground int = %E\n", background_int);
		printf("\tsignal + background int = %E\n", signal_int + background_int);
		printf("\tbackground/signal = %E\n", background_int / signal_int);
		*/
	}
	
	printf("* fitting p_x_vs_y_L_F\n");
	p_x_vs_y_L_F->Fit("pol1", "", "", -15., +15.);
	printf("* fitting p_x_vs_y_L_N\n");
	p_x_vs_y_L_N->Fit("pol1", "", "", -15., +15.);
	printf("* fitting p_x_vs_y_R_N\n");
	p_x_vs_y_R_N->Fit("pol1", "", "", -15., +15.);
	printf("* fitting p_x_vs_y_R_F\n");
	p_x_vs_y_R_F->Fit("pol1", "", "", -15., +15.);

	th_y_diffLR_safe->Fit("gaus");
	th_x_diffLR_safe->Fit("gaus");

	printf("* fitting p_ta_th_x_vs_th_x\n");
	p_ta_th_x_vs_th_x->Fit("pol1", "", "", -100E-6, +100E-6);
	
	printf("* fitting p_ta_th_y_vs_th_y\n");
	p_ta_th_y_vs_th_y->Fit("pol1", "", "", th_y_low_bound, th_y_high_bound);

	// apply unfolding correction
	map<unsigned int, TH1D *>  bh_t_normalized_unsmeared;
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		bh_t_normalized_unsmeared[bi] = new TH1D(* bh_t_normalized[bi]);
		bh_t_normalized_unsmeared[bi]->SetName("h_t_normalized_unsmeared");

		map<unsigned int, TH1D *>::iterator it = map_unsmearing_correction.find(bi);
		TH1D *h_corr = (it == map_unsmearing_correction.end()) ? NULL : it->second;

		for (int bin = 1; bin <= bh_t_normalized_unsmeared[bi]->GetNbinsX(); ++bin)
		{
			//double c = bh_t_normalized_unsmeared[bi]->GetBinCenter(bin);
			double v = bh_t_normalized_unsmeared[bi]->GetBinContent(bin);
			double v_u = bh_t_normalized_unsmeared[bi]->GetBinError(bin);

			double corr = (h_corr == NULL) ? 0. : h_corr->GetBinContent(bin);

			bh_t_normalized_unsmeared[bi]->SetBinContent(bin, v * corr);
			bh_t_normalized_unsmeared[bi]->SetBinError(bin, v_u * corr);
		}
	}

	// derived plots
	map<unsigned int, TGraphErrors *> bh_t_normalized_rel_diff, bh_t_normalized_unsmeared_rel_diff;
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		bh_t_normalized_rel_diff[bi] = MakeRelDiff(bh_t_normalized[bi]);
		bh_t_normalized_rel_diff[bi]->SetName("h_t_eb_normalized_rel_diff");
	
		bh_t_normalized_unsmeared_rel_diff[bi] = MakeRelDiff(bh_t_normalized_unsmeared[bi]);
		bh_t_normalized_unsmeared_rel_diff[bi]->SetName("h_t_eb_normalized_unsmeared_rel_diff");
	}

	// save histograms
	TCanvas *c;
	
	gDirectory = outF->mkdir("metadata");
	if (detailsLevel >= 2)
	{
		h_timestamp_dgn->Write();
		h_timestamp_B0->SetLineColor(4);
		h_timestamp_B0->Write();
		h_timestamp_sel->SetLineColor(2);
		h_timestamp_sel->Write();

		c = new TCanvas("rate cmp");
		h_timestamp_dgn->Draw();
		h_timestamp_B0->Draw("sames");
		h_timestamp_sel->Draw("sames");
		c->Write();
	
		//g_timestamp_vs_ev_idx_dgn->Write();
		g_timestamp_vs_ev_idx_sel->Write();

		g_run_vs_timestamp->Write();
		//g_ev_num_vs_timestamp->Write();
		//g_tr_num_vs_timestamp->Write();
		g_bunch_num_vs_timestamp->Write();
	}
	
	TDirectory *hitDistDir = outF->mkdir("hit distributions");
	gDirectory = hitDistDir->mkdir("vertical, aligned, before selection");
	h_y_L_F_vs_x_L_F_al_nosel->Write();
	h_y_L_N_vs_x_L_N_al_nosel->Write();
	h_y_R_N_vs_x_R_N_al_nosel->Write();
	h_y_R_F_vs_x_R_F_al_nosel->Write();
	
	gDirectory = hitDistDir->mkdir("vertical, not aligned, after selection");
	h_y_L_F_vs_x_L_F_noal_sel->Write();
	h_y_L_N_vs_x_L_N_noal_sel->Write();
	h_y_R_N_vs_x_R_N_noal_sel->Write();
	h_y_R_F_vs_x_R_F_noal_sel->Write();
	
	gDirectory = hitDistDir->mkdir("vertical, aligned, after selection");
	h_y_L_F_vs_x_L_F_al_sel->Write();
	h_y_L_N_vs_x_L_N_al_sel->Write();
	h_y_R_N_vs_x_R_N_al_sel->Write();
	h_y_R_F_vs_x_R_F_al_sel->Write();

	g_y_L_F_vs_x_L_F_al_sel->Write();
	g_y_L_N_vs_x_L_N_al_sel->Write();
	g_y_R_N_vs_x_R_N_al_sel->Write();
	g_y_R_F_vs_x_R_F_al_sel->Write();
	
	gDirectory = hitDistDir->mkdir("horizontal, aligned, before selection");
	h_y_L_FH_vs_x_L_FH_al_nosel->Write();
	h_y_L_NH_vs_x_L_NH_al_nosel->Write();
	h_y_R_NH_vs_x_R_NH_al_nosel->Write();
	h_y_R_FH_vs_x_R_FH_al_nosel->Write();
	
	gDirectory = hitDistDir->mkdir("horizontal, aligned, after selection");
	h_y_L_FH_vs_x_L_FH_al_sel->Write();
	h_y_L_NH_vs_x_L_NH_al_sel->Write();
	h_y_R_NH_vs_x_R_NH_al_sel->Write();
	h_y_R_FH_vs_x_R_FH_al_sel->Write();

	TDirectory *alDir = outF->mkdir("alignment");

	TF1 *ff = new TF1("ff", "[0] + [1]*x");

	TGraphErrors* g_ext_diffLR_th_x_vs_time = new TGraphErrors(); g_ext_diffLR_th_x_vs_time->SetName("g_ext_diffLR_th_x_vs_time");
	TGraphErrors* g_ext_diffLR_th_y_vs_time = new TGraphErrors(); g_ext_diffLR_th_y_vs_time->SetName("g_ext_diffLR_th_y_vs_time");

	TGraphErrors* g_L_L_F_vs_time = new TGraphErrors(); g_L_L_F_vs_time->SetName("g_L_L_F_vs_time");
	TGraphErrors* g_L_R_F_vs_time = new TGraphErrors(); g_L_R_F_vs_time->SetName("g_L_R_F_vs_time");

	for (map<signed int, TGraph *>::iterator pid = g_y_L_N_vs_x_L_N_sel.begin(); pid != g_y_L_N_vs_x_L_N_sel.end(); ++pid)
	{
		signed int period = pid->first;

		char buf[100];
		sprintf(buf, "%i", period);
		gDirectory = alDir->mkdir(buf);

		g_y_L_N_vs_x_L_N_sel[period]->SetName("g_y_L_N_vs_x_L_N_sel"); g_y_L_N_vs_x_L_N_sel[period]->Write();
		g_y_L_F_vs_x_L_F_sel[period]->SetName("g_y_L_F_vs_x_L_F_sel"); g_y_L_F_vs_x_L_F_sel[period]->Write();
		g_y_R_N_vs_x_R_N_sel[period]->SetName("g_y_R_N_vs_x_R_N_sel"); g_y_R_N_vs_x_R_N_sel[period]->Write();
		g_y_R_F_vs_x_R_F_sel[period]->SetName("g_y_R_F_vs_x_R_F_sel"); g_y_R_F_vs_x_R_F_sel[period]->Write();
		g_w_vs_timestamp_sel[period]->SetName("g_w_vs_timestamp_sel"); g_w_vs_timestamp_sel[period]->Write();

		tm_h_th_x_L[period]->Write("tm_h_th_x_L");
		tm_h_th_x_R[period]->Write("tm_h_th_x_R");

		TProfile *p_th_x = tm_p_diffLR_th_x[period];
		p_th_x->SetName("p_diffLR_th_x");
		unsigned int reasonableBins_x = SuppressLowStatisticsBins(p_th_x, 5);
		p_th_x->Fit(ff, "Q");
		p_th_x->Write();
		double v_x = ff->GetParameter(0), u_x = ff->GetParError(0);

		TProfile *p_th_y = tm_p_diffLR_th_y[period];
		p_th_y->SetName("p_diffLR_th_y");
		unsigned int reasonableBins_y = SuppressLowStatisticsBins(p_th_y, 5);
		p_th_y->Fit(ff, "Q");
		p_th_y->Write();
		double v_y = ff->GetParameter(0), u_y = ff->GetParError(0);

		double time = anal.alignment_t0 + (period + 0.5) * anal.alignment_ts;
		double time_beg = anal.alignment_t0 + (period + 0.0) * anal.alignment_ts;
		double time_end = anal.alignment_t0 + (period + 1.0) * anal.alignment_ts;
		printf("period %u: from %.0f to %0.f\n", period, time_beg, time_end);
		
		if (reasonableBins_x > 9)
		{
			int idx = g_ext_diffLR_th_x_vs_time->GetN();
			g_ext_diffLR_th_x_vs_time->SetPoint(idx, time, v_x);
			g_ext_diffLR_th_x_vs_time->SetPointError(idx, 0., u_x);
		}
		
		if (reasonableBins_y > 9)
		{
			int idx = g_ext_diffLR_th_y_vs_time->GetN();
			g_ext_diffLR_th_y_vs_time->SetPoint(idx, time, v_y);
			g_ext_diffLR_th_y_vs_time->SetPointError(idx, 0., u_y);
		}


		TProfile *p_L_L = tm_p_x_L_F_vs_th_x_L[period];
		p_L_L->SetName("p_x_L_F_vs_th_x_L");
		unsigned int reasonableBins_L_L = SuppressLowStatisticsBins(p_L_L, 5);
		p_L_L->Fit(ff, "Q");
		p_L_L->Write();
		double L_L = ff->GetParameter(1), u_L_L = ff->GetParError(1);

		if (reasonableBins_L_L > 9)
		{
			int idx = g_L_L_F_vs_time->GetN();
			g_L_L_F_vs_time->SetPoint(idx, time, -L_L);
			g_L_L_F_vs_time->SetPointError(idx, 0., u_L_L);
		}

		TProfile *p_L_R = tm_p_x_R_F_vs_th_x_R[period];
		p_L_R->SetName("p_x_R_F_vs_th_x_R");
		unsigned int reasonableBins_L_R = SuppressLowStatisticsBins(p_L_R, 5);
		p_L_R->Fit(ff, "Q");
		p_L_R->Write();
		double L_R = ff->GetParameter(1), u_L_R = ff->GetParError(1);

		if (reasonableBins_L_R > 9)
		{
			int idx = g_L_R_F_vs_time->GetN();
			g_L_R_F_vs_time->SetPoint(idx, time, +L_R);
			g_L_R_F_vs_time->SetPointError(idx, 0., u_L_R);
		}
	}

	TDirectory *cutDir = outF->mkdir("elastic cuts");
	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
	{
		char buf[100];
		sprintf(buf, "cut %u", ci);
		gDirectory = cutDir->mkdir(buf);

		printf("* RMS of cut distribution %u\n", ci);
		printf("%E\n", h_cq[ci]->GetRMS());

		printf("* fitting cut distribution %u\n", ci);
		h_cq[ci]->Fit("gaus", "", "", -3., +3.);

		h_cq[ci]->Write();
		h2_cq[ci]->Write();
		h2_cq_full[ci]->Write();
		g_cq[ci]->Write();
		
		p_cq[ci]->Write();

		p_cq_time[ci]->Write();

		TGraphErrors *g_cq_time = new TGraphErrors();
		sprintf(buf, "g_cq_RMS%u", ci); g_cq_time->SetName(buf);
		sprintf(buf, ";time   (s);RMS of cq%u", ci); g_cq_time->SetTitle(buf);
		ProfileToRMSGraph(p_cq_time[ci], g_cq_time);
		g_cq_time->Write();

		h2_cq_full[ci]->Draw();

		sprintf(buf, "plot_before_cq%u", ci);
		c = new TCanvas(buf);
		c->SetLogz(1);
		h2_cq_full[ci]->Draw("colz");
		double lim = h2_cq_full[ci]->GetXaxis()->GetXmax();

		double qa[2] = {-lim, +lim};
		double qbp[2]= {(+anal.n_si*anal.csi[ci] - anal.cca[ci]*qa[0] - anal.ccc[ci])/anal.ccb[ci],
			(+anal.n_si*anal.csi[ci] - anal.cca[ci]*qa[1] - anal.ccc[ci])/anal.ccb[ci]};
		double qbm[2]= {(-anal.n_si*anal.csi[ci] - anal.cca[ci]*qa[0] - anal.ccc[ci])/anal.ccb[ci],
			(-anal.n_si*anal.csi[ci] - anal.cca[ci]*qa[1] - anal.ccc[ci])/anal.ccb[ci]};
		TGraph *gP = new TGraph(2, qa, qbp); gP->Draw("l");
		TGraph *gM = new TGraph(2, qa, qbm); gM->Draw("l");
		c->Write();
		
		sprintf(buf, "plot_after_cq%u", ci);
		c = new TCanvas(buf);
		c->SetLogz(1);
		h2_cq[ci]->Draw("colz");
		gP = new TGraph(2, qa, qbp); gP->Draw("l");
		gM = new TGraph(2, qa, qbm); gM->Draw("l");
		gP->Draw("l");
		gM->Draw("l");
		c->Write();
	}

	gDirectory = outF->mkdir("selected - hits");
	p_x_vs_y_L_F->Write();
	p_x_vs_y_L_N->Write();
	p_x_vs_y_R_N->Write();
	p_x_vs_y_R_F->Write();
	
	p_x_vs_y_L_F_noal->Write();
	p_x_vs_y_L_N_noal->Write();
	p_x_vs_y_R_N_noal->Write();
	p_x_vs_y_R_F_noal->Write();

	h_y_L_diffFN_vs_y_L_N->Write();
	h_y_R_diffFN_vs_y_R_N->Write();

	h_y_L_ratioFN_vs_y_L_N->Write();
	h_y_R_ratioFN_vs_y_R_N->Write();
	
	gDirectory = outF->mkdir("selected - angles");
	th_x_diffLR->Write();
	th_y_diffLR->Write();

	th_x_diffLF->Write();
	th_x_diffRF->Write();
	
	h_th_x_diffLR_vs_th_x->Write();
	h_th_y_diffLR_vs_th_y->Write();

	p_th_x_diffLR_vs_th_x->Write();
	p_th_y_diffLR_vs_th_y->Write();
	p_th_y_L_diffNF_vs_th_y_L->Write();
	p_th_y_R_diffNF_vs_th_y_R->Write();

	th_x_sigmaLR_vs_th_x->Write();
	th_y_sigmaLR_vs_th_y->Write();
	
	th_x_diffLR_safe->Write();
	th_y_diffLR_safe->Write();

	p_th_x_L_vs_th_y_L->Write();
	p_th_x_R_vs_th_y_R->Write();

	h_th_y_L_vs_th_x_L->Write();
	h_th_y_R_vs_th_x_R->Write();
	h_th_y_vs_th_x->Write();

	g_th_y_L_vs_th_x_L->Write();
	g_th_y_R_vs_th_x_R->Write();
	g_th_y_vs_th_x->Write();
	
	h_th_y_L_vs_th_y_R->Write();
	g_th_y_L_vs_th_y_R->Write();
	
	h_th_x_raw->Write();
	h_th_y_raw->Write();
	h_th_y_raw_flipped->Write();
	
	h_th_x_L->Write();
	h_th_x_R->Write();
	
	h_th_y_L->Write();
	h_th_y_R->Write();

	{
		double x[] = {0, 1, 2, 3};
		double y[] = {anal.th_y_lcut_L, anal.th_y_hcut_L, anal.th_y_lcut_R, anal.th_y_hcut_R};
		TGraph *g = new TGraph(4, x, y);
		g->SetName("g_th_y_cuts");
		g->Write();
	}

	gDirectory = outF->mkdir("selected - angles, alternative");

	h_ta_th_x->Write();
	p_ta_th_x_vs_th_x->Write();

	h_ta_th_y_R->Write();
	h_ta_th_y_L->Write();
	h_ta_th_y->Write();
	h_ta_y_diffLR->Write();
	p_ta_th_y_vs_th_y->Write();
	
	gDirectory = outF->mkdir("selected - vertex");

	h_vtx_x_L->Write();
	h_vtx_x_R->Write();
	h_vtx_x->Write();
	
	p_x_L_F_vs_th_x->Write();
	p_x_L_N_vs_th_x->Write();
	p_x_R_F_vs_th_x->Write();
	p_x_R_N_vs_th_x->Write();
	
	p_x_L_F_vs_vtx_x->Write();
	p_x_L_N_vs_vtx_x->Write();
	p_x_R_F_vs_vtx_x->Write();
	p_x_R_N_vs_vtx_x->Write();

	p_vtx_x_L_vs_th_x_L->Write();
	p_vtx_x_L_vs_th_x->Write();
	p_vtx_x_R_vs_th_x_R->Write();
	p_vtx_x_R_vs_th_x->Write();
	
	p_vtx_x_L_vs_th_y_L->Write();
	p_vtx_x_R_vs_th_y_R->Write();

	h_vtx_x_L_vs_th_x_L->Write();
	h_vtx_x_R_vs_th_x_R->Write();

	h_vtx_x_diffLR->Write();
	h_vtx_x_L_vs_vtx_x_R->Write();
	p_vtx_x_diffLR_vs_th_y->Write();
	p_vtx_x_diffLR_vs_th_x->Write();
	h_vtx_x_diffLR_vs_th_y->Write();
	h_vtx_x_diffLR_vs_th_x->Write();

	h_vtx_y_L_vs_vtx_y_R->Write();

	gDirectory = outF->mkdir("optics");
	h_x_L_F_vs_th_x_L->Write();
	h_x_R_F_vs_th_x_R->Write();

	p_x_L_N_vs_th_x_L->Write();
	p_x_L_F_vs_th_x_L->Write();
	p_x_R_N_vs_th_x_R->Write();
	p_x_R_F_vs_th_x_R->Write();

	h_ratio_ta_x_to_th_x->Write();
	h_rel_diff_ta_x_th_x->Write();
	
	gDirectory = outF->mkdir("binning");
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		TGraph *g = new TGraph();
		g->SetName(("g_binning_"+binnings[bi]).c_str());
		g->SetTitle(";bin center;bin width");

		TH1D *h = bh_t_Nev_before[bi];
		for (int bin = 1; bin <= h->GetNbinsX(); bin++)
			g->SetPoint(g->GetN(), h->GetBinCenter(bin), h->GetBinWidth(bin));

		g->Write();
	}

	TDirectory *timeDepDir = outF->mkdir("time dependences");
	gDirectory = timeDepDir;

	p_diffLR_th_x_vs_time->Write();
	ProfileToRMSGraph(p_diffLR_th_x_vs_time, gRMS_diffLR_th_x_vs_time);
	gRMS_diffLR_th_x_vs_time->Write();

	p_diffLR_th_y_vs_time->Write();
	ProfileToRMSGraph(p_diffLR_th_y_vs_time, gRMS_diffLR_th_y_vs_time);
	gRMS_diffLR_th_y_vs_time->Write();

	p_diffNF_th_y_L_vs_time->Write();
	ProfileToRMSGraph(p_diffNF_th_y_L_vs_time, gRMS_diffNF_th_y_L_vs_time);
	gRMS_diffNF_th_y_L_vs_time->Write();

	p_diffNF_th_y_R_vs_time->Write();
	ProfileToRMSGraph(p_diffNF_th_y_R_vs_time, gRMS_diffNF_th_y_R_vs_time);
	gRMS_diffNF_th_y_R_vs_time->Write();

	p_vtx_x_vs_time->Write();
	ProfileToRMSGraph(p_vtx_x_vs_time, gRMS_vtx_x_vs_time);
	gRMS_vtx_x_vs_time->Write();

	TGraphErrors *g_beam_div_x_vs_time = new TGraphErrors;
	TGraphErrors *g_sensor_res_x_vs_time = new TGraphErrors;

	DecoupleResolution(gRMS_vtx_x_vs_time, gRMS_diffLR_th_x_vs_time, g_beam_div_x_vs_time, g_sensor_res_x_vs_time);

	g_beam_div_x_vs_time->Write();
	g_sensor_res_x_vs_time->Write();

	p_th_x_R_vs_time->Write();
	p_th_y_R_vs_time->Write();
	p_th_x_L_vs_time->Write();
	p_th_y_L_vs_time->Write();

	g_ext_diffLR_th_x_vs_time->Write();
	g_ext_diffLR_th_y_vs_time->Write();

	p_input_beam_div_x_vs_time->Write();
	p_input_beam_div_y_vs_time->Write();

	g_L_L_F_vs_time->Write();
	g_L_R_F_vs_time->Write();

	for (map<unsigned int, TProfile *>::iterator it = p_bunch_diffLR_th_x_vs_time.begin(); it != p_bunch_diffLR_th_x_vs_time.end(); ++it)
	{
		char buf[20];
		sprintf(buf, "bunch %u", it->first);
		gDirectory = timeDepDir->mkdir(buf);
		
		p_bunch_diffLR_th_x_vs_time[it->first]->Write();
		TGraphErrors *gRMS_diffLR_th_x_vs_time = new TGraphErrors;
		gRMS_diffLR_th_x_vs_time->SetName("gRMS_diffLR_th_x_vs_time");
		gRMS_diffLR_th_x_vs_time->SetTitle(";timestamp;RMS of #Delta^{R-L}#theta_{x}");
		ProfileToRMSGraph(p_bunch_diffLR_th_x_vs_time[it->first], gRMS_diffLR_th_x_vs_time);
		gRMS_diffLR_th_x_vs_time->Write();

		p_bunch_diffLR_th_y_vs_time[it->first]->Write();
		TGraphErrors *gRMS_diffLR_th_y_vs_time = new TGraphErrors;
		gRMS_diffLR_th_y_vs_time->SetName("gRMS_diffLR_th_y_vs_time");
		gRMS_diffLR_th_y_vs_time->SetTitle(";timestamp;RMS of #Delta^{R-L}#theta_{y}");
		ProfileToRMSGraph(p_bunch_diffLR_th_y_vs_time[it->first], gRMS_diffLR_th_y_vs_time);
		gRMS_diffLR_th_y_vs_time->Write();

		p_bunch_vtx_x_vs_time[it->first]->Write();
		TGraphErrors *gRMS_vtx_x_vs_time = new TGraphErrors;
		gRMS_vtx_x_vs_time->SetName("gRMS_vtx_x_vs_time");
		gRMS_vtx_x_vs_time->SetTitle(";timestamp;RMS of x^{*}");
		ProfileToRMSGraph(p_bunch_vtx_x_vs_time[it->first], gRMS_vtx_x_vs_time);
		gRMS_vtx_x_vs_time->Write();
	
		TGraphErrors *g_beam_div_x_vs_time = new TGraphErrors;
		TGraphErrors *g_sensor_res_x_vs_time = new TGraphErrors;

		DecoupleResolution(gRMS_vtx_x_vs_time, gRMS_diffLR_th_x_vs_time, g_beam_div_x_vs_time, g_sensor_res_x_vs_time);

		g_beam_div_x_vs_time->Write();
		g_sensor_res_x_vs_time->Write();
	}

	TDirectory *accDir = outF->mkdir("acceptance correction");
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		gDirectory = accDir->mkdir(binnings[bi].c_str());
		bh_t_Nev_before[bi]->Write();
		bh_t_Nev_after_no_corr[bi]->Write();
		bh_t_before[bi]->Write();
		bh_t_after_no_corr[bi]->Write();
		bh_t_after[bi]->Write();
		bp_t_phi_corr[bi]->Write();
		bp_t_full_corr[bi]->Write();
		bp_t_full_acc[bi]->Write();
	
		c = new TCanvas("t cmp");
		c->SetLogy(1);
		bh_t_after[bi]->Draw("");
		bh_t_before[bi]->Draw("same");
		c->Write();
	}
		
	gDirectory = accDir;
	
	p_t_ub_div_corr->Write();
	
	h_th_y_vs_th_x_before->Write();
	h_th_y_vs_th_x_after->Write();
	h_th_vs_phi_after->Write();
	
	if (detailsLevel >= 1)
		g_weight_vs_th_y->Write();
	
	TDirectory *normDir = outF->mkdir("normalization");
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		gDirectory = normDir->mkdir(binnings[bi].c_str());
		bh_t_normalized[bi]->Write();
		bh_t_normalized_rel_diff[bi]->Write();
	}

	gDirectory = normDir;
	p_norm_corr->Write();
	p_3outof4_corr->Write();

	h_th_y_vs_th_x_normalized->Write();

	g_norm_corr_vs_div_corr->Write();
	
	TDirectory *normUnfDir = outF->mkdir("normalization+unfolding");
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		gDirectory = normUnfDir->mkdir(binnings[bi].c_str());
		
		bh_t_normalized_unsmeared[bi]->Write();
		bh_t_normalized_unsmeared_rel_diff[bi]->Write();
	}

	gDirectory = outF->mkdir("background");
	for (map<unsigned int, TH1D *>::iterator it = hb_cq.begin(); it != hb_cq.end(); ++it)
		it->second->Write();
	
	if (detailsLevel >= 2)
	{
		hb_th_y_L->SetLineColor(2);
		hb_th_y_R->SetLineColor(4);
		
		hb_th_y_L->Write();
		hb_th_y_R->Write();
	
		c = new TCanvas("th_y LR cmp");
		c->SetLogy(1);
		hb_th_y_L->Draw("");
		hb_th_y_R->Draw("same");
		c->Write();
		
		hb_th_x_L->SetLineColor(2);
		hb_th_x_R->SetLineColor(4);
		
		hb_th_x_L->Write();
		hb_th_x_R->Write();
		
		c = new TCanvas("th_x LR cmp");
		c->SetLogy(1);
		hb_th_x_L->Draw("");
		hb_th_x_R->Draw("same");
		c->Write();
	
		hb_th_y_L_vs_th_x_L->Write();
		hb_th_y_R_vs_th_x_R->Write();
		
		hb_th_y_diffLR->Write();
		hb_th_x_diffLR->Write();
		
		hb_th_y_diffLR_vs_th_y->Write();
		hb_th_x_diffLR_vs_th_y->Write();
		
		hb_th_y_L_vs_th_y_R->Write();
		hb_th_x_L_vs_th_x_R->Write();
	
		hb_th_y_6cut->Write();
		hb_th_y_6cut_cut7fail->Write();
		TH1D *hb_th_y_6cut_cut7failRatio = new TH1D(*hb_th_y_6cut_cut7fail);
		hb_th_y_6cut_cut7failRatio->SetName("hb_th_y_6cut_cut7failRatio");
		hb_th_y_6cut_cut7failRatio->Divide(hb_th_y_6cut);
		hb_th_y_6cut_cut7failRatio->Write();
	}

	// print counters
	for (map<unsigned int, unsigned long>::iterator it = n_ev_cut.begin(); it != n_ev_cut.end(); ++it)
		printf("\tcut %u: %lu\n", it->first, it->second);

	return 0;
}
