#include "input_files.h"

#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"
#include "stat.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void AnalyzeModelDependence(const vector<TH1D *> &hists)
{
	printf("\n>> AnalyzeModelDependence\n");

	if (hists.size() < 1)
		return;

	TH1D *h_ref = hists[0];

	int bins = h_ref->GetNbinsX();

	Stat st(bins);

	for (unsigned int hi = 0; hi < hists.size(); hi++)
	{
		vector<double> corr(bins);
		for (int bi = 1; bi <= bins; bi++)
			corr[bi - 1] = hists[hi]->GetBinContent(bi);

		st.Fill(corr);
	}

	// histogram showing mean correction and its std. deviation
	TH1D *h_mean_stddev = new TH1D(*h_ref);
	h_mean_stddev->SetName("h_mean_stddev");
	TH1D *h_unc = new TH1D(*h_ref);
	h_unc->SetName("h_unc");
	for (int bi = 1; bi <= bins; bi++)
	{
		int i = bi - 1;
		h_mean_stddev->SetBinContent(bi, st.GetMean(i));
		h_mean_stddev->SetBinError(bi, st.GetStdDev(i));

		h_unc->SetBinContent(bi, st.GetStdDev(i));
		h_unc->SetBinError(bi, 0.);
	}

	h_mean_stddev->Write();
	h_unc->Write();

	// correlation matrix
	const double* binEdges = h_ref->GetXaxis()->GetXbins()->GetArray();
	TH2D *h2_corr_mat = new TH2D("h_corr_mat", "", bins, binEdges, bins, binEdges);
	
	for (int bi = 1; bi <= bins; bi++)
	{
		for (int bj = 1; bj <= bins; bj++)
		{
			h2_corr_mat->SetBinContent(bi, bj, st.GetCorrelation(bi-1, bj-1));
		}
	}

	h2_corr_mat->Write();
}

//----------------------------------------------------------------------------------------------------

void ExtractDifference(const vector<TH1D *> &hists)
{
	printf("\n>> ExtractDifference\n");

	if (hists.size() != 2)
		return;

	TH1D *h1 = hists[0];
	TH1D *h2 = hists[1];

	TH1D *h_unc = new TH1D(*h1);

	for (int bi = 1; bi <= h1->GetNbinsX(); ++bi)
	{
		const double v1 = h1->GetBinContent(bi);
		const double v2 = h2->GetBinContent(bi);
		h_unc->SetBinContent(bi, v2 - v1);
	}

	h_unc->Write("h_unc");
}

//----------------------------------------------------------------------------------------------------

void AddHistToList(TFile *f, const string &path, vector<TH1D *> &lst)
{
	TH1D *h = (TH1D *) f->Get(path.c_str());
	if (!h)
	{
		printf("ERROR: can't load object '%s' from file '%s'.\n", path.c_str(), f->GetName());
		abort();
	}

	lst.push_back(h);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	string dgn_str = argv[1];

	Init(dgn_str);
	if (diagonal == dCombined || diagonal == ad45b_56b || diagonal == ad45t_56t)
		return rcIncompatibleDiagonal;

	// binnings
	vector<string> binnings;
	binnings.push_back("ob-1-30-0.10");
	binnings.push_back("bt1");

	// get input
	TFile *f_in_cf = new TFile(("unfolding_cf_" + dgn_str + ".root").c_str());
	TFile *f_in_gr = new TFile(("unfolding_gr_" + dgn_str + ".root").c_str());

	// prepare output
	TFile *f_out = new TFile(("unfolding_summarize_" + dgn_str + ".root").c_str(), "recreate");

	// process data
	for (const auto &binning : binnings)
	{
		TDirectory *binningDir = f_out->mkdir(binning.c_str());

		// model uncertainty
		gDirectory = binningDir->mkdir("model");

		vector<TH1D *> histList;
		AddHistToList(f_in_cf, binning+"/exp3+exp4/+0,+0/corr_final", histList);
		//AddHistToList(f_in_cf, binning+"/exp5+erf*exp2/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binning+"/p1*exp3+p1*exp1/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binning+"/p1*exp3+p2*exp2/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binning+"/exp3-intf-exp1/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binning+"/(exp3-intf-exp1)*expG/+0,+0/corr_final", histList);
		
		AddHistToList(f_in_gr, binning+"/smearing_matrix_mc_" + dgn_str + ".root,p1*exp3+p2*exp2," + binning + "/alpha=1.00E+00/h_corr", histList);
		AddHistToList(f_in_gr, binning+"/smearing_matrix_mc_" + dgn_str + ".root,exp3-intf-exp1," + binning + "/alpha=1.00E+00/h_corr", histList);
		AddHistToList(f_in_gr, binning+"/smearing_matrix_mc_" + dgn_str + ".root,p1*exp3+p2*exp2," + binning + "/alpha=1.00E-01/h_corr", histList);
		AddHistToList(f_in_gr, binning+"/smearing_matrix_mc_" + dgn_str + ".root,exp3-intf-exp1," + binning + "/alpha=1.00E-01/h_corr", histList);

		AnalyzeModelDependence(histList);

		// sigma (x and y) uncertainty
		if (binning == "bt1")
		{
			const string model = "exp3-intf-exp1";

			gDirectory = binningDir->mkdir("sigma x");
			histList.clear();
			AddHistToList(f_in_cf, binning+"/" + model + "/+0,+0/corr_final", histList);
			AddHistToList(f_in_cf, binning+"/" + model + "/+1,+0/corr_final", histList);
			ExtractDifference(histList);

			gDirectory = binningDir->mkdir("sigma y");
			histList.clear();
			AddHistToList(f_in_cf, binning+"/" + model + "/+0,+0/corr_final", histList);
			AddHistToList(f_in_cf, binning+"/" + model + "/+0,+1/corr_final", histList);
			ExtractDifference(histList);
		}
	}

	// clean up
	delete f_in_cf;
	delete f_in_gr;
	delete f_out;

	return 0;
}
