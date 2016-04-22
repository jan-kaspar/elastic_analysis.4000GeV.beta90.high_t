#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"
#include "../stat.h"

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
	for (int bi = 1; bi <= bins; bi++)
	{
		int i = bi - 1;
		h_mean_stddev->SetBinContent(bi, st.GetMean(i));
		h_mean_stddev->SetBinError(bi, st.GetStdDev(i));
	}

	h_mean_stddev->Write();

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
	//binnings.push_back("ob-3-10-0.10");

	// get input
	TFile *f_in_cf = new TFile(("unfolding_cf_" + dgn_str + ".root").c_str());
	TFile *f_in_gr = new TFile(("unfolding_gr_" + dgn_str + ".root").c_str());

	// prepare output
	TFile *f_out = new TFile(("unfolding_summarize_" + dgn_str + ".root").c_str(), "recreate");

	// process data
	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		TDirectory *binningDir = f_out->mkdir(binnings[bi].c_str());

		// model uncertainty
		gDirectory = binningDir->mkdir("model");

		vector<TH1D *> histList;
		AddHistToList(f_in_cf, binnings[bi]+"/exp3+exp4/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binnings[bi]+"/exp5+erf*exp2/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binnings[bi]+"/p1*exp3+p1*exp1/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binnings[bi]+"/p1*exp3+p2*exp2/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binnings[bi]+"/exp3-intf-exp1/+0,+0/corr_final", histList);
		AddHistToList(f_in_cf, binnings[bi]+"/(exp3-intf-exp1)*expG/+0,+0/corr_final", histList);
		
		// TODO: resolve binning and diagonal dependence
		AddHistToList(f_in_gr, binnings[bi]+"/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E+00/h_corr", histList);
		AddHistToList(f_in_gr, binnings[bi]+"/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr", histList);
		AddHistToList(f_in_gr, binnings[bi]+"/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E-01/h_corr", histList);
		AddHistToList(f_in_gr, binnings[bi]+"/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E-01/h_corr", histList);

		AnalyzeModelDependence(histList);

		// sigma (x and y) uncertainty
		// TODO


		// TODO: careful about histogram "jumps" at low and high-|t| edges !!
	}

	delete f_in_cf;
	delete f_in_gr;
	delete f_out;
	return 0;
}
