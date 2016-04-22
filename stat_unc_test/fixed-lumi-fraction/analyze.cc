#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <vector>
#include <string>

#include "analyze_common.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

vector<string> datasets;

void AnalyzeOneHistogramOneDiagonal(const string &obj_path, const string &diagonal, TDirectory *d_out)
{
	vector<string> files;
	
	// build list of input files
	for (unsigned int dsi = 0; dsi < datasets.size(); dsi++)
	{
		string f = datasets[dsi] + "/distributions_" + diagonal.c_str() + ".root";
		files.push_back(f);
	}

	AnalyzeOne(files, obj_path, d_out);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: analyze <base directory> <number of subsamples>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// parse command line
	if (argc != 3)
	{
		PrintUsage();
		return 1;
	}

	string baseDir = argv[1];
	unsigned int nSamples = atoi(argv[2]);

	// input
	for (unsigned int i = 0; i < nSamples; i++)
	{
		char buf[200];
		sprintf(buf, "%s/%u-%u", baseDir.c_str(), nSamples, i);
		datasets.push_back(buf);
	}
	
	// diagonals
	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	// binnings
	vector<string> binnings = {
		"ub",
		"ob-1-30-0.10",
		"ob-2-20-0.10",
		"ob-3-10-0.10",
		"ob-3-10-0.20",
		"ob-3-10-0.30"
	};

	// run processing
	char buf[100];
	sprintf(buf, "analyze_%u.root", nSamples);
	TFile *f_out = new TFile(buf, "recreate");

	for (unsigned int bi = 0; bi < binnings.size(); bi++)
	{
		printf("* %s\n", binnings[bi].c_str());

		TDirectory *binDir = f_out->mkdir(binnings[bi].c_str());

		for (unsigned int dgni = 0; dgni < diagonals.size(); dgni++)
		{
			printf("\t%s\n", diagonals[dgni].c_str());

			TDirectory *dgnDir = binDir->mkdir(diagonals[dgni].c_str());

			AnalyzeOneHistogramOneDiagonal("acceptance correction/"+binnings[bi]+"/h_t_Nev_before", diagonals[dgni], dgnDir->mkdir("h_t_Nev_before"));

			//AnalyzeOneHistogram("acceptance correction/"+binnings[bi]+"/h_t_Nev_before", "h_t_Nev_before", dgn);
			//AnalyzeOneHistogram("acceptance correction/ub/h_t_Nev_after_no_corr", "ub-h_t_Nev_after_no_corr", f_out);
			//AnalyzeOneHistogram("acceptance correction/ub/h_t_after", "ub-h_t_after", f_out);
		}
	}

	delete f_out;
	return 0;
}
