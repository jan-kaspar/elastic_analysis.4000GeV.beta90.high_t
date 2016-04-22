#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;
	
	// get input
	TFile *f_in = new TFile("eff3outof4_details.root");
	
	// prepare output
	TFile *f_out = new TFile("eff3outof4_details_fit.root", "recreate");

	// combinations
	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	vector<string> units;
	units.push_back("L_F");
	units.push_back("L_N");
	units.push_back("R_N");
	units.push_back("R_F");

	// fit function
	TF1 *ff = new TF1("ff", "[0] + [1]*x");

	for (unsigned int dgni = 0; dgni < diagonals.size(); dgni++)
	{
		printf("\n\n------------------------------ %s ------------------------------\n", diagonals[dgni].c_str());

		Init(diagonals[dgni]);

		printf("\n\n");

		double th_y_min = max(anal.th_y_lcut_L, anal.th_y_lcut_R) + 8E-6;
		double th_y_max = min(anal.th_y_hcut_L, anal.th_y_hcut_R) - 2E-6;

		th_y_min *= th_y_sign;
		th_y_max *= th_y_sign;
		if (th_y_max < th_y_min)
			swap(th_y_max, th_y_min);

		printf("th_y_min = %E\n", th_y_min);
		printf("th_y_max = %E\n", th_y_max);

		TDirectory *dgnDir = f_out->mkdir(diagonals[dgni].c_str());

		for (unsigned int ui = 0; ui < units.size(); ui++)
		{
			printf("\n\n>> %s\n", units[ui].c_str());

			TDirectory *unitDir = dgnDir->mkdir(units[ui].c_str());
			gDirectory = unitDir;

			TH1D *h_eff = (TH1D *) f_in->Get((diagonals[dgni]+"/"+units[ui]+"/track_compatible/th_y : rel").c_str());

			h_eff->Fit(ff, "", "", th_y_min*1E6, th_y_max*1E6);
			h_eff->Write();

			ff->Write("fit");
		}
	}

	delete f_out;
	return 0;
}
