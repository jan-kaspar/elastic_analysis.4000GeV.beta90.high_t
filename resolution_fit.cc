#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TF1.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void SetRange(TGraph *g, TF1 *f)
{
	double x_min = 1E100, x_max = -1E100;
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, y;
		g->GetPoint(i, x, y);

		if (y > 0.)
		{
			x_min = min(x_min, x);
			x_max = max(x_max, x);
		}
	}

	printf(">> SetRange: %E, %E\n", x_min, x_max);

	f->SetRange(x_min, x_max);
}

//----------------------------------------------------------------------------------------------------

void ProcessOne(TDirectory *dir)
{
	TF1 *ff = new TF1("ff", "[0] + [1]*x");

	TGraphErrors *g_x = (TGraphErrors *) dir->Get("gRMS_diffLR_th_x_vs_time");
	g_x->Fit(ff);
	g_x->Write();

	SetRange(g_x, ff);
	ff->Write("fit_x");

	TGraphErrors *g_y = (TGraphErrors *) dir->Get("gRMS_diffLR_th_y_vs_time");
	g_y->Fit(ff);
	g_y->Write();
	SetRange(g_y, ff);
	ff->Write("fit_y");
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	Init(argv[1]);
	if (diagonal == dCombined || diagonal == ad45b_56b || diagonal == ad45t_56t)
		return rcIncompatibleDiagonal;

	// default parameters
	string outputDir = ".";
	string inputDir = ".";

	// get input
	TFile *f_in = new TFile((inputDir + "/distributions_" + argv[1] + ".root").c_str());
	
	// prepare output
	TFile *f_out = new TFile((outputDir + "/resolution_fit_" + argv[1] + ".root").c_str(), "recreate");

	// process distributions for all bunches
	TDirectory *baseDir = (TDirectory *) f_in->Get("time dependences");
	TIter next(baseDir->GetListOfKeys());
	TObject *o;
	while ((o = next()))
	{
		TKey *k = (TKey *) o;
		if (k->IsFolder())
		{
			TDirectory *bunchDir = (TDirectory *) k->ReadObj();

			gDirectory = f_out->mkdir(k->GetName());
			ProcessOne(bunchDir);
		}
	}

	delete f_out;
	delete f_in;
	return 0;
}
