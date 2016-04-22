#include "TGraph.h"
#include "TSpline.h"

#include <string>

using namespace std;


//----------------------------------------------------------------------------------------------------

TGraph *input_dist_t_data_fit1 = NULL;
TSpline *spline_dist_t_data_fit1 = NULL;

//----------------------------------------------------------------------------------------------------

TGraph* CropTDistribution(TGraph *input, double cut_off = 2.1)
{
	printf(">> CropTDistribution > cut_off = %.2f\n", cut_off);

	TGraph *g = new TGraph();
	for (int i = 0; i < input->GetN(); i++)
	{
		double t, v;
		input->GetPoint(i, t, v);

		if (t > cut_off)
			break;

		int idx = g->GetN();
		g->SetPoint(idx, t, v);
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

TSpline* BuildSpline(TGraph *g)
{
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	s->SetName(g->GetName());
	return s;
}

//----------------------------------------------------------------------------------------------------

TObject* GetObject(const string &file, const string &path)
{
	TFile *f = TFile::Open(file.c_str());
	if (!f)
	{
		printf("ERROR: can't load file `%s'.\n", file.c_str());
		return NULL;
	}

	TObject *o = f->Get(path.c_str());
	if (!o)
	{
		printf("ERROR: can't load object `%s' from file `%s'.\n", path.c_str(), file.c_str());
		return NULL;
	}

	return o;
}

//----------------------------------------------------------------------------------------------------

int LoadTDistributions()
{
	input_dist_t_data_fit1 = (TGraph *) GetObject("../unfolding/models/exp3-intf-exp1.root", "g_dsdt");

	if (!input_dist_t_data_fit1)
		return 1;

	input_dist_t_data_fit1 = CropTDistribution(input_dist_t_data_fit1);
	input_dist_t_data_fit1->SetName("input_dist_t_data_fit1");
	spline_dist_t_data_fit1 = BuildSpline(input_dist_t_data_fit1);

	return 0;
}

//----------------------------------------------------------------------------------------------------

void WriteTDistributions()
{
	input_dist_t_data_fit1->Write();
	spline_dist_t_data_fit1->Write();
}
