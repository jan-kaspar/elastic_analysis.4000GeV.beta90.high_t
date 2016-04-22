#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"

#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void SampleModel(const string &fn_in, const string &tag, const string &fn_out)
{
	// get input
	TFile *f_in = TFile::Open(fn_in.c_str());
	TF1 *ff = (TF1 *) f_in->Get(("ob-1-30-0.10/" + tag + "/+0,+0/iteration 2/ff").c_str());

	// make output
	TFile *f_out = TFile::Open(fn_out.c_str(), "recreate");

	TGraph *g_dsdt = new TGraph();
	g_dsdt->SetName("g_dsdt");

	for (double t = 0; t <= 2.1; t += 0.01)
	{
		double v = ff->Eval(t);

		g_dsdt->SetPoint(g_dsdt->GetN(), t, v);
	}

	g_dsdt->Write();

	delete f_out;
	delete f_in;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	string fn;
   
	fn = "../DS4/unfolding_cf_45b_56t.root";
	SampleModel(fn, "exp3-intf-exp1", "exp3-intf-exp1.root");

	fn = "../DS4/unfolding_cf_45t_56b.root";
	SampleModel(fn, "p1*exp3+p2*exp2", "p1*exp3+p2*exp2.root");

	return 0;
}
