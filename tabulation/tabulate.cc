#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"

#include <vector>
#include <cmath>

#include "../stat.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

vector<TF1 *> fit_graphs;

void AddFit(const string &fn, const string &on)
{
	TFile *f = new TFile(fn.c_str());
	fit_graphs.push_back( (TF1 *) f->Get(on.c_str()) );
}

//----------------------------------------------------------------------------------------------------

double GetBinRepresentativePoint(TF1 *g_fit, double l, double r)
{
	double w = r - l;

	// calculate integral mean
	const unsigned int n_div = 1000;
	double V = 0.;
	for (unsigned int i = 0; i < n_div; i++)
	{
		double x = l + (0.5 + i) * w / n_div;
		V += g_fit->Eval(x) / n_div;
	}

	bool decreasing = (g_fit->Eval(r) < g_fit->Eval(l));


	// find representative point
	double xr;
	while (r - l > w/100000)
	{
		xr = (r+l)/2.;

		const double f_at_xr = g_fit->Eval(xr);

		//printf("  %.4f, %.4f | %.8f, V = %.8f, cond = %u\n", l, r, f_at_xr, V, (f_at_xr > V));

		bool cond = (f_at_xr < V);
		if (!decreasing)
			cond = !cond;

		if (cond)
		{
			r -= (r-l)/2.;
		} else {
			l += (r-l)/2.;
		}
	}
	xr = (r+l)/2.;

	return xr;
}

//----------------------------------------------------------------------------------------------------

string Format(double v, bool sign, int N_bef, int N_aft, int N_sig)
{
	string output = "";

	if (sign)
	{
		if (v < 0)
			output += "-";
		else
			output += "+";
	}

	signed int order = floor(log10(fabs(v)));

	signed int pos_left = order;
	signed int pos_right = order - N_sig + 1;
	signed int pos_start = N_bef - 1;
	signed int pos_end = -N_aft;

	signed int pos_round = max(pos_right, pos_end);

	int v_rounded = floor(fabs(v) * pow(10, -pos_round) + 0.5);
	char buf[20];
	sprintf(buf, "%i", v_rounded);

	/*
	printf("\n");
	printf("v = %E\n", v);
	printf("buf = %s\n", buf);
	printf("order = %i\n", order);
	printf("pos_left = %i\n", pos_left);
	printf("pos_right = %i\n", pos_right);
	printf("pos_start = %i\n", pos_start);
	printf("pos_end = %i\n", pos_end);
	*/

	string spaceStr = "@";

	for (signed int pos = pos_start; pos >= pos_end; pos--)
	{
		//printf("pos = %i\n", pos);

		// add decimal point
		if (pos == -1)
			output += ".";

		if (pos > pos_left)
		{
			if (pos > 0)
				output += spaceStr;
			else
				output += "0";
		}

		if (pos < pos_right)
		{
			if (pos < 0)
				output += spaceStr;
			else
				output += "0";
		}

		if (pos >= pos_right && pos <= pos_left)
			output += buf[pos_left - pos];

		//printf("\t<%s>\n", output.c_str());
	}

	return output;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	string dataset = "DS4-sc";
	string diagonal = "combined";
	string binning = "bt1";

	// get data
	TFile *f_data_in = new TFile("../DS-merged/merged.root");
	TH1D *h_data = (TH1D *) f_data_in->Get((binning + "/" + dataset + "/" + diagonal + "/h_dsdt").c_str());
	printf("h_data = %p\n", h_data);

	// get uncertainties
	string f_unc_name = "../DS4/systematics_matrix.root";
	TFile *f_unc_in = new TFile(f_unc_name.c_str());

	TH1D *h_unc_all = (TH1D *) f_unc_in->Get(("matrices/all/"+diagonal+"/"+binning+"/h_stddev").c_str());
	//TH1D *h_unc_all_anal = (TH1D *) f_unc_in->Get(("matrices/all-anal/"+binning+"/h_stddev").c_str());
	
	printf("h_unc_all = %p\n", h_unc_all);
	//printf("h_unc_all_anal = %p\n", h_unc_all_anal);

	// fits
	AddFit("make_fits.root", "exp3-intf-exp1");
	AddFit("make_fits.root", "p1*exp3+p2*exp2");

	if (fit_graphs.empty())
	{
		printf("ERROR: fit_graphs empty.\n");
		return 1;
	}

	TF1 *ref_fit = fit_graphs[0];

	// output
	TFile *f_out = new TFile("tabulate.root", "recreate");

	TH1D *h_dsdt_val_exp_unc_stat = new TH1D(*h_data);
	TH1D *h_dsdt_val_exp_unc_syst = new TH1D(*h_data);
	TH1D *h_dsdt_val_ref_unc_stat = new TH1D(*h_data);
	TH1D *h_dsdt_val_ref_unc_syst = new TH1D(*h_data);

	h_dsdt_val_exp_unc_stat->Reset();
	h_dsdt_val_exp_unc_syst->Reset();
	h_dsdt_val_ref_unc_stat->Reset();
	h_dsdt_val_ref_unc_syst->Reset();

	TGraphErrors *g_dsdt_val_exp_unc_stat = new TGraphErrors();

	// process all bins
	for (int bi = 1; bi <= h_data->GetNbinsX(); bi++)
	{
		const double b_left = h_data->GetBinLowEdge(bi);
		const double b_right = b_left + h_data->GetBinWidth(bi);
		const double b_cen = (b_left + b_right) / 2.;

		const double v = h_data->GetBinContent(bi);
		const double v_unc_stat = h_data->GetBinError(bi);

		if (v == 0.)
			continue;

		// set range
		if (b_left < 0.19 || b_left >= 1.86)
			continue;

		// determine representative point
		/*
		Stat st_x_rep(1);
		for (unsigned int fi = 0; fi < fit_graphs.size(); fi++)
		{
			const double b_rep = GetBinRepresentativePoint(fit_graphs[fi], b_left, b_right);
			st_x_rep.Fill(b_rep);
		}

		double b_rep_mean = st_x_rep.GetMean(0);
		double b_rep_mean_unc = st_x_rep.GetMeanUnc(0);

		if (b_rep_mean <= b_left || b_rep_mean >= b_right)
		{
			printf("ERROR: problem in calculation of representative point, bin %u.\n", bi);
		}

		if (fabs((b_left + b_right)/2. - b_rep_mean) > 0.10 * (b_right - b_left))
			printf("WARNING: strange b_rep_mean = %.4f, while b_cen = %.4f\n", b_rep_mean, (b_left + b_right)/2.);
		*/

		// check matching of binning in data and systematics
		const double b_cen_data = h_data->GetBinCenter(bi);
		const double b_cen_syst = h_unc_all->GetBinCenter(bi);

		if (fabs(b_cen_data - b_cen_syst) > 1E-4 * b_cen_data)
			printf("ERROR: binning in data and systematics don't match\n");

		// reference bin content
		//const double v_ref = ref_fit->Eval(b_rep_mean);
		const double v_ref = ref_fit->Eval(b_cen);

		// systematics
		const double v_unc_syst_all = h_unc_all->GetBinContent(bi) * v_ref;
		//const double v_unc_syst_all_anal = h_unc_all_anal->GetBinContent(bi) * v_ref;

		// print formatted line
		if (false)
			printf("%.4f, %.4f | %.5f, %.5f, %.5f\n",
				b_left, b_right,
				v, v_unc_stat, v_unc_syst_all
			);

		// TeX simple
		if (false)
			printf("$%.5f$ & $%.5f$ & $%.5f$ & $%.5f$ & $%.5f$\\cr\n",
				b_left, b_right,
				v, v_unc_stat, v_unc_syst_all
			);

		// TeX nice 
		if (true)
			printf("$%.5f$ & $%.5f$ & $%s$ & $%s$ & $%s$\\cr\n",
				b_left, b_right,
				Format(v, false, 2, 7, 5).c_str(),
				Format(v_unc_stat, false, 2, 6, 3).c_str(),
				Format(v_unc_syst_all, false, 2, 6, 3).c_str()
			);

		// fill output
		h_dsdt_val_exp_unc_stat->SetBinContent(bi, v);
		h_dsdt_val_exp_unc_stat->SetBinError(bi, v_unc_stat);

		h_dsdt_val_exp_unc_syst->SetBinContent(bi, v);
		h_dsdt_val_exp_unc_syst->SetBinError(bi, v_unc_syst_all);

		h_dsdt_val_ref_unc_stat->SetBinContent(bi, v_ref);
		h_dsdt_val_ref_unc_stat->SetBinError(bi, v_unc_stat);

		h_dsdt_val_ref_unc_syst->SetBinContent(bi, v_ref);
		h_dsdt_val_ref_unc_syst->SetBinError(bi, v_unc_syst_all);

		/*
		int idx = g_dsdt_val_exp_unc_stat->GetN();
		g_dsdt_val_exp_unc_stat->SetPoint(idx, b_rep_mean, v);
		g_dsdt_val_exp_unc_stat->SetPointError(idx, b_rep_mean, v_unc_stat);
		*/
	}

	// save output
	h_dsdt_val_exp_unc_stat->Write("h_dsdt_val_exp_unc_stat");
	h_dsdt_val_exp_unc_syst->Write("h_dsdt_val_exp_unc_syst");
	h_dsdt_val_ref_unc_stat->Write("h_dsdt_val_ref_unc_stat");
	h_dsdt_val_ref_unc_syst->Write("h_dsdt_val_ref_unc_syst");

	g_dsdt_val_exp_unc_stat->Write("g_dsdt_val_exp_unc_stat");

	delete f_out;

	return 0;
}
