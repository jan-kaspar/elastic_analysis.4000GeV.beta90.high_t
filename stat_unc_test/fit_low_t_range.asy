import root;
import pad_layout;
include "../plots/t_distributions/common_code.asy";

A_ref = 519.5;
B_ref = 19.38;
ref_str = MakeRefStr("");

xSizeDef = 10cm;
ySizeDef = 8cm;

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b--56t", "45t--56b" };

string binnings[] = {
	"ob-1-30-0.10",
	"ob-2-20-0.20",
	"ob-3-10-0.30"
};

//----------------------------------------------------------------------------------------------------

void DrawRelDiff2(transform t = identity(), rObject o, pen p, marker m = mCi+2pt+black, string label="")
{
	if (o.InheritsFrom("TH1"))
	{
		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xl = o.rExec("GetBinLowEdge", i);
			real xw = o.rExec("GetBinWidth", i);
			real xr = xl + xw;
			//real xc = xl + xw/2;
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
	
			real y_ref_l = A_ref * exp(-B_ref * xl);
			real y_ref_r = A_ref * exp(-B_ref * xr);
	
			real y_rel_l = (y - y_ref_l) / y_ref_l;
			real y_rel_r = (y - y_ref_r) / y_ref_r;


			//real y_rel_unc = y_unc / y_ref;
	
			draw(t * ((xl, y_rel_l)--(xr, y_rel_r)), p);


			//draw(t * ((xc, y_rel-y_rel_unc)--(xc, y_rel+y_rel_unc)), p+0.1pt);	
			//draw(t * (xc, y_rel), mCi+0.001pt+p);
		}
	
		if (label != "")
			AddToLegend(label, mPl+4pt+p);
	}
}


//----------------------------------------------------------------------------------------------------

for (int dgni : diagonals.keys)
{
	NewRow();

	NewPad(false);
	label(replace(diagonals[dgni], "_", "--"));

	for (int bi : binnings.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "${\d\si/\d t - \hbox{ref}\over\hbox{ref}}\ ,\quad\hbox{ref} = "+ref_str+"$");

		rObject obj_h = rGetObj("fit_low_t_range.root", diagonals[dgni]+"/"+binnings[bi]+"/h_dsdt");
		rObject obj_g = rGetObj("fit_low_t_range.root", diagonals[dgni]+"/"+binnings[bi]+"/g_dsdt");
		rObject obj_f = rGetObj("fit_low_t_range.root", diagonals[dgni]+"/"+binnings[bi]+"/ff");

		DrawRelDiff(obj_g, black);
		DrawRelDiff(obj_f, red, obj_f.sExec("GetTitle"));
	
		limits((0, -0.02), (0.2, +0.08), Crop);
		AttachLegend(binnings[bi]);
	}
}
