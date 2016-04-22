import root;
import pad_layout;

string diagonals[] = {
	"45b_56t",
	"45t_56b"
};

string binnings[] = {
	"ob-1-30-0.10",
//	"ob-2-20-0.20",
//	"ob-3-10-0.30",
};

TGraph_errorBar = None;

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

void DrawOneSet(string diagonal, string binning)
{
	string f = "fit_test.root";
	
	NewPad("$|t|\ung{GeV^2}$", "$\ch^2$");
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	currentpad.yTicks = RightTicks(1., 0.2);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_3,exp1"), "p", black, mCi+2pt+black, "3 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_4,exp2"), "p", red, mCi+2pt+red, "4 bins, exp2");
	limits((0, 0), (2, 8), Crop);
	xaxis(YEquals(1, false), dashed);
	yaxis(XEquals(0.5, false), dotted);
	AttachLegend("ndf=1");

	NewPad("$|t|\ung{GeV^2}$", "$\ch^2$");
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	currentpad.yTicks = RightTicks(1., 0.2);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_4,exp1"), "p", black, mCi+2pt+black, "4 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_5,exp2"), "p", red, mCi+2pt+red, "5 bins, exp2");
	limits((0, 0), (2, 8), Crop);
	xaxis(YEquals(2, false), dashed);
	yaxis(XEquals(0.5, false), dotted);
	AttachLegend("ndf=2");

	NewPad("$|t|\ung{GeV^2}$", "$\ch^2$");
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	currentpad.yTicks = RightTicks(1., 0.2);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_5,exp1"), "p", black, mCi+2pt+black, "5 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_6,exp2"), "p", red, mCi+2pt+red, "6 bins, exp2");
	limits((0, 0), (2, 8), Crop);
	xaxis(YEquals(3, false), dashed);
	yaxis(XEquals(0.5, false), dotted);
	AttachLegend("ndf=3");

	NewRow();

	NewPad("$\ch^2$", "");
	draw(rGetObj(f, "g_chisq_ndf1"), "l", heavygreen);
	draw(rGetObj(f, "h_chisq_ndf1"), "vl", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/h_3,exp1"), "vl", black, "3 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/h_4,exp2"), "vl", red, "4 bins, exp2");
	limits((0, 0), (10, 1.), Crop);
	AttachLegend("ndf=1");

	NewPad("$\ch^2$", "");
	draw(rGetObj(f, "g_chisq_ndf2"), "l", heavygreen);
	draw(rGetObj(f, "h_chisq_ndf2"), "vl", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/h_4,exp1"), "vl", black, "4 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/h_5,exp2"), "vl", red, "5 bins, exp2");
	limits((0, 0), (10, 1.), Crop);
	AttachLegend("ndf=2");

	NewPad("$\ch^2$", "");
	draw(rGetObj(f, "g_chisq_ndf3"), "l", heavygreen);
	draw(rGetObj(f, "h_chisq_ndf3"), "vl", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/h_5,exp1"), "vl", black, "5 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/h_6,exp2"), "vl", red, "6 bins, exp2");
	limits((0, 0), (10, 1.), Crop);
	AttachLegend("ndf=3");
	
	NewRow();

	NewPad("$\ch^2$", "cumulative distribution");
	currentpad.xTicks = LeftTicks(1., 0.5);
	draw(rGetObj(f, "g_cum_chisq_ndf1"), "l", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_3,exp1"), "l", black, "3 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_4,exp2"), "l", red, "4 bins, exp2");
	limits((-0.5, 0), (10, 1.1), Crop);
	xaxis(YEquals(1, false), dashed);
	AttachLegend("ndf=1", E, E);

	NewPad("$\ch^2$", "cumulative distribution");
	currentpad.xTicks = LeftTicks(1., 0.5);
	draw(rGetObj(f, "g_cum_chisq_ndf2"), "l", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_4,exp1"), "l", black, "4 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_5,exp2"), "l", red, "5 bins, exp2");
	limits((-0.5, 0), (10, 1.1), Crop);
	xaxis(YEquals(1, false), dashed);
	AttachLegend("ndf=2", E, E);

	NewPad("$\ch^2$", "cumulative distribution");
	currentpad.xTicks = LeftTicks(1., 0.5);
	draw(rGetObj(f, "g_cum_chisq_ndf3"), "l", heavygreen);
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_5,exp1"), "l", black, "5 bins, exp1");
	draw(rGetObj(f, diagonal+"/"+binning+"/g_kol_6,exp2"), "l", red, "6 bins, exp2");
	limits((-0.5, 0), (10, 1.1), Crop);
	xaxis(YEquals(1, false), dashed);
	AttachLegend("ndf=3", E, E);
}

//----------------------------------------------------------------------------------------------------

for (int dgni : diagonals.keys)
{
	for (int bi : binnings.keys)
	{
		NewPage();

		NewPad(false, -1, yGridHintDef);
		label("{\SetFontSizesXX " + replace(diagonals[dgni], "_", "--") + ", " + binnings[bi] + "}");

		DrawOneSet(diagonals[dgni], binnings[bi]);
	}
}
