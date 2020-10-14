import root;
import pad_layout;

xSizeDef = 10cm;
ySizeDef = 8cm;

string f_matrix = "../../DS4/systematics_matrix.root";

real y_min = -0.04, y_max = 0.04;
yTicksDef = RightTicks(0.02, 0.01);

xTicksDef = LeftTicks(0.5, 0.1);

//----------------------------------------------------------------------------------------------------

void DrawOne(string tag)
{
	NewPad("$|t|\ung{GeV^2}$", "relative effect on $\d\sigma/\d t$");

	draw(RootGetObject(f_matrix, "contributions/" + tag + "/g_eff_45b"), red);
	draw(RootGetObject(f_matrix, "contributions/" + tag + "/g_eff_45t"), blue);
	draw(RootGetObject(f_matrix, "contributions/" + tag + "/g_eff_comb1"), heavygreen+2pt);

	RootObject obj = RootGetObject(f_matrix, "contributions/" + tag + "/g_eff_comb2", error=false);
	if (obj.valid)
		draw(obj, magenta+2pt);

	limits((0, y_min), (2., y_max), Crop);
	AttachLegend(tag, NW, NW);
	
	yaxis(XEquals(1.86, false), dashed);
}

//----------------------------------------------------------------------------------------------------

y_min = -0.02; y_max = 0.05;
yTicksDef = RightTicks(0.01, 0.002);
NewRow();
DrawOne("alig-sh-thx");
DrawOne("alig-sh-thy");
DrawOne("tilt-thx-thy");

NewRow();
DrawOne("opt-m1");
DrawOne("opt-m2");

//y_min = -0.002; y_max = 0.002;
//yTicksDef = RightTicks(0.001, 0.0002);
NewRow();
DrawOne("acc-corr-sigma-unc");
DrawOne("acc-corr-sigma-asym");
DrawOne("acc-corr-non-gauss");

NewRow();
DrawOne("eff-slp");
DrawOne("beam-mom");

NewRow();
DrawOne("unsm-sigma-x");
DrawOne("unsm-sigma-y");
DrawOne("unsm-model");

//y_min = -0.00; y_max = 0.08;
//yTicksDef = RightTicks(0.01, 0.002);
NewRow();
DrawOne("norm");

NewPad(false);
AddToLegend("diagonal 45b -- 56t", red);
AddToLegend("diagonal 45t -- 56b", blue);
AddToLegend("combination 1", heavygreen + 2pt);
AddToLegend("combination 2 (only if uncorrelated)", magenta + 2pt);
AttachLegend();

GShipout(vSkip=1mm);
