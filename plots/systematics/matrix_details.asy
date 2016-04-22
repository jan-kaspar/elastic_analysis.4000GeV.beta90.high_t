import root;
import pad_layout;

xSizeDef = 10cm;
ySizeDef = 8cm;


string f_matrix = "../../DS4/systematics_matrix.root";

real y_min = -0.04, y_max = 0.04;
yTicksDef = RightTicks(0.02, 0.01);

//----------------------------------------------------------------------------------------------------

void DrawOne(string tag)
{
	NewPad();

	draw(rGetObj(f_matrix, "contributions/" + tag + "/g_eff_45b"), red);
	draw(rGetObj(f_matrix, "contributions/" + tag + "/g_eff_45t"), blue);
	draw(rGetObj(f_matrix, "contributions/" + tag + "/g_eff_comb1"), heavygreen+1pt);

	rObject obj = rGetObj(f_matrix, "contributions/" + tag + "/g_eff_comb2", error=false);
	if (obj.valid)
		draw(obj, magenta+1pt);

	limits((0, y_min), (2., y_max), Crop);
	AttachLegend(tag, NW, NW);
	
	yaxis(XEquals(0.16, false), dotted);
}

//----------------------------------------------------------------------------------------------------

y_min = -0.04; y_max = 0.04;
yTicksDef = RightTicks(0.01, 0.002);
NewRow();
DrawOne("alig-sh-thx");
DrawOne("alig-sh-thy");
DrawOne("tilt-thx-thy");

NewRow();
DrawOne("opt-m1");
DrawOne("opt-m2");

y_min = -0.002; y_max = 0.002;
yTicksDef = RightTicks(0.001, 0.0002);
NewRow();
DrawOne("acc-corr-sigma-unc");
DrawOne("acc-corr-sigma-asym");
DrawOne("acc-corr-non-gauss");

NewRow();
DrawOne("eff-slp");

y_min = -0.04; y_max = 0.04;
yTicksDef = RightTicks(0.01, 0.002);
NewRow();
DrawOne("beam-mom");

NewRow();
// unfolding
//DrawOne("");

y_min = -0.00; y_max = 0.08;
yTicksDef = RightTicks(0.01, 0.002);
NewRow();
DrawOne("norm");
