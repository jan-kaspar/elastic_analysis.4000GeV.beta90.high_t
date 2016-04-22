import root;
import pad_layout;

xSizeDef = 10cm;
ySizeDef = 8cm;

yTicksDef = RightTicks(0.02, 0.01);

string f_45b = "../../DS4/systematics_ni_45b_56t.root";
string f_45t = "../../DS4/systematics_ni_45t_56b.root";

//string f_matrix = "../../DS4/systematics_matrix.root";

string models[] = {
	"exp3-intf-exp1",
	"p1*exp3+p2*exp2"
};

real y_min = 0.98, y_max = 1.02;

//yTicksDef = RightTicks(0.001, 0.0005);

//----------------------------------------------------------------------------------------------------

void DrawOne(string tag)
{
	NewPad();

	for (int mi : models.keys)
	{
		pen p = solid;
		if (mi == 1) p = dashed;
		if (mi == 2) p = longdashed;
		
		draw(rGetObj(f_45b, tag + "/" + models[mi] + "/g_r"), blue+p);
		draw(rGetObj(f_45t, tag + "/" + models[mi] + "/g_r"), red+p);

		//draw(rGetObj(f_matrix, "contributions/" + tag + "/g_eff_comb1"), heavygreen+p);
	}

	//limits((0, y_min), (2., y_max), Crop);
	limits((0, 0.97), (2., 1.05), Crop);
	AttachLegend(tag, NW, NW);
	
	//yaxis(XEquals(0.16, false), dotted);
}

//----------------------------------------------------------------------------------------------------

y_min = 0.96; y_max = 1.04;
NewRow();
DrawOne("alig-sh-thx");
DrawOne("alig-sh-thy");
DrawOne("tilt-thx-thy");

NewRow();
DrawOne("opt-m1");
DrawOne("opt-m2");

y_min = 0.995; y_max = 1.005;
NewRow();
DrawOne("acc-corr-sigma-unc");
DrawOne("acc-corr-sigma-asym");
DrawOne("acc-corr-non-gauss");

NewRow();
DrawOne("eff-slp");

y_min = 0.96; y_max = 1.04;
NewRow();
DrawOne("beam-mom");

NewRow();
// unfolding
//DrawOne("");

NewRow();
DrawOne("norm");
