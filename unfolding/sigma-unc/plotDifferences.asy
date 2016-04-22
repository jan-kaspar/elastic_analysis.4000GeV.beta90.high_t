import root;
import pad_layout;

string f = "unfolding_cf_45b_56t.root";

string binning = "ob-1-30-0.10";
string iteration = "iteration 1";

string models[] = {
	"exp3+exp4",
	"exp5+erf*exp2",
	"p1*exp3+p1*exp1",
	"p1*exp3+p2*exp2",
	"exp3-intf-exp1",
	"(exp3-intf-exp1)*expG",
};

xSizeDef = 10cm;
ySizeDef = 8cm;

drawGridDef = true;

xTicksDef = LeftTicks(0.5, 0.1);

//----------------------------------------------------------------------------------------------------

void DrawDifference(rObject o1, rObject o2, pen p, string label = "")
{
	guide g;

	int N = o1.iExec("GetN");
	for (int i = 0; i < N; ++i)
	{
		real ax[] = {0.};
		real ay1[] = {0.};
		real ay2[] = {0.};

		o1.vExec("GetPoint", i, ax, ay1);
		o2.vExec("GetPoint", i, ax, ay2);

		g = g--(ax[0], ay2[0] - ay1[0]);
	}

	draw(g, p, label);
}

//----------------------------------------------------------------------------------------------------

void DrawDifferenceSum(rObject ob, rObject o1, rObject o2, pen p, string label = "")
{
	guide g;

	int N = o1.iExec("GetN");
	for (int i = 0; i < N; ++i)
	{
		real ax[] = {0.};
		real ayb[] = {0.};
		real ay1[] = {0.};
		real ay2[] = {0.};

		ob.vExec("GetPoint", i, ax, ayb);
		o1.vExec("GetPoint", i, ax, ay1);
		o2.vExec("GetPoint", i, ax, ay2);

		g = g--(ax[0], (-ay2[0] + ayb[0]) + (-ay1[0] + ayb[0]));
	}

	draw(g, p, label);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "difference in unfolding correction");
for (int mi : models.keys)
	DrawDifference(
		rGetObj(f, binning+"/"+models[mi]+"/+1,+0/"+iteration+"/g_corr"),
		rGetObj(f, binning+"/"+models[mi]+"/+0,+0/"+iteration+"/g_corr"),
		StdPen(mi), models[mi]
	);
limits((0, -0.002), (2.0, 0.01), Crop);
AttachLegend("effect due to uncertainty in $\si[\th_x]$", NW, NE);
yaxis(XEquals(1.9, false), dashed);

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad("$|t|\ung{GeV^2}$", "difference in unfolding correction");
for (int mi : models.keys)
	DrawDifference(
		rGetObj(f, binning+"/"+models[mi]+"/+0,+1/"+iteration+"/g_corr"),
		rGetObj(f, binning+"/"+models[mi]+"/+0,+0/"+iteration+"/g_corr"),
		StdPen(mi), models[mi]
	);
limits((0, -0.002), (2.0, 0.01), Crop);
AttachLegend("effect due to uncertainty in $\si[\th_y]$", NW, NE);
yaxis(XEquals(1.9, false), dashed);

//----------------------------------------------------------------------------------------------------
/*
NewRow();

NewPad("$|t|\ung{GeV^2}$", "difference in unfolding correction");
for (int mi : models.keys)
{
	DrawDifference(
		rGetObj(f, binning+"/"+models[mi]+"/+1,+1/"+iteration+"/g_corr"),
		rGetObj(f, binning+"/"+models[mi]+"/+0,+0/"+iteration+"/g_corr"),
		StdPen(mi), models[mi]
	);

	DrawDifferenceSum(
		rGetObj(f, binning+"/"+models[mi]+"/+0,+0/"+iteration+"/g_corr"),
		rGetObj(f, binning+"/"+models[mi]+"/+0,+1/"+iteration+"/g_corr"),
		rGetObj(f, binning+"/"+models[mi]+"/+1,+0/"+iteration+"/g_corr"),
		StdPen(mi)+dashed, models[mi]
	);

}
limits((0, -0.002), (2.0, 0.01), Crop);
AttachLegend("effect due to uncertainty in $\si[\th_x]$ and $\si[\th_y]$", NW, NE);
yaxis(XEquals(1.9, false), dashed);
*/
