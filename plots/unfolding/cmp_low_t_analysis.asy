import root;
import pad_layout;

string dataset = "DS4";

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };

string topDir = "../../";
string topDir_old = "../../../4000GeV,beta90/";

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "multiplicative correction",
		xTicks=LeftTicks(0.05, 0.01), yTicks=RightTicks(0.005, 0.001));

for (int dgi : diagonals.keys)
{
	AddToLegend(dgn_labels[dgi] + ":");

	string f_old = topDir_old + dataset + "/unfolding_"+diagonals[dgi]+".root";
	//string ff = topDir_old + dataset + "/unfolding_fit_"+diagonals[dgi]+".root";
	draw(RootGetObject(f_old, "cf,eb/exp3/corr_final"), "d0,eb", StdPen(dgi));
	AddToLegend("low-$|t|$ analysis", mPl+5pt+StdPen(dgi));
	
	string f_this = topDir + dataset + "/unfolding_cf_"+diagonals[dgi]+".root";
	draw(RootGetObject(f_this, "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", StdPen(dgi), "this analysis");
}

limits((0, 0.97), (0.2, 1.01), Crop);
for (real x=0.0; x <= 0.2; x += 0.05)
	yaxis(XEquals(x, false), dotted);
for (real y=0.97; y <= 1.01; y += 0.005)
	xaxis(YEquals(y, false), dotted);
	
AttachLegend(SW, SW);
