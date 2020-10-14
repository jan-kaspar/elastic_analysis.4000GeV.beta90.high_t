import root;
import pad_layout;

string topDir = "../../";

xSizeDef = 12cm;
ySizeDef = 8cm;

string label;

void NewPlot(string _label)
{
	label = _label;

	NewPad("$|t|\ung{GeV^2}$", "multiplicative unsmearing correction");
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	currentpad.yTicks = RightTicks(0.05, 0.01);
}

void FinalizePlot()
{
	limits((0, 0.80), (2, 1.05), Crop);
	//limits((0, 0.50), (0.4, 1.5), Crop);
	//xlimits(0, 0.4, Crop);

	AttachLegend(label, NW, NE);

	xaxis(YEquals(1, false), dashed);
}

string diagonals[], dgn_labels[];
diagonals.push("45b_56t"); dgn_labels.push("45b -- 56t");
diagonals.push("45t_56b"); dgn_labels.push("45t -- 56b");

//----------------------------------------------------------------------------------------------------

// ---------- CF method ----------

NewPlot("diagonal cmp");
AddToLegend("<bt1, exp3-intf-exp1, iteration 0");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "bt1/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", black, "45 bot -- 56 top");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "bt1/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", red, "45 top -- 56 bot");
FinalizePlot();

NewRow();

for (int dgni : diagonals.keys)
{
	NewPlot("iteration cmp");
	AddToLegend("<bt1, exp3-intf-exp1, " + dgn_labels[dgni]);
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", black, "1");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", red, "2");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3-intf-exp1/+0,+0/iteration 2/g_corr"), "l", blue, "3");
	FinalizePlot();
}

NewRow();

for (int dgni : diagonals.keys)
{
	NewPlot("parametrisation cmp");
	AddToLegend("<bt1, iteration 2, " + dgn_labels[dgni]);
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3+exp4/+0,+0/iteration 2/g_corr"), "l", black, "exp3+exp4");
	//draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp5+erf*exp2/+0,+0/iteration 2/g_corr"), "l", red, "exp5+erf*exp2");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/p1*exp3+p1*exp1/+0,+0/iteration 2/g_corr"), "l", blue, "p1*exp3+p1*exp1");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/p1*exp3+p2*exp2/+0,+0/iteration 2/g_corr"), "l", cyan, "p1*exp3+p2*exp2");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3-intf-exp1/+0,+0/iteration 2/g_corr"), "l", heavygreen, "exp3-intf-exp1");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/(exp3-intf-exp1)*expG/+0,+0/iteration 2/g_corr"), "l", magenta, "(exp3-intf-exp1)*expG");
	FinalizePlot();
}

NewRow();

for (int dgni : diagonals.keys)
{
	NewPlot("binning cmp, " + dgn_labels[dgni]);
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", black, "bt1");
	draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", red, "ob-1-30-0.10");
	//draw(RootGetObject(topDir + "DS4/unfolding_cf_" + diagonals[dgni]+ ".root", "ob-3-10-0.30/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", blue, "ob-1-30-0.30");
	FinalizePlot();
}

GShipout("comparison_CF", margin=1mm);

NewPlot("validation of the NI calculation (wrt. MC calculation)");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", black, "NI: graph");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/corr"), "vl", black, "NI: histogram");

draw(RootGetObject(topDir + "DS4/smearing_matrix_mc_45b_56t.root", "exp3-intf-exp1/bt1/C"), "vl", red, "MC");
FinalizePlot();

GShipout("comparison_CF_validation", margin=1mm);

// ---------- GR method ----------

NewPlot("smoothing-level comparison");
AddToLegend("<bt1, 45t -- 56b, exp3-intf-exp1");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E-03/h_corr"), "vl", heavygreen, "$\al = 1\cdot10^{-3}$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E-02/h_corr"), "vl", blue, "$\al = 1\cdot10^{-2}$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E-01/h_corr"), "vl", red, "$\al = 1\cdot10^{-1}$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", black+1pt, "$\al = 1\cdot10^{0}$");
FinalizePlot();

NewRow();

NewPlot("diagonal and binning comparison");
AddToLegend("<45b -- 56t");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black, "ob-1-30-0.10, $\al = 1$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", red, "bt1, $\al = 1$");

AddToLegend("<45t -- 56b");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", blue, "ob-1-30-0.10, $\al = 1$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", heavygreen, "bt1, $\al = 1$");
FinalizePlot();

/*
NewPlot("uncertainty estimate");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,bt1/alpha=1.00E-03/h_corr"), "vl,eb", black, "45 bot -- 56 top, bt1");
FinalizePlot();
*/

NewRow();

NewPlot("model and binning comparison");
AddToLegend("<exp3-intf-exp1");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", blue, "ob-1-30-0.10, $\al = 1$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", magenta, "bt1, $\al = 1$");

AddToLegend("<p1*exp3+p2*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black, "ob-1-30-0.10, $\al = 1$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,bt1/alpha=1.00E+00/h_corr"), "vl", red, "bt1, $\al = 1$");
FinalizePlot();

GShipout("comparison_GR", margin=1mm);

// ---------- CF vs. GR method ----------

NewPlot("45 bot -- 56 top, bt1");
AddToLegend("<CF:");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "bt1/p1*exp3+p2*exp2/+0,+0/iteration 1/g_corr"), "l", blue, "p1*exp3+p2*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", heavygreen, "exp3-intf-exp1");

AddToLegend("<GR:");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,bt1/alpha=1.00E+00/h_corr"), "vl", black, "p1*exp3+p2*exp2, $\al=1\cdot10^{0}$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "bt1/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", red, "exp3-intf-exp1, $\al=1\cdot10^{0}$");
FinalizePlot();

NewPlot("45 top -- 56 bot, bt1");
AddToLegend("<CF:");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "bt1/p1*exp3+p2*exp2/+0,+0/iteration 1/g_corr"), "l", blue, "p1*exp3+p2*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "bt1/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", heavygreen, "exp3-intf-exp1");

AddToLegend("<GR:");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,p1*exp3+p2*exp2,bt1/alpha=1.00E+00/h_corr"), "vl", black, "p1*exp3+p2*exp2, $\al=1\cdot10^{0}$");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "bt1/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,bt1/alpha=1.00E+00/h_corr"), "vl", red, "exp3-intf-exp1, $\al=1\cdot10^{0}$");
FinalizePlot();

GShipout("comparison_CF_GR", margin=1mm);
