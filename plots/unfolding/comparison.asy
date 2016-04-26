import root;
import pad_layout;

string topDir = "../../";

xSizeDef = 12cm;
ySizeDef = 8cm;

string label;

void NewPlot(string _label)
{
	label = _label;

	NewRow();
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

	for (real x = 0; x <= 2; x += 0.2)
		yaxis(XEquals(x, false), dotted);
	
	for (real y = 0.80; y <= 1.05; y += 0.05)
		xaxis(YEquals(y, false), (fabs(y - 1) < 1e-4) ? dashed : dotted);
}

//----------------------------------------------------------------------------------------------------

// ---------- CF method ----------

NewPlot("diagonal cmp");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", black, "45 bot -- 56 top");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", red, "45 top -- 56 bot");
FinalizePlot();

NewPlot("iteration cmp");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 0/g_corr"), "l", black, "1");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", red, "2");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 2/g_corr"), "l", blue, "3");
FinalizePlot();

NewPlot("parametrisation cmp");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp3+exp4/+0,+0/iteration 1/g_corr"), "l", black, "exp3+exp4");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp5+erf*exp2/+0,+0/iteration 1/g_corr"), "l", red, "exp5+erf*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/p1*exp3+p1*exp1/+0,+0/iteration 1/g_corr"), "l", blue, "p1*exp3+p1*exp1");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/p1*exp3+p2*exp2/+0,+0/iteration 1/g_corr"), "l", cyan, "p1*exp3+p2*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", heavygreen, "exp3-intf-exp1");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/(exp3-intf-exp1)*expG/+0,+0/iteration 1/g_corr"), "l", magenta, "(exp3-intf-exp1)*expG");
FinalizePlot();

NewPlot("binning cmp");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", black, "ob-1-30-0.10");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-2-20-0.20/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", red, "ob-2-20-0.20");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-3-10-0.30/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", blue, "ob-1-30-0.30");
FinalizePlot();

NewPlot("validation of the NI calculation (wrt. MC calculation)");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", black, "NI: graph");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/corr"), "vl", black, "NI: histogram");

draw(RootGetObject(topDir + "DS4/smearing_matrix_mc_45b_56t.root", "exp3-intf-exp1/ob-1-30-0.10/C"), "vl", red, "MC");
FinalizePlot();

GShipout("comparison_CF", margin=1mm);

// ---------- GR method ----------

NewPlot("diagonal and binning comparison");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black, "45 bot -- 56 top, ob-1-30-0.10");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-3-10-0.30/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-3-10-0.30/alpha=5.00E-01/h_corr"), "vl", red, "45 bot -- 56 top, ob-3-10-0.30");

draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", blue, "45 top -- 56 bot, ob-1-30-0.10");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-3-10-0.30/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-3-10-0.30/alpha=5.00E-01/h_corr"), "vl", heavygreen, "45 top -- 56 bot, ob-3-10-0.30");
FinalizePlot();

NewPlot("uncertainty estimate");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl,eb", black, "45 bot -- 56 top, ob-1-30-0.10");
FinalizePlot();

NewPlot("model and binning comparison");
AddToLegend("<exp3-intf-exp1");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", blue, "ob-1-30-0.10");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-3-10-0.30/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-3-10-0.30/alpha=5.00E-01/h_corr"), "vl", magenta, "ob-3-10-0.30");

AddToLegend("<p1*exp3+p2*exp2");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black, "ob-1-30-0.10");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-3-10-0.30/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-3-10-0.30/alpha=5.00E-01/h_corr"), "vl", red, "ob-3-10-0.30");
FinalizePlot();

NewPlot("smoothing-level comparison");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=5.00E-02/h_corr"), "vl", black, "5.00E-02");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E-01/h_corr"), "vl", red, "1.00E-01");
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", blue+1pt, "1.00E+00");
FinalizePlot();

GShipout("comparison_GR", margin=1mm);

// ---------- CF vs. GR method ----------

NewPlot("CF vs.~GR, 45 bot -- 56 top");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/p1*exp3+p2*exp2/+0,+0/iteration 1/g_corr"), "l", blue);
draw(RootGetObject(topDir + "DS4/unfolding_cf_45b_56t.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", heavygreen);

draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black);
draw(RootGetObject(topDir + "DS4/unfolding_gr_45b_56t.root", "ob-1-30-0.10/smearing_matrix_mc_45b_56t.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", red);
FinalizePlot();

NewPlot("CF vs.~GR, 45 top -- 56 bot");
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/p1*exp3+p2*exp2/+0,+0/iteration 1/g_corr"), "l", blue);
draw(RootGetObject(topDir + "DS4/unfolding_cf_45t_56b.root", "ob-1-30-0.10/exp3-intf-exp1/+0,+0/iteration 1/g_corr"), "l", heavygreen);

draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,p1*exp3+p2*exp2,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", black);
draw(RootGetObject(topDir + "DS4/unfolding_gr_45t_56b.root", "ob-1-30-0.10/smearing_matrix_mc_45t_56b.root,exp3-intf-exp1,ob-1-30-0.10/alpha=1.00E+00/h_corr"), "vl", red);
FinalizePlot();

GShipout("comparison_CF_GR", margin=1mm);
