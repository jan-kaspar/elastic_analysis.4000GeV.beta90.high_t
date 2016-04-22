import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string binning = "ob-1-30-0.10";
//string binning = "ob-2-20-0.20";

drawGridDef = true;

TH1_x_max = 1.9;

xSizeDef = 10cm;

pen p_unc_full = yellow;
pen p_unc_anal = red;

//----------------------------------------------------------------------------------------------------

void DrawUncertaintyHist(real yScale=100, rObject o, pen p=black)
{
	int N = o.iExec("GetNbinsX");
	for (int bi = 1; bi <= N; ++bi)
	{
		real c = o.rExec("GetBinCenter", bi);
		real w = o.rExec("GetBinWidth", bi);
		real v = o.rExec("GetBinContent", bi);
		real u = o.rExec("GetBinError", bi);

		if (c > TH1_x_max)
			continue;

		real u_rel = (v > 0) ? u/v : 0.;

		if (v > 0)
		{
			draw((c-w/2, 0.)--(c+w/2, 0.), p);
			draw((c, -u_rel*yScale)--(c, +u_rel*yScale), p);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void DrawBandGraph(real yScale=100, rObject o, pen p=lightblue)
{
	guide g_u, g_b;
	
	int N = o.iExec("GetN");
	for (int i = 0; i < N; ++i)
	{
		real xa[] = {0.};
		real ya[] = {0.};
		o.vExec("GetPoint", i, xa, ya);
		real x = xa[0];
		real y = ya[0];

		//if (x < rel_diff_low_limit || x > rel_diff_high_limit)
		//	continue;
	
		g_u = g_u -- (x, +y*yScale);
		g_b = g_b -- (x, -y*yScale);
	}

	g_b = reverse(g_b);
	filldraw(g_u--g_b--cycle, p, nullpen);
}

//----------------------------------------------------------------------------------------------------

void DrawBandAroundFit(rObject of, rObject ou, pen p=lightblue)
{
	guide g_u, g_b;
	
	int N = of.iExec("GetN");
	for (int i = 0; i < N; ++i)
	{
		real xa[] = {0.};
		real ya[] = {0.};
		of.vExec("GetPoint", i, xa, ya);
		real x = xa[0];
		real y = ya[0];

		real u = ou.rExec("Eval", x);
	
		g_u = g_u -- Scale((x, y * (1.+u)));
		g_b = g_b -- Scale((x, y * (1.-u)));
	}

	g_b = reverse(g_b);
	filldraw(g_u--g_b--cycle, p, nullpen);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

// TODO: replace the temp objects
rObject o_dsdt = rGetObj(topDir+"DS-merged/merged.root", binning+"/"+dataset+"/combined/h_dsdt");
rObject o_unc_anal = rGetObj(topDir+dataset+"/systematics_matrix.root", "matrices/all-temp-anal/combined/g_envelope");
rObject o_unc_full = rGetObj(topDir+dataset+"/systematics_matrix.root", "matrices/all-temp/combined/g_envelope");
rObject o_fit = rGetObj(topDir+"models/exp3-intf-exp1.root", "g_dsdt");

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si_{\rm el}/\d t \ung{mb/GeV^2}$", ySize = 8cm);
currentpad.xTicks = LeftTicks(0.2, 0.1);
scale(Linear, Log);

DrawBandAroundFit(o_fit, o_unc_full, p_unc_full);
DrawBandAroundFit(o_fit, o_unc_anal, p_unc_anal);

draw(o_dsdt, "eb", black);

limits((0, 1e-4), (2., 1e3), Crop);

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad("$|t|\ung{GeV^2}$", "relative uncertainty $\ung{\%}$", ySize = 3.5cm);
currentpad.xTicks = LeftTicks(0.2, 0.1);
currentpad.yTicks = RightTicks(20., 10.);

DrawBandGraph(o_unc_full, p_unc_full);
DrawBandGraph(o_unc_anal, p_unc_anal);
DrawUncertaintyHist(o_dsdt);

limits((0, -60), (2., +60), Crop);

//----------------------------------------------------------------------------------------------------
/*
NewRow();

NewPad("$|t|\ung{GeV^2}$", "relative uncertainty $\ung{\%}$", ySize = 3.5cm);
currentpad.xTicks = LeftTicks(0.2, 0.1);
currentpad.yTicks = RightTicks(2., 1.);

DrawBandGraph(o_unc_full, p_unc_full);
DrawBandGraph(o_unc_anal, p_unc_anal);
DrawUncertaintyHist(o_dsdt);

limits((0, -8), (2., +8), Crop);
*/

GShipout(vSkip=1mm, margin=1mm);
