import root;
import pad_layout;

string datasets[] = { "DS4" };

string topDir = "../../";

string binnings[] = {
	"ob-1-30-0.10",
	"ob-2-20-0.20",
	"ob-3-10-0.30"
};

//----------------------------------------------------------------------------------------------------

void PlotRatio(RootObject num, RootObject den, pen p = black, string label = "")
{
	int bins = num.iExec("GetNbinsX");

	for (int bi = 1; bi <= bins; ++bi)
	{
		real c = num.rExec("GetBinCenter", bi);
		real w = num.rExec("GetBinWidth", bi);

		real v_n = num.rExec("GetBinContent", bi);
		real u_n = num.rExec("GetBinError", bi);

		real v_d = den.rExec("GetBinContent", bi);
		real u_d = den.rExec("GetBinError", bi);

		real r = (v_d > 0) ? v_n / v_d : 0.;
		real r_u = (v_n > 0 && v_d > 0) ? r * sqrt((u_n/v_n)^2 + (u_d/v_d)^2) : 0.;

		if (r != 0)
		{
			draw((c-w/2, r)--(c+w/2, r), p);
			draw((c, r-r_u)--(c, r+r_u), p);
		}
	}

	if (label != "")
		AddToLegend(label, mPl+5pt+p);
}

//----------------------------------------------------------------------------------------------------

void PlotPulls(RootObject num, RootObject den, pen pe = black, string label = "")
{
	int bins = num.iExec("GetNbinsX");

	for (int bi = 1; bi <= bins; ++bi)
	{
		real c = num.rExec("GetBinCenter", bi);
		real w = num.rExec("GetBinWidth", bi);

		real v_n = num.rExec("GetBinContent", bi);
		real u_n = num.rExec("GetBinError", bi);

		real v_d = den.rExec("GetBinContent", bi);
		real u_d = den.rExec("GetBinError", bi);

		real r = (v_d > 0) ? v_n / v_d : 0.;
		real r_u = (v_n > 0 && v_d > 0) ? r * sqrt((u_n/v_n)^2 + (u_d/v_d)^2) : 0.;

		real p = (r_u > 0) ? (r - 1.) / r_u : 0.;
		real p_u = (r_u > 0.) ? 1. : 0.;

		if (p != 0 && p_u != 0)
		{
			draw((c-w/2, p)--(c+w/2, p), pe);
			draw((c, p-p_u)--(c, p+p_u), pe);
		}
	}

	if (label != "")
		AddToLegend(label, mPl+5pt+pe);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewPad("$|t|\ung{GeV^2}$", "$\displaystyle R = {\hbox{45 bot -- 56 top} \over \hbox{45 top -- 56 bot}}$", 12cm, 8cm);
	currentpad.xTicks = LeftTicks(0.5, 0.1);
	currentpad.yTicks = RightTicks(0.5, 0.1);
	
	for (int bi : binnings.keys)
	{
		RootObject o_45b = RootGetObject(topDir+datasets[dsi]+"/distributions_45b_56t.root", "normalization/"+binnings[bi]+"/h_t_normalized");
		RootObject o_45t = RootGetObject(topDir+datasets[dsi]+"/distributions_45t_56b.root", "normalization/"+binnings[bi]+"/h_t_normalized");

		PlotRatio(o_45b, o_45t, StdPen(bi+1), binnings[bi]);
	}
	
	limits((0, -0.5), (1.9, 2.5), Crop);
	for (real y = -0.5; y <= 2.5; y += 0.5)
		xaxis(YEquals(y, false), (fabs(y - 1) < 1e-3) ? dashed : dotted);

	AttachLegend(datasets[dsi]);

	//--------------------

	NewPad("$|t|\ung{GeV^2}$", "$\displaystyle R = {\hbox{45 bot -- 56 top} \over \hbox{45 top -- 56 bot}}$", 12cm, 8cm);
	currentpad.xTicks = LeftTicks(0.1, 0.02);
	currentpad.yTicks = RightTicks(0.05, 0.01);
	
	for (int bi : binnings.keys)
	{
		RootObject o_45b = RootGetObject(topDir+datasets[dsi]+"/distributions_45b_56t.root", "normalization/"+binnings[bi]+"/h_t_normalized");
		RootObject o_45t = RootGetObject(topDir+datasets[dsi]+"/distributions_45t_56b.root", "normalization/"+binnings[bi]+"/h_t_normalized");

		PlotRatio(o_45b, o_45t, StdPen(bi+1), binnings[bi]);
	}
	
	limits((0, 0.9), (0.4, 1.1), Crop);
	for (real y = 0.9; y <= 1.1; y += 0.05)
		xaxis(YEquals(y, false), (fabs(y - 1) < 1e-3) ? dashed : dotted);

	AttachLegend(datasets[dsi]);
}

//----------------------------------------------------------------------------------------------------

NewRow();

for (int dsi : datasets.keys)
{
	NewPad("$|t|\ung{GeV^2}$", "$\displaystyle {R - 1\over \si[R]}\ ,\quad R = {\hbox{45 bot -- 56 top} \over \hbox{45 top -- 56 bot}}$", 12cm, 8cm);
	currentpad.xTicks = LeftTicks(0.5, 0.1);
	currentpad.yTicks = RightTicks(1., 0.2);
	
	for (int bi : binnings.keys)
	{
		RootObject o_45b = RootGetObject(topDir+datasets[dsi]+"/distributions_45b_56t.root", "normalization/"+binnings[bi]+"/h_t_normalized");
		RootObject o_45t = RootGetObject(topDir+datasets[dsi]+"/distributions_45t_56b.root", "normalization/"+binnings[bi]+"/h_t_normalized");

		PlotPulls(o_45b, o_45t, StdPen(bi+1), binnings[bi]);
	}
	
	limits((0, -5), (1.9, +5), Crop);
	for (real y = -5; y <= 5; y += 1)
		xaxis(YEquals(y, false), (fabs(y - 0) < 1e-3) ? dashed : dotted);

	AttachLegend(datasets[dsi]);
}
