import root;
import pad_layout;
include "../run_info.asy";

string dataset = "DS4";

string bunches[] = { "648", "2990", "26" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b--56t", "45t--56b" };

string topDir = "../../";

xSizeDef = 10cm;
ySizeDef = 6cm;

xTicksDef = LeftTicks(1., 0.5);

TGraph_errorBar = None;

real t_min = 19.5, t_max = 31.5;

drawGridDef = true;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

for (int dgni : diagonals.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX "+dgn_labels[dgni]+"}");
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "std.~dev.~of $\De^{R-L}\th_x^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.2, 0.1);
	DrawRunBands(-5, +10);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/resolution_fit_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "bunch "+bunches[bi]+"/gRMS_diffLR_th_x_vs_time");
		draw(swToHours*scale(1, 1e6), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);

		RootObject fit = RootGetObject(f, "bunch "+bunches[bi]+"/fit_x");
		draw(shift(0,  0.0)*swToHours*scale(1, 1e6), fit, "def", p+1pt);
	}
		
	limits((t_min, 8.2), (t_max, 9.8), Crop);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "std.~dev.~of $\De^{R-L}\th_y^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.1, 0.05);
	DrawRunBands(3, 3.8);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/resolution_fit_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "bunch "+bunches[bi]+"/gRMS_diffLR_th_y_vs_time");
		draw(swToHours*scale(1, 1e6), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);

		RootObject fit = RootGetObject(f, "bunch "+bunches[bi]+"/fit_y");
		draw(shift(0,  0.0)*swToHours*scale(1, 1e6), fit, "def", p+1pt);
	}
	
	limits((t_min, 3), (t_max, 3.8), Crop);
}

frame f_leg = BuildLegend();

NewPad(false);
attach(f_leg);

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "std.~dev.~of reconstructed vertex $x^*\ung{\mu m}$");
	currentpad.yTicks = RightTicks(10., 2);
	DrawRunBands(140, 200);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/distributions_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "time dependences/bunch "+bunches[bi]+"/gRMS_vtx_x_vs_time");
		draw(swToHours*scale(1, 1e3), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);
	}
	
	limits((t_min, 140), (t_max, 200), Crop);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "\vbox{\hbox{beam divergence in x $\ung{\mu m}$}\hbox{$\displaystyle \si^{\rm bd}_x = \sqrt 2\ \si[x^*]\over \be^*$}}");
	currentpad.yTicks = RightTicks(0.1, 0.05);
	DrawRunBands(2.2, 3.2);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/distributions_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "time dependences/bunch "+bunches[bi]+"/g_beam_div_x_vs_time");
		draw(swToHours*scale(1, 1e6), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);
	}
	
	limits((t_min, 2.2), (t_max, 3.2), Crop);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "\vbox{\hbox{sensor contribution to resolution $\ung{\mu rad}$}\hbox{$\displaystyle\sqrt{ {\si^2[\De^{R-L} \th_x^*]\over 2} - {\si^{\rm bd}_x}^2 }$}}");
	currentpad.yTicks = RightTicks(0.2, 0.1);
	DrawRunBands(5.2, 6.6);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/distributions_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "time dependences/bunch "+bunches[bi]+"/g_sensor_res_x_vs_time");
		draw(swToHours*scale(1, 1e6), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);
	}

	limits((t_min, 5.2), (t_max, 6.6), Crop);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "final (two-arm) resolution in $\th_x^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.1, 0.02);
	DrawRunBands(4, 5);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/resolution_fit_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "bunch "+bunches[bi]+"/gRMS_diffLR_th_x_vs_time");
		draw(swToHours*scale(1, 1e6)*scale(1, 0.5), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);
	}
		
	limits((t_min, 4), (t_max, 5), Crop);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "final (two-arm) resolution in $\th_y^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.05, 0.01);
	DrawRunBands(1.5, 1.9);

	for (int bi : bunches.keys)
	{
		pen p = StdPen(bi);
		string f = topDir+dataset+"/resolution_fit_"+diagonals[dgni]+".root";

		RootObject data = RootGetObject(f, "bunch "+bunches[bi]+"/gRMS_diffLR_th_y_vs_time");
		draw(swToHours*scale(1, 1e6)*scale(1, 0.5), data, "p,eb,d0", p, mCi+1pt+p, bunches[bi]);
	}
	
	limits((t_min, 1.5), (t_max, 1.9), Crop);
}


GShipout(vSkip=0mm);
