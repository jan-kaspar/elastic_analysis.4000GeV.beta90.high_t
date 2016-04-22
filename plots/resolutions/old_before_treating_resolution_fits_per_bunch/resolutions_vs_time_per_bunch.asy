import root;
import pad_layout;
include "../run_info.asy";

//string datasets[] = { "DS4-b26" };
string datasets[] = { /*"DS4",*/ "DS4-b648", "DS4-b2990", "DS4-b26" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b-56t", "45t-56b" };

string topDir = "../../";

xSizeDef = 10cm;
ySizeDef = 6cm;

xTicksDef = LeftTicks(1., 0.5);

TGraph_errorBar = None;

real t_min = 19.5, t_max = 31.5;

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

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/gRMS_diffLR_th_x_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);

		rObject fit = rGetObj(topDir+"resolutions/fit_"+datasets[dsi]+".root", diagonals[dgni]+"/x/g_fit");
		//draw(shift(0, +0.2)*swToHours*scale(1, 1e6), fit, "l", black+dashed);
		draw(shift(0,  0.0)*swToHours*scale(1, 1e6), fit, "l", p+1pt);
		//draw(shift(0, -0.2)*swToHours*scale(1, 1e6), fit, "l", black+dashed);
	}
		
	limits((t_min, 8), (t_max, 10), Crop);
	for (real y=8; y <= 10; y += 0.2)
		xaxis(YEquals(y, false), dotted);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "std.~dev.~of $\De^{R-L}\th_y^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.1, 0.05);
	DrawRunBands(3, 3.8);

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/gRMS_diffLR_th_y_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);

		rObject fit = rGetObj(topDir+"resolutions/fit_"+datasets[dsi]+".root", diagonals[dgni]+"/y/g_fit");
		//draw(shift(0, +0.05)*swToHours*scale(1, 1e6), fit, "l", black+dashed);
		draw(shift(0,  0.00)*swToHours*scale(1, 1e6), fit, "l", p+1pt);
		//draw(shift(0, -0.05)*swToHours*scale(1, 1e6), fit, "l", black+dashed);
	}
	
	limits((t_min, 3), (t_max, 3.8), Crop);
	for (real y=3; y <= 3.8; y += 0.1)
		xaxis(YEquals(y, false), dotted);
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

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e3),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/gRMS_vtx_x_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);
	}
	
	limits((t_min, 140), (t_max, 200), Crop);
	for (real y = 140; y <= 200; y += 10)
		xaxis(YEquals(y, false), dotted);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "\vbox{\hbox{beam divergence in x $\ung{\mu m}$}\hbox{$\displaystyle \si^{\rm bd}_x = \sqrt 2\ \si[x^*]\over \be^*$}}");
	currentpad.yTicks = RightTicks(0.1, 0.05);
	DrawRunBands(2.2, 3.2);

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/g_beam_div_x_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);
	}
	
	limits((t_min, 2.2), (t_max, 3.2), Crop);
	for (real y=2.2; y <= 3.2; y += 0.1)
		xaxis(YEquals(y, false), dotted);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "\vbox{\hbox{sensor contribution to resolution $\ung{\mu m}$}\hbox{$\displaystyle\sqrt{ {\si^2[\De^{R-L} \th_x^*]\over 2} - {\si^{\rm bd}_x}^2 }$}}");
	currentpad.yTicks = RightTicks(0.2, 0.1);
	DrawRunBands(5.2, 6.6);

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/g_sensor_res_x_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);
	}

	limits((t_min, 5.2), (t_max, 6.6), Crop);
	for (real y=5.2; y <= 6.6; y += 0.2)
		xaxis(YEquals(y, false), dotted);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "final (two-arm) resolution in $\th_x^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.1, 0.02);
	DrawRunBands(4, 5);

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6)*scale(1, 0.5),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/gRMS_diffLR_th_x_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);
	}
		
	limits((t_min, 4), (t_max, 5), Crop);
	for (real y = 4; y <= 5; y += 0.1)
		xaxis(YEquals(y, false), dotted);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dgni : diagonals.keys)
{
	NewPad("time $\ung{h}$", "final (two-arm) resolution in $\th_y^*\ung{\mu rad}$");
	currentpad.yTicks = RightTicks(0.05, 0.01);
	DrawRunBands(1.5, 1.9);

	for (int dsi : datasets.keys)
	{
		pen p = StdPen(dsi);

		draw(swToHours*scale(1, 1e6)*scale(1, 0.5),
			rGetObj(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "time dependences/gRMS_diffLR_th_y_vs_time"), "p,eb,d0",
			p, mCi+1pt+p, datasets[dsi]);
	}
	
	limits((t_min, 1.5), (t_max, 1.9), Crop);
	for (real y = 1.5; y <= 1.9; y += 0.05)
		xaxis(YEquals(y, false), dotted);
}



GShipout(vSkip=0mm);
