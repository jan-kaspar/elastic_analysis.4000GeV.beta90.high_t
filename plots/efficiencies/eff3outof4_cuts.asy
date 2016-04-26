import root;
import pad_layout;

string datasets[] = { "DS4" };

string diagonals[] = { "45b_56t", "45t_56b" };

string RPs[] = { "L_F", "L_N", "R_N", "R_F" };
string RP_labels[] = { "left far", "left near", "right near", "right far" };

string topDir = "../../";

TH2_palette = Gradient(blue, cyan, heavygreen, yellow, magenta);
//TH2_palette.push(black);

xSizeDef = 10cm;
ySizeDef = 5cm;

TH2_paletteBarSpacing = 0.02;
TH2_paletteBarWidth = 0.05;

int gx=0, gy=0;

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int dgi : diagonals.keys)
	{
		string f = topDir + datasets[dsi] + "/eff3outof4_"+diagonals[dgi]+".root";
		string fmc = topDir + datasets[dsi] + "/eff3outof4_more_cuts_"+diagonals[dgi]+".root";

		real sgn = (diagonals[dgi] == "45b_56t") ? +1 : +1;
		string opt = "vl,eb";
	
		++gy; gx = 0;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad(false, gx, gy);
			label(RP_labels[rpi]);
		}
		
		NewPad(false, -1, gy);
		label(replace("\vbox{\SetFontSizesXX\hbox{dataset: "+datasets[dsi]+"}\hbox{diagonal: "+diagonals[dgi]+"}}", "_", "\_"));

		frame fLegend;
		
		// th_y dependence
		xTicksDef = LeftTicks(20., 10);
		yTicksDef = RightTicks(1., 0.2);
		drawGridDef = true;

		++gy; gx = 0;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad("$\th_y^*\ung{\mu rad}$", "$N_{\rm ev}(\hbox{4 tracks}) / N_{\rm ev}(\hbox{3 tracks})\ung{\%}$", gx, gy);
			draw(scale(sgn*1e6, 100), RootGetObject(f, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_y dependence/h_refined_ratio.th_y"), opt, black, "only $\th_y^*$ L-R cut, $3\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 1.0/th_y dependence/h_refined_ratio"), opt, red, "all cuts, $1\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 5.0/th_y dependence/h_refined_ratio"), opt, heavygreen, "all cuts, $5\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_y dependence/h_refined_ratio"), opt, blue+1pt, "all cuts, $3\un{\si}$");

			limits((20, 95), (120, 100), Crop);
			fLegend = BuildLegend();
		}
		
		// th_x dependence
		xTicksDef = LeftTicks(100., 50.);
		yTicksDef = RightTicks(2., 1.);
		drawGridDef = true;

		++gy; gx = 0;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad("$\th_x^*\ung{\mu rad}$", "\ung{\%}", gx, gy);
			draw(scale(sgn*1e6, 100), RootGetObject(f, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_x dependence/h_refined_ratio.th_x"), opt, black, "only $\th_x^*$ L-R cut, $3\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 1.0/th_x dependence/h_refined_ratio"), opt, red, "all cuts, $1\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 5.0/th_x dependence/h_refined_ratio"), opt, heavygreen, "all cuts, $5\un{\si}$");
			draw(scale(sgn*1e6, 100), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_x dependence/h_refined_ratio"), opt, blue+1pt, "all cuts, $3\un{\si}$");

			limits((-320, 90), (320, 100), Crop);
			fLegend = BuildLegend();
		}

		NewPad(false, -1, gy);
		add(fLegend);

		// th_x, th_y dependence
		drawGridDef = false;
		xTicksDef = LeftTicks(100., 50);
		yTicksDef = RightTicks(20., 10.);
		TH2_paletteTicks = PaletteTicks(0.01, 0.002);

		++gy; gx = 0;
		TH2_y_min = 10e-6;
		TH2_y_max = 110e-6;
		TH2_z_min = 0.96;
		TH2_z_max = 1;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$", gx, gy);
			draw(scale(1e6, 1e6), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_x, th_y dependence/h_simple_ratio"),
				"p,bar", blue, "all cuts, $3\un{\si}$");

			limits((-320, 20), (+320, 120), Crop);
			fLegend = BuildLegend();
		}

		// th_x, th_y dependence, large z scale
		drawGridDef = false;
		xTicksDef = LeftTicks(100., 50);
		yTicksDef = RightTicks(20., 10.);
		TH2_paletteTicks = PaletteTicks(0.05, 0.01);

		++gy; gx = 0;
		TH2_y_min = 10e-6;
		TH2_y_max = 110e-6;
		TH2_z_min = 0.7;
		TH2_z_max = 1;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$", gx, gy);
			draw(scale(1e6, 1e6), RootGetObject(fmc, "excluded RPs "+RPs[rpi]+"/n_si 3.0/th_x, th_y dependence, coarse/hc_simple_ratio"),
				"p,bar", blue, "all cuts, $3\un{\si}$");

			limits((-320, 20), (+320, 120), Crop);
			fLegend = BuildLegend();
		}

	}
}

GShipout(vSkip=0pt);
