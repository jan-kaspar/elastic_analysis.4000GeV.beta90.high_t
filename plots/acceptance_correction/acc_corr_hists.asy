import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS4" };
real t_min[] = { 0.027, 0.02 };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };

string binning = "ob-1-30-0.10";

for (int dsi : datasets.keys)
{
	string dataset = datasets[dsi];

	for (int dgi : dgns.keys)
	{
		string dgn = dgns[dgi];
		string f = topDir+dataset+"/distributions_" + dgn + ".root";
		
		NewRow();

		drawGridDef = true;
		NewPad(false);
		label("\vbox{\SetFontSizesXX\hbox{"+dataset+"}\hbox{"+dgn_labels[dgi]+"}}");

		NewPad("$|t_y|\ung{GeV^2}$", "correction factor", xTicks=LeftTicks(0.1, 0.05), yTicks=RightTicks(0.2, 0.1));
		scale(Linear, Linear);
		TH1_x_min = 3e-4;
		draw(RootGetObject(f, "acceptance correction/p_t_ub_div_corr"), "eb", blue+1pt, "divergence");
		limits((0, 0.9), (0.2, 2), Crop);
		AttachLegend(NE, NE);
		
		NewPad("$|t|\ung{GeV^2}$", "\vbox{\hbox{correction factor}\hbox{(mean $\pm$ std.~dev.~per bin)}}", xSize=10cm, xTicks=LeftTicks(0.5, 0.1));
		scale(Linear, Log);
		draw(RootGetObject(f, "acceptance correction/" + binning + "/p_t_phi_corr"), "d0,eb", heavygreen, "phi");
		draw(RootGetObject(f, "acceptance correction/" + binning + "/p_t_full_corr"), "d0,eb", red, "full = divergence * phi");
		limits((0, 1e0), (1.9, 1e2), Crop);
		AttachLegend(NW, NW);
		
		drawGridDef = false;
		NewPad("$|t|\ung{GeV^2}$", "\vbox{\hbox{correction factor}\hbox{(mean $\pm$ std.~dev.~per bin)}}", xTicks=LeftTicks(0.02, 0.005), yTicks=RightTicks(1., 0.5));
		scale(Linear, Linear);
		draw(RootGetObject(f, "acceptance correction/" + binning + "/p_t_phi_corr"), "d0,eb", heavygreen, "phi");
		draw(RootGetObject(f, "acceptance correction/" + binning + "/p_t_full_corr"), "d0,eb", red, "full");
		limits((0, 2), (0.1, 10), Crop);
		AttachLegend(NE, NE);
		xaxis(YEquals(5, false), dotted);
		yaxis(XEquals(t_min[dsi], false), dotted);
		
		/*
		NewPad("$|t|\ung{GeV^2}$", "acceptance");
		//scale(Linear, Log);
		draw(RootGetObject(f, "acceptance correction/eb/p_t_full_acc"), "d0,eb,vl", magenta, "acceptance");
		limits((0, 0), (0.25, 1), Crop);
		AttachLegend(NE, NE);
		*/
	}
}
