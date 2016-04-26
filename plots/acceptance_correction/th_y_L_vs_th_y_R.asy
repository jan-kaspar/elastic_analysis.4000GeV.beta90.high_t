import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS4" };
real z_scale[] = { 2e4 };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };
real dgn_scale[] = { +1, -1 };

TH2_palette = Gradient(blue, cyan, heavygreen, yellow, red, magenta);
//TH2_palette = new pen[] { paleblue, blue, cyan, heavygreen, yellow, orange, red, magenta, black };

xTicksDef = LeftTicks(5., 1.);
yTicksDef = RightTicks(5., 1.);

for (int di : datasets.keys)
{
	write(">> " + datasets[di]);

	NewPad(false);
	label("{\SetFontSizesXX "+datasets[di]+"}");

	for (int dgni : diagonals.keys)
	{
		real s = dgn_scale[dgni];

		TH2_z_max = log10(z_scale[di]);
	
		NewPad("$\th_y^{*R}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$");
		scale(Linear, Linear, Log);
		TH2_x_min = TH2_y_min = min(s*15e-6, s*45e-6);
		TH2_x_max = TH2_y_max = max(s*15e-6, s*45e-6);
		draw(scale(s*1e6, s*1e6), RootGetObject(topDir+datasets[di]+"/distributions_"+diagonals[dgni]+".root", "selected - angles/h_th_y_L_vs_th_y_R"), "p,bar");
		limits((25, 25), (45, 45), Crop);

		for (real x = 26; x <= 45; x += 2)
			yaxis(XEquals(x, false), dotted, above=true);

		for (real y = 26; y <= 45; y += 2)
			xaxis(YEquals(y, false), dotted, above=true);
		
		AttachLegend(dgn_labels[dgni], NW, NW);
	}

	NewRow();
}

write(">> done");
