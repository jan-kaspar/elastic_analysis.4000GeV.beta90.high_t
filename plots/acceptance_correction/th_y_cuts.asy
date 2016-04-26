import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS4" };
real y_lim[] = { 3e4 };
real y_edge[] = { 30 };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };
real x_scale[] = { +1e6, -1e6 };


for (int dsi : datasets.keys)
{
	for (int dgni : diagonals.keys)
	{
		string f = topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root";

		// get th_y cuts
		RootObject obj = RootGetObject(f, "selected - angles/g_th_y_cuts");
		real x[] = {0};
		real y[] = {0};
		obj.vExec("GetPoint", 0, x, y); real th_y_lcut_L = y[0]*1e6; 
		obj.vExec("GetPoint", 1, x, y); real th_y_hcut_L = y[0]*1e6; 
		obj.vExec("GetPoint", 2, x, y); real th_y_lcut_R = y[0]*1e6; 
		obj.vExec("GetPoint", 3, x, y); real th_y_hcut_R = y[0]*1e6; 
	
		NewPad(false);
		label(replace("\vbox{\SetFontSizesXX\hbox{dataset: "+datasets[dsi]+"}\hbox{diagonal: "+dgn_labels[dgni]+"}}", "_", "\_"));

		NewPad("$\th_y^*\ung{\mu rad}$");
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_L"), "vl,eb", blue, "left arm");
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_R"), "vl,eb", red, "right arm");
		limits((0, 0), (120, y_lim[dsi]), Crop);
		AttachLegend();

		NewPad("$\th_y^*\ung{\mu rad}$", xTicks=LeftTicks(2., 1.));
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_L"), "vl,eb", blue, "left arm");
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_R"), "vl,eb", red, "right arm");
		limits((y_edge[dsi]-2, 0), (y_edge[dsi]+13, y_lim[dsi]), Crop);
		yaxis(XEquals(th_y_lcut_L, false), dashed+blue);
		yaxis(XEquals(th_y_lcut_R, false), dashed+red);
		AttachLegend(SE, SE);

		NewPad("$\th_y^*\ung{\mu rad}$", xTicks=LeftTicks(5., 1.));
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_L"), "vl,eb", blue, "left arm");
		draw(scale(x_scale[dgni], 1), RootGetObject(f, "selected - angles/h_th_y_R"), "vl,eb", red, "right arm");
		limits((95, 0), (115, y_lim[dsi]/12), Crop);
		yaxis(XEquals(th_y_hcut_L, false), dashed+blue);
		yaxis(XEquals(th_y_hcut_R, false), dashed+red);
		AttachLegend(NE, NE);

		NewRow();
	}
}
