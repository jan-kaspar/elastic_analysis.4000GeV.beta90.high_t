import root;
import pad_layout;

string datasets[] = { "DS4" };

string diagonals[] = { "45b_56t", "45t_56b", "anti_45b_56b", "anti_45t_56t" };
string dgn_labels[] = { "45 bot -- 56 top", "45 top -- 56 bot", "45 bot -- 56 bot", "45 top -- 56 top" };

string topDir = "../../";

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string dataset;

void MakeComparison(string cut_str, string quantity, real xscale, string unit, string obj, real xlimit, real sigma,
	real xStep, real xstep)
{
	NewPad(quantity+"$\ung{"+unit+"}$", "", xTicks=LeftTicks(xStep, xstep));
	scale(Linear, Log(true));
	
	string dir = topDir+dataset+"/background_study/"+cut_str;

	for (int dgni : diagonals.keys)
		draw(scale(xscale, 1), RootGetObject(dir + "/distributions_"+diagonals[dgni]+".root", "elastic cuts/"+obj),
			"vl", StdPen(dgni+1), dgn_labels[dgni]);
	
	// TODO
	//real n_si = 4;
	real n_si = 5;
	yaxis(XEquals(-n_si * sigma, false), dashed);
	yaxis(XEquals(+n_si * sigma, false), dashed);

	xlimits(-xlimit, +xlimit, Crop);
	AttachLegend(quantity, NW, NE);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewRow();

	dataset = datasets[dsi];

	NewPad(false);
	label(replace("{\SetFontSizesXX\vbox{\hbox{"+dataset+"}}}", "_", "\_"));

	MakeComparison("cuts:2,7,5,6", "discriminator 1: $\De\th_x^*$", 1e6, "\mu rad", "cut 1/h_cq1", 200, 9.4, 100, 20);
	
	MakeComparison("cuts:1,7,5,6", "discriminator 2: $\De\th_y^*$", 1e6, "\mu rad", "cut 2/h_cq2", 50, 3.3, 50, 10);
		
	MakeComparison("cuts:1,2,5,6", "discriminator 7: $\De x^*$", 1e3, "\mu m", "cut 7/h_cq7", 200, 8.7, 100, 20);
	
	MakeComparison("cuts:1,2,7,6", "discriminator 5: $\De^{F-N} y^{R}\hbox{ vs. }y^{RN}$", 1e3, "\mu m", "cut 5/h_cq5", 300, 18., 100, 20);
		
	MakeComparison("cuts:1,2,7,5", "discriminator 6: $\De^{F-N} y^{L}\hbox{ vs. }y^{LN}$", 1e3, "\mu m", "cut 6/h_cq6", 300, 18., 100, 20);
}
