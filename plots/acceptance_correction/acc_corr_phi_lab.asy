import root;
import pad_layout;

string top_dir = "../../";

string dataSets[] = { "DS4" };

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

real cut_th_x_low_top = -1000, cut_th_x_high_top = 1000;
real cut_th_y_low_top = 33.8, cut_th_y_high_top = +100;

real cut_th_x_low_bot = -1000, cut_th_x_high_bot = 1000;
real cut_th_y_low_bot = -33.8, cut_th_y_high_bot = -105;

void DrawAcceptedArcs(real th)
{
	real dphi = 360. / 1000;
	bool inside = false;
	real phi_start;
	for (real phi = 0; phi < 360; phi += dphi)
	{
		real x = th * Cos(phi);
		real y = th * Sin(phi);

		bool p_in = false;
		if (x > cut_th_x_low_top && x < cut_th_x_high_top && y > cut_th_y_low_top && y < cut_th_y_high_top)
			p_in = true;
		if (x > cut_th_x_low_bot && x < cut_th_x_high_bot && y < cut_th_y_low_bot && y > cut_th_y_high_bot)
			p_in = true;

		if (!inside && p_in)
		{
			phi_start = phi;
			inside = true;
		}

		if (inside && !p_in)
		{
			inside = false;
			real phi_end = phi - dphi;
			draw(arc((0, 0), th, phi_start, phi_end), black+1pt);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void DrawFullArc(real th)
{
	draw(scale(th)*unitcircle, dotted);
	label(rotate(-90)*Label(format("\SmallerFonts $%.0f$", th)), (th, 0), 0.5E, Fill(white+opacity(0.8)));
}

//----------------------------------------------------------------------------------------------------

real arc_th_y[] = { 50, 100, 150, 200 };
//real arc_th_y[] = { 177, 250, 345 };

real th_x_max = 350;

for (int dsi : dataSets.keys)
{
	real ySize = 6cm;

	NewPad("$\th_x^{*}\ung{\mu rad}$", "$\th_y^{*}\ung{\mu rad}$", ySize/200*th_x_max, ySize);
	scale(Linear, Linear, Log);
	TH2_zLabel = "(corrected) events per bin";
	TH2_paletteBarWidth = 0.05;
	
	label("$\th^*\!=$", (50, 0), 0.5W, Fill(white+opacity(0.8)));
	for (real th_y : arc_th_y)
		DrawFullArc(th_y);
	label(rotate(-90)*Label("\SmallerFonts $\rm\mu rad$"), (215, 0), 0.5E, Fill(white+opacity(0.8)));

	// z scale
	//TH2_z_min = 5.5;
	//TH2_z_max = 3.75;

	// 45 bottom - 56 top
	draw(scale(1e6, 1e6), RootGetObject(top_dir+"/"+dataSets[dsi]+"/distributions_45b_56t.root", "normalization/h_th_y_vs_th_x_normalized"), "def");
	
	// 45 top - 56 bottom
	draw(scale(1e6, 1e6), RootGetObject(top_dir+"/"+dataSets[dsi]+"/distributions_45t_56b.root", "normalization/h_th_y_vs_th_x_normalized"), "p");
	
	draw((-th_x_max, cut_th_y_low_top)--(+th_x_max, cut_th_y_low_top), magenta+1pt);
	draw((-th_x_max, cut_th_y_low_bot)--(+th_x_max, cut_th_y_low_bot), magenta+1pt);
	
	draw((-th_x_max, cut_th_y_high_top)--(+th_x_max, cut_th_y_high_top), red+1pt);
	draw((-th_x_max, cut_th_y_high_bot)--(+th_x_max, cut_th_y_high_bot), red+1pt);

	/*
	draw((cut_th_x_low_top , 0)--(cut_th_x_low_top , +200), magenta+1pt);
	draw((cut_th_x_low_bot , -200)--(cut_th_x_low_bot , 0), magenta+1pt);
	draw((cut_th_x_high_top, 0)--(cut_th_x_high_top, +200), magenta+1pt);
	draw((cut_th_x_high_bot, -200)--(cut_th_x_high_bot, 0), magenta+1pt);
	*/
	
	for (real th_y : arc_th_y)
		DrawAcceptedArcs(th_y);

	label("\vbox{\hbox{detector}\hbox{edges}}", (-340, -130), SE, magenta, Fill(white));
	draw((-300, -135)--(-310, cut_th_y_low_bot), magenta, EndArrow());
	draw((-300, -135)--(-290, cut_th_y_low_top), magenta, EndArrow());

	label("\vbox{\hbox{LHC}\hbox{appertures}}", (-240, 190), S, red, Fill(white));
	draw((-230, 130)--(-210, cut_th_y_high_top), red, EndArrow);
	draw((-230, 130)--(-240, cut_th_y_high_bot), red, EndArrow);

	/*
	label("\vbox{\hbox{horiz.}\hbox{RPs}}", (200, -150), W, magenta, Fill(white));
	draw((130, -150)--(cut_th_x_high_bot, -140), magenta, EndArrow);
	draw((130, -150)--(cut_th_x_low_bot, -160), magenta, EndArrow);
	*/
	
	limits((-th_x_max, -200), (th_x_max, 200), Crop);
	AttachLegend(dataSets[dsi]);
}

GShipout(margin=0mm);
