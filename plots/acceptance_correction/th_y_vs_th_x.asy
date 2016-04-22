import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS4" };
real z_scale[] = { 1e3 };

string diagonals[] = { "45b_56t", "45t_56b" };

TH2_palette = new pen[] { paleblue, blue, cyan, heavygreen, yellow, orange, red, magenta, black };

//----------------------------------------------------------------------------------------------------

for (int di : datasets.keys)
{
	write(">> " + datasets[di]);

	NewPad(false);
	label("{\SetFontSizesXX "+datasets[di]+"}");

	NewRow();
	TH2_z_max = log10(z_scale[di]);

	NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$");
	scale(Linear, Linear, Log);
	TH2_x_min = -150e-6; TH2_x_max = +150e-6;
	TH2_y_min = -50e-6; TH2_y_max = -20e-6;
	draw(scale(1e6, -1e6), rGetObj(topDir+datasets[di]+"/distributions_45t_56b.root", "selected - angles/h_th_y_L_vs_th_x_L"), "p,bar");
	limits((-150, 20), (150, 50), Crop);
	AttachLegend("left top", N, N);

	NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$");
	scale(Linear, Linear, Log);
	TH2_x_min = -150e-6; TH2_x_max = +150e-6;
	TH2_y_min = 20e-6; TH2_y_max = 50e-6;
	draw(scale(1e6, +1e6), rGetObj(topDir+datasets[di]+"/distributions_45b_56t.root", "selected - angles/h_th_y_R_vs_th_x_R"), "p,bar");
	limits((-150, 20), (150, 50), Crop);
	AttachLegend("right top", N, N);
	
	NewRow();

	NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$");
	scale(Linear, Linear, Log);
	TH2_x_min = -150e-6; TH2_x_max = +150e-6;
	TH2_y_min = 20e-6; TH2_y_max = 50e-6;
	draw(scale(1e6, -1e6), rGetObj(topDir+datasets[di]+"/distributions_45b_56t.root", "selected - angles/h_th_y_L_vs_th_x_L"), "p,bar");
	limits((-150, -50), (150, -20), Crop);
	AttachLegend("left bottom", N, N);

	NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$");
	scale(Linear, Linear, Log);
	TH2_x_min = -150e-6; TH2_x_max = +150e-6;
	TH2_y_min = -50e-6; TH2_y_max = -20e-6;
	draw(scale(1e6, +1e6), rGetObj(topDir+datasets[di]+"/distributions_45t_56b.root", "selected - angles/h_th_y_R_vs_th_x_R"), "p,bar");
	limits((-150, -50), (150, -20), Crop);
	AttachLegend("right bottom", N, N);

	NewPage();
}
