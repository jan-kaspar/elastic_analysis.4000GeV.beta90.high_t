import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string file_45b = topDir+"DS4/distributions_45b_56t.root";
string file_45t = topDir+"DS4/distributions_45t_56b.root";

string rps[] = { "L_F", "L_N", "R_N", "R_F" };
string rp_labels[] = { "left far", "left near", "right near", "right far" };
real sh_top[] = { 9.5, 9.7, 8.5, 8.1 };
real sh_bot[] = { -7.7, -7.95, -9.7, -9.5 };
real sh_hor[] = { 3.6, 5.1, 7.8, 5.2 };

TH2_x_min = -11; TH2_x_max = +11;
TH2_y_min = -30; TH2_y_max = +30;

//----------------------------------------------------------------------------------------------------

real edge = 3.6101;
real cutEdge = 2.22721 / sqrt(2);
int strips = 11;
real margin_v_e = 0.2;
real margin_v_b = 0.4;
real margin_u = 0.1;

path det_shape = (cutEdge, 0)--(edge, 0)--(edge, edge)--(0, edge)--(0, cutEdge)--cycle;
det_shape = scale(10) * rotate(45) * det_shape;

path hor_det_shape = shift(0, -cutEdge/sqrt(2)*10) * det_shape;

//----------------------------------------------------------------------------------------------------

//string rps[] = { "L_F" };
//string rp_labels[] = { "left near", "left near", "bla" };

xTicksDef = LeftTicks(Step=10, step=5);
yTicksDef = RightTicks(Step=10, step=5);
for (int ri : rps.keys)
{
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 10cm);
	scale(Linear, Linear, Log);

	//TH2_z_max = log10(200);
	
	draw(rGetObj(file_45b, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "o", lightblue);
	draw(rGetObj(file_45t, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "o", lightblue);
	
	draw(rGetObj(file_45b, "hit distributions/horizontal, aligned, after selection/h_y_"+rps[ri]+"H_vs_x_"+rps[ri]+"H_al_sel"), "o", red);
	draw(rGetObj(file_45t, "hit distributions/horizontal, aligned, after selection/h_y_"+rps[ri]+"H_vs_x_"+rps[ri]+"H_al_sel"), "o", heavygreen);

	draw(shift(0, -sh_top[ri])*det_shape);
	draw(shift(0, -sh_bot[ri])*scale(1, -1)*det_shape);
	draw(shift(sh_hor[ri], 0)*rotate(-90)*hor_det_shape);

	//limits((0, -20), (+10, +20), Crop);
	limits((-30, -50), (+30, +50), Crop);

	for (real x = -30; x <= +30; x += 5)
		yaxis(XEquals(x, false), dotted+gray);

	for (real y = -50; y <= +50; y += 5)
		xaxis(YEquals(y, false), dotted+gray);

	AttachLegend(replace(rp_labels[ri], "_", "\_"));
}

GShipout("hit_distributions_with_horizontals");

//----------------------------------------------------------------------------------------------------

xTicksDef = LeftTicks(Step=10, step=2);
yTicksDef = RightTicks(Step=10, step=2);

//string rps[] = { "L_F", "L_N", "R_N" };
//string rp_labels[] = { "left near", "left near", "bla" };

for (int ri : rps.keys)
{
	write(rps[ri]);

	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 10cm);
	scale(Linear, Linear, Log);

	//TH2_z_max = log10(200);
	
	draw(rGetObj(file_45b, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p,bar");
	draw(rGetObj(file_45t, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p");
	
	draw(shift(0, -sh_top[ri])*det_shape);
	draw(shift(0, -sh_bot[ri])*scale(1, -1)*det_shape);
	draw(shift(sh_hor[ri], 0)*rotate(-90)*hor_det_shape);

	//limits((0, -20), (+10, +20), Crop);
	limits((-30, -50), (+30, +50), Crop);

	AttachLegend(replace(rp_labels[ri], "_", "\_") + " : before");
}

NewRow();

for (int ri : rps.keys)
{
	write(rps[ri]);

	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 10cm);
	scale(Linear, Linear, Log);

	//TH2_z_max = log10(200);
	
	draw(rGetObj(file_45b, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "p,bar");
	draw(rGetObj(file_45t, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "p");
	
	draw(shift(0, -sh_top[ri])*det_shape);
	draw(shift(0, -sh_bot[ri])*scale(1, -1)*det_shape);
	draw(shift(sh_hor[ri], 0)*rotate(-90)*hor_det_shape);

	//limits((0, -20), (+10, +20), Crop);
	limits((-30, -50), (+30, +50), Crop);

	AttachLegend(replace(rp_labels[ri], "_", "\_") + " : after");
}

GShipout("hit_distributions_before_after");
