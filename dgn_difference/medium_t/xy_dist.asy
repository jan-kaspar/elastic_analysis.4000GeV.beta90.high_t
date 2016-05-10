import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string diagonal = "45b_56t";

string rps[] = { "L_N", "L_F", "R_N", "R_F" };

xSizeDef = 12cm;
ySizeDef = 8cm;

dotfactor = 5;

//----------------------------------------------------------------------------------------------------

for (int rpi : rps.keys)
{
	if (rpi == 2)
		NewRow();

	NewPad("$x\ung{mm}$", "$y\ung{mm}$");
	
	string f = topDir+dataset+"/distributions_" + diagonal + ".root";
	string dir = "hit distributions/vertical, aligned, after selection, t selection/";
	draw(RootGetObject(f, dir + "g_y_"+rps[rpi]+"_vs_x_"+rps[rpi]+"_al_sel_tsel"), "d", red);

	if (rpi < 2)
		limits((-1.4, -30), (+1.4, -5), Crop);
	else
		limits((-1.4, 5), (+1.4, 30), Crop);

	AttachLegend(replace(dataset + ", " + diagonal + ", " + rps[rpi], "_", "\_"));
}

GShipout(margin=1mm);
