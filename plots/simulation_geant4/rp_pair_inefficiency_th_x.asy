import root;
import pad_layout;

string topDir = "../../simulation_geant4/";

//string f = topDir + "analysis_full.root";
string f = topDir + "analysis.root";

string arms[] = {
	"L",
	"R",
};

string pairs[] = {
	"T",
	"B",
};

xSizeDef = 10cm;

drawGridDef = true;

//----------------------------------------------------------------------------------------------------

frame f_legend;

for (int ai : arms.keys)
{
	for (int pi : pairs.keys)
	{
		NewPad("$\th_x^*\ung{\mu rad}$", "inefficiency\ung{\%}", ai, pi);

		AddToLegend("<track reconstructed in:");
		
		string obj_n = replace("N_#.prot & F_#.prot & N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot/h_th_x", "#", pairs[pi]);
		rObject obj = rGetObj(f, arms[ai] + "/" + obj_n, error=false);
		if (obj.valid)
			draw(scale(1e6, 1e2), obj, "eb", black, "near \& not far ($\equiv$ 3/4 ineff.~far)");

		string obj_n = replace("N_#.prot & F_#.prot & !N_#.tr_val & F_#.tr_val OVER N_#.prot & N_F.prot/h_th_x", "#", pairs[pi]);
		rObject obj = rGetObj(f, arms[ai] + "/" + obj_n, error=false);
		if (obj.valid)
			draw(scale(1e6, 1e2), obj, "eb", red, "not near \& far ($\equiv$ 3/4 ineff.~near)");

		string obj_n = replace("N_#.prot & F_#.prot & !N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot/h_th_x", "#", pairs[pi]);
		rObject obj = rGetObj(f, arms[ai] + "/" + obj_n, error=false);
		if (obj.valid)
			draw(scale(1e6, 1e2), obj, "eb", blue, "not near \& not far ($\equiv$ 2/4 ineff.)");

		limits((-350, 0), (+350, 3.5), Crop);
		
		if (obj.valid)
			f_legend = BuildLegend();
		currentpicture.legend.delete();

		AttachLegend(arms[ai] + ", " + pairs[pi]);
	}
}

NewPad(false, 2, 0);
add(f_legend);
