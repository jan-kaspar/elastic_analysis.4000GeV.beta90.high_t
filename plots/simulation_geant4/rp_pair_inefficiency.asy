import root;
import pad_layout;

string topDir = "../../simulation_geant4/";

//string f = topDir + "analysis_full.root";
//string f = topDir + "analysis.root";
string f = topDir + "analysis_elastic.root";

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

void DrawOneSet(string desc, pen p, string label)
{
	rObject obj = rGetObj(f, desc, error=false);
	if (obj.valid)
	{
		draw(scale(1e6, 1e2), obj, "eb", p, label);
		draw(scale(1e6, 1e2), rGetObj(f, desc + "|pol1"), p+1pt);
	}
}

//----------------------------------------------------------------------------------------------------

frame f_legend;

for (int ai : arms.keys)
{
	for (int pi : pairs.keys)
	{
		NewPad("$\th_y^*\ung{\mu rad}$", "inefficiency\ung{\%}", ai, pi);

		AddToLegend("<track reconstructed in:");
		
		string obj = replace("N_#.prot & F_#.prot & N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot/h_th_y", "#", pairs[pi]);
		DrawOneSet(arms[ai] + "/" + obj, black, "near \& not far ($\equiv$ 3/4 ineff.~far)");

		string obj = replace("N_#.prot & F_#.prot & !N_#.tr_val & F_#.tr_val OVER N_#.prot & N_F.prot/h_th_y", "#", pairs[pi]);
		DrawOneSet(arms[ai] + "/" + obj, red, "not near \& far ($\equiv$ 3/4 ineff.~near)");

		string obj = replace("N_#.prot & F_#.prot & !N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot/h_th_y", "#", pairs[pi]);
		DrawOneSet(arms[ai] + "/" + obj, blue, "not near \& not far ($\equiv$ 2/4 ineff.)");

		real x_min, x_max;
		if ((ai + pi) % 2 == 1)
		{
			x_min = 20;
			x_max = 110;
		} else	{
			x_min = -110;
			x_max = -20;
		}

		limits((x_min, 0), (x_max, 3.5), Crop);
		
		f_legend = BuildLegend();
		currentpicture.legend.delete();

		AttachLegend(arms[ai] + ", " + pairs[pi]);

	}
}

NewPad(false, 2, 0);
add(f_legend);
