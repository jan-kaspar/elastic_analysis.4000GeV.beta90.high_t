import root;
import pad_layout;

string topDir = "../../simulation_geant4/";

//string f = topDir + "analysis_full.root";
string f = topDir + "analysis.root";

string dataDir = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta90/DS4/";

string arms[] = {
	"L",
	"R",
};

string pairs[] = {
	"T",
	"B",
};

xSizeDef = 10cm;

//----------------------------------------------------------------------------------------------------

frame f_legend;

for (int ai : arms.keys)
{
	for (int pi : pairs.keys)
	{
		NewPad("$\th_y^*\ung{\mu rad}$", "inefficiency\ung{\%}", ai, pi);

		// Monte-Carlo
		AddToLegend("<Monte-Carlo:");

		string obj = replace("N_#.prot & F_#.prot & N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot/h_th_y", "#", pairs[pi]);
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj), "eb", blue);
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj + "|pol1"), blue+dashed, "3/4 ineff.~far");

		string obj = replace("N_#.prot & F_#.prot & !N_#.tr_val & F_#.tr_val OVER N_#.prot & N_F.prot/h_th_y", "#", pairs[pi]);
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj), "eb", red);
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj + "|pol1"), red+dashed, "3/4 ineff.~near");

		/*
		string obj = replace("N_#.prot & F_#.prot & !N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot", "#", pairs[pi]);
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj), "eb", blue, "not near \& not far ($\equiv$ 2/4 ineff.)");
		draw(scale(1e6, 1e2), rGetObj(f, arms[ai] + "/" + obj + "|pol1"), blue+1pt);
		*/

		// horizontal scale
		real x_min, x_max;
		if ((ai + pi) % 2 == 1)
		{
			x_min = 30;
			x_max = 100;
		} else	{
			x_min = -100;
			x_max = -30;
		}

		// data fits
		AddToLegend("<data fits:");

		string fp = arms[ai] + "_" + pairs[pi];
		string diagonal = (fp == "L_B" || fp == "R_T") ? "45b_56t" : "45t_56b";
		string f = dataDir + "eff3outof4_details_fit.root";

		rObject fit = rGetObj(f, diagonal + "/" + arms[ai] + "_N/th_y : rel|ff");
		real y_min = fit.rExec("Eval", x_min);
		real y_max = fit.rExec("Eval", x_max);
		draw((x_min, 100*(1 - y_min))--(x_max, 100*(1 - y_max)), red+1pt, "3/4 ineff.~near");

		rObject fit = rGetObj(f, diagonal + "/" + arms[ai] + "_F/th_y : rel|ff");
		real y_min = fit.rExec("Eval", x_min);
		real y_max = fit.rExec("Eval", x_max);
		draw((x_min, 100*(1 - y_min))--(x_max, 100*(1 - y_max)), blue+1pt, "3/4 ineff.~far");

		limits((x_min, 0), (x_max, 3.5), Crop);
		
		f_legend = BuildLegend();
		currentpicture.legend.delete();

		AttachLegend(arms[ai] + ", " + pairs[pi]);

	}
}

NewPad(false, 2, 0);
add(f_legend);

//----------------------------------------------------------------------------------------------------


/*
// x scale
real x_min, x_max;
x_min = 30;
x_max = 100;

// data
string f = dataDir + "eff3outof4_details_fit.root";

AddToLegend("<data:");
AddToLegend("3/4 ineff.~far", heavygreen);
AddToLegend("3/4 ineff.~near", magenta);

for (int dgi : diagonals.keys)
{
	for (int rpi : RPs.keys)
	{
		string d = diagonals[dgi] + "/" + RPs[rpi];
		rObject fit = rGetObj(f, d+"/th_y : rel|ff");
		
		real y_min = fit.rExec("Eval", x_min);
		real y_max = fit.rExec("Eval", x_max);

		bool near = (find(RPs[rpi], "_N") != -1);

		draw((x_min, 100*(1 - y_min))--(x_max, 100*(1 - y_max)), (near) ? magenta : heavygreen);
	}
}

limits((x_min, 0), (x_max, 3.5), Crop);
*/
