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

string signatures[] = {
	"N_#.prot & F_#.prot & N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot",
	"N_#.prot & F_#.prot & !N_#.tr_val & F_#.tr_val OVER N_#.prot & N_F.prot",
	"N_#.prot & F_#.prot & !N_#.tr_val & !F_#.tr_val OVER N_#.prot & N_F.prot"
};

string signature_labels[] = {
	"3/4 far",
	"3/4 near",
	"2/4"
};

xSizeDef = 10cm;

TH2_z_min = 0;
TH2_z_max = 0.04;

//----------------------------------------------------------------------------------------------------
	
NewPad(false);
for (int ai : arms.keys)
{
	for (int pi : pairs.keys)
	{
		NewPad(false);
		label(arms[ai]+", " + pairs[pi]);
	}
}

for (int si : signatures.keys)
{
	NewRow();

	NewPad(false);
	label(signature_labels[si]);

	for (int ai : arms.keys)
	{
		for (int pi : pairs.keys)
		{
			NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$", axesAbove=true);

			RootObject obj = RootGetObject(f, arms[ai] + "/" + replace(signatures[si], "#", pairs[pi]) + "/h_th_x_th_y", error=false);
			if (obj.valid)
				draw(scale(1e6, 1e6), obj, "p,bar");
		}
	}
}

GShipout();
