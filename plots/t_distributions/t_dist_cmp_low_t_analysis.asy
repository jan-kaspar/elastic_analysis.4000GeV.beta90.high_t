import root;
import pad_layout;
include common_code;

A_ref = 519.5;
B_ref = 19.38;
ref_str = MakeRefStr("");

xSizeDef = 10cm;
ySizeDef = 8cm;

drawGridDef = true;

string dgns[] = { "45b_56t", "45t_56b", "combined" };
string dgn_labels[] = { "45b--56t", "45t--56b", "combined" };

string datasets[] = { "DS4-sc" };
//string datasets[] = { "DS4-sc", "DS4-b648", "DS4-b2990", "DS4-b26" };

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewRow();

	for (int dgni : dgns.keys)
	{
		RootObject o_this = RootGetObject("../../DS-merged/merged.root", "ub/"+datasets[dsi]+"/"+dgns[dgni]+"/h_dsdt");
		RootObject o_old = RootGetObject("../../../4000GeV,beta90/DS-merged/merged.root", "ub/DS4-sc/"+dgns[dgni]+"/h_dsdt");
	
		NewPad("$|t|\ung{GeV^2}$", "${\d\si/\d t - \hbox{ref}\over\hbox{ref}}\ ,\quad\hbox{ref} = "+ref_str+"$");
		DrawRelDiff(o_this, blue, "this analysis, " + datasets[dsi]);
		DrawRelDiff(o_old, red, "low-$|t|$ analysis, DS4-sc");
		
		limits((0, -0.05), (0.3, 0.06), Crop);
		
		currentpad.xTicks = LeftTicks(0.05, 0.01);
		currentpad.yTicks = RightTicks(0.01, 0.002);
		
		AttachLegend(dgn_labels[dgni], NW, NW);
	}
}
