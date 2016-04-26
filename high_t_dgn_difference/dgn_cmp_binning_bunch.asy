import root;
import pad_layout;

string datasets[] = { /*"DS4",*/ "DS4-r8369", "DS4-r8371", "DS4-r8372"};

string topDir = "../";

string binnings[] = {
	"ob-1-30-0.10",
	"ob-2-20-0.20",
	"ob-3-10-0.30"
};

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b--56t", "45t--56b" };

xSizeDef = 10cm;
ySizeDef = 8cm;

drawGridDef = true;

xTicksDef = LeftTicks(0.5, 0.1);

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int dsi : datasets.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX "+datasets[dsi]+"}");
}

for (int bi : binnings.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX "+binnings[bi]+"}");

	for (int dsi : datasets.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
		scale(Linear, Log);
		for (int dgni : diagonals.keys)
		{
			RootObject obj = RootGetObject(topDir+datasets[dsi]+"/distributions_"+diagonals[dgni]+".root", "normalization/"+binnings[bi]+"/h_t_normalized");
			draw(obj, "vl,eb", StdPen(dgni+1), dgn_labels[dgni]);
		}
	
		limits((0.4, 1e-4), (1.9, 1e-1), Crop);
		//yaxis(XEquals(1.9, false), dashed);
	
		AttachLegend("before unfolding");
	}
}

GShipout(margin=1mm);
