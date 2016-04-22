import root;
import pad_layout;

string topDir = "../";

string dataset = "DS4";

string binnings[] = {
	"ob-1-30-0.10",
	"ob-2-20-0.20",
	"ob-3-10-0.30"
};

string n_si[] = {
	"4",
	"5",
	"6",
};

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b--56t", "45t--56b" };

xSizeDef = 10cm;
ySizeDef = 8cm;

drawGridDef = true;

xTicksDef = LeftTicks(0.5, 0.1);

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int ni : n_si.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX $"+n_si[ni]+"\un{\si}$}");
}

for (int bi : binnings.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX "+binnings[bi]+"}");

	for (int ni : n_si.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
		scale(Linear, Log);

		for (int dgni : diagonals.keys)
		{
			rObject obj = rGetObj(topDir+dataset+"/background_study/n_si="+n_si[ni]
				+"/distributions_"+diagonals[dgni]+".root", "normalization/"+binnings[bi]+"/h_t_normalized");
			draw(obj, "vl,eb", StdPen(dgni+1), dgn_labels[dgni]);
		}
	
		limits((0, 1e-4), (2, 1e3), Crop);
		yaxis(XEquals(1.9, false), dashed);
	
		AttachLegend("before unfolding");
	}
}

