import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string suffixes[] = { "_L", "_R" };
string suffixLabels[] = { "L", "R" };

string diagonals[] = { "45b_56t", "45t_56b" };

xSizeDef = 12cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

for (int di : diagonals.keys)
{
	NewRow();

	for (int si : suffixes.keys)
	{
		NewPad("$\ung{\mu rad}$");
		scale(Linear, Log);

		string f = topDir+dataset+"/distributions_" + diagonals[di] + ".root";

		RootGetObject(f, "selected - angles/h_th_x" + suffixes[si]);
		robj.vExec("Rebin", 2);

		draw(scale(+1e6, 1), robj, "eb", blue, "pos");
		draw(scale(-1e6, 1), robj, "eb", red, "neg");

		limits((150, 1), (300, 100), Crop);

		AttachLegend(replace(diagonals[di] + ", " + suffixLabels[si], "_", "\_"));
	}
}



/*
for (int bi : binnings.keys)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si_{\rm el}/\d t \ung{mb/GeV^2}$", 10cm, 8cm);
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	scale(Linear, Log);
	
	for (int di : datasets.keys)
	{
		draw(RootGetObject(topDir+"DS-merged/merged.root", binnings[bi]+"/"+datasets[di]+"/combined/h_dsdt" ),
			"eb", StdPen(di + 1), datasets[di]);
	}
	
	limits((0, 1e-4), (2., 1e3), Crop);

	AttachLegend(binnings[bi]);
}
*/
