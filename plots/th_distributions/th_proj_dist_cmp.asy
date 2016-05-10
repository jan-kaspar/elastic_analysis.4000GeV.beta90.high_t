import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string quantities[] = { "th_x", "th_y", "th" };

string suffixes[] = { "_R" };
string suffixLabels[] = { "L", "R", "L+R" };

string diagonals[] = { "45b_56t", "45t_56b" };

xSizeDef = 12cm;

//----------------------------------------------------------------------------------------------------

for (int qi : quantities.keys)
{
	NewRow();

	for (int di : diagonals.keys)
	{
		NewPad("$\ung{\mu rad}$");
		scale(Linear, Log);

		for (int si : suffixes.keys)
		{
			pen p = StdPen(si);
			string f = topDir+dataset+"/distributions_" + diagonals[di] + ".root";
			RootGetObject(f, "selected - angles/h_" + quantities[qi] + suffixes[si]);
			robj.vExec("Rebin", 5);
			draw(scale(1e6, 1), robj, "eb", p, suffixLabels[si]);
		}

		AttachLegend(replace(quantities[qi] + ", " + diagonals[di], "_", "\_"));
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
