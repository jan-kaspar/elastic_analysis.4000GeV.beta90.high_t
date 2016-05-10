import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string suffixes[] = { "_L", "_R" };
string suffixLabels[] = { "L", "R" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45 bot -- 56 top", "45 top -- 56 bot" };

xSizeDef = 12cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

for (int di : diagonals.keys)
{
	NewRow();

	for (int si : suffixes.keys)
	{
		NewPad("$\th^*\ung{\mu rad}$");
		scale(Linear, Log);

		string f = topDir+dataset+"/distributions_" + diagonals[di] + ".root";

		RootGetObject(f, "selected - angles/h_th" + suffixes[si]);
		robj.vExec("Rebin", 2);

		draw(scale(+1e6, 1), robj, "eb", red);

		limits((150, 1), (300, 100), Crop);

		AttachLegend(replace(dgn_labels[di] + ", " + suffixLabels[si], "_", "\_"));
	}
}
