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
		NewPad("$\th_x^*\ung{\mu rad}$");
		scale(Linear, Log);

		string f = topDir+dataset+"/distributions_" + diagonals[di] + ".root";

		RootGetObject(f, "selected - angles/h_th_x" + suffixes[si]);
		robj.vExec("Rebin", 2);

		draw(scale(+1e6, 1), robj, "cl,eb", blue, "pos");
		draw(scale(-1e6, 1), robj, "cl,eb", red, "neg");

		limits((150, 1), (300, 100), Crop);

		AttachLegend(replace(diagonals[di] + ", " + suffixLabels[si], "_", "\_"));
	}
}
