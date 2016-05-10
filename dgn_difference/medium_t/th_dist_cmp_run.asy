import root;
import pad_layout;

string datasets[] = { "DS4-r8369", "DS4-r8371", "DS4-r8372" };
real dataset_sh[] = { 0, 1., -1.1 };

string topDir = "../../";

string suffixes[] = { "_L", "_R" };
string suffixLabels[] = { "Left", "Right" };

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

		for (int dsi : datasets.keys)
		{
			pen p = StdPen(dsi);

			string f = topDir+datasets[dsi]+"/distributions_" + diagonals[di] + ".root";
			RootGetObject(f, "selected - angles/h_th" + suffixes[si]);
			robj.vExec("Rebin", 2);
			draw(scale(+1e6, 1) * shift(0, dataset_sh[dsi]), robj, "eb", p, datasets[dsi]);
		}

		limits((150, 1e-1), (300, 1e3), Crop);

		yaxis(XEquals(181.7, false), heavygreen+1.5pt);
		yaxis(XEquals(205, false), magenta+1.5pt);
		yaxis(XEquals(226.2, false), orange+1.5pt);
		yaxis(XEquals(234.8, false), orange+1.5pt);

		AttachLegend(replace(dgn_labels[di] + ", " + suffixLabels[si], "_", "\_"));
	}
}
