import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

xSizeDef = 8cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

void MakeOne(string model, string fits[], string ranges[][])
{
	for (int fi : fits.keys)
	{
		NewRow();

		NewPadLabel(replace(fits[fi], "_", "\_"));

		for (int padi : ranges.keys)
		{
			NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");
			scale(Linear, Log);
			
			draw(RootGetObject(f, "h_dsdt"), "eb", black);

			for (int crvi : ranges[padi].keys)
			{
				string range = ranges[padi][crvi];

				RootObject fit = RootGetObject(f, model + "/fit range/" + range + "/" + fits[fi]);

				RootObject disabledObj = RootGetObject(f, model + "/fit range/" + range + "/disabled", error=false);
				

				string l = range;
				l += format(", $\ch^2/ndf = %#.2f", fit.rExec("GetChisquare")) + format(" / %i", fit.iExec("GetNDF"))
					+ format(" = %#.2f$", fit.rExec("GetChisquare") / fit.iExec("GetNDF"));

				pen p = StdPen(crvi+1);
				if (disabledObj.valid)
					p += dashed;

				draw(fit, "l", p, l);
			}

			limits((0.4, 0.5e-2), (1.2, 0.20), Crop);
			AttachLegend();
		}
	}

	GShipout("plots_fits_ranges_" + model, hSkip=0mm, vSkip=0mm);
}

//----------------------------------------------------------------------------------------------------

string ranges[][] = {
	new string[] { "LH,LH", "LC,LC", "LL,LL"},
	new string[] { "CH,CH", "CC,CC", "CL,CL"},
	new string[] { "HH,HH", "HC,HC", "HL,HL"},
};
MakeOne("local", new string[] {"f_dip", "f_bump"}, ranges);

string ranges[][] = {
	new string[] { "LH", "LC", "LL"},
	new string[] { "CH", "CC", "CL"},
	new string[] { "HH", "HC", "HL"},
};
MakeOne("exp2+exp3", new string[] {"f_global"}, ranges);
