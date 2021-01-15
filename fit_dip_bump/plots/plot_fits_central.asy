import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

xSizeDef = 8cm;
ySizeDef = 8cm;

string models[];
string m_fits[][];
models.push("local"); m_fits.push(new string[] {"f_dip", "f_bump"});
models.push("exp2+exp3"); m_fits.push(new string[] {"f_global"});

//----------------------------------------------------------------------------------------------------

for (int mi : models.keys)
{
	string model = models[mi];

	NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);

	draw(RootGetObject(f, "h_dsdt"), "eb", black);

	for (int fi : m_fits[mi].keys)
	{
		pen p = StdPen(fi + 1);

		RootObject fit = RootGetObject(f, model + "/central/" + m_fits[mi][fi]);

		string l = replace(m_fits[mi][fi], "_", "\_");
		l += format(", $\ch^2/ndf = %#.2f", fit.rExec("GetChisquare")) + format(" / %i", fit.iExec("GetNDF"))
			+ format(" = %#.2f$", fit.rExec("GetChisquare") / fit.iExec("GetNDF"));


		draw(fit, "l", p, l);
	}

	limits((0.40, 0.008), (1.1, 0.20), Crop);

	AttachLegend(model);
}
