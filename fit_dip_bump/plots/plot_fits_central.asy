import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

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
		string l = replace(m_fits[mi][fi], "_", "\_");
		draw(RootGetObject(f, model + "/central/" + m_fits[mi][fi]), "l", p, l);
	}

	limits((0.25, 0.008), (1.1, 0.20), Crop);

	AttachLegend(model);
}
