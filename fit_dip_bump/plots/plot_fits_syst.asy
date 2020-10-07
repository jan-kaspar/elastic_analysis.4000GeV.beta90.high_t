import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

string model = "local";

string modes[][];
modes.push(new string[] {"alig-sh-thx:1", "alig-sh-thy:1", "tilt-thx-thy:1"});
modes.push(new string[] {"opt-m1:1", "opt-m2:1"});
modes.push(new string[] {"eff-slp:1", "eff-slp:2", "beam-mom:1", "norm:1"});

for (int ri : modes.keys)
{
	NewRow();

	for (int ci : modes[ri].keys)
	{
		string mode = modes[ri][ci];

		NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");

		draw(RootGetObject(f, "h_dsdt"), "eb", heavygreen);
		draw(RootGetObject(f, model + "/central/f_dip"), "l", blue+dashed);
		draw(RootGetObject(f, model + "/central/f_bump"), "l", red+dashed);

		draw(RootGetObject(f, model + "/" + mode + "/h_dsdt_mod"), "eb", black);
		draw(RootGetObject(f, model + "/" + mode + "/f_dip"), "l", blue);
		draw(RootGetObject(f, model + "/" + mode + "/f_bump"), "l", red);

		limits((0.3, 1e-2), (1., 0.05), Crop);
		AttachLegend(mode);
	}
}

GShipout(hSkip=0mm, vSkip=0mm);
