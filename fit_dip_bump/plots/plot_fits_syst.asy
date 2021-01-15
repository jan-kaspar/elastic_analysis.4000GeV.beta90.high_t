import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

string models[];
string m_fits[][];
models.push("local"); m_fits.push(new string[] {"f_dip", "f_bump"});
models.push("exp2+exp3"); m_fits.push(new string[] {"f_global"});

string modes[][];
modes.push(new string[] {"alig-sh-thx:1", "alig-sh-thy:1", "tilt-thx-thy:1"});
modes.push(new string[] {"opt-m1:1", "opt-m2:1"});
modes.push(new string[] {"acc-corr-sigma-unc:1", "acc-corr-sigma-asym:1", "acc-corr-non-gauss:1"});
modes.push(new string[] {"eff-slp:1", "eff-slp:2", "beam-mom:1"});
modes.push(new string[] {"unsm-sigma-x:1", "unsm-sigma-y:1", "unsm-model:1"});
modes.push(new string[] {"norm:1"});

//----------------------------------------------------------------------------------------------------

for (int mi : models.keys)
{
	string model = models[mi];

	for (int ri : modes.keys)
	{
		NewRow();

		for (int ci : modes[ri].keys)
		{
			string mode = modes[ri][ci];

			NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");

			draw(RootGetObject(f, "h_dsdt"), "eb", heavygreen);
			for (int fi : m_fits[mi].keys)
				draw(RootGetObject(f, model + "/central/" + m_fits[mi][fi]), "l", StdPen(fi+1)+dashed);

			draw(RootGetObject(f, model + "/systematics/" + mode + "/h_dsdt_mod"), "eb", black);
			for (int fi : m_fits[mi].keys)
				draw(RootGetObject(f, model + "/systematics/" + mode + "/" + m_fits[mi][fi]), "l", StdPen(fi+1));

			limits((0.3, 1e-2), (1., 0.05), Crop);
			AttachLegend(mode);
		}
	}

	NewPad(false);
	AddToLegend("original histogram", mPl+5pt+heavygreen);
	AddToLegend("original dip fit", blue+dashed);
	AddToLegend("original bump fit", red+dashed);
	AddToLegend("modified histogram", mPl+5pt+black);
	AddToLegend("modified dip fit", blue);
	AddToLegend("modified bump fit", red);
	AttachLegend();

	GShipout("plots_fits_syst_" + model, hSkip=0mm, vSkip=0mm);
}
