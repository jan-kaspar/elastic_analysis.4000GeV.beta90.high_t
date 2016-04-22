import root;
import pad_layout;

string datasets[] = { "DS4-b648", "DS4-b2990", "DS4-b26" };

string topDir = "../../";

string binnings[] = {
	"ob-1-30-0.10",
	"ob-2-20-0.20",
	"ob-3-10-0.30"
};

//----------------------------------------------------------------------------------------------------

TH1_x_max = 1.85;

for (int bi : binnings.keys)
{
	NewPad("$|t|\ung{GeV^2}$", "$\d\si_{\rm el}/\d t \ung{mb/GeV^2}$", 10cm, 8cm);
	currentpad.xTicks = LeftTicks(0.2, 0.1);
	scale(Linear, Log);
	
	for (int di : datasets.keys)
	{
		draw(rGetObj(topDir+"DS-merged/merged.root", binnings[bi]+"/"+datasets[di]+"/combined/h_dsdt" ),
			"eb", StdPen(di + 1), datasets[di]);
	}
	
	limits((0, 1e-4), (2., 1e3), Crop);

	AttachLegend(binnings[bi]);
}

