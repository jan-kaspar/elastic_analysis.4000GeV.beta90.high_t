import root;
import pad_layout;

string topDir = "../";

string f = topDir + "make_fits.root";

string model = "local";

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");
//scale(Linear, Log);

draw(RootGetObject(f, "h_dsdt"), "eb", black);

draw(RootGetObject(f, model + "/central/f_dip"), "l", blue);
draw(RootGetObject(f, model + "/central/f_bump"), "l", red);

limits((0.3, 1e-2), (1., 0.05), Crop);
