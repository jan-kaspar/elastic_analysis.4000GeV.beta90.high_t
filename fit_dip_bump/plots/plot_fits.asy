import root;
import pad_layout;

string f = "make_fits.root";

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");
//scale(Linear, Log);

draw(RootGetObject(f, "h_dsdt"), "eb", black);

//draw(RootGetObject(f, "h_dsdt|ff0"), "l", red);
//draw(RootGetObject(f, "h_dsdt|ff1"), "l", heavygreen);

draw(RootGetObject(f, "h_dsdt|ff2"), "l", blue);
draw(RootGetObject(f, "h_dsdt|ff3"), "l", red);

limits((0.3, 1e-2), (1., 0.05), Crop);
