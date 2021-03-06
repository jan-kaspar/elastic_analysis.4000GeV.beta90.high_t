import root;
import pad_layout;

string topDir = "../";

string iteration = "iteration 1";

string data_fit = "exp3-intf-exp1";

//string binning = "ob-1-30-0.10";
string binning = "bt1";

TF1_x_min = 0;
TF1_x_max = 2;

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$", LeftTicks(0.5, 0.1));
scale(Linear, Log);

string fn = topDir + "DS4/unfolding_cf_45b_56t.root";
draw(RootGetObject(fn, binning + "/" + data_fit + "/+0,+0/" + iteration + "/input_corr"), "eb", black, "data 45 bot -- 56 top");

string fn = topDir + "DS4/unfolding_cf_45t_56b.root";
draw(RootGetObject(fn, binning + "/" + data_fit + "/+0,+0/" + iteration + "/input_corr"), "eb", red, "data 45 top -- 56 bot");

string fn = topDir + "DS4/unfolding_cf_45b_56t.root";
draw(RootGetObject(fn, binning + "/exp3-intf-exp1/+0,+0/" + iteration + "/ff"), blue, "fit exp3-intf-exp1");

string fn = topDir + "DS4/unfolding_cf_45t_56b.root";
draw(RootGetObject(fn, binning + "/p1*exp3+p2*exp2/+0,+0/" + iteration + "/ff"), green, "fit p1*exp3+p2*exp2");

limits((0, 1e-4), (1.9, 1e3), Crop);

AttachLegend();
