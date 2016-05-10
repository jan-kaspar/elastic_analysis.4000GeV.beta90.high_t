import root;
import pad_layout;

string dataset = "DS4";

string topDir = "../../";

string diagonal = "45b_56t";

string binning = "ob-2-20-0.20";

xSizeDef = 12cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$");

scale(Linear, Log);
string f = topDir+dataset+"/distributions_" + diagonal + ".root";
draw(RootGetObject(f, "acceptance correction/"+binning+"/h_t_before"), "eb", blue, "before acc.~corr.");
draw(RootGetObject(f, "acceptance correction/"+binning+"/h_t_after"), "eb", red, "after acc.~corr.");

limits((0.4, 1e2), (1.3, 1e5), Crop);

yaxis(XEquals(0.528, false), heavygreen+1.5pt);
yaxis(XEquals(0.672, false), magenta+1.5pt);
yaxis(XEquals(0.819, false), orange+1.5pt);
yaxis(XEquals(0.882, false), orange+1.5pt);

AttachLegend(replace(dataset + ", " + diagonal + ", " + binning, "_", "\_"));

GShipout(margin=1mm);
