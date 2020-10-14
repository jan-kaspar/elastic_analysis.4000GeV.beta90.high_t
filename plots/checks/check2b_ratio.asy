
import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "DS-merged/merged.root";
string o = "bt1/DS4-sc/combined/h_dsdt";
RootObject hist = RootGetObject(f, o);

xSizeDef = 8cm;

yTicksDef = RightTicks(0.005, 0.001);

//----------------------------------------------------------------------------------------------------

void DrawPoint(real x, real x_unc, real y, real y_unc)
{
	pen p = magenta;

	real x_left = x - x_unc;
	real x_right = x + x_unc;
	
	int bi = hist.iExec("FindBin", x);
	real y_this_anal = hist.rExec("GetBinContent", bi);

	real ratio = y_this_anal / y;

	draw(Scale((x_left, ratio))--Scale((x_right, ratio)), p);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "this analysis / D0-TOTEM table (not rounded)");
TH1_x_min = 0.3;
TH1_x_max = 1.2;

DrawPoint(0.38333382, 0.00628798, 0.22269673, 0.00753097);
DrawPoint(0.39599909, 0.00637729, 0.15783754, 0.00632121);
DrawPoint(0.40875975, 0.00638337, 0.11992759, 0.00552650);
DrawPoint(0.42160380, 0.00646068, 0.09007623, 0.00478101);
DrawPoint(0.43454648, 0.00648200, 0.05503046, 0.00373155);
DrawPoint(0.44753336, 0.00650488, 0.04161669, 0.00324487);
DrawPoint(0.46060757, 0.00656934, 0.03450659, 0.00296052);
DrawPoint(0.47376724, 0.00659033, 0.02323922, 0.00244567);
DrawPoint(0.48697963, 0.00662206, 0.02262528, 0.00244716);
DrawPoint(0.50027076, 0.00666906, 0.01719063, 0.00217241);
DrawPoint(0.51366690, 0.00672708, 0.01416986, 0.00201507);
DrawPoint(0.52715199, 0.00675801, 0.01515624, 0.00212531);
DrawPoint(0.54741000, 0.01350000, 0.01823910, 0.00170942);
DrawPoint(0.57441000, 0.01350000, 0.02335870, 0.00198143);
DrawPoint(0.61398503, 0.02607503, 0.02866802, 0.00162735);
DrawPoint(0.66508795, 0.02502789, 0.02616203, 0.00162001);
DrawPoint(0.71500955, 0.02489371, 0.03098456, 0.00180209);
DrawPoint(0.76714893, 0.02724567, 0.02855684, 0.00168114);
DrawPoint(0.82644768, 0.03205309, 0.02225333, 0.00139488);
DrawPoint(0.87902149, 0.02052071, 0.01831780, 0.00160556);
DrawPoint(1.03366284, 0.03102641, 0.00855712, 0.00092874);

limits((0.3, 0.99), (1.3, 1.01), Crop);
