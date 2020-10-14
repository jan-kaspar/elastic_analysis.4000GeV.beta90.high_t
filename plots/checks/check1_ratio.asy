import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "DS-merged/merged.root";
string o = "NPB/DS4-sc/combined/h_dsdt";
RootObject hist = RootGetObject(f, o);

xSizeDef = 8cm;
xTicksDef = LeftTicks(0.05, 0.01);

//----------------------------------------------------------------------------------------------------

void DrawPoint(real x_left, real x_right, real y, real y_unc)
{
	pen p = heavygreen;
	
	real x = (x_left + x_right)/2.;

	int bi = hist.iExec("FindBin", x);
	real y_this_anal = hist.rExec("GetBinContent", bi);

	real ratio = y_this_anal / y;

	draw(Scale((x_left, ratio))--Scale((x_right, ratio)), p);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "this analysis / NPB");
//scale(Linear, Log);

DrawPoint(0.02697, 0.03005, 305.09 , 0.527 );
DrawPoint(0.03005, 0.03325, 287.95 , 0.478 );
DrawPoint(0.03325, 0.03658, 269.24 , 0.436 );
DrawPoint(0.03658, 0.04005, 251.31 , 0.401 );
DrawPoint(0.04005, 0.04365, 235.15 , 0.371 );
DrawPoint(0.04365, 0.04740, 218.32 , 0.343 );
DrawPoint(0.04740, 0.05129, 202.64 , 0.318 );
DrawPoint(0.05129, 0.05534, 187.10 , 0.295 );
DrawPoint(0.05534, 0.05956, 173.06 , 0.274 );
DrawPoint(0.05956, 0.06394, 158.77 , 0.255 );
DrawPoint(0.06394, 0.06850, 144.93 , 0.236 );
DrawPoint(0.06850, 0.07324, 133.12 , 0.219 );
DrawPoint(0.07324, 0.07817, 121.24 , 0.203 );
DrawPoint(0.07817, 0.08329, 109.77 , 0.188 );
DrawPoint(0.08329, 0.08862,  99.077, 0.174 );
DrawPoint(0.08862, 0.09417,  89.126, 0.161 );
DrawPoint(0.09417, 0.09994,  79.951, 0.148 );
DrawPoint(0.09994, 0.10593,  71.614, 0.137 );
DrawPoint(0.10593, 0.11217,  63.340, 0.125 );
DrawPoint(0.11217, 0.11866,  56.218, 0.115 );
DrawPoint(0.11866, 0.12540,  49.404, 0.105 );
DrawPoint(0.12540, 0.13242,  43.300, 0.0961);
DrawPoint(0.13242, 0.13972,  37.790, 0.0876);
DrawPoint(0.13972, 0.14730,  32.650, 0.0795);
DrawPoint(0.14730, 0.15520,  28.113, 0.0720);
DrawPoint(0.15520, 0.16340,  24.155, 0.0659);
DrawPoint(0.16340, 0.17194,  20.645, 0.0616);
DrawPoint(0.17194, 0.18082,  17.486, 0.0574);
DrawPoint(0.18082, 0.19005,  14.679, 0.0543);
DrawPoint(0.19005, 0.19965,  12.291, 0.0504);

limits((0, 0.999), (0.20, 1.001), Crop);

AttachLegend();
