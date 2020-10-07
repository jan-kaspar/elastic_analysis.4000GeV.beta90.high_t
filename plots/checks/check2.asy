
import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "DS-merged/merged.root";
string o = "bt1/DS4-sc/combined/h_dsdt";

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

void DrawPoint(real x, real x_unc, real y, real y_unc)
{
	pen p = blue;

	real x_left = x - x_unc;
	real x_right = x + x_unc;
	draw(Scale((x_left, y))--Scale((x_right, y)), p);
	draw(Scale((x, y-y_unc))--Scale((x, y+y_unc)), p);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\un{mb/GeV^{2}}$");
scale(Linear, Log);

TH1_x_min = 0.3;
TH1_x_max = 1.2;

draw(RootGetObject(f, o), "vl", red, "this analysis");

// D0-TOTEM note, Table XXIII
AddToLegend("D0-TOTEM note, Table XXIII", blue);
DrawPoint(0.3833, 0.0063, 0.2227, 0.0157 );
DrawPoint(0.396 , 0.0064, 0.1579, 0.0119 );
DrawPoint(0.4088, 0.0064, 0.1199, 0.0091 );
DrawPoint(0.4216, 0.0065, 0.0901, 0.0071 );
DrawPoint(0.4345, 0.0065, 0.055 , 0.0053 );
DrawPoint(0.4475, 0.0065, 0.0416, 0.0042 );
DrawPoint(0.4606, 0.0066, 0.0345, 0.0036 );
DrawPoint(0.4738, 0.0066, 0.0232, 0.0029 );
DrawPoint(0.487 , 0.0066, 0.0226, 0.0027 );
DrawPoint(0.5003, 0.0067, 0.0172, 0.0024 );
DrawPoint(0.5137, 0.0067, 0.0142, 0.0023 );
DrawPoint(0.5272, 0.0068, 0.0151, 0.0023 );
DrawPoint(0.5474, 0.0135, 0.0182, 0.002  );
DrawPoint(0.5744, 0.0135, 0.0234, 0.0023 );
DrawPoint(0.614 , 0.0261, 0.0287, 0.0023 );
DrawPoint(0.6651, 0.025 , 0.0262, 0.0025 );
DrawPoint(0.715 , 0.0249, 0.031 , 0.0026 );
DrawPoint(0.7671, 0.0272, 0.0286, 0.0024 );
DrawPoint(0.8264, 0.0321, 0.0223, 0.0021 );
DrawPoint(0.879 , 0.0205, 0.0183, 0.0021 );
DrawPoint(1.0337, 0.031 , 0.0086, 0.001  );

AttachLegend();
