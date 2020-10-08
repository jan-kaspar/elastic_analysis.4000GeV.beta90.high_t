import root;
import pad_layout;

string topDir = "../../";

xSizeDef = 12cm;
ySizeDef = 8cm;

string datasets[] = { "DS4" };

string diagonals[] = { "45b_56t", "45t_56b" };

//string binning = "ob-1-30-0.10";
string binning = "bt1";

string iteration = "iteration 2";

//----------------------------------------------------------------------------------------------------

void DrawOne(string fits[])
{
	for (int dsi : datasets.keys)
	{
		NewRow();

		NewPad(false);
	
		for (int dgni : diagonals.keys)
		{
			NewPad("$|t|\ung{GeV^2}$");
			scale(Linear, Log);
			currentpad.xTicks = LeftTicks(0.2, 0.1);
			
			string fn = topDir + datasets[dsi] + "/unfolding_cf_" + diagonals[dgni] + ".root";
			draw(RootGetObject(fn, binning + "/" + fits[0] + "/+0,+0/" + iteration + "/input_corr"), "eb");
	
			for (int fi : fits.keys)
			{
				write("* " + fits[fi]);

				RootObject obj = RootGetObject(fn, binning + "/" + fits[fi] + "/+0,+0/" + iteration + "/ff", error=false);
				if (!obj.valid)
					continue;
	
				TF1_x_min = 0;
				TF1_x_max = 2;
				draw(obj, StdPen(fi+1), fits[fi]);
			}
	
			limits((0, 1e-4), (2, 1e3), Crop);
			
			for (real x = 0; x <= 2; x += 0.2)
				yaxis(XEquals(x, false), dotted);
			
			for (real y = -4; y <= 3; y += 1)
				xaxis(YEquals(10^y, false), dotted);
	
			AttachLegend(replace(datasets[dsi] + ", " + diagonals[dgni], "_", "\_"));
		}
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

NewPad(false);
AddToLegend("<binning: " + binning);
AddToLegend("<iteration: " + iteration);
AttachLegend();

for (int dgni : diagonals.keys)
	NewPadLabel(replace(diagonals[dgni], "_", "\_"));


string fits[] = {
	"exp3+exp4",
	//"exp5+erf*exp2",
};

NewRow();
DrawOne(fits);

//--------------------

string fits[] = {
	"p1*exp3+p1*exp1",
	"p1*exp3+p2*exp2",
};

NewRow();
DrawOne(fits);

//--------------------

string fits[] = {
	"exp3-intf-exp1",
	"(exp3-intf-exp1)*expG",
};

NewRow();
DrawOne(fits);

//--------------------

GShipout(vSkip=0mm);
