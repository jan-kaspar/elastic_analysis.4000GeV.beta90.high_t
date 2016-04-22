import root;
import pad_layout;

string topDir = "../../";

xSizeDef = 10cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string binnings[] = {
	//"ub",
	"ob-1-30-0.10",
	"ob-2-20-0.10",
	"ob-3-10-0.10",
};

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$t$ bin size$\ung{GeV^2}$");
currentpad.xTicks = LeftTicks(0.3, 0.1);

AddToLegend("<{\it resolution considerations:}");
draw(scale(1, 1), rGetObj(topDir+"binning/generators.root", "g_rms_t"), red+longdashed, "1 smearing sigma");
draw(scale(1, 2), rGetObj(topDir+"binning/generators.root", "g_rms_t"), blue+longdashed, "2 smearing sigma");
draw(scale(1, 3), rGetObj(topDir+"binning/generators.root", "g_rms_t"), heavygreen+longdashed, "3 smearing sigma");

AddToLegend("<{\it fixed statistical uncertainty, DS4:}");

draw(rGetObj(topDir+"binning/generators.root", "g_bs_stat_unc_30"), red+dashed, "$30\un{\%}$");
draw(rGetObj(topDir+"binning/generators.root", "g_bs_stat_unc_20"), blue+dashed, "$20\un{\%}$");
draw(rGetObj(topDir+"binning/generators.root", "g_bs_stat_unc_10"), heavygreen+dashed, "$10\un{\%}$");

AddToLegend("<{\it binnings in analysis:}");

for (int bi : binnings.keys)
{
	pen p = StdPen(bi);
	draw(rGetObj(topDir+"DS4/distributions_45b_56t.root", "binning/g_binning_"+binnings[bi]), "p,l,d0", p+0.2pt, mCi+0.8pt+p, binnings[bi]);
}

limits((0, 0), (2.0, 0.12), Crop);

for (real x = 0; x <= 2.0; x += 0.3)
	yaxis(XEquals(x, false), dotted);

for (real y = 0; y <= 0.12; y += 0.01)
	xaxis(YEquals(y, false), dotted);

frame f_legend = BuildLegend();

NewPad(false);
add(shift(0, 100)*f_legend);
