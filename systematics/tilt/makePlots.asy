import pad_layout;
import root;


NewPad("$|t|\ung{GeV^2}$", "systematic effect", 12cm, 8cm);
currentpad.yTicks = RightTicks(0.01, 0.005);

draw(rGetObj("simu,exp,N_ev=1.0E+08,ga=+0.020,thysign=+1.root", "h_t_ratio"), "vl", black);
draw(rGetObj("simu,mod,N_ev=1.0E+08,ga=+0.020,thysign=-1.root", "h_t_ratio"), "vl", red);
draw(rGetObj("simu,mod,N_ev=1.0E+08,ga=+0.020,thysign=+1.root", "h_t_ratio"), "vl", blue);

limits((0, 0.98), (2.0, 1.04), Crop);
AttachLegend(NW, NW);

yaxis(XEquals(0.16, false), dotted);
yaxis(XEquals(0.2, false), dashed);
