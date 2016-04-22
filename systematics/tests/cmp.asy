import root;
import pad_layout;

string f = "test.root";

NewPad("$|t|\ung{GeV^2}$", "relative error");
draw(scale(1, 3) * shift(0, -1), rGetObj(f, "opt-m1/data fit 1/g_r"), black, "$3\cdot 1\un{\si}$");
draw(scale(1, 1) * shift(0, -1), rGetObj(f, "opt-m1*3/data fit 1/g_r"), red, "$1\cdot 3\un{\si}$");
AttachLegend("opt-m1", SE, SE);

NewPad("$|t|\ung{GeV^2}$", "relative error");
draw(scale(1, 3) * shift(0, -1), rGetObj(f, "opt-m2/data fit 1/g_r"), black, "$3\cdot 1\un{\si}$");
draw(scale(1, 1) * shift(0, -1), rGetObj(f, "opt-m2*3/data fit 1/g_r"), red, "$1\cdot 3\un{\si}$");
AttachLegend("opt-m2", SE, SE);

NewPad("$|t|\ung{GeV^2}$", "relative error");
draw(scale(1, 3) * shift(0, -1), rGetObj(f, "beam-mom1/data fit 1/g_r"), black, "$3\cdot 1\un{\si}$");
draw(scale(1, 1) * shift(0, -1), rGetObj(f, "beam-mom3/data fit 1/g_r"), red  , "$1\cdot 3\un{\si}$");
AttachLegend("beam mom.", SE, SE);

GShipout(margin=0mm);
