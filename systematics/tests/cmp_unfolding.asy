import root;
import pad_layout;

string f = "unfolding_test.root";

xSizeDef = 10cm;

NewPad("$|t|\ung{GeV^2}$", "unfolding correction");
draw(rGetObj(f, "none/g_corr"), black, "nominal");
draw(rGetObj(f, "beam-mom3/g_corr"), red, "beam mom, $3\un{\si}$");
draw(rGetObj(f, "opt-m1*3/g_corr"), blue, "opt-m1, $3\un{\si}$");
draw(rGetObj(f, "opt-m2*3/g_corr"), heavygreen, "opt-m2, $3\un{\si}$");

limits((0, 0.8), (2.1, 1.05), Crop);
AttachLegend(SE, SE);

GShipout(margin=0mm);
