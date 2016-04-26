import root;
import pad_layout;

xSizeDef = 10cm;
ySizeDef = 8cm;

yTicksDef = RightTicks(0.005, 0.001);

string f_mlt = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta90/" + "systematics/DS4/matrix_direct_45b_56t.root";
string f_lt = "/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta90/" + "systematics/DS4/simu_direct_45b_56t.root";

string f_ht = "../../DS4/systematics_ni_45b_56t.root";

TH1_x_max = 0.3;

real y_min = 0.98, y_max = 1.02;

drawGridDef = true;

//----------------------------------------------------------------------------------------------------

NewPad();
//draw(RootGetObject(f_mlt, "alignment-shy/contribution: de_th_y"), "eb", red);
draw(RootGetObject(f_lt, "de_th_x/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "alig-sh-thx/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("alig-sh-thx");

NewPad();
draw(RootGetObject(f_lt, "de_th_y/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "alig-sh-thy/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("alig-sh-thy");

NewPad();
draw(RootGetObject(f_lt, "tilt/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "tilt-thx-thy/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("tilt-thx-thy");


//----------------------------------------------------------------------------------------------------
NewRow();

NewPad();
draw(RootGetObject(f_lt, "scale_mode1/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "opt-m1/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("opt-m1");

NewPad();
draw(RootGetObject(f_lt, "scale_mode2/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "opt-m2/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("opt-m2");

//----------------------------------------------------------------------------------------------------
NewRow();

y_min = 0.995; y_max = 1.005;
yTicksDef = RightTicks(0.001, 0.0005);

NewPad();
draw(RootGetObject(f_lt, "de_si_th_y/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "acc-corr-sigma-unc/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("beam divergence: RMS uncertainty");

NewPad();
draw(RootGetObject(f_lt, "sm_asym/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "acc-corr-sigma-asym/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("beam divergence: left-right assymetry");

NewPad();
draw(RootGetObject(f_lt, "sm_non_gauss/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "acc-corr-non-gauss/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("beam divergence: non-Gaussianity");

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad();
draw(RootGetObject(f_lt, "eff_slp/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "eff-slp/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("efficiency: slope");

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad();
draw(RootGetObject(f_lt, "de_p/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "beam-mom/exp3-intf-exp1/g_r"), blue);
limits((0, y_min), (0.4, y_max), Crop);
AttachLegend("beam-mom");

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad(false);
AttachLegend("unfolding");

//----------------------------------------------------------------------------------------------------
NewRow();

yTicksDef = RightTicks(0.01, 0.002);

NewPad();
draw(RootGetObject(f_lt, "norm/h_eff_syst"), "eb", red);
draw(RootGetObject(f_ht, "norm/exp3-intf-exp1/g_r"), blue);
limits((0, 1.), (0.4, 1.05), Crop);
currentpad.yTicks = RightTicks();
AttachLegend("norm");
