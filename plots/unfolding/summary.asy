import root;
import pad_layout;

string topDir = "../../";

string diagonal = "45b_56t";
string binning = "ob-1-30-0.10";

xSizeDef = 10cm;

drawGridDef = true;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "unfolding correction");


draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/exp3+exp4/+0,+0/corr_final"), "vl", black+1pt);
draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/exp5+erf*exp2/+0,+0/corr_final"), "vl", black+1pt);
draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/p1*exp3+p1*exp1/+0,+0/corr_final"), "vl", black+1pt);
draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/p1*exp3+p2*exp2/+0,+0/corr_final"), "vl", black+1pt);
draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/exp3-intf-exp1/+0,+0/corr_final"), "vl", black+1pt);
draw(rGetObj(topDir + "DS4/unfolding_cf_"+diagonal+".root", binning+"/(exp3-intf-exp1)*expG/+0,+0/corr_final"), "vl", black+1pt);

draw(rGetObj(topDir + "DS4/unfolding_gr_"+diagonal+".root", binning+"/smearing_matrix_mc_"+diagonal+".root,p1*exp3+p2*exp2,"+binning+"/alpha=1.00E+00/h_corr"), "vl", blue+1pt);
draw(rGetObj(topDir + "DS4/unfolding_gr_"+diagonal+".root", binning+"/smearing_matrix_mc_"+diagonal+".root,p1*exp3+p2*exp2,"+binning+"/alpha=1.00E-01/h_corr"), "vl", blue+1pt);
draw(rGetObj(topDir + "DS4/unfolding_gr_"+diagonal+".root", binning+"/smearing_matrix_mc_"+diagonal+".root,exp3-intf-exp1,"+binning+"/alpha=1.00E+00/h_corr"), "vl", blue+1pt);
draw(rGetObj(topDir + "DS4/unfolding_gr_"+diagonal+".root", binning+"/smearing_matrix_mc_"+diagonal+".root,exp3-intf-exp1,"+binning+"/alpha=1.00E-01/h_corr"), "vl", blue+1pt);

draw(rGetObj(topDir + "DS4/unfolding_summarize_"+diagonal+".root", binning+"/model/h_mean_stddev"), "eb", red);

limits((0, 0.8), (2., 1.10), Crop);
currentpad.xTicks = LeftTicks(0.2, 0.1);
currentpad.yTicks = RightTicks(0.05, 0.01);
