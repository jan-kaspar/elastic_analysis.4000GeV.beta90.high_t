#include <string>
#include <vector>
#include <map>
#include <cmath>

double timestamp0 = 1342044000;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

void Init_base()
{
	// selection of bunches
	bunchMap[8369].push_back(648); bunchMap[8369].push_back(2990);
	bunchMap[8371].push_back(648); bunchMap[8371].push_back(2990);
	bunchMap[8372].push_back(648); bunchMap[8372].push_back(2990); bunchMap[8372].push_back(26);

	// alignment settings
	AlignmentSource alSrc1;
	alSrc1.SetAlignmentA(atNone);
	alSrc1.SetAlignmentB(atTimeDependent, "alignment_fine_corrections.root");
	alSrc1.SetAlignmentC(atNone);
	alignmentSources.push_back(alSrc1);

	// environment settings
	env.InitNominal();
	env.UseMatchedOptics();
	env.si_tilt = 0.2E-3;	// mrad

	// binning
	anal.t_min = 0.02; anal.t_max = 1.9;
	anal.t_min_full = 0.; anal.t_max_full = 2.1;

	// analysis settings
	anal.th_y_lcut = 0E-6; anal.th_y_lcut = 0E-6;
	anal.th_y_lcut_L = 0E-6; anal.th_y_lcut_R = 0E-6;
	anal.th_y_hcut_L = 0E-6; anal.th_y_hcut_R = 0E-6;
	
	anal.th_x_lcut = -1.;
	anal.th_x_hcut = +1.;
	
	anal.use_time_dependent_resolutions = true;
	
	// approximate (time independent) resolutions
	anal.si_th_y_1arm = 3.3E-6 / sqrt(2.);
	anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
	anal.si_th_x_1arm = anal.si_th_x_2arm = 0.;	// diagonal dependent

	anal.si_th_y_1arm_unc = 0.05E-6 / sqrt(2.);
	anal.si_th_x_1arm_unc = 0.2E-6 / sqrt(2.);
	
	anal.use_3outof4_efficiency_fits = true;
	anal.use_pileup_efficiency_fits = false;
	anal.inefficiency_3outof4 = 0.;				// diagonal dependent
	anal.inefficiency_shower_near = 0.03;
	anal.inefficiency_pile_up = 0.;				// diagonal dependent
	anal.inefficiency_trigger = 0.;
	anal.inefficiency_DAQ = 0.020;

	anal.bckg_corr = 1.;
	
	anal.L_int = 734.494E3;		// scaled to match with DS2 (equal dsigma/dt integral over "ub" bins 5 to 12)

	anal.alignment_t0 = 67300;	// beginning of the first time-slice
	anal.alignment_ts = 10*60.;	// time-slice in s
	
	anal.eff_th_y_min = 35E-6;

	anal.t_min_fit = 0.027;

	cmsZeroBiasSelection = true;

	anal.alignmentYRanges["L_F"] = Analysis::AlignmentYRange(-26., -8.8, 9.2, 27.);
	anal.alignmentYRanges["L_N"] = Analysis::AlignmentYRange(-23., -7.8, 8.0, 24.);
	anal.alignmentYRanges["R_N"] = Analysis::AlignmentYRange(-23., -7.8, 8.0, 23.);
	anal.alignmentYRanges["R_F"] = Analysis::AlignmentYRange(-26., -8.8, 9.0, 26.);

	unsmearing_file = "";	// diagonal dependent
	unsmearing_object = "<binning>/exp3-intf-exp1/+0,+0/corr_final";

	luminosity_data_file = "../fill_2836_lumiCalc2.py_V04-02-09_lumibylsXing.csv";
}

//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	AlignmentSource alSrc2;
	alSrc2.SetAlignmentA(atConstant);
	alSrc2.SetAlignmentB(atConstant);
	alSrc2.SetAlignmentC(atConstant);
	alSrc2.cnst.a_L_F = +0.14E-3; alSrc2.cnst.b_L_F = +2.1E-3; alSrc2.cnst.c_L_F = 0E-3;
	alSrc2.cnst.a_L_N = +0.52E-3; alSrc2.cnst.b_L_N = +8.1E-3; alSrc2.cnst.c_L_N = 0E-3;
	alSrc2.cnst.a_R_N = -0.04E-3; alSrc2.cnst.b_R_N = +1.7E-3; alSrc2.cnst.c_R_N = 0E-3;
	alSrc2.cnst.a_R_F = +0.17E-3; alSrc2.cnst.b_R_F = -2.0E-3; alSrc2.cnst.c_R_F = 0E-3;
	alignmentSources.push_back(alSrc2);

	anal.cut1_a = 1.0; anal.cut1_c = 0E-6; anal.cut1_si = 9.6E-6;
	anal.cut2_a = 1.0; anal.cut2_c = 0E-6; anal.cut2_si = 3.3E-6;

	anal.cut3_a = 0.; anal.cut3_b = 0.; anal.cut3_si = 0.17;
	anal.cut4_a = 0.; anal.cut4_b = 0.; anal.cut4_si = 0.17;

	anal.cut5_a = 0.10711; anal.cut5_b = -4.4E-3; anal.cut5_si = 18E-3;
	anal.cut6_a = 0.10716; anal.cut6_b = +2.2E-3; anal.cut6_si = 18E-3;

	anal.cut7_al = 0.; anal.cut7_c = 0E-3; anal.cut7_si = 8.9E-3;

	anal.th_y_lcut_L = 32.3E-6; anal.th_y_lcut_R = 31.0E-6;	anal.th_y_lcut = 33.8E-6;
	anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;
	
	anal.si_th_x_1arm = 9.3E-6 / sqrt(2.);
	anal.si_th_x_2arm = anal.si_th_x_1arm / sqrt(2.);
	
	unsmearing_file = "unfolding_cf_45b_56t.root";

	anal.inefficiency_3outof4 = 0.080;
	anal.inefficiency_pile_up = 0.01;
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	AlignmentSource alSrc2;
	alSrc2.SetAlignmentA(atConstant);
	alSrc2.SetAlignmentB(atConstant);
	alSrc2.SetAlignmentC(atConstant);
	alSrc2.cnst.a_L_F = -0.17E-3; alSrc2.cnst.b_L_F = +2.8E-3; alSrc2.cnst.c_L_F = 0E-3;
	alSrc2.cnst.a_L_N = -0.29E-3; alSrc2.cnst.b_L_N = +3.7E-3; alSrc2.cnst.c_L_N = 0E-3;
	alSrc2.cnst.a_R_N = +0.19E-3; alSrc2.cnst.b_R_N = +1.7E-3; alSrc2.cnst.c_R_N = 0E-3;
	alSrc2.cnst.a_R_F = +0.16E-3; alSrc2.cnst.b_R_F = +2.6E-3; alSrc2.cnst.c_R_F = 0E-3;
	alignmentSources.push_back(alSrc2);

	anal.cut1_a = 1.0; anal.cut1_c = 0E-6; anal.cut1_si = 8.9E-6;
	anal.cut2_a = 1.0; anal.cut2_c = 0E-6; anal.cut2_si = 3.3E-6;

	anal.cut3_a = 0.; anal.cut3_b = 0.; anal.cut3_si = 0.17;
	anal.cut4_a = 0.; anal.cut4_b = 0.; anal.cut4_si = 0.17;

	anal.cut5_a = 0.10730; anal.cut5_b = -2.4E-3; anal.cut5_si = 17E-3;
	anal.cut6_a = 0.10698; anal.cut6_b = -1.0E-3; anal.cut6_si = 17E-3;

	anal.cut7_al = 0.; anal.cut7_c = 0E-3; anal.cut7_si = 8.5E-3;

	anal.th_y_lcut_L = 31.4E-6; anal.th_y_lcut_R = 32.3E-6; anal.th_y_lcut = 33.8E-6;
	anal.th_y_hcut_L = 107E-6; anal.th_y_hcut_R = 107E-6; anal.th_y_hcut = 105E-6;
	
	anal.si_th_x_1arm = 8.7E-6 / sqrt(2.);
	anal.si_th_x_2arm = anal.si_th_x_1arm / sqrt(2.);
	
	unsmearing_file = "unfolding_cf_45t_56b.root";

	anal.inefficiency_3outof4 = 0.079;
	anal.inefficiency_pile_up = 0.01;
}
