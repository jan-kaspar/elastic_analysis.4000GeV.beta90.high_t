void Environment::InitNominal()
{
	// beam momentum (GeV)
	p = p_L = p_R = 4000.;

	// momentum uncertainty
	si_de_p = 1E-3 * p;

	// angular (one-side) beam smearing (rad)
	si_th_x = 2.6E-6;
	si_th_y = 2.3E-6;

	// vertex smearing (mm)
	si_vtx_x = si_vtx_y = 168E-3;

	// pitch-induced error
	si_de_P_L = si_de_P_R = 12E-3;

	// optics
	L_x_L_N = 2.878659768E3;		L_x_R_N = 2.878659768E3;		// mm
	L_x_L_F = -0.000005011626066E3;	L_x_R_F = -0.000005011626066E3;	// mm

	v_x_L_N = -2.164426367;			v_x_R_N = -2.164426367;			// 1
	v_x_L_F = -1.866140123;			v_x_R_F = -1.866140123;			// 1

	// dispersion, xi of diffractive protons is treated as negative
	// D_x_N = -0.056814277	// m
	// D_x_F = -0.043283791	// m

	L_y_L_N = 237.6888478E3;		L_y_R_N = 237.6888478E3;		// mm
	L_y_L_F = 263.1640863E3;		L_y_R_F = 263.1640863E3;		// mm

	v_y_L_N = 0.02036775253;		v_y_R_N = 0.02036775253;		// 1
	v_y_L_F = -0.00005020156816;	v_y_R_F = -0.00005020156816;	// 1

	// optics imperfections
	double opt_cov_data[] = {
		1.6415982538E-05,	8.3528909654E-04,	-9.2615663054E-05,	-3.1171124845E-03,	1.4312822946E-05,	7.0270087260E-04,	-1.0282458307E-04,	-3.4474855920E-03,		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		8.3528909654E-04,	7.9460907454E-02,	-4.0524897063E-03,	-1.9403987680E-01,	6.0020502912E-04,	6.8081967288E-02,	-4.5022866642E-03,	-2.1752327565E-01,      0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		-9.2615663054E-05,	-4.0524897063E-03,	9.8198811839E-04,	2.9916636217E-02,	-8.4213447878E-05,	-3.3562650855E-03,	1.0903026143E-03,	3.3088576215E-02,       0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		-3.1171124845E-03,	-1.9403987680E-01,	2.9916636217E-02,	1.0140738647E+00,	-2.5898010263E-03,	-1.6373502865E-01,	3.3220868493E-02,	1.1266606006E+00,       0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		1.4312822946E-05,	6.0020502912E-04,	-8.4213447878E-05,	-2.5898010263E-03,	1.3386771445E-05,	4.9929455397E-04,	-9.3487135895E-05,	-2.8454523500E-03,      0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		7.0270087260E-04,	6.8081967288E-02,	-3.3562650855E-03,	-1.6373502865E-01,	4.9929455397E-04,	5.8361304436E-02,	-3.7289177419E-03,	-1.8366186756E-01,      0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		-1.0282458307E-04,	-4.5022866642E-03,	1.0903026143E-03,	3.3220868493E-02,	-9.3487135895E-05,	-3.7289177419E-03,	1.2105646284E-03,	3.6743376699E-02,       0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,
		-3.4474855920E-03,	-2.1752327565E-01,	3.3088576215E-02,	1.1266606006E+00,	-2.8454523500E-03,	-1.8366186756E-01,	3.6743376699E-02,	1.2521286529E+00,       0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,

		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,		1.6415982538E-05,	8.3528909654E-04,	-9.2615663054E-05,	-3.1171124845E-03,	1.4312822946E-05,	7.0270087260E-04,	-1.0282458307E-04,	-3.4474855920E-03,
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    8.3528909654E-04,	7.9460907454E-02,	-4.0524897063E-03,	-1.9403987680E-01,	6.0020502912E-04,	6.8081967288E-02,	-4.5022866642E-03,	-2.1752327565E-01,
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    -9.2615663054E-05,	-4.0524897063E-03,	9.8198811839E-04,	2.9916636217E-02,	-8.4213447878E-05,	-3.3562650855E-03,	1.0903026143E-03,	3.3088576215E-02, 
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    -3.1171124845E-03,	-1.9403987680E-01,	2.9916636217E-02,	1.0140738647E+00,	-2.5898010263E-03,	-1.6373502865E-01,	3.3220868493E-02,	1.1266606006E+00, 
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    1.4312822946E-05,	6.0020502912E-04,	-8.4213447878E-05,	-2.5898010263E-03,	1.3386771445E-05,	4.9929455397E-04,	-9.3487135895E-05,	-2.8454523500E-03,
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    7.0270087260E-04,	6.8081967288E-02,	-3.3562650855E-03,	-1.6373502865E-01,	4.9929455397E-04,	5.8361304436E-02,	-3.7289177419E-03,	-1.8366186756E-01,
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    -1.0282458307E-04,	-4.5022866642E-03,	1.0903026143E-03,	3.3220868493E-02,	-9.3487135895E-05,	-3.7289177419E-03,	1.2105646284E-03,	3.6743376699E-02, 
		0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	0.0000000000E+00,	    -3.4474855920E-03,	-2.1752327565E-01,	3.3088576215E-02,	1.1266606006E+00,	-2.8454523500E-03,	-1.8366186756E-01,	3.6743376699E-02,	1.2521286529E+00
	};
	opt_cov.SetMatrixArray(opt_cov_data);

	// IMPORTANT: uncrease uncertainties by factor 3 - experimental observation
	opt_cov *= 3.*3.;

	TMatrixDSymEigen eig_decomp(opt_cov);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(16);
	for (unsigned int i = 0; i < 16; i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	opt_per_gen = eig_decomp.GetEigenVectors() * S;

	// alignment uncertainties, per RP (uncorrelated)
	si_de_x = 2E-3;
	si_de_y = 100E-3;
	si_tilt = 0.5E-3;
}

//----------------------------------------------------------------------------------------------------

void Environment::ApplyOpticsCorrections()
{
	printf(">> Environment::ApplyOpticsCorrections\n");

	// leading optics-perturbations mode
	vector<double> opt_per_mode1(8);
	opt_per_mode1[0] = 0.00366995;	// v_x_N
	opt_per_mode1[1] = 0.220136E3;	// L_x_N, mm
	opt_per_mode1[2] = -0.02826;	// v_y_N
	opt_per_mode1[3] = -0.967594E3;	// L_y_N, mm
	opt_per_mode1[4] = 0.00310878;	// v_x_F
	opt_per_mode1[5] = 0.185767E3;	// L_x_F, mm
	opt_per_mode1[6] = -0.0313802;	// v_y_F
	opt_per_mode1[7] = -1.07369E3;	// L_y_F, mm
	
	// experimental L_x_F
	double L_x_L_F_exp = -600.;	// mm
	double L_x_R_F_exp = -370.;	// mm
	
	// apply optics corrections
	printf(">> applying optics corrections (leading mode, x and y)\n");
	double f_L = L_x_L_F_exp / opt_per_mode1[5];
	double f_R = L_x_R_F_exp / opt_per_mode1[5];
	printf("\tf_L = %.3f, f_R = %.3f\n", f_L, f_R);

	v_x_L_N += f_L * opt_per_mode1[0]; v_x_R_N += f_R * opt_per_mode1[0];
	L_x_L_N += f_L * opt_per_mode1[1]; L_x_R_N += f_R * opt_per_mode1[1];
	v_y_L_N += f_L * opt_per_mode1[2]; v_y_R_N += f_R * opt_per_mode1[2];
	L_y_L_N += f_L * opt_per_mode1[3]; L_y_R_N += f_R * opt_per_mode1[3];

	v_x_L_F += f_L * opt_per_mode1[4]; v_x_R_F += f_R * opt_per_mode1[4];
	L_x_L_F += f_L * opt_per_mode1[5]; L_x_R_F += f_R * opt_per_mode1[5];
	v_y_L_F += f_L * opt_per_mode1[6]; v_y_R_F += f_R * opt_per_mode1[6];
	L_y_L_F += f_L * opt_per_mode1[7]; L_y_R_F += f_R * opt_per_mode1[7];

	// updated (reduced optics uncertainties)
	double opt_cov_data[] = {
		+1.03836275E-04, +2.73215249E-04, -8.60780959E-04, -2.03064601E-02, +1.08682917E-04, +4.93891542E-05, -9.55297956E-04, -2.20794305E-02,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		+2.73215249E-04, +4.04227854E-03, -2.60251412E-03, -6.60187321E-02, +2.88581966E-04, +2.96252833E-03, -2.88885011E-03, -7.20746055E-02,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		-8.60780959E-04, -2.60251412E-03, +1.18609444E-02, +3.05914979E-01, -9.12373424E-04, -3.80315126E-04, +1.31666956E-02, +3.35780897E-01,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		-2.03064601E-02, -6.60187321E-02, +3.05914979E-01, +8.09347800E+00, -2.12435697E-02, -1.30216077E-02, +3.39603518E-01, +8.90096140E+00,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		+1.08682917E-04, +2.88581966E-04, -9.12373424E-04, -2.12435697E-02, +1.17879280E-04, +4.18370862E-05, -1.01257562E-03, -2.30208234E-02,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		+4.93891542E-05, +2.96252833E-03, -3.80315126E-04, -1.30216077E-02, +4.18370862E-05, +2.50000000E-03, -4.22305899E-04, -1.44494178E-02,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		-9.55297956E-04, -2.88885011E-03, +1.31666956E-02, +3.39603518E-01, -1.01257562E-03, -4.22305899E-04, +1.46161963E-02, +3.72759400E-01,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,
		-2.20794305E-02, -7.20746055E-02, +3.35780897E-01, +8.90096140E+00, -2.30208234E-02, -1.44494178E-02, +3.72759400E-01, +9.79166758E+00,    0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,

		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    +1.03836275E-04, +2.73215249E-04, -8.60780959E-04, -2.03064601E-02, +1.08682917E-04, +4.93891542E-05, -9.55297956E-04, -2.20794305E-02,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    +2.73215249E-04, +4.04227854E-03, -2.60251412E-03, -6.60187321E-02, +2.88581966E-04, +2.96252833E-03, -2.88885011E-03, -7.20746055E-02,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    -8.60780959E-04, -2.60251412E-03, +1.18609444E-02, +3.05914979E-01, -9.12373424E-04, -3.80315126E-04, +1.31666956E-02, +3.35780897E-01,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    -2.03064601E-02, -6.60187321E-02, +3.05914979E-01, +8.09347800E+00, -2.12435697E-02, -1.30216077E-02, +3.39603518E-01, +8.90096140E+00,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    +1.08682917E-04, +2.88581966E-04, -9.12373424E-04, -2.12435697E-02, +1.17879280E-04, +4.18370862E-05, -1.01257562E-03, -2.30208234E-02,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    +4.93891542E-05, +2.96252833E-03, -3.80315126E-04, -1.30216077E-02, +4.18370862E-05, +2.50000000E-03, -4.22305899E-04, -1.44494178E-02,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    -9.55297956E-04, -2.88885011E-03, +1.31666956E-02, +3.39603518E-01, -1.01257562E-03, -4.22305899E-04, +1.46161963E-02, +3.72759400E-01,
		0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00,    -2.20794305E-02, -7.20746055E-02, +3.35780897E-01, +8.90096140E+00, -2.30208234E-02, -1.44494178E-02, +3.72759400E-01, +9.79166758E+00
	};
	opt_cov.SetMatrixArray(opt_cov_data);

	TMatrixDSymEigen eig_decomp(opt_cov);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(16);
	for (unsigned int i = 0; i < 16; i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	opt_per_gen = eig_decomp.GetEigenVectors() * S;
}

//----------------------------------------------------------------------------------------------------

void Environment::UseMatchedOptics()
{
	printf(">> Environment::UseMatchedOptics\n");
	printf("\tmatched optics sent by Frici on 17 December 2013\n");

	// optics
	L_x_L_N = 2.20334790903323E3;		L_x_R_N = 2.45111090171735E3;		// mm
	L_x_L_F = -0.575998338144321E3;		L_x_R_F = -0.366458055428883E3;		// mm

	v_x_L_N = -2.17331142381815;		v_x_R_N = -2.1661935965758;			// 1
	v_x_L_F = -1.87288728332818;		v_x_R_F = -1.86665604150867;		// 1

	L_y_L_N = 240.529357858677E3;		L_y_R_N = 238.509126939124E3;		// mm
	L_y_L_F = 266.332295363533E3;		L_y_R_F = 264.09748697684E3;		// mm

	v_y_L_N = 0.0928326363392919;		v_y_R_N = 0.0396581314265547;		// 1
	v_y_L_F = 0.0804304419601075;		v_y_R_F = 0.0214013552003497;		// 1

	printf("\tcovariance matrix sent by Frici on 31 January 2014, scaled by (1.35)^2\n");

	// L's in m, v's in 1
	double opt_cov_data[] = {
		7.79924295E-06, 1.43015767E-05, -5.41575923E-05, -1.21807517E-03, 7.81863435E-06, -2.49881153E-07, -6.00960263E-05, -1.32345029E-03, 6.97902683E-06, 1.55856920E-05, -5.42228378E-05, -1.22680312E-03, 7.64961570E-06, 2.29722480E-07, -6.01754873E-05, -1.32406994E-03,
		1.43015767E-05, 3.60175208E-05, -1.14158302E-04, -2.50043355E-03, 1.68904562E-05, 4.19360895E-07, -1.26698560E-04, -2.67563048E-03, 1.54708562E-05, 3.41879130E-05, -1.14140988E-04, -2.48731155E-03, 1.70479019E-05, -1.20749009E-07, -1.26666848E-04, -2.67546645E-03,
		-5.41575923E-05, -1.14158302E-04, 5.23042920E-04, 1.32510877E-02, -5.89824788E-05, 1.29170234E-06, 5.80537328E-04, 1.44268553E-02, -5.40553500E-05, -1.14818594E-04, 5.19917333E-04, 1.32383849E-02, -5.80832573E-05, -1.35582154E-06, 5.77065465E-04, 1.44279853E-02,
		-1.21807517E-03, -2.50043355E-03, 1.32510877E-02, 3.53845665E-01, -1.31024081E-03, 2.09744235E-05, 1.47085592E-02, 3.86787353E-01, -1.22389076E-03, -2.50644780E-03, 1.32023176E-02, 3.53503035E-01, -1.29613649E-03, -2.38536090E-05, 1.46542669E-02, 3.86756370E-01,
		7.81863435E-06, 1.68904562E-05, -5.89824788E-05, -1.31024081E-03, 8.59511048E-06, -2.58375825E-07, -6.54567278E-05, -1.41163378E-03, 7.64259908E-06, 1.71995157E-05, -5.85642150E-05, -1.31013693E-03, 8.40945240E-06, 2.57382563E-07, -6.49929015E-05, -1.41235184E-03,
		-2.49881153E-07, 4.19360895E-07, 1.29170234E-06, 2.09744235E-05, -2.58375825E-07, 7.25050643E-07, 1.43232280E-06, 2.21852925E-05, 1.74435667E-07, -3.24556268E-07, 1.04949209E-07, 2.08275300E-05, 1.62823790E-07, -5.19667650E-07, 1.18107113E-07, 2.35915335E-05,
		-6.00960263E-05, -1.26698560E-04, 5.80537328E-04, 1.47085592E-02, -6.54567278E-05, 1.43232280E-06, 6.44352165E-04, 1.60136150E-02, -5.99901390E-05, -1.27419723E-04, 5.77069110E-04, 1.46943619E-02, -6.44596380E-05, -1.50385410E-06, 6.40497578E-04, 1.60148725E-02,
		-1.32345029E-03, -2.67563048E-03, 1.44268553E-02, 3.86787353E-01, -1.41163378E-03, 2.21852925E-05, 1.60136150E-02, 4.23060570E-01, -1.32098627E-03, -2.69564153E-03, 1.43815298E-02, 3.86539493E-01, -1.39692256E-03, -2.51065778E-05, 1.59632228E-02, 4.23022298E-01,
		6.97902683E-06, 1.54708562E-05, -5.40553500E-05, -1.22389076E-03, 7.64259908E-06, 1.74435667E-07, -5.99901390E-05, -1.32098627E-03, 7.77166853E-06, 1.43271281E-05, -5.41067445E-05, -1.21509538E-03, 7.78748783E-06, -2.02430543E-07, -6.00389820E-05, -1.32029737E-03,
		1.55856920E-05, 3.41879130E-05, -1.14818594E-04, -2.50644780E-03, 1.71995157E-05, -3.24556268E-07, -1.27419723E-04, -2.69564153E-03, 1.43271281E-05, 3.63594218E-05, -1.14906803E-04, -2.51842163E-03, 1.69218396E-05, 5.76629888E-07, -1.27528526E-04, -2.69534993E-03,
		-5.42228378E-05, -1.14140988E-04, 5.19917333E-04, 1.32023176E-02, -5.85642150E-05, 1.04949209E-07, 5.77069110E-04, 1.43815298E-02, -5.41067445E-05, -1.14906803E-04, 5.20793955E-04, 1.32014975E-02, -5.85609345E-05, -1.35814887E-07, 5.78036858E-04, 1.43787778E-02,
		-1.22680312E-03, -2.48731155E-03, 1.32383849E-02, 3.53503035E-01, -1.31013693E-03, 2.08275300E-05, 1.46943619E-02, 3.86539493E-01, -1.21509538E-03, -2.51842163E-03, 1.32014975E-02, 3.53417378E-01, -1.29578110E-03, -2.37579278E-05, 1.46534468E-02, 3.86534025E-01,
		7.64961570E-06, 1.70479019E-05, -5.80832573E-05, -1.29613649E-03, 8.40945240E-06, 1.62823790E-07, -6.44596380E-05, -1.39692256E-03, 7.78748783E-06, 1.69218396E-05, -5.85609345E-05, -1.29578110E-03, 8.54767080E-06, -1.80909369E-07, -6.49881630E-05, -1.39597121E-03,
		2.29722480E-07, -1.20749009E-07, -1.35582154E-06, -2.38536090E-05, 2.57382563E-07, -5.19667650E-07, -1.50385410E-06, -2.51065778E-05, -2.02430543E-07, 5.76629888E-07, -1.35814887E-07, -2.37579278E-05, -1.80909369E-07, 7.22590268E-07, -1.52653147E-07, -2.67716138E-05,
		-6.01754873E-05, -1.26666848E-04, 5.77065465E-04, 1.46542669E-02, -6.49929015E-05, 1.18107113E-07, 6.40497578E-04, 1.59632228E-02, -6.00389820E-05, -1.27528526E-04, 5.78036858E-04, 1.46534468E-02, -6.49881630E-05, -1.52653147E-07, 6.41571030E-04, 1.59601793E-02,
		-1.32406994E-03, -2.67546645E-03, 1.44279853E-02, 3.86756370E-01, -1.41235184E-03, 2.35915335E-05, 1.60148725E-02, 4.23022298E-01, -1.32029737E-03, -2.69534993E-03, 1.43787778E-02, 3.86534025E-01, -1.39597121E-03, -2.67716138E-05, 1.59601793E-02, 4.23024120E-01
	};
	opt_cov.SetMatrixArray(opt_cov_data);

	TMatrixDSymEigen eig_decomp(opt_cov);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(16);
	for (unsigned int i = 0; i < 16; i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	opt_per_gen = eig_decomp.GetEigenVectors() * S;
}
