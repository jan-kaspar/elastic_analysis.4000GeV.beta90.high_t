#include "../DS4/parameters.h"

#define USE_INIT_ADDITIONAL 1

void Init_additional()
{
	// select subset of bunches
	bunchMap.clear();
	bunchMap[8369].push_back(2990);
	bunchMap[8371].push_back(2990);
	bunchMap[8372].push_back(2990);

	// constant-alignment parameters
	AlignmentSource alSrc3;
	alSrc3.SetAlignmentA(atConstant);
	alSrc3.SetAlignmentB(atConstant);
	alSrc3.SetAlignmentC(atConstant);
	alSrc3.cnst.a_L_F = 0E-3; alSrc3.cnst.b_L_F = +1.0E-3; alSrc3.cnst.c_L_F = 0E-3;
	alSrc3.cnst.a_L_N = 0E-3; alSrc3.cnst.b_L_N = +2.0E-3; alSrc3.cnst.c_L_N = 0E-3;
	alSrc3.cnst.a_R_N = 0E-3; alSrc3.cnst.b_R_N = +1.4E-3; alSrc3.cnst.c_R_N = 0E-3;
	alSrc3.cnst.a_R_F = 0E-3; alSrc3.cnst.b_R_F = +1.5E-3; alSrc3.cnst.c_R_F = 0E-3;
	alignmentSources.push_back(alSrc3);

	anal.L_int = 383.924E3;	// scaled to match with DS2 (equal dsigma/dt integral over "ub" bins 5 to 12)

	if (diagonal == d45b_56t)
	{
		anal.cut1_si = 9.6E-6;
		anal.cut2_si = 3.2E-6;

		anal.cut7_si = 8.9E-3;
	}

	if (diagonal == d45t_56b)
	{
		anal.cut1_si = 8.9E-6;
		anal.cut2_si = 3.2E-6;

		anal.cut7_si = 8.5E-3;
	}
}
