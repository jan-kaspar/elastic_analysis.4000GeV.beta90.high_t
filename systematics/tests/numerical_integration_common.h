#include "common.h"

//----------------------------------------------------------------------------------------------------

double DoubleInt(double (*fcn)(double*, double*), double a, double b, double *params, double epsilon)
{	
	const double Z1 = 1;
	const double HF = Z1/2;
	const double CST = 5*Z1/1000;
	
	double x[12] = { 0.96028985649753623,	0.79666647741362674,
										 0.52553240991632899,	0.18343464249564980,
										 0.98940093499164993,	0.94457502307323258,
										 0.86563120238783174,	0.75540440835500303,
										 0.61787624440264375,	0.45801677765722739,
										 0.28160355077925891,	0.09501250983763744};
	
	double w[12] = { 0.10122853629037626,	0.22238103445337447,
										 0.31370664587788729,	0.36268378337836198,
										 0.02715245941175409,	0.06225352393864789,
										 0.09515851168249278,	0.12462897125553387,
										 0.14959598881657673,	0.16915651939500254,
										 0.18260341504492359,	0.18945061045506850};
	
	double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
	double xx[1];
	int i;
	
	h = 0;
	if (b == a) return h;
	
	aconst = CST/TMath::Abs(b-a);
	bb = a;
	
		CASE1:
	aa = bb;
	bb = b;
	
		CASE2:
	c1 = HF*(bb+aa);
	c2 = HF*(bb-aa);
	
	s8 = 0;
	for (i=0; i<4; i++) {
		u		 = c2*x[i];
		xx[0] = c1+u;
		f1		= (*fcn)(xx, params);
		xx[0] = c1-u;
		f2		= (*fcn)(xx, params);
		s8	 += w[i]*(f1 + f2);
		}
	
	s16 = 0;
	for (i=4; i<12; i++) {
		u		 = c2*x[i];
		xx[0] = c1+u;
		f1		= (*fcn)(xx, params);
		xx[0] = c1-u;
		f2		= (*fcn)(xx, params);
		s16	+= w[i]*(f1 + f2);
		}
	
	s16 = c2*s16;
	
	if (TMath::Abs(s16-c2*s8) <= epsilon*(1. + TMath::Abs(s16))) {
			h += s16;
			if(bb != b) goto CASE1;
		} else {
			bb = c1;
			if(1. + aconst*TMath::Abs(c2) != 1) goto CASE2;
			h = s8;	//this is a crude approximation (cernlib function returned 0 !)
			printf("WARNING in DoubleInt > Crude approximation.\n");
		}
		
	return h;
}

//----------------------------------------------------------------------------------------------------

string TDistTypeName(unsigned int t)
{
	if (t == 0)
		return "exp1";

	if (t == 1)
		return "exp2";

	if (t == 2)
		return "data fit 1";
	
	return "unknown";
}

//----------------------------------------------------------------------------------------------------

unsigned int t_dist_type;

double dist_true_t(double t)
{
	if (t_dist_type == 0)
	{
		double B = 20.;
		return exp(-B * t);
	}

	if (t_dist_type == 1)
	{
		double B1 = -19.89, B2 = 2.65;
		return exp(B1*t + B2*t*t);
	}

	if (t_dist_type == 2)
		return spline_dist_t_data_fit1->Eval(t);

	return 0.;
}
