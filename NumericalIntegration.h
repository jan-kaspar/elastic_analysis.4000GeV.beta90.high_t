#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#ifndef _numerical_integration_h_
#define _numerical_integration_h_

//----------------------------------------------------------------------------------------------------

typedef double (* RealFunction)(double x, double *par, const void *obj);

//----------------------------------------------------------------------------------------------------

struct RealIntegPar
{
	RealFunction fcn;
	double *parameters;
	const void *object;

	RealIntegPar(RealFunction _f, double *_p, const void *_o) : fcn(_f), parameters(_p), object(_o) {}
};

//----------------------------------------------------------------------------------------------------

double RealIntegFcn(double x, void *vpar)
{
	RealIntegPar *par = (RealIntegPar *) vpar;
	return par->fcn(x, par->parameters, par->object);
}

//----------------------------------------------------------------------------------------------------

double RealIntegrate(RealFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char *errorLabel)
{
	// prepare structures
	RealIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = RealIntegFcn;
  	F.params = &ocip;

	// real part
	double result, unc;
	int status = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS15, work_space, &result, &unc);
	if (status != 0) 
	{
		printf("WARNING in %s > Integration failed: %s.\n", errorLabel, gsl_strerror(status));
	}

	return result;
}

#endif
