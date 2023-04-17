//
//  MellinConvolution.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/22/22.
//

#include <gsl/gsl_integration.h>
#include "gsl_RAII_wrapper.h"


#ifndef MellinConvolution_h
#define MellinConvolution_h

#define NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE 1E-5
#define NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE 1E-5
#define NUMERICAL_INTEGRATION_MAX_INTERVALS 100

// Generic definition with static functions
double MellinConvolution (double x, double(*P)(double),double(*f)(double))
{
    double result;
    double abserr;
    double p[2] = {x, 1};
    IntegrationWorkspace wsp(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    
    auto integrand = make_gsl_function( [&](double chi) {return P(chi)*f(x/chi)/chi;});
    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp, &result, &abserr);
    //gsl_integration_qag(integrand, x, 1, epsabs, epsrel, limit, 1, wsp, &result, &abserr);
    
    return result;
}




#endif /* MellinConvolution_h */
