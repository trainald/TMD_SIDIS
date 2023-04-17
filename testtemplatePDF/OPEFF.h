//
//  OPEFF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/27/22.
//

#ifndef OPEFF_h
#define OPEFF_h

#include "CollinearFF.h"
#include "MellinConvolution.h"
#include "QCDSplittingKernels.h"

// small bT OPE expression
template<int t_HadronFlavour, int t_PartonFlavour>
double EvaluateOPE_bT_1(CollinearFF<t_HadronFlavour, t_PartonFlavour> collFF,
                      double z,
                      double mubstar,
                      gsl_integration_workspace* wsp)
{
    double result, abserr;

    auto integrand = make_gsl_function( [&](double xi) {return (4.0*(2.0/(1.0-xi)+1.0/(xi*xi)+1.0/xi)*log(xi) + 2.0/(xi*xi) - 2.0/xi)*collFF.Evaluate(z/xi)/xi;});
    gsl_integration_qag(integrand, z, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_ITERATIONS, 6, wsp, &result, &abserr);


    return CF*( result - M_PI*M_PI/6.0*collFF.Evaluate(z) )/(4*M_PI);// to be multiplied by alphaS
}

#endif /* OPEFF_h */
