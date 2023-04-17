//
//  OPEPDF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/27/22.
//

#ifndef OPEPDF_h
#define OPEPDF_h

#include "CollinearPDF.h"
#include "MellinConvolution.h"
#include "QCDSplittingKernels.h"

// OPE ----------------------------------------------------------------
template<int t_PartonFlavour, int t_HadronFlavour>
double EvaluateOPE_bT_1(CollinearPDF<t_PartonFlavour, t_HadronFlavour> collPDF,
                      double x,
                      double mubstar,
                      gsl_integration_workspace* wsp)
{
    double result, abserr;

    auto integrand = make_gsl_function( [&](double xi) {return 2*(1-xi)*collPDF.Evaluate(x/xi)/xi;});
    gsl_integration_qag(integrand, x, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6, wsp, &result, &abserr);


    return CF*( result - M_PI*M_PI/6.0*collPDF.Evaluate(x) )/(4*M_PI);// to be multiplied by alphaS
}


#endif /* OPEPDF_h */
