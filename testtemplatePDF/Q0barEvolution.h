//
//  Q0barEvolution.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 2/3/23.
//

#ifndef Q0barEvolution_h
#define Q0barEvolution_h

//#include "TMDPDF.h"
//#include "TMDFF.h"
#include "QCDParameters.h"
#


// Define Q0 bar
double Q0bar(double bT)
{
    double result = INPUT_HARD_SCALE*( 1.0 - ( 1.0 - 2.0*exp(-M_EULER)/(INPUT_HARD_SCALE*bT))*exp(-pow(bT*INPUT_HARD_SCALE,2)) );
//    double result = INPUT_HARD_SCALE*( 1.0 + 2.0*exp(-M_EULER)*exp(-bT*INPUT_HARD_SCALE)/(INPUT_HARD_SCALE*bT) );
    
    return result;
}


//// Define evolution from Q0 bar to Q0
//double E_Q0bar2Q0(double bT)
//{
//    double Q0_bar = Q0bar(bT);
//
//    double integral_result;
//    double abserr;
//    double p[2] = {Q0_bar, INPUT_HARD_SCALE};
//    IntegrationWorkspace wsp(NUMERICAL_INTEGRATION_MAX_INTERVALS);
//
//    auto integrand = make_gsl_function( [&](double mu) {return (gamma2Parton(evolve_alphas(mu*mu, 3, 0)) -log(INPUT_HARD_SCALE/mu)*gammaCSKernel(evolve_alphas(mu*mu, 3, 0)))/mu;});
////    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp, &integral_result, &abserr);
//    gsl_integration_qag(integrand, Q0_bar, INPUT_HARD_SCALE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 1, wsp, &integral_result, &abserr);
//
//
//    double output = integral_result + log(INPUT_HARD_SCALE/Q0_bar)*CSKernel(bT, Q0_bar, evolve_alphas(pow(Q0_bar,2), 3, 0), true);
//
//    return exp(output);
//}



#endif /* Q0barEvolution_h */
