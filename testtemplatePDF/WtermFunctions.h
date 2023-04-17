//
//  WtermFunctions.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/24/22.
//

#ifndef WtermFunctions_h
#define WtermFunctions_h

#include "TMDPDF.h"
#include "TMDFF.h"
#include "Evolution.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

// Momentum space -------------------------------------------------------------------------------------------------
// Integrand of exact W term in kT space
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double WtermSingleFlavour_integrand(class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf,
                                    class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff,
                                    double qT, double kT, double phi)
{
    double x = InputTMDpdf.Get_xBjorken();
    double z = InputTMDff.Get_z();
    
    double sqrt_p = sqrt(kT*kT + qT*qT/4.0 + kT*qT*cos(phi));
    double sqrt_m = sqrt(kT*kT + qT*qT/4.0 - kT*qT*cos(phi));
    
    return InputTMDpdf.Evaluate(x, sqrt_m) * InputTMDff.Evaluate(z, sqrt_p);
}

// Compute integral of exact W term over kT
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_exact(class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf,
                                    class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff,
                                    double qT,
                                    gsl_integration_workspace * wsp1,
                                    gsl_integration_workspace * wsp2)
{
    double inner_result, inner_abserr;
    double result, abserr;
    
    auto outer = make_gsl_function( [&](double phi) {
    auto inner = make_gsl_function( [&](double kT) {return kT*WtermSingleFlavour_integrand(InputTMDpdf, InputTMDff, qT, kT, phi);});
    gsl_integration_qagiu(inner, 0.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp1,
                         &inner_result, &inner_abserr);
//    gsl_integration_qag(inner, 0.0, 1E+6 ,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 4,wsp1, &inner_result, &inner_abserr);
    return inner_result;
      } );
    gsl_integration_qag(outer, 0.0, 2.0*M_PI ,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 4,wsp2, &result, &abserr);

    
    return result;
}

// Compute integral of W term LO Asymptotic
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_AsyLO(class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf,
                                    class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf_pert,
                                    class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff,
                                    class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff_pert,
                                    double qT)
{
    double z      = InputTMDff.Get_z();
    double x      = InputTMDpdf.Get_xBjorken();
    double kc     = InputTMDff.Get_RenormalizationScale();// kc = mu
    
    return InputTMDpdf_pert.Evaluate(x, -qT) * InputTMDff.Evaluate_CutoffIntegral(z, kc)/(z*z) + InputTMDff_pert.Evaluate(z, qT) * InputTMDpdf.Evaluate_CutoffIntegral(x, kc);
    
}


// Structure containing parameters and functions for VEGAS inputs parameters
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
struct my_Vegas_parameters {class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDPDF_pert;
                            class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDFF_pert;
                            double qT;};


// --------------------------------
double HeavisideTheta(double x){
    if (x > 0)
        return 1.0;
    else
        return 0.0;
    }

//----------------------------------

template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double WtermSingleFlavour_Correction_integrand(double *k, size_t dim, void * VEGAS_params)
{
    my_Vegas_parameters<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> * fV = (my_Vegas_parameters<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> *)VEGAS_params;
    double qT = fV->qT;
    
    class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDPDF_pert = fV->InputTMDPDF_pert;
    class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDFF_pert = fV->InputTMDFF_pert;
    
    double kT  = (1-k[0])/k[0];// transformation from (0,inf) to (0,1)
    double phi = k[1];

    double muQ = InputTMDFF_pert.Get_RenormalizationScale();
    double x   = InputTMDPDF_pert.Get_xBjorken();
    double z   = InputTMDFF_pert.Get_z();

    double sqrt_p = sqrt(kT*kT + qT*qT/4.0 + kT*qT*cos(phi));
    double sqrt_m = sqrt(kT*kT + qT*qT/4.0 - kT*qT*cos(phi));
    
    double first  = InputTMDPDF_pert.Evaluate(x,sqrt_m)*InputTMDFF_pert.Evaluate(z,sqrt_p);
    double second = InputTMDPDF_pert.Evaluate(x,sqrt_m)*InputTMDFF_pert.Evaluate(z,qT) *HeavisideTheta(muQ - sqrt_m);
    double third  = InputTMDFF_pert.Evaluate(z,sqrt_p)*InputTMDPDF_pert.Evaluate(x,-qT)*HeavisideTheta(muQ - sqrt_p);
    
    return kT*(first - second - third)/pow(k[0],2);
}


// Correction integral
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_Correction(const gsl_rng_type * T,
                                         gsl_rng * r,
                                         my_Vegas_parameters<t_QuarkFlavour,t_BosonFlavour ,t_HadronFlavour> VEGAS_params,
                                         gsl_monte_vegas_state *s)
{
    double VEGAS_res;
    double res, err;
    gsl_monte_function F;
    F.f = &WtermSingleFlavour_Correction_integrand<t_QuarkFlavour,t_BosonFlavour ,t_HadronFlavour>;//integrand
    F.dim = 2; // dimension
    F.params = &VEGAS_params;// parameters

    double xl[2] = { 0.0, 0.0 };// lower limits of integration
    double xu[2] = { 1.0, 2*M_PI };// upper limits of integration

    size_t calls = 100000;// Monte Carlo calls

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    //do vegas integration
    gsl_monte_vegas_integrate (&F, xl, xu, 2, 10000, r, s, &res, &err);
    //printf ("converging...\n");
    do
      {
        gsl_monte_vegas_integrate (&F, xl, xu, 2, calls/5, r, s,
                                   &res, &err);
//            printf ("result = % .6f sigma = % .6f "
//                    "chisq/dof = %.1f\n", res, err, s->chisq);
          VEGAS_res = res;
      }
    while (fabs (s->chisq - 1.0) > 0.5);
    
    return VEGAS_res;
}


// Coordinate space -------------------------------------------------------------------------------------------------

// Integrand of W term in bT space
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double WtermSingleFlavour_integrand_bT(class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf,
                                       class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff,
                                       double bT)
{
    double z = InputTMDff.Get_z();
    double x = InputTMDpdf.Get_xBjorken();
    
    return InputTMDff.Evaluate_inbTspace(z, bT) * InputTMDpdf.Evaluate_inbTspace(x, bT);
}

// W term from bT space expression
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_exact_bT(class TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> InputTMDpdf,
                                       class TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> InputTMDff,
                                       double qT,
                                       gsl_integration_workspace * wsp)
{
    double result, abserr;
    
    auto integrand = make_gsl_function( [&](double bT) {return
        gsl_sf_bessel_J0(qT*bT) * bT * WtermSingleFlavour_integrand_bT(InputTMDpdf, InputTMDff, bT)/(2.0*M_PI);});
    // The domain should be (0,+inf)
    gsl_integration_qag(integrand, 0.0, 3500.0,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6,wsp,
                         &result, &abserr);
    
    return result;
    
}




// old OPE way of calculating W term

// integrand
template<int t_PartonFlavour, int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_OPEbstar_integrand(double bT,
                                                 double mustr,
                                       class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                                       class TMDPDF<t_PartonFlavour,t_BosonFlavour,t_HadronFlavour> tmdPDF_pert,
                                       class CollinearFF<t_HadronFlavour, t_PartonFlavour> colFF,
                                       class TMDFF<t_HadronFlavour, t_PartonFlavour, t_BosonFlavour> tmdFF_pert,
                                       gsl_integration_workspace * wsp)
{
    double muQ = tmdPDF_pert.Get_RenormalizationScale();
    
    double x = tmdPDF_pert.Get_xBjorken();
    double z = tmdFF_pert.Get_z();
    
    double output = tmdPDF_pert.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp) * tmdFF_pert.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp);
    
    return output;
}



// integral
template<int t_PartonFlavour, int t_BosonFlavour, int t_HadronFlavour>
double Get_WtermSingleFlavour_OPEbstar(double qT,
                                       class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                                       class TMDPDF<t_PartonFlavour,t_BosonFlavour,t_HadronFlavour> tmdPDF_pert,
                                       class CollinearFF<t_HadronFlavour, t_PartonFlavour> colFF,
                                       class TMDFF<t_HadronFlavour, t_PartonFlavour, t_BosonFlavour> tmdFF_pert,
                                       gsl_integration_workspace * wsp1,
                                       gsl_integration_workspace * wsp2)
{
    double result, abserr;

    auto integrand = make_gsl_function( [&](double bT) {
        
//        double mustar = mubstar(bstar(bT));
        return gsl_sf_bessel_J0(qT*bT) * bT * Get_WtermSingleFlavour_OPEbstar_integrand(bT,mubstar(bstar(bT)), colPDF, tmdPDF_pert, colFF, tmdFF_pert, wsp1)/(2.0*M_PI);});
     // The domain should be (0,+inf)
//    gsl_integration_qag(integrand, 0.0, 3600.0,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6,wsp2,
//                         &result, &abserr);
    gsl_integration_qagiu(integrand, 0.0,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS,wsp2,
                         &result, &abserr);
    
    return result;
}



#endif /* WtermFunctions_h */
