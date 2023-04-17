//
//  YtermFunctions.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 11/15/22.
//

#include "CollinearPDF.h"
#include "CollinearFF.h"

#ifndef YtermFunctions_h
#define YtermFunctions_h

// Compute Fixed Order
//template<int t_PartonFlavour, int t_BosonFlavour, int t_HadronFlavour>
//double FixedOrder(class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
//                  class CollinearFF<t_HadronFlavour, t_PartonFlavour> colFF,
//                  )


// Nadolsky et all notation

// Unpolarized beam and target --------------
double A1(double y) // cosh(psi) = 2/y - 1
{
    return 2.0*(1.0 + 2.0*(1-y)/(y*y));
}

double A2()
{
    return -2.0;
}

// From https://arxiv.org/pdf/hep-ph/0602188.pdf
// e + p --> e + pi + X

// [f•D•sigma_k_hat] = sum_qqbar

// qq sector
// 1
double sigma1_hat_qq(double x_hat, double z_hat, double Q, double qT)
{
    return 2.0*CF*x_hat*z_hat*( 6.0 + ( pow(Q,4)/(pow(x_hat,2)*pow(z_hat,2)) + pow(Q*Q - qT*qT,2) )/(pow(Q,2)*pow(qT,2)) );
}
// 2
double sigma2_hat_qq(double x_hat, double z_hat)
{
    return 8.0*CF*x_hat*z_hat;
}
// 3 POLARIZED CASE
double sigma3_hat_qq(double x_hat, double z_hat, double Q, double qT)
{
    return 4.0*CF*x_hat*z_hat*(Q*Q + qT*qT)/(Q*qT);
}
// 4 POLARIZED CASE
double sigma4_hat_qq(double x_hat, double z_hat)
{
    return 4.0*CF*x_hat*z_hat;
}

// qg sector
// 1
double sigma1_hat_qg(double x_hat, double z_hat, double Q, double qT)
{
    return x_hat*(1.0 - x_hat)*( pow(Q/qT,2)*( 1.0/pow(x_hat*z_hat,2) - 2.0/(x_hat*z_hat) + 2.0 ) + 10.0 - 2.0/x_hat - 2.0/z_hat );
}
// 2
double sigma2_hat_qg(double x_hat)
{
    return 8.0*x_hat*(1.0 - x_hat);
}
// 3 POLARIZED CASE
double sigma3_hat_qg(double x_hat, double z_hat, double Q, double qT)
{
    return x_hat*(1.0 - x_hat)*2.0/(Q*qT)*( 2.0*(Q*Q + qT*qT) - Q*Q/(x_hat*z_hat) );
}
// 4 POLARIZED CASE
double sigma4_hat_qg(double x_hat)
{
    return 4.0*x_hat*(1.0 - x_hat);
}

// gq sector
// 1
double sigma1_hat_gq(double x_hat, double z_hat, double Q, double qT)
{
    return 2.0*CF*x_hat*(1.0 - z_hat)*( 1.0/pow(Q*qT,2)*( pow(Q,4)/pow(x_hat*z_hat,2) + pow(1.0 - z_hat,2)/pow(z_hat,2)*pow( Q*Q - pow(z_hat*qT/(1.0 - z_hat),2) ,2) ) + 6.0 );
}
// 2
double sigma2_hat_gq(double x_hat, double z_hat)
{
    return 8.0*CF*x_hat*(1.0 - z_hat);
}
// 3 POLARIZED CASE
double sigma3_hat_gq(double x_hat, double z_hat, double Q, double qT)
{
    return -4.0*CF*x_hat*pow(1.0 - z_hat,2)/(z_hat*Q*qT)*( Q*Q + pow(z_hat*qT/(1.0 - z_hat),2) );
    
}
// 4 POLARIZED CASE
double sigma4_hat_gq(double x_hat, double z_hat)
{
    return 4.0*CF*x_hat*(1.0 - z_hat);
}

// ---------------------------------------------------------------------------------------------------
// Standard procedure for FIXED ORDER

// 1
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double colPDF_colFF_Partonic_Xsec1(class CollinearPDF<t_QuarkFlavour, t_HadronFlavour> QuarkColPDF,
                                   class CollinearPDF<t_BosonFlavour, t_HadronFlavour> BosonColPDF,
                                   class CollinearFF<t_HadronFlavour, t_QuarkFlavour> QuarkColFF,
                                   class CollinearFF<t_HadronFlavour, t_BosonFlavour> BosonColFF,
                                   double x, double xi, double z, double zeta, double qT)
{
    double QuarkCharge = QuarkColPDF.Get_PartonCharge();
    
    double qq_term = QuarkColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma1_hat_qq(x/xi, z/zeta, QuarkColPDF.Get_HardScale(), qT);
    double qg_term = BosonColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma1_hat_qg(x/xi, z/zeta, BosonColPDF.Get_HardScale(), qT);
    double gq_term = QuarkColPDF.Evaluate(xi)*BosonColFF.Evaluate(zeta)*sigma1_hat_gq(x/xi, z/zeta, QuarkColPDF.Get_HardScale(), qT);
    
    return pow(QuarkCharge,2)*( qq_term + qg_term + gq_term);
    
}

// 2
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double colPDF_colFF_Partonic_Xsec2(class CollinearPDF<t_QuarkFlavour, t_HadronFlavour> QuarkColPDF,
                                   class CollinearPDF<t_BosonFlavour, t_HadronFlavour> BosonColPDF,
                                   class CollinearFF<t_HadronFlavour, t_QuarkFlavour> QuarkColFF,
                                   class CollinearFF<t_HadronFlavour, t_BosonFlavour> BosonColFF,
                                   double x, double xi, double z, double zeta, double qT)
{
    double QuarkCharge = QuarkColPDF.Get_PartonCharge();
    
    double qq_term = QuarkColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma2_hat_qq(x/xi, z/zeta);
    double qg_term = BosonColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma2_hat_qg(x/xi);
    double gq_term = QuarkColPDF.Evaluate(xi)*BosonColFF.Evaluate(zeta)*sigma2_hat_gq(x/xi, z/zeta);
    
    return pow(QuarkCharge,2)*( qq_term + qg_term + gq_term);
    
}

// Solution to the delta
double zeta_bar(double x, double z, double xi, double Q, double qT)
{
    return z*( 1.0 - qT*qT*x/(Q*Q*(x - xi)) );
}
double xi_bar(double x, double z, double zeta, double Q, double qT)
{
    return x*( 1.0 - qT*qT*z/(Q*Q*(z - zeta)) );
}

// Equations A.4 From https://arxiv.org/pdf/hep-ph/0602188.pdf
double x_min(double x, double z, double Q, double qT)
{
    return x*(1.0 + z/(1.0 - z)*pow(qT/Q,2));
}
double z_min(double x, double z, double Q, double qT)
{
    return z*(1.0 + x/(1.0 - x)*pow(qT/Q,2));
}

// Fixed order
// Standard procedure
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double Get_TMDFixedOrder_unpolarized(class CollinearPDF<t_QuarkFlavour, t_HadronFlavour> QuarkColPDF,
                                 class CollinearPDF<t_BosonFlavour, t_HadronFlavour> BosonColPDF,
                                 class CollinearFF<t_HadronFlavour, t_QuarkFlavour> QuarkColFF,
                                 class CollinearFF<t_HadronFlavour, t_BosonFlavour> BosonColFF,
                                 double x, double z, double qT, double y,
                                 gsl_integration_workspace * wsp)// to be multiplied by alphaS*alphaEM^2
{
    double result, abserr;
    
    double Q = QuarkColPDF.Get_HardScale();
    
    double p[2] = {x_min(x, z, Q, qT), 1.0};
    
    auto integrand = make_gsl_function( [&](double xi) {
        return x*z*y/(xi-x)/(xi*zeta_bar(x, z, xi, Q, qT)) *
        ( A1(y) * colPDF_colFF_Partonic_Xsec1(QuarkColPDF, BosonColPDF, QuarkColFF, BosonColFF, x, xi, z, zeta_bar(x, z, xi, Q, qT), qT) +
         A2() * colPDF_colFF_Partonic_Xsec2(QuarkColPDF, BosonColPDF, QuarkColFF, BosonColFF, x, xi, z, zeta_bar(x, z, xi, Q, qT), qT) );
        
    });
    gsl_integration_qagp(integrand, p, 2,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS,wsp,
                         &result, &abserr);
    
    return result/(4.0*pow(Q,4));
}

//// Hadronic Structure Oriented (HSO) procedure for FIXED ORDER
//
//// 1
//template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
//double colPDF_colFF_Partonic_Xsec1_HSO(class TMDPDF<t_QuarkFlavour, t_BosonFlavour ,t_HadronFlavour> QuarkColPDF,
//                                       class CollinearPDF<t_BosonFlavour, t_HadronFlavour> BosonColPDF,
//                                       class CollinearFF<t_HadronFlavour, t_QuarkFlavour> QuarkColFF,
//                                       class CollinearFF<t_HadronFlavour, t_BosonFlavour> BosonColFF,
//                                       double x, double xi, double z, double zeta, double qT)
//{
//    double QuarkCharge = QuarkColPDF.Get_PartonCharge();
//    
//    double qq_term = QuarkColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma1_hat_qq(x/xi, z/zeta, QuarkColPDF.Get_HardScale(), qT);
//    double qg_term = BosonColPDF.Evaluate(xi)*QuarkColFF.Evaluate(zeta)*sigma1_hat_qg(x/xi, z/zeta, BosonColPDF.Get_HardScale(), qT);
//    double gq_term = QuarkColPDF.Evaluate(xi)*BosonColFF.Evaluate(zeta)*sigma1_hat_gq(x/xi, z/zeta, QuarkColPDF.Get_HardScale(), qT);
//    
//    return QuarkCharge*QuarkCharge*( qq_term + qg_term + gq_term);
//    
//}





#endif /* YtermFunctions_h */
