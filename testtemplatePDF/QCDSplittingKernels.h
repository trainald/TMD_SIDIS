//
//  QCDSplittingKernels.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/22/22.
//

#include "QCDParameters.h"

#ifndef QCDSplittingKernels_h
#define QCDSplittingKernels_h


double Pg_q(double z)
{
    return CF * ( 1.0 + pow((1.0 - z),2) )/z;
}

double Pq_g(double x)
{
    return TF * ( x*x + pow((1.0 - x),2) );
}


// MSbar to Cutoff
double CDelta_iH_tilde(double x)
{
    return CF * (1.0 - x);// only non delta part
}

double CDelta_gH(double x)
{
    return 2.0 * TF * x * (1.0 - x);
}

double CDelta_Hi_tilde(double z)
{
    return CF*( 2.0*(1.0 + z*z)*log(z) + pow((1.0 - z),2) )/(1.0 - z);
}

double CDelta_Hg(double z)
{
    return 2.0*Pg_q(z)*log(z) + CF*z;
}

#endif /* QCD_SplittingKernels_h */

