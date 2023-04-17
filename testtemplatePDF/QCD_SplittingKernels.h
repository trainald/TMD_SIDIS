//
//  QCD_SplittingKernels.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/22/22.
//

#ifndef QCD_SplittingKernels_h
#define QCD_SplittingKernels_h

const int Nc = 3;

const double CF = (Nc*Nc - 1.)/(2. * Nc);
const double CA = Nc;
const double TF = 0.5;

double Pg_q(double z)
{
    return CF * ( 1. + pow((1. - z),2) )/z;
}

double Pq_g(double x)
{
    return TF * ( x*x + pow((1. - x),2) );
}

double CDelta_i_p_tilde(double x)
{
    return CF * (1.-x);// only non delta part
}

double CDelta_g_p(double x)
{
    return 2. * TF * x * (1.-x);
}

#endif /* QCD_SplittingKernels_h */
