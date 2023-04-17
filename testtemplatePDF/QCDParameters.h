//
//  QCDParameters.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/28/22.
//

//#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include "bstar.h"
#include "MassParameters.h"

#ifndef QCDParameters_h
#define QCDParameters_h

#define INPUT_HARD_SCALE 4.0

const int Nc = 3;// Number of colors

const double CF = (Nc*Nc - 1.)/(2. * Nc);// SU(Nc) Casimir - Fundamental representation
const double CA = Nc;// SU(Nc) Casimir - Adjoint representation
const double TF = 0.5;// SU(Nc) Dynkin index

//
double gamma2Parton(double alphaS)
{
    return 6.0*CF*alphaS/(4.0*M_PI);// this is only the fist order
}

double gammaCSKernel(double alphaS)
{
    return 8.0*CF*alphaS/(4.0*M_PI);// this is only the fist order
}

double CSKernel(double bT, double muQ, double alphaS, bool NonPerturbative)
{
    double output;
    
    if(NonPerturbative)
    {
        // HSO
        output = alphaS*2.0*CF/M_PI*( gsl_sf_bessel_K0(bT*mass_CSKernel) + log(mass_CSKernel/muQ));
    }
    else
    {
        // OPE or perturbative HSO
        output = -alphaS * 8.0 * CF * log( bT * muQ * exp(M_EULER)/2.0 );// Choose this for OPE
    }
    return output;
}

double gCSKernel(double bT, double muQ, double alphaS, bool NonPerturbative, bool isAnsatz)
{
    if(!isAnsatz)
    {
        // HSO parametrization
         return CSKernel(bstar(bT), muQ, alphaS, NonPerturbative) - CSKernel(bT, muQ, alphaS, NonPerturbative);
    }
    else
    {
        // OPE parametrization
        double M_K = 0.1;
        double g2 = (2.0*pow(M_K,2.0));
        
        return g2/(2.0*pow(M_K,2.0)) * log( 1.0 + pow(M_K*bT,2.0) );
    }
}



#endif /* QCDParameters_h */
