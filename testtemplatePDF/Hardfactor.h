//
//  Hardfactor.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 11/9/22.
//

#ifndef Hardfactor_h
#define Hardfactor_h

double Hardfactor(double Q, double muQ, double alphaS)
{
    return ( 1 + alphaS*CF/(4*M_PI)*(-16.0 + M_PI*M_PI/3.0 + 6.0*log(Q*Q/(muQ*muQ)) -2.0*pow(log(pow(Q,2)/pow(muQ,2)),2 )) );
}


#endif /* Hardfactor_h */
