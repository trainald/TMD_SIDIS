//
//  bstar.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 11/4/22.
//

#ifndef bstar_h
#define bstar_h

// bstar definition
double bstar(double bT)
{
    double bmax = 1.0;// GeV-1 // default is 1.0 
    return bT/sqrt(1.0 + pow(bT/bmax,2.0));
}

// mubstar
double mubstar(double bstar)
{
    double b0 = 2.0*exp(-M_EULER);
    
    return b0/bstar;
}


#endif /* bstar_h */
