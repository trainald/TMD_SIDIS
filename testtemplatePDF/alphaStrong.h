//
//  alphaStrong.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/23/22.
//

#ifndef alphaStrong_h
#define alphaStrong_h


#include <math.h>

// orders of the QCD beta function
double get_order_beta(int Nf, int o)
{
    if(o == 0)
        return 11 - 2*Nf/3.0;
    else if (o == 1)
        return 102 - 38*Nf/3.0;
    else
        return 2857/2.0 - 5033/18.0*Nf + 325/54.0*Nf*Nf;
}

// QCD beta function
double get_beta(double a, int Nf, int o)
{
    if (o >= 2) // no more than order 2
        o = 2;
    double beta = 0;
    for (int i = 0; i <= o; i++){
        beta += - get_order_beta(Nf, i);
    }
    return beta*a*a;
}

//double evolve_alphas (double Q20, double alphas_IN, double Q2, int Nf, int order)
//{
////    double XK0, XK1, XK2, XK3;
////    double iteration = 20;
////    double LR = log(Q2/Q20)/iteration;
////    for(int i=0; i<= iteration; i++)
////    {
////        XK0 = LR * get_beta(alphas,Nf,order);
////        XK1 = LR * get_beta(alphas + 0.5 * XK0,Nf,order);
////        XK2 = LR * get_beta(alphas + 0.5 * XK1,Nf,order);
////        XK3 = LR * get_beta(alphas + XK2,Nf,order);
////
////        alphas += (XK0 + 2.0 * XK1 + 2.0 * XK2 + XK3)/6.0;
////    }
//    double alphas = 0;
//    double LR = log(Q2/Q20);
//
//    alphas = alphas_IN/(1 + alphas_IN * get_beta(alphas_IN,Nf,0)*LR/(4*M_PI));
//    std::cout << "alphas = " << alphas << std::endl;
//
//    return alphas;
//}



double beta_0(int Nf)
{
    return  (11.0 - 2.0/3.0*Nf)/4.0/M_PI;
}

double evolve_alphas(double muQ2, int Nf, int order){
    double output=0.0;
    double mZ  = 91.1876;
    double mZ2 = mZ*mZ;
    double alphaSMZ = 0.118;
    double beta0val=beta_0(Nf);
    
    output = alphaSMZ/(1.0 + alphaSMZ * beta0val*log(muQ2/(mZ2)));
    
    return output;

    
    
}


#endif /* alphaStrong_h */
