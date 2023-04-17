//
//  Evolution.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/28/22.
//

#include "QCDParameters.h"
#include "CollinearPDF.h"
#include "CollinearFF.h"
#include "bstar.h"

#ifndef Evolution_h
#define Evolution_h


// OPE RG evolution
template<int t_PartonFlavour, int t_BosonFlavour, int t_HadronFlavour>
double OPE_RG_Evolution_bstar(double bT,
                              double mustr,
                              double mu_f,
                              double Q,
                              class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                              class TMDPDF<t_PartonFlavour,t_BosonFlavour,t_HadronFlavour> tmdPDF,
                              class CollinearFF<t_HadronFlavour, t_PartonFlavour> colFF,
                              class TMDFF<t_HadronFlavour, t_PartonFlavour, t_BosonFlavour> tmdFF,
                              gsl_integration_workspace * wsp)
{
    double bstr = bstar(bT);
    
    return colPDF.RG_Evolution_bT(bstr, mustr, mu_f, Q, colPDF, wsp) * exp(-tmdPDF.Evaluate_g_QuarkHadron(bT) - tmdFF.Evaluate_g_HadronQuark(bT) - 2*log(Q/INPUT_HARD_SCALE)*gCSKernel(bT, mu_f, colPDF.Get_AlphaStrong(mu_f), tmdPDF.Get_NonPerturbative()));
//    return colPDF.RG_Evolution_bT(bstr, mustr, mu_f, Q, colPDF, wsp) * exp(-tmdPDF.Evaluate_g_QuarkHadron(bT) - tmdFF.Evaluate_g_HadronQuark(bT) );// add non input scale term
}



#endif /* Evolution_h */
