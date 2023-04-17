//
//  CollinearPDF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/21/22.
//
#include <math.h>
#include <LHAPDF/LHAPDF.h>

#include "Flavour.h"
#include "PDFBase.h"
#include "alphaStrong.h"
#include "MellinConvolution.h"
#include "QCDParameters.h"



#ifndef CollinearPDF_h
#define CollinearPDF_h

template<int t_PartonFlavour, int t_HadronFlavour>
class CollinearPDF : public PDFBase
{
public:
    
    // Constructors -----------------------------------------------------------
    CollinearPDF(double Q, double muQ, bool isLHAPDF)
    :   PDFBase(Q, muQ, t_PartonFlavour, t_HadronFlavour),
        m_isLHAPDF(isLHAPDF),
        m_HardScale(Q),
        m_RenormalizationScale(muQ)
    {
        // Import LHAPDF pdf
        m_lhapdf_pdf = LHAPDF::mkPDF(m_setname, m_setmember);
    };
    
    // Methods ----------------------------------------------------------------
    double Evaluate(double x) const;
    
    // Evolution exponent for PDF
    double RG_Evolution_gammas_exponent(double mu_i,
                                        double mu_f,
                                        double Q,
                                        class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                                        gsl_integration_workspace * wsp);
    // Evolution exponential
    double RG_Evolution_bT(double bT,
                           double mu_i,
                           double mu_f,
                           double Q,
                           class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                           gsl_integration_workspace * wsp,
                           bool NonPerturbative);
    
    // Mellin Convolution type
    double MellinConvolutionPDF(double x, double(*P)(double),
                                double(CollinearPDF<t_PartonFlavour,t_HadronFlavour>::*Evaluate)(double) const);
    
    // Get ---------------------------------------
    double Get_AlphaStrong(double mu);// to be defined from LHAPDF data
    
    double Get_PartonCharge();
    
    // Set ---------------------------------------
    void Set_LHAPDF_Setname_and_member(const std::string LHAPDFSET_PDF, int member)
    {
        m_setname = LHAPDFSET_PDF;
        m_setmember = member;
    }
    
    
private:
    
    /* \brief Hard Scale*/
    double m_HardScale;
    
    /* \brief Renormalization Scale*/
    double m_RenormalizationScale;
    
    /* \brief t_PartonFlavour*/
    int m_PartonFlavour = t_PartonFlavour;
    
    /* \brief Charge of t_PartonFlavour*/
    double m_PartonCharge = Get_PartonCharge();
    
    /* \brief Flag for LHAPDF data */
    bool m_isLHAPDF;
    
    /* \brief LHAPDF set name */
    const std::string m_setname = "cteq66";//default (can be changed)
//    const std::string m_setname = "MMHT2014nlo68cl";//default (can be changed)
    /* \brief LHAPDF set member */
    int m_setmember = 0;//default (can be changed)
    
    /* LHAPDF pointer to PDF set*/
    LHAPDF::PDF* m_lhapdf_pdf;

};// END of Collinear PDF class



// Collinear PDF for quarks and antiquarks
// Default collinear PDF
template<>
double CollinearPDF<UNDEFINED_PARTON,UNDEFINED_HADRON>::Evaluate(double x) const
{
    return x * (1.0 - x);
}


// UP
template<>
double CollinearPDF<UP, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(UP,x,Q2)/x;
    }
    else
    {
        return pow(x,2) * (1.0 - x);
    }
}

// DOWN
template<>
double CollinearPDF<Flavour::DOWN, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(DOWN,x,Q2)/x;
    }
    else
    {
        return x*pow((1.0 - x), 2);
    }
    
}


// STRANGE
template<>
double CollinearPDF<Flavour::STRANGE, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(STRANGE,x,Q2)/x;
    }
    else
    {
        return pow(x,3)*pow((1.0 - x), 2);
    }
}

// Antiquarks
// UBAR
template<>
double CollinearPDF<UBAR, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(UBAR,x,Q2)/x;
    }
    else
    {
        return pow(x,2) * (1.0 - x);
    }
}

// DBAR
template<>
double CollinearPDF<Flavour::DBAR, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(DBAR,x,Q2)/x;
    }
    else
    {
        return x*pow((1.0 - x), 2);
    }
    
}


// STRANGE
template<>
double CollinearPDF<Flavour::SBAR, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(SBAR,x,Q2)/x;
    }
    else
    {
        return pow(x,3)*pow((1.0 - x), 2);
    }
}


// Collinear PDF for bosons

template<>
double CollinearPDF<Flavour::GLUON, Flavour::UNDEFINED_HADRON>::Evaluate(double x) const
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_pdf->xfxQ2(GLUON,x,Q2)/x;
    }
    else
    {
        return pow(x,2)*pow(1.0 - x, 2);
    }
}

// --------------------------------------------------------------------------------------------

// Evolution exponent for PDF
template<int t_PartonFlavour, int t_HadronFlavour>
double CollinearPDF<t_PartonFlavour, t_HadronFlavour>::RG_Evolution_gammas_exponent(double mu_i,
                             double mu_f,
                             double Q,
                             class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                             gsl_integration_workspace * wsp)
{
    double result, abserr;
    
    auto integrand = make_gsl_function( [&](double mu) {return
        (gamma2Parton(colPDF.Get_AlphaStrong(mu)) - log(Q/mu)*gammaCSKernel(colPDF.Get_AlphaStrong(mu)))/mu;});
    gsl_integration_qag(integrand, mu_i, mu_f,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6,wsp,
                         &result, &abserr);
    
    return result;
}

// RG evolution
template<int t_PartonFlavour, int t_HadronFlavour>
double CollinearPDF<t_PartonFlavour, t_HadronFlavour>::RG_Evolution_bT(double bT,
                                                                    double mu_i,
                             double mu_f,
                             double Q,
                             class CollinearPDF<t_PartonFlavour, t_HadronFlavour> colPDF,
                             gsl_integration_workspace * wsp,
                             bool NonPerturbative)
{
    return exp( RG_Evolution_gammas_exponent(mu_i, mu_f, Q, colPDF, wsp) + log(Q/mu_i)*CSKernel(bT, mu_i, colPDF.Get_AlphaStrong(mu_i), NonPerturbative) );// old
//    return exp( RG_Evolution_gammas_exponent(mu_i, mu_f, Q, colPDF, wsp) );
}


// Mellin convolution for PDF type
template<int t_PartonFlavour, int t_HadronFlavour>
double CollinearPDF<t_PartonFlavour, t_HadronFlavour>::MellinConvolutionPDF(double x, double(*P)(double),
                                                double(CollinearPDF<t_PartonFlavour,t_HadronFlavour>::*Evaluate)(double) const)
{
    double result;
    double abserr;
    double p[2] = {x, 1};
    IntegrationWorkspace wsp(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    
    auto integrand = make_gsl_function( [&](double chi) {return P(chi)*(this->*Evaluate)(x/chi)/chi;});
    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp, &result, &abserr);
    
    return result;
}


// Get alpha Strong
template<int t_PartonFlavour, int t_HadronFlavour>
double CollinearPDF<t_PartonFlavour, t_HadronFlavour>::Get_AlphaStrong(double mu)
{
    int Nf = 3;
    if(!m_isLHAPDF)
    {
        return evolve_alphas(mu*mu, Nf, 0);
    }
    else// to change if LHAPDF works
    {
        double q2 = pow(mu,2);// or HardScale ???
        return m_lhapdf_pdf->alphasQ2(q2);
    }
    
    return 0;
}

// Get Parton's charge
template<int t_PartonFlavour, int t_HadronFlavour>
double CollinearPDF<t_PartonFlavour, t_HadronFlavour>::Get_PartonCharge()
{
    if(t_PartonFlavour == UP || t_PartonFlavour == CHARM || t_PartonFlavour == TOP)
    {
        return 2.0/3.0;
    }
    else if (t_PartonFlavour == DOWN || t_PartonFlavour == STRANGE || t_PartonFlavour == BOTTOM)
    {
        return -1.0/3.0;
    }
    else if (t_PartonFlavour == UBAR || t_PartonFlavour == CBAR || t_PartonFlavour == TBAR)
    {
        return -2.0/3.0;
    }
    else if (t_PartonFlavour == DBAR || t_PartonFlavour == SBAR || t_PartonFlavour == BBAR)
    {
        return 1.0/3.0;
    }
    else if (t_PartonFlavour == GLUON || t_PartonFlavour == PHOTON)
    {
        return 0.0;
    }
    else if (t_PartonFlavour == PROTON)
    {
        return 1.0;
    }
    else
    {
        std::cout << "No flavour found" << std::endl;
        exit(1);
    }
        
}


#endif /* CollinearPDF_h */
