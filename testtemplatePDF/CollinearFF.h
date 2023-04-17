//
//  CollinearFF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/24/22.
//

#include <math.h>

#include "MellinConvolution.h"
#include "Flavour.h"
#include "FFBase.h"
#include "QCDParameters.h"

#ifndef CollinearFF_h
#define CollinearFF_h

template<int t_HadronFlavour, int t_PartonFlavour>
class CollinearFF : public FFBase
{
public:
    
    // Constructors -----------------------------------------------------------
    CollinearFF(double Q, double muQ, bool isLHAPDF)
    :   FFBase(Q, muQ, t_HadronFlavour, t_PartonFlavour),
        m_isLHAPDF(isLHAPDF),
        m_HardScale(Q),
        m_RenormalizationScale(muQ)
        
    {
        // Import LHAPDF ff
        m_lhapdf_ff = LHAPDF::mkPDF(m_setname, m_setmember);
    };
    
    // Methods ----------------------------------------------------------------
    double Evaluate(double z);
    
    // Evolution exponent for FF
    double RG_Evolution_gammas_exponent(double mu_i,
                                        double mu_f,
                                        double Q,
                                        class CollinearFF<t_HadronFlavour, t_PartonFlavour> colFF,
                                        gsl_integration_workspace * wsp);
    // RG Evolution
    double RG_Evolution_bT(double bT,
                           double mu_i,
                           double mu_f,
                           double Q,
                           class CollinearFF<t_HadronFlavour,t_PartonFlavour> colFF,
                           gsl_integration_workspace * wsp,
                           bool NonPerturbative);
    
    // Mellin Convolution type
    double MellinConvolutionFF(double z, double(*P)(double),
                                double(CollinearFF<t_HadronFlavour,t_PartonFlavour>::*Evaluate)(double));
    
    // Get ---------------------------------------
    double Get_AlphaStrong(double mu);
    
    double Get_PartonCharge();
    
    // Set ---------------------------------------
    void Set_LHAPDF_Setname_and_member(const std::string LHAPDFSET_FF, int member)
    {
        m_setname = LHAPDFSET_FF;
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
    const std::string m_setname = "MAPFF10NLOPIp2";//default (can be changed)
    /* \brief LHAPDF set member */
    int m_setmember = 0;//default (can be changed)
    
    /* LHAPDF pointer to PDF set*/
    LHAPDF::PDF* m_lhapdf_ff;

};


// Collinear FF for quarks and antiquarks
// Default collinear FF
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::UNDEFINED_PARTON>::Evaluate(double z)
{
    return z * (1.0 - z);
}


// UP
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::UP>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(UP,z,Q2)/z;
    }
    else
    {
        return pow(z,2) * (1.0 - z);
    }
}

// DOWN
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::DOWN>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(DOWN,z,Q2)/z;
    }
    else
    {
        return z*pow((1.0 - z), 2);
    }
}


// STRANGE
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::STRANGE>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(STRANGE,z,Q2)/z;
    }
    else
    {
        return pow(z,3)*pow((1.0 - z), 2);
    }
}


// Antiquarks
// UBAR
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::UBAR>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(UBAR,z,Q2)/z;
    }
    else
    {
        return pow(z,2) * (1.0 - z);
    }
}

// DOWN
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::DBAR>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(DBAR,z,Q2)/z;
    }
    else
    {
        return z*pow((1.0 - z), 2);
    }
}


// STRANGE
template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::SBAR>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(SBAR,z,Q2)/z;
    }
    else
    {
        return pow(z,3)*pow((1.0 - z), 2);
    }
}

// Collinear FF for bosons


template<>
double CollinearFF<Flavour::UNDEFINED_HADRON, Flavour::GLUON>::Evaluate(double z)
{
    if(m_isLHAPDF)
    {
        double Q2 = pow(Get_RenormalizationScale(),2);
        return m_lhapdf_ff->xfxQ2(GLUON,z,Q2)/z;
    }
    else
    {
        return pow(z,2)*pow((1.0 - z), 2);
    }
}

// --------------------------------------------------------------------------------------------

// Evolution exponent for FF
template<int t_HadronFlavour,int t_PartonFlavour>
double CollinearFF<t_HadronFlavour,t_PartonFlavour>::RG_Evolution_gammas_exponent(double mu_i,
                             double mu_f,
                             double Q,
                             class CollinearFF<t_HadronFlavour,t_PartonFlavour> colFF,
                             gsl_integration_workspace * wsp)
{
    double result, abserr;
    
    auto integrand = make_gsl_function( [&](double mu) {return
        (gamma2Parton(colFF.Get_AlphaStrong(mu)) - log(Q/mu)*gammaCSKernel(colFF.Get_AlphaStrong(mu)))/mu;});
    gsl_integration_qag(integrand, mu_i, mu_f,NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6,wsp,
                         &result, &abserr);
    
    return result;
}

// RG evolution
template<int t_HadronFlavour, int t_PartonFlavour>
double CollinearFF<t_HadronFlavour,t_PartonFlavour>::RG_Evolution_bT(double bT,
                             double mu_i,
                             double mu_f,
                             double Q,
                             class CollinearFF<t_HadronFlavour,t_PartonFlavour> colFF,
                             gsl_integration_workspace * wsp,
                             bool NonPerturbative)
{
    return exp( RG_Evolution_gammas_exponent(mu_i, mu_f, Q, colFF, wsp) + log(Q/mu_i)*CSKernel(bT, mu_i, colFF.Get_AlphaStrong(mu_i), NonPerturbative) );
}


// Mellin convolution for FF type
template<int t_HadronFlavour, int t_PartonFlavour>
double CollinearFF<t_HadronFlavour,t_PartonFlavour>::MellinConvolutionFF(double z, double(*P)(double),
                                            double(CollinearFF<t_HadronFlavour, t_PartonFlavour>::*Evaluate)(double))
{
    double result;
    double abserr;
    double p[2] = {z, 1};
    IntegrationWorkspace wsp(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    
    auto integrand = make_gsl_function( [&](double chi) {return P(chi)*(this->*Evaluate)(z/chi)/chi;});
    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp, &result, &abserr);
    //gsl_integration_qag(integrand, x, 1, epsabs, epsrel, limit, 1, wsp, &result, &abserr);
    
    return result;
}


// Get alpha Strong
template<int t_HadronFlavour, int t_PartonFlavour>
double CollinearFF<t_HadronFlavour, t_PartonFlavour>::Get_AlphaStrong(double mu)
{
    int Nf = 3;
    if(!m_isLHAPDF)
    {
        return evolve_alphas(mu*mu, Nf, 0);
    }
    else// to change if LHAPDF works
    {
        double q2 = pow(mu,2);// or HardScale ???
        return m_lhapdf_ff->alphasQ2(q2);
    }
    
    return 0;
}

// Get Parton's charge
template<int t_HadronFlavour, int t_PartonFlavour>
double CollinearFF<t_HadronFlavour, t_PartonFlavour>::Get_PartonCharge()
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

#endif /* CollinearFF_h */
