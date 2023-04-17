//
//  TMDFF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/24/22.
//
#include <gsl/gsl_errno.h>

#include "CollinearFF.h"
#include "QCDSplittingKernels.h"
#include "alphaStrong.h"
#include "MellinConvolution.h"
#include "bstar.h"

#ifndef TMDFF_h
#define TMDFF_h

// Struct for mass parameters
struct InputTMDFF_massParameters {double m_Hi; double m_Hg; double M_Hi;};


template<int t_HadronFlavour, int t_QuarkFlavour,int t_BosonFlavour>
class TMDFF
{
    // FlavourQuark can never be equal to FlavourBoson
    //static_assert( std::is_same<t_QuarkFlavour, t_BosonFlavour> );
    
    using QuarkFF_T = class CollinearFF<t_HadronFlavour, t_QuarkFlavour>;
    using BosonFF_T = class CollinearFF<t_HadronFlavour, t_BosonFlavour>;

public:
    
    TMDFF(QuarkFF_T quarkFF, BosonFF_T bosonFF, double z, InputTMDFF_massParameters massParam, bool NonPerturbative)
    :   m_collinearQuark(quarkFF),
        m_collinearBoson(bosonFF),
        m_z(z),
        m_HardScale(quarkFF.Get_HardScale()),
        m_RenormalizationScale(quarkFF.Get_RenormalizationScale()),
        m_NonPerturbative(NonPerturbative),
        m_mass_FFQuarkHadron(massParam.m_Hi),
        m_mass_FFBosonHadron(massParam.m_Hg),
        m_Mass_FFQuarkHadron(massParam.M_Hi)
    {
        // Set mass parameters 
        if(!m_NonPerturbative)// default for perturbative to avoid divergences
        {
            Set_mass_FFQuarkHadron(1e-12);
            Set_mass_FFBosonHadron(1e-12);
            Set_Mass_FFQuarkHadron(1e-12);
        }
        
        
    };
    
    // Methods -----------------------------------------------------------------------------------------
    // Evaluate TMDFF in momentum space
    double Evaluate(double z, double kT);
    // Evaluate cutoff integral over kT space up to kc
    double Evaluate_CutoffIntegral(double z, double kc);
    // Evaluate TMDFF in cooerdinate space
    double Evaluate_inbTspace(double z, double bT);
    // Evaluate improved TMDFF in coordinate space
    double Evaluate_inbTspace_improved(double x, double bT, gsl_integration_workspace* wsp);
    // Evaluate g function
    double Evaluate_g_HadronQuark(double bT);
    // Evaluate order alphaS OPE bstar contribution to TMDFF
    double EvaluateOPE_qq_bTstar_1(double mustar, double z, gsl_integration_workspace* wsp);
    double EvaluateOPE_gq_bTstar_1(double mustar, double z, gsl_integration_workspace* wsp);
    // OPE evolution for TMD FF in bT space
    double OPE_RG_Evolution_bstar_TMDFF(double bT, gsl_integration_workspace* wsp);
    // Evaluate TMD in bT space with conventional bstar method
    double Evaluate_Quark_inbTspace_bstarMethod(double z, double bT, gsl_integration_workspace* wsp);
    // Get average of squared transverse momentum
    double Get_mean_squared_transverse_momentum(double max_kT2);// input is maximum squared transverse momentum
    
    // Declaration of Input TMD FF Coefficients
    double AD_Hi(double z, double muQ);
    double BD_Hi(double z, double muQ);
    double AD_Hg(double z, double muQ);
    // Non perturbative small Transverse Momentum Coefficient
    double CD_Hi_0(double z, double muQ);
    double CD_Hi_1(double z, double muQ, double alphaS);
    // Difference order alphaS from cutoff and MSbar
    double Delta_MSbar_Hi_1(double z, double muQ, double alphaS);

    
    // Get -----------------------------------------------------------------------------------------
    // Get coefficients
    double Get_A_FFQuarkHadron() const {return m_A_FFQuarkHadron;}
    double Get_B_FFQuarkHadron() const {return m_B_FFQuarkHadron;}
    double Get_A_FFBosonHadron() const {return m_A_FFBosonHadron;}
    double Get_C_FFQuarkHadron() const {return m_C_FFQuarkHadron;}
    
    // Get z
    double Get_z() {return m_z;}
    
    // Get Hard scale and Renormalization scale
    double Get_RenormalizationScale() {return m_collinearQuark.Get_RenormalizationScale();}
    double Get_HardScale() {return m_collinearQuark.Get_HardScale();}
    
    // Non Perturbative flag
    bool Get_NonPerturbative() const {return m_NonPerturbative;}
    
    // Set --------------------------------------------------------------------
    
    void Set_mass_FFQuarkHadron(double m_Hi) {m_mass_FFQuarkHadron = m_Hi;}
    void Set_mass_FFBosonHadron(double m_Hg) {m_mass_FFBosonHadron = m_Hg;}
    void Set_Mass_FFQuarkHadron(double M_Hi) {m_Mass_FFQuarkHadron = M_Hi;}
    
    // Set Hard Scale
    void Set_HardScale(double Q) {m_HardScale = Q;}
    // Set Renormalization Scale
    void Set_RenormalizationScale(double muQ) {m_RenormalizationScale = muQ;}
    
    // Set alpha Strong
    void Set_alphaS(double mu) {m_alphaS = m_collinearQuark.Get_AlphaStrong(mu);}
    
    void Set_z(double z) {m_z = z;}
    // Set bool for ansatz of g functions
    void Set_gAnsatz(bool gAnsatz) {m_gAnsatz = gAnsatz;}
    // Set bool for NonPerturbative
    void Set_NonPerturbative(bool NonPerturbative) {m_NonPerturbative = NonPerturbative;}
    
    // Set coefficients
    void Set_A_FFQuarkHadron() {m_A_FFQuarkHadron = m_alphaS*AD_Hi(m_z, m_RenormalizationScale);}
    void Set_B_FFQuarkHadron() {m_B_FFQuarkHadron = m_alphaS*BD_Hi(m_z, m_RenormalizationScale);}
    void Set_A_FFBosonHadron() {m_A_FFBosonHadron = m_alphaS*AD_Hg(m_z, m_RenormalizationScale);}
    void Set_C_FFQuarkHadron() {m_C_FFQuarkHadron = CD_Hi_0(m_z, m_RenormalizationScale) + CD_Hi_1(m_z, m_RenormalizationScale, m_alphaS);}
    
private:
    
    // Member variables -------------------------------------------------------
    
    /* \brief Collinear quark ff*/
    QuarkFF_T m_collinearQuark;
    /* \brief Collinear boson ff*/
    BosonFF_T m_collinearBoson;
    
    /* \brief z*/
    double m_z;
    
    /* \brief Hard Scale*/
    double m_HardScale;
    
    /* \brief Renormalization Scale*/
    double m_RenormalizationScale;
    
    /* \brief Input TMD FF mass parameters */
    double m_mass_FFQuarkHadron;
    double m_mass_FFBosonHadron;
    double m_Mass_FFQuarkHadron;
    
    /* \brief Input TMD FF coefficients */
    double m_alphaS = m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale);// inherited from CollinearFF
    // Assign coefficient parameters to class members
    double m_A_FFQuarkHadron = m_alphaS*AD_Hi(m_z, m_RenormalizationScale);
    double m_B_FFQuarkHadron = m_alphaS*BD_Hi(m_z, m_RenormalizationScale);
    double m_A_FFBosonHadron = m_alphaS*AD_Hg(m_z, m_RenormalizationScale);
    double m_C_FFQuarkHadron = CD_Hi_0(m_z, m_RenormalizationScale) + CD_Hi_1(m_z, m_RenormalizationScale, m_alphaS);

    
    /* \brief Choose if masses are small or not*/
    bool m_NonPerturbative;
    
    /* \brief  Choose if ansatz is used for g functions*/
    bool m_gAnsatz = true;// Default for now
    
    
};// end of TMDFF class


// Definition of functions ---------------------------------------------------------------------
// Evaluate TMDFF in momentum space
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate(double z, double kT)
{
    double Q  = m_collinearQuark.Get_HardScale();
    double mi, mg, Mi;
    
    double AB_i;
    double A_g;
    double C_Mi;
    
    mi = m_mass_FFQuarkHadron;
    mg = m_mass_FFBosonHadron;
    Mi = m_Mass_FFQuarkHadron;
    
    AB_i = 1.0/(2.0*M_PI*z*z)/(kT*kT + mi*mi) *( m_A_FFQuarkHadron + m_B_FFQuarkHadron*log(Q*Q/(kT*kT + mi*mi)) );
    A_g  = 1.0/(2.0*M_PI*z*z)/(kT*kT + mg*mg) *  m_A_FFBosonHadron;
    C_Mi = m_C_FFQuarkHadron/(M_PI*Mi*Mi)*exp(-z*z*kT*kT/(Mi*Mi));
    
    return AB_i + A_g + C_Mi;
}

// Evaluate cutoff integral over kT space up to Renormalization Scale
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate_CutoffIntegral(double z, double kc)
{
    double Q   = m_HardScale;
    double muQ = kc;
    
    double output;
    
    if(m_NonPerturbative)
    {
        output = m_C_FFQuarkHadron*(1 - exp(-z*z*muQ*muQ/pow(m_Mass_FFQuarkHadron,2))) + 0.5*m_A_FFBosonHadron*log(1 + muQ*muQ/pow(m_mass_FFBosonHadron,2)) + 0.5*m_A_FFQuarkHadron*log(1 + muQ*muQ/pow(m_mass_FFQuarkHadron,2)) + 0.25*m_B_FFQuarkHadron*( pow(log(pow(m_mass_FFQuarkHadron/Q,2)),2) - pow(log((muQ*muQ + pow(m_mass_FFQuarkHadron,2))/pow(Q,2)),2) );
    }
    else
    {
        output = 0.5*m_A_FFBosonHadron*log(1 + muQ*muQ/pow(m_mass_FFBosonHadron,2)) + 0.5*m_A_FFQuarkHadron*log(1 + muQ*muQ/pow(m_mass_FFQuarkHadron,2)) + 0.25*m_B_FFQuarkHadron*( pow(log(pow(m_mass_FFQuarkHadron/Q,2)),2) - pow(log((muQ*muQ + pow(m_mass_FFQuarkHadron,2))/pow(Q,2)),2) );
    }
    
    
    
    return output;
}

// Evaluate TMDFF in cooerdinate space
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate_inbTspace(double z, double bT)
{
    if (m_NonPerturbative)
    {
        gsl_set_error_handler_off();// comment when NOT strictly necessary
        double output = gsl_sf_bessel_K0(m_mass_FFQuarkHadron*bT) * ( m_A_FFQuarkHadron + m_B_FFQuarkHadron*log(m_HardScale*m_HardScale*bT*exp(M_EULER)/(2*m_mass_FFQuarkHadron)) ) + gsl_sf_bessel_K0(m_mass_FFBosonHadron*bT) * m_A_FFBosonHadron + m_C_FFQuarkHadron*exp(-bT*bT*m_Mass_FFQuarkHadron*m_Mass_FFQuarkHadron/(4.0*z*z));
        return (output)/(z*z);
    }
    else
        return -log(bT*m_RenormalizationScale*exp(M_EULER)/2.0)*(m_A_FFQuarkHadron + m_B_FFQuarkHadron*log(bT*m_RenormalizationScale*exp(M_EULER)*m_HardScale*m_HardScale/(2.0*m_RenormalizationScale*m_RenormalizationScale))) + m_C_FFQuarkHadron/(z*z);
}


// Evaluate improved TMDPDF in coordinate space with Q0 bar method
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate_inbTspace_improved(double z, double bT, gsl_integration_workspace* wsp)
{
    double in_HardScale = m_collinearQuark.Get_HardScale();
    double in_RenormalizationScale = m_collinearQuark.Get_RenormalizationScale();
    
    double new_HardScale = Q0bar(bT);// RG improvement scale
    double new_RenormalizationScale = new_HardScale;
    
    m_collinearQuark.Set_RenormalizationScale(new_RenormalizationScale);
    m_collinearQuark.Set_HardScale(new_HardScale);
    Set_RenormalizationScale(new_RenormalizationScale);
    Set_HardScale(new_HardScale);
    
    Set_alphaS(new_RenormalizationScale);
    Set_A_FFQuarkHadron();
    Set_B_FFQuarkHadron();
    Set_A_FFBosonHadron();
    Set_C_FFQuarkHadron();

    
    
    
//    double E_Q0bar_to_Q0 = m_collinearQuark.RG_Evolution_bT(bT, new_RenormalizationScale, in_RenormalizationScale, in_HardScale, m_collinearQuark, wsp) * exp( log(INPUT_HARD_SCALE/new_HardScale) * CSKernel(bT, new_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(new_RenormalizationScale), true));//old
    double E_Q0bar_to_Q0 = m_collinearQuark.RG_Evolution_bT(bT, new_RenormalizationScale, in_RenormalizationScale, in_HardScale, m_collinearQuark, wsp, true);
    
    double output = Evaluate_inbTspace(z, bT)*E_Q0bar_to_Q0;
    
    // reset to initial values
    m_collinearQuark.Set_RenormalizationScale(in_RenormalizationScale);
    m_collinearQuark.Set_HardScale(in_HardScale);
    m_collinearBoson.Set_RenormalizationScale(in_RenormalizationScale);
    m_collinearBoson.Set_HardScale(in_HardScale);
    Set_RenormalizationScale(in_RenormalizationScale);
    Set_HardScale(in_HardScale);
    Set_alphaS(new_RenormalizationScale);
    Set_A_FFQuarkHadron();
    Set_B_FFQuarkHadron();
    Set_A_FFBosonHadron();
    Set_C_FFQuarkHadron();
    
    return output;
    
}


// Evaluate g function
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate_g_HadronQuark(double bT)
{
    if(m_gAnsatz)
    {
//        std::cout << " FF m_gAnsatz = " << std::endl;
        double MgD = 0.3*m_z;
        return pow(MgD*bT/(2*m_z),2);
    }
    else
    {
        bool in_NonPerturbative = m_NonPerturbative;

        Set_NonPerturbative(false);
        double P_tmdFF = Evaluate_inbTspace(m_z, bstar(bT));
        Set_NonPerturbative(true);
        double NP_tmdFF = Evaluate_inbTspace(m_z, bT);

        Set_NonPerturbative(in_NonPerturbative);

        double g_Hq = -log(NP_tmdFF/P_tmdFF);

        return g_Hq;
    }
}


// Evaluate order alphaS OPE bstar contribution to TMDFF
template< int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::EvaluateOPE_qq_bTstar_1(double mustar,
                                                                                  double z,                               gsl_integration_workspace* wsp)
{
    double result, abserr;
    
    double mu = m_RenormalizationScale;

    m_collinearQuark.Set_RenormalizationScale(mustar);

    auto integrand = make_gsl_function( [&](double xi) {return xi*(4.0*(2.0/(1.0-xi)+1.0/(xi*xi)+1.0/xi)*log(xi) + 2.0/(xi*xi) - 2.0/xi)*m_collinearQuark.Evaluate(z/xi)/(m_z*m_z);});
//    auto integrand = make_gsl_function( [&](double xi) {return m_z*(4.0*(2.0/(1.0-xi)+1.0/(xi*xi)+1.0/xi)*log(xi) + 2.0/(xi*xi) - 2.0/xi)*m_collinearQuark.Evaluate(z/xi)/(xi*xi);});
    
    gsl_integration_qag(integrand, z, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6, wsp, &result, &abserr);


    double output =  CF*( result  - M_PI*M_PI/6.0*m_collinearQuark.Evaluate(z)/(z*z) );
    
    m_collinearQuark.Set_RenormalizationScale(mu);
    
    return output/(4*M_PI);// to be multiplied by alphaS

}

template< int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::EvaluateOPE_gq_bTstar_1(double mustar,
                                                                                  double z,                               gsl_integration_workspace* wsp)
{
    double result, abserr;
    double mu = m_RenormalizationScale;

    m_collinearBoson.Set_RenormalizationScale(mustar);

    auto integrand = make_gsl_function( [&](double xi) {return xi/(m_z*m_z)*2.0*( 2.0*(1.0 + pow(1-xi,2))*log(xi) + xi*xi)/(pow(xi,3))*m_collinearBoson.Evaluate(z/xi);});
//    auto integrand = make_gsl_function( [&](double xi) {return m_collinearBoson.Evaluate(z/xi)*m_z/(xi*xi)*(2.0*( 2.0*(1.0 + pow(1-xi,2))*log(xi) + xi*xi)/(pow(xi,3)));});
    
    gsl_integration_qag(integrand, z, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6, wsp, &result, &abserr);

    double output =  CF * result ;
    
    m_collinearBoson.Set_RenormalizationScale(mu);// reset value
    
    return output/(4*M_PI);// to be multiplied by alphaS


}

// OPE evolution factor for TMD FF in bt space
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::OPE_RG_Evolution_bstar_TMDFF(double bT, gsl_integration_workspace* wsp)
{
    double bstr = bstar(bT);
    
    double output = m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp, false) * exp(-Evaluate_g_HadronQuark(bT) - log(m_HardScale/INPUT_HARD_SCALE)*gCSKernel(bT, m_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale), false,true));
    
    // without gK
//    double output = m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp) * exp(-Evaluate_g_HadronQuark(bT) );
    
//     without the g function
//    double output =  m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp) * exp( - log(m_HardScale/INPUT_HARD_SCALE)*gCSKernel(bT, m_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale), false,true));
    
    
    return output;

}

// Conventional TMD ff in bT space with bstar method
// Quark TMD ff
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Evaluate_Quark_inbTspace_bstarMethod(double z, double bT, gsl_integration_workspace* wsp)
{
    double output;
    double mustr = mubstar(bstar(bT));
    double in_RenormalizationScale = m_collinearQuark.Get_RenormalizationScale();
    
    double alphaS_mustr = m_collinearQuark.Get_AlphaStrong(mustr);
    m_collinearQuark.Set_RenormalizationScale(mustr);// set to mubstar
    m_collinearBoson.Set_RenormalizationScale(mustr);// set to Renormalization Scale to mubstar
    
    bool in_NonPerturbative = Get_NonPerturbative();
    Set_NonPerturbative(false);// nothing should depend on the model
    
    double D_OPE = m_collinearQuark.Evaluate(z)/pow(z,2) + alphaS_mustr * ( EvaluateOPE_qq_bTstar_1(mustr, z, wsp) + EvaluateOPE_gq_bTstar_1(mustr, z, wsp) );
    
    Set_NonPerturbative(in_NonPerturbative);
    
    output = D_OPE * OPE_RG_Evolution_bstar_TMDFF(bT, wsp);
    
    m_collinearQuark.Set_RenormalizationScale(in_RenormalizationScale);// reset to input value
    m_collinearBoson.Set_RenormalizationScale(in_RenormalizationScale);// reset to input value
    
    return output;
}



// Get average of squared transverse momentum
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Get_mean_squared_transverse_momentum(double max_kT2)// input is maximum squared transverse momentum
{
    double kc2 = max_kT2;
    double Q   = m_HardScale;
    double z   = m_z;
    double mi  = m_mass_FFQuarkHadron;
    double mg  = m_mass_FFBosonHadron;
    double Mi  = m_Mass_FFQuarkHadron;
    
    return (0.5*m_A_FFQuarkHadron*( kc2 - mi*mi*log( 1.0 + kc2/(mi*mi) ) )
          + 0.5*m_B_FFQuarkHadron*( kc2 + 0.5*mi*mi*pow(log(Q*Q/(kc2 + mi*mi)),2) + (kc2 + mi*mi)*log(Q*Q/(kc2 + mi*mi)) -0.5*mi*mi*log(Q*Q/(mi*mi))*( log(Q*Q/(mi*mi))+ 2.0) )
          + 0.5*m_A_FFBosonHadron*( kc2 - mg*mg*log( 1.0 + kc2/(mg*mg) ) )
          + m_C_FFQuarkHadron/(z*z)*(Mi*Mi - (kc2 + Mi*Mi)*exp(-kc2/(Mi*Mi)) ))/(z*z);
}


// Coefficients for input TMD FF

// Coefficient A_FFQuarkHadron
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::AD_Hi (double z, double muQ)
{
    double result;
    double abserr;
    double p[2] = {z, 1};
    IntegrationWorkspace wsp1(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    auto integrand = make_gsl_function( [&](double xi) {return ( (1+xi*xi)*m_collinearQuark.Evaluate(z/xi)/xi - 2*m_collinearQuark.Evaluate(z) )/(1-xi);});
    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp1, &result, &abserr);
    
    return CF/M_PI * (result + 2*m_collinearQuark.Evaluate(z)*log(1-z));
}

// Coefficient B_FFQuarkHadron
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::BD_Hi (double z, double muQ)
{
    return CF/M_PI * m_collinearQuark.Evaluate(z);
}

// Coefficient A_FFBosonHadron
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::AD_Hg (double z, double muQ)
{
//    return 1.0/M_PI * MellinConvolution(z,Pg_q, m_collinearBoson.Evaluate);
    return 1.0/M_PI * m_collinearBoson.MellinConvolutionFF(z,Pg_q, &CollinearFF<t_HadronFlavour,t_BosonFlavour>::Evaluate);
}

// Terms of C coefficient -------------------------------------------------------------------------
// Coefficient C_FFQuarkHadron order one in alpha strong times alpha strong ----------------------
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::CD_Hi_1 (double z, double muQ, double alphaS)
{
    double mi_D = m_mass_FFQuarkHadron;
    double mg_D = m_mass_FFBosonHadron;
    double Q  = m_collinearQuark.Get_HardScale();
    
    double output = - m_A_FFQuarkHadron*log(muQ/mi_D) - m_B_FFQuarkHadron*log(muQ/mi_D)*log(Q*Q/(muQ*mi_D)) - m_A_FFBosonHadron*log(muQ/mg_D) + alphaS/(2.0*M_PI)*(m_collinearQuark.MellinConvolutionFF(z,CDelta_Hi_tilde,&CollinearFF<t_HadronFlavour,t_QuarkFlavour>::Evaluate) - CF*M_PI*M_PI/12.0*m_collinearQuark.Evaluate(z) + m_collinearBoson.MellinConvolutionFF(z,CDelta_Hg,&CollinearFF<t_HadronFlavour,t_BosonFlavour>::Evaluate));
    
    return output;
}
// Coefficient C_FFQuarkHadron order zero in alpha strong ------------------------------------------
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::CD_Hi_0 (double z, double muQ)
{
    return m_collinearQuark.Evaluate(z);
}

// Delta Cutoff MSbar
template<int t_HadronFlavour,int t_QuarkFlavour,int t_BosonFlavour>
double TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour>::Delta_MSbar_Hi_1(double z, double muQ,double alphaS)
{
    double output;
    // with Dirac Delta contribution = Standard MSbar
    output = - alphaS*CF*M_PI/24.0*m_collinearQuark.Evaluate(z) + alphaS/(2.0*M_PI)*(m_collinearQuark.MellinConvolutionFF(z,CDelta_Hi_tilde,&CollinearFF<t_HadronFlavour,t_QuarkFlavour>::Evaluate)  + m_collinearBoson.MellinConvolutionFF(z,CDelta_Hg,&CollinearFF<t_HadronFlavour,t_BosonFlavour>::Evaluate));
    
    // without Dirac Delta contribution = Collins MSbar
//    output = alphaS/(2.0*M_PI)*(m_collinearQuark.MellinConvolutionFF(z,CDelta_Hi_tilde,&CollinearFF<t_HadronFlavour,t_QuarkFlavour>::Evaluate)  + m_collinearBoson.MellinConvolutionFF(z,CDelta_Hg,&CollinearFF<t_HadronFlavour,t_BosonFlavour>::Evaluate));
//
    // Delta with k_c = mu_F = mu/z
    
    return output;
    
}

#endif /* TMDFF_h */
