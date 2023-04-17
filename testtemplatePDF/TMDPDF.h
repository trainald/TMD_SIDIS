//
//  TMDPDF.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/21/22.
//

#include <gsl/gsl_sf_bessel.h>

#include "CollinearPDF.h"
#include "QCDSplittingKernels.h"
#include "alphaStrong.h"
#include "MellinConvolution.h"
#include "bstar.h"
#include "Q0barEvolution.h"



#ifndef TMDPDF_h
#define TMDPDF_h
 
// Struct for mass parameters
struct InputTMDPDF_massParameters {double m_iH; double m_gH; double M_iH;};



template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
class TMDPDF
{
    // FlavourQuark can never be equal to FlavourBoson
    //static_assert( std::is_same<t_QuarkFlavour, t_BosonFlavour> );
    
    using QuarkPDF_T = class CollinearPDF<t_QuarkFlavour, t_HadronFlavour>;
    using BosonPDF_T = class CollinearPDF<t_BosonFlavour, t_HadronFlavour>;

public:
    
    TMDPDF(QuarkPDF_T quarkPDF, BosonPDF_T bosonPDF, double x, InputTMDPDF_massParameters massParam, bool NonPerturbative)
    :   m_collinearQuark(quarkPDF),
        m_collinearBoson(bosonPDF),
        m_xBjorken(x),
        m_HardScale(quarkPDF.Get_HardScale()),
        m_RenormalizationScale(quarkPDF.Get_RenormalizationScale()),
        m_NonPerturbative(NonPerturbative),
        m_mass_PDFQuarkHadron(massParam.m_iH),
        m_mass_PDFBosonHadron(massParam.m_gH),
        m_Mass_PDFQuarkHadron(massParam.M_iH)
    {
        // Set mass parameters
        if(!m_NonPerturbative)// default for perturbative to avoid divergences
        {
            Set_mass_PDFQuarkHadron(1e-12);
            Set_mass_PDFBosonHadron(1e-12);
            Set_Mass_PDFQuarkHadron(1e-12);
        }
        
        // Import TMD PDF from LHAPDF if needed
        m_lhapdf_TMDpdf = LHAPDF::mkPDF(m_setname, m_setmember);
        
    };
    
    
    
    // Methods -----------------------------------------------------------------------------------------
    // Evaluate TMDPDF in momentum space
    double Evaluate(double x, double kT) const;
    // Evaluate cutoff integral over kT space up to kc
    double Evaluate_CutoffIntegral(double x, double kc);
    // Evaluate TMDPDF in cooerdinate space
    double Evaluate_inbTspace(double x, double bT);
    // Evaluate improved TMDPDF in coordinate space
    double Evaluate_inbTspace_improved(double x, double bT, gsl_integration_workspace* wsp);
    // Evaluate g function with or without ansatz
    double Evaluate_g_QuarkHadron(double bT);
    // Evaluate order alphaS OPE bstar contribution to TMDPDF
    double EvaluateOPE_qq_bTstar_1(double mustar, double x, gsl_integration_workspace* wsp);
    double EvaluateOPE_qg_bTstar_1(double mustar, double x, gsl_integration_workspace* wsp);
    // OPE evolution for TMD PDF in bT space
    double OPE_RG_Evolution_bstar_TMDPDF(double bT, gsl_integration_workspace* wsp);
    // Evaluate TMD in bT space with conventional bstar method
    double Evaluate_Quark_inbTspace_bstarMethod(double x, double bT, gsl_integration_workspace* wsp);
    // Get average of squared transverse momentum
    double Get_mean_squared_transverse_momentum(double max_kT2);// input is maximum squared transverse momentum
    
    // Declaration of Input TMD PDF Coefficients
    double Af_iH(double x, double muQ);
    double Bf_iH(double x, double muQ);
    double Af_gH(double x, double muQ);
    // Non perturbative small Transverse Momentum Coefficient
    double Cf_iH_0(double x, double muQ);
    double Cf_iH_1(double x, double muQ, double alphaS);
    // Difference order alphaS from cutoff and MSbar
    double Delta_MSbar_iH_1(double x, double muQ, double alphaS);
    
    // Get -----------------------------------------------------------------------------------------
    
    // HSO approach
    
    // Get coefficients
    double Get_A_PDFQuarkHadron() const {return m_A_PDFQuarkHadron;}
    double Get_B_PDFQuarkHadron() const {return m_B_PDFQuarkHadron;}
    double Get_A_PDFBosonHadron() const {return m_A_PDFBosonHadron;}
    double Get_C_PDFQuarkHadron() const {return m_C_PDFQuarkHadron;}
    
    // Get x Bjorken
    double Get_xBjorken() const {return m_xBjorken;}
    // Get Hard scale and Renormalization scale
    double Get_RenormalizationScale() const {return m_collinearQuark.Get_RenormalizationScale();}
    double Get_HardScale() const {return m_collinearQuark.Get_HardScale();}
    
    // Non Perturbative flag
    bool Get_NonPerturbative() const {return m_NonPerturbative;}
    
    // Set --------------------------------------------------------------------
    
    void Set_mass_PDFQuarkHadron(double m_iH) {m_mass_PDFQuarkHadron = m_iH;}
    void Set_mass_PDFBosonHadron(double m_gH) {m_mass_PDFBosonHadron = m_gH;}
    void Set_Mass_PDFQuarkHadron(double M_iH) {m_Mass_PDFQuarkHadron = M_iH;}
    
    // Set Hard Scale
    void Set_HardScale(double Q) {m_HardScale = Q;}
    // Set Renormalization Scale
    void Set_RenormalizationScale(double muQ) {m_RenormalizationScale = muQ;}
    
    // Set alpha Strong
    void Set_alphaS(double mu) {m_alphaS = m_collinearQuark.Get_AlphaStrong(mu);}
    
    void Set_xBjorken(double x) {m_xBjorken = x;}
    // Set bool for ansatz of g functions
    void Set_gAnsatz(bool gAnsatz) {m_gAnsatz = gAnsatz;}
    // Set bool for NonPerturbative
    void Set_NonPerturbative(bool NonPerturbative) {m_NonPerturbative = NonPerturbative;}
    
    // Set coefficients
    void Set_A_PDFQuarkHadron() {m_A_PDFQuarkHadron = m_alphaS*Af_iH(m_xBjorken, m_RenormalizationScale);}
    void Set_B_PDFQuarkHadron() {m_B_PDFQuarkHadron = m_alphaS*Bf_iH(m_xBjorken, m_RenormalizationScale);}
    void Set_A_PDFBosonHadron() {m_A_PDFBosonHadron = m_alphaS*Af_gH(m_xBjorken, m_RenormalizationScale);}
    void Set_C_PDFQuarkHadron() {m_C_PDFQuarkHadron = Cf_iH_0(m_xBjorken, m_RenormalizationScale) + Cf_iH_1(m_xBjorken, m_RenormalizationScale, m_alphaS);}
    
    // TMD directly from LHAPDF
    void Set_LHAPDF_Setname_and_member(const std::string LHAPDFSET_PDF, int member)
    {
        m_setname   = LHAPDFSET_PDF;
        m_setmember = member;
    }
    
private:
    
    // Member variables -------------------------------------------------------
    
    /* \brief Collinear quark pdf*/
    QuarkPDF_T m_collinearQuark;
    /* \brief Collinear boson pdf*/
    BosonPDF_T m_collinearBoson;
    
    /* \brief x Bjorken*/
    double m_xBjorken;
    
    /* \brief Hard Scale*/
    double m_HardScale;
    
    /* \brief Renormalization Scale*/
    double m_RenormalizationScale;
    
    /* \brief Input TMD PDF mass parameters */
    double m_mass_PDFQuarkHadron;
    double m_mass_PDFBosonHadron;
    double m_Mass_PDFQuarkHadron;
    
    /* \brief Input TMD PDF coefficients */
    double m_alphaS = m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale);// inherited from CollinearPDF
    // Assign coefficient parameters to class members
    double m_A_PDFQuarkHadron = m_alphaS*Af_iH(m_xBjorken, m_RenormalizationScale);
    double m_B_PDFQuarkHadron = m_alphaS*Bf_iH(m_xBjorken, m_RenormalizationScale);
    double m_A_PDFBosonHadron = m_alphaS*Af_gH(m_xBjorken, m_RenormalizationScale);
    double m_C_PDFQuarkHadron = Cf_iH_0(m_xBjorken, m_RenormalizationScale) + Cf_iH_1(m_xBjorken, m_RenormalizationScale, m_alphaS);
    
    /* \brief Choose if masses are small or not*/
    bool m_NonPerturbative;
    
    /* \brief  Choose if ansatz is used for g functions*/
    bool m_gAnsatz = true;// Default for now
    
    // Diractly from LHAPDF
    /* \brief LHAPDF set name */
    const std::string m_setname = "JAM21PionPDFnlonll_double_Mellin";//default (can be changed)
    /* \brief LHAPDF set member */
    int m_setmember = 0;//default (can be changed)
    
    /* LHAPDF pointer to PDF set*/
    LHAPDF::PDF* m_lhapdf_TMDpdf;
    
    
};// end of TMDPDF class



// Definition of functions ------------------------------------------------------------------

// Evaluate TMDPDF in momentum space
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate(double x, double kT) const
{    
    double AB_i;
    double A_g;
    double C_Mi;

//    mi = m_mass_PDFQuarkHadron;
//    mg = m_mass_PDFBosonHadron;
//    Mi = m_Mass_PDFQuarkHadron;
    
    AB_i = 1.0/(2.0*M_PI)/(kT*kT + m_mass_PDFQuarkHadron*m_mass_PDFQuarkHadron) *( m_A_PDFQuarkHadron + m_B_PDFQuarkHadron*log(m_HardScale*m_HardScale/(kT*kT + m_mass_PDFQuarkHadron*m_mass_PDFQuarkHadron)) );
    A_g  = 1.0/(2.0*M_PI)/(kT*kT + m_mass_PDFBosonHadron*m_mass_PDFBosonHadron) *  m_A_PDFBosonHadron;
    C_Mi = m_C_PDFQuarkHadron/(M_PI*m_Mass_PDFQuarkHadron*m_Mass_PDFQuarkHadron)*exp(-kT*kT/(m_Mass_PDFQuarkHadron*m_Mass_PDFQuarkHadron));
    
    return AB_i + A_g + C_Mi;
}

// Evaluate cutoff integral over kT space up to Renormalization Scale
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate_CutoffIntegral(double x, double kc)
{
//    double Q   = m_HardScale;
//    double muQ = m_RenormalizationScale;
    
    double output;
    
    if(m_NonPerturbative)
    {
//    output = m_C_PDFQuarkHadron*(1 - exp(-m_RenormalizationScale*m_RenormalizationScale/pow(m_Mass_PDFQuarkHadron,2))) + 0.5*m_A_PDFBosonHadron*log(1 + m_RenormalizationScale*m_RenormalizationScale/pow(m_mass_PDFBosonHadron,2)) + 0.5*m_A_PDFQuarkHadron*log(1 + m_RenormalizationScale*m_RenormalizationScale/pow(m_mass_PDFQuarkHadron,2)) + 0.25*m_B_PDFQuarkHadron*( pow(log(pow(m_mass_PDFQuarkHadron/m_HardScale,2)),2) - pow(log((m_RenormalizationScale*m_RenormalizationScale + pow(m_mass_PDFQuarkHadron,2))/pow(m_HardScale,2)),2) );
    output = m_C_PDFQuarkHadron*(1 - exp(-kc*kc/pow(m_Mass_PDFQuarkHadron,2))) + 0.5*m_A_PDFBosonHadron*log(1 + kc*kc/pow(m_mass_PDFBosonHadron,2)) + 0.5*m_A_PDFQuarkHadron*log(1 + kc*kc/pow(m_mass_PDFQuarkHadron,2)) + 0.25*m_B_PDFQuarkHadron*( pow(log(pow(m_mass_PDFQuarkHadron/m_HardScale,2)),2) - pow(log((kc*kc + pow(m_mass_PDFQuarkHadron,2))/pow(m_HardScale,2)),2) );
    
    }
    else
    {
//        output =  0.5*m_A_PDFBosonHadron*log(1 + m_RenormalizationScale*m_RenormalizationScale/pow(m_mass_PDFBosonHadron,2)) + 0.5*m_A_PDFQuarkHadron*log(1 + m_RenormalizationScale*m_RenormalizationScale/pow(m_mass_PDFQuarkHadron,2)) + 0.25*m_B_PDFQuarkHadron*( pow(log(pow(m_mass_PDFQuarkHadron/m_HardScale,2)),2) - pow(log((m_RenormalizationScale*m_RenormalizationScale + pow(m_mass_PDFQuarkHadron,2))/pow(m_HardScale,2)),2) );
        output =  0.5*m_A_PDFBosonHadron*log(1 + kc*kc/pow(m_mass_PDFBosonHadron,2)) + 0.5*m_A_PDFQuarkHadron*log(1 + kc*kc/pow(m_mass_PDFQuarkHadron,2)) + 0.25*m_B_PDFQuarkHadron*( pow(log(pow(m_mass_PDFQuarkHadron/m_HardScale,2)),2) - pow(log((kc*kc + pow(m_mass_PDFQuarkHadron,2))/pow(m_HardScale,2)),2) );
    }
    
    return output;
}

// Evaluate input TMDPDF in cooerdinate space
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate_inbTspace(double x, double bT)
{
//    double Q   = m_HardScale;
//    double muQ = m_RenormalizationScale;
    bool NonPerturbative = m_NonPerturbative;
//    double mi, mg, Mi;
//    mi = m_mass_PDFQuarkHadron;
//    mg = m_mass_PDFBosonHadron;
//    Mi = m_Mass_PDFQuarkHadron;
    double output;
    
    if (NonPerturbative)
    {
        output =  gsl_sf_bessel_K0(m_mass_PDFQuarkHadron*bT) * ( m_A_PDFQuarkHadron + m_B_PDFQuarkHadron*log(m_HardScale*m_HardScale*bT*exp(M_EULER)/(2*m_mass_PDFQuarkHadron)) )  + gsl_sf_bessel_K0(m_mass_PDFBosonHadron*bT) *   m_A_PDFBosonHadron + m_C_PDFQuarkHadron*exp(-bT*bT*m_Mass_PDFQuarkHadron*m_Mass_PDFQuarkHadron/4.0);
        
    }
    else
    {
        output =  (-log(bT*m_RenormalizationScale*exp(M_EULER)/2.0)*(m_A_PDFQuarkHadron + m_B_PDFQuarkHadron*log(bT*m_RenormalizationScale*exp(M_EULER)*m_HardScale*m_HardScale/(2.0*m_RenormalizationScale*m_RenormalizationScale))) + m_C_PDFQuarkHadron);
    }
    

    return output;
        
}

// Evaluate improved TMDPDF in coordinate space with Q0 bar method
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate_inbTspace_improved(double x, double bT, gsl_integration_workspace* wsp)
{
    double in_HardScale = m_collinearQuark.Get_HardScale();
    double in_RenormalizationScale = m_collinearQuark.Get_RenormalizationScale();
    
    double new_HardScale = Q0bar(bT);
    double new_RenormalizationScale = new_HardScale;// could be different
    
    m_collinearQuark.Set_RenormalizationScale(new_RenormalizationScale);
    m_collinearQuark.Set_HardScale(new_HardScale);
    Set_RenormalizationScale(new_RenormalizationScale);
    Set_HardScale(new_HardScale);
    
    Set_alphaS(new_RenormalizationScale);
    Set_A_PDFQuarkHadron();
    Set_B_PDFQuarkHadron();
    Set_A_PDFBosonHadron();
    Set_C_PDFQuarkHadron();
    
    
//    double E_Q0bar_to_Q0 = m_collinearQuark.RG_Evolution_bT(bT, new_RenormalizationScale, in_RenormalizationScale, in_HardScale, m_collinearQuark, wsp) * exp( log(INPUT_HARD_SCALE/new_HardScale) * CSKernel(bT, new_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(new_RenormalizationScale), true));
    double E_Q0bar_to_Q0 = m_collinearQuark.RG_Evolution_bT(bT, new_RenormalizationScale, in_RenormalizationScale, in_HardScale, m_collinearQuark, wsp, true);

    
    double output = Evaluate_inbTspace(x, bT)*E_Q0bar_to_Q0;

    
    
    // reset to initial values
    m_collinearQuark.Set_RenormalizationScale(in_RenormalizationScale);
    m_collinearQuark.Set_HardScale(in_HardScale);
    m_collinearBoson.Set_RenormalizationScale(in_RenormalizationScale);
    m_collinearBoson.Set_HardScale(in_HardScale);
    Set_RenormalizationScale(in_RenormalizationScale);
    Set_HardScale(in_HardScale);
    Set_alphaS(in_RenormalizationScale);
    Set_A_PDFQuarkHadron();
    Set_B_PDFQuarkHadron();
    Set_A_PDFBosonHadron();
    Set_C_PDFQuarkHadron();
    
    return output;
    
}


// Evaluate g function
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate_g_QuarkHadron(double bT)
{
    if(m_gAnsatz)
    {
//        std::cout << "PDF m_gAnsatz = " << std::endl;
        double MgF = 0.4;
        return pow(MgF*bT,2)/4.0;
    }
    else
    {
        bool in_NonPerturbative = m_NonPerturbative;

        Set_NonPerturbative(false);
        double P_tmdPDF = Evaluate_inbTspace(m_xBjorken, bstar(bT));
        Set_NonPerturbative(true);
        double NP_tmdPDF = Evaluate_inbTspace(m_xBjorken, bT);

        Set_NonPerturbative(in_NonPerturbative);

        double g_qH = -log(NP_tmdPDF/P_tmdPDF);

        return g_qH;
    }

}


// Evaluate order alphaS OPE bstar contribution to TMDPDF
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::EvaluateOPE_qq_bTstar_1(double mustar,
                                                                                   double x,                               gsl_integration_workspace* wsp)
{
    double result, abserr;
    
    double mu = m_RenormalizationScale;
  
    m_collinearQuark.Set_RenormalizationScale(mustar);
    
    auto integrand = make_gsl_function( [&](double xi) {return 2.0*(1.0-xi)*m_collinearQuark.Evaluate(x/xi)/xi;});
    gsl_integration_qag(integrand, x, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6, wsp, &result, &abserr);

    double output =  CF*( result - M_PI*M_PI/6.0*m_collinearQuark.Evaluate(x) );
    
//    std::cout << "m_collinearQuark.Evaluate(x) inside = " << m_collinearQuark.Evaluate(x) << std::endl;
    m_collinearQuark.Set_RenormalizationScale(mu);

    
    return output/(4*M_PI);// to be multiplied by alphaS
}

template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::EvaluateOPE_qg_bTstar_1(double mustar,
                                                                                   double x,                               gsl_integration_workspace* wsp)
{
    double result, abserr;
    double mu = m_RenormalizationScale;

    m_collinearBoson.Set_RenormalizationScale(mustar);

    auto integrand = make_gsl_function( [&](double xi) {return 2.0*(1.0-xi)*m_collinearBoson.Evaluate(x/xi);});
    gsl_integration_qag(integrand, x, 1.0, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, 6, wsp, &result, &abserr);


    double output =  2.0 * TF *  result;
    m_collinearBoson.Set_RenormalizationScale(mu);
    
    return output/(4*M_PI);// to be multiplied by alphaS
}

// OPE evolution factor for TMD PDF in bt space
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::OPE_RG_Evolution_bstar_TMDPDF(double bT, gsl_integration_workspace* wsp)
{
    double bstr = bstar(bT);
    
    double output = m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp, false) * exp(-Evaluate_g_QuarkHadron(bT) - log(m_HardScale/INPUT_HARD_SCALE)*gCSKernel(bT, m_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale), false,true));
    
    //without gK
//    double output = m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp) * exp(-Evaluate_g_QuarkHadron(bT) );
    
    // without g function
//    double output = m_collinearQuark.RG_Evolution_bT(bstr, mubstar(bstr), m_RenormalizationScale, m_HardScale, m_collinearQuark, wsp) * exp( - log(m_HardScale/INPUT_HARD_SCALE)*gCSKernel(bT, m_RenormalizationScale, m_collinearQuark.Get_AlphaStrong(m_RenormalizationScale), false,true));
    
    
    return output;
}

// Conventional TMD pdf in bT space with bstar method
// Quark TMD pdf
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Evaluate_Quark_inbTspace_bstarMethod(double x, double bT, gsl_integration_workspace* wsp)
{
    double output;
    double mustr = mubstar(bstar(bT)); //correct one
    double in_RenormalizationScale = m_collinearQuark.Get_RenormalizationScale();
    
    double alphaS_mustr = m_collinearQuark.Get_AlphaStrong(mustr);
    m_collinearQuark.Set_RenormalizationScale(mustr);// set to Renormalization Scale to mubstar
    m_collinearBoson.Set_RenormalizationScale(mustr);// set to Renormalization Scale to mubstar
    
    bool in_NonPerturbative = Get_NonPerturbative();
    Set_NonPerturbative(false); // model independence to be safe

    double f_OPE = m_collinearQuark.Evaluate(x) + alphaS_mustr*( EvaluateOPE_qq_bTstar_1(mustr, x, wsp) + EvaluateOPE_qg_bTstar_1(mustr, x, wsp) );
    
    Set_NonPerturbative(in_NonPerturbative);// reset to initial value
    
    output = f_OPE * OPE_RG_Evolution_bstar_TMDPDF(bT, wsp);
    
    m_collinearQuark.Set_RenormalizationScale(in_RenormalizationScale);// reset to input value
    m_collinearBoson.Set_RenormalizationScale(in_RenormalizationScale);// reset to input value
    
    return output;
}


// Get average of squared transverse momentum
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Get_mean_squared_transverse_momentum(double max_kT2)// input is maximum squared transverse momentum
{
    double kc2 = max_kT2;
    double Q   = m_HardScale;
    double mi  = m_mass_PDFQuarkHadron;
    double mg  = m_mass_PDFBosonHadron;
    double Mi  = m_Mass_PDFQuarkHadron;
    
    return 0.5*m_A_PDFQuarkHadron*( kc2 - mi*mi*log( 1.0 + kc2/(mi*mi) ) )
         + 0.5*m_B_PDFQuarkHadron*( kc2 + 0.5*mi*mi*pow(log(Q*Q/(kc2 + mi*mi)),2) + (kc2 + mi*mi)*log(Q*Q/(kc2 + mi*mi)) -0.5*mi*mi*log(Q*Q/(mi*mi))*( log(Q*Q/(mi*mi))+ 2.0) )
         + 0.5*m_A_PDFBosonHadron*( kc2 - mg*mg*log( 1.0 + kc2/(mg*mg) ) )
         + m_C_PDFQuarkHadron*(Mi*Mi - (kc2 + Mi*Mi)*exp(-kc2/(Mi*Mi)) );
}


// Coefficients for input TMD PDF ---------------------------------------------------
// Coefficient A_PDFQuarkHadron
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Af_iH (double x, double muQ)
{
    double result;
    double abserr;
    double p[2] = {x, 1.0};
    IntegrationWorkspace wsp1(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    auto integrand = make_gsl_function( [&](double xi)
                                       {
        return ( (1+xi*xi)*m_collinearQuark.Evaluate(x/xi)/xi - 2.0*m_collinearQuark.Evaluate(x) )/(1.0-xi);
                                        });
    gsl_integration_qagp(integrand, p, 2, NUMERICAL_INTEGRATION_ABSOLUTE_TOLERANCE, NUMERICAL_INTEGRATION_RELATIVE_TOLERANCE, NUMERICAL_INTEGRATION_MAX_INTERVALS, wsp1, &result, &abserr);
    
    
    return CF/M_PI * (result + 2.0*m_collinearQuark.Evaluate(x)*log(1.0-x));
}

// Coefficient B_PDFQuarkHadron
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Bf_iH (double x, double muQ)
{
    return CF/M_PI * m_collinearQuark.Evaluate(x);
}

// Coefficient A_PDFBosonHadron
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Af_gH (double x, double muQ)
{
    return 1.0/M_PI *m_collinearBoson.MellinConvolutionPDF(x,Pq_g, &CollinearPDF<t_BosonFlavour,t_HadronFlavour>::Evaluate);
}

// Terms of C coefficient -------------------------------------------------------------------------
// Coefficient C_PDFQuarkHadron order one in alpha strong times alpha strong ----------------------
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Cf_iH_1 (double x, double muQ, double alphaS)
{
    double mi = m_mass_PDFQuarkHadron;
    double mg = m_mass_PDFBosonHadron;
    double Q  = m_collinearQuark.Get_HardScale();
    
//    double output =  -alphaS*M_PI/24.0*CF*m_collinearQuark.Evaluate(x) - m_A_PDFQuarkHadron*log(muQ/mi) - m_B_PDFQuarkHadron*log(muQ/mi)*log(Q*Q/(muQ*mi)) - m_A_PDFBosonHadron*log(muQ/mg) + alphaS/(2*M_PI)*(MellinConvolution(x, CDelta_iH_tilde, m_collinearQuark.Evaluate) + MellinConvolution(x,CDelta_gH, m_collinearBoson.Evaluate));
    double output =  -alphaS*M_PI/24.0*CF*m_collinearQuark.Evaluate(x) - m_A_PDFQuarkHadron*log(muQ/mi) - m_B_PDFQuarkHadron*log(muQ/mi)*log(Q*Q/(muQ*mi)) - m_A_PDFBosonHadron*log(muQ/mg) + alphaS/(2*M_PI)*(m_collinearQuark.MellinConvolutionPDF(x, CDelta_iH_tilde, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate) + m_collinearQuark.MellinConvolutionPDF(x,CDelta_gH, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate));
    
    return output ;
}
// Coefficient C_PDFQuarkHadron order zero in alpha strong ------------------------------------------
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Cf_iH_0 (double x, double muQ)
{
    return m_collinearQuark.Evaluate(x);
}

// Delta Cutoff MSbar
template<int t_QuarkFlavour,int t_BosonFlavour, int t_HadronFlavour>
double TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour>::Delta_MSbar_iH_1(double x, double muQ,double alphaS)
{
    // with Dirac Delta contribution = Standard MSbar
    return -alphaS*M_PI/24.0*CF*m_collinearQuark.Evaluate(x) + alphaS/(2*M_PI)*(m_collinearQuark.MellinConvolutionPDF(x, CDelta_iH_tilde, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate) + m_collinearQuark.MellinConvolutionPDF(x,CDelta_gH, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate));
    
    // without Dirac Delta contribution = Collins MSbar
//    return  alphaS/(2*M_PI)*(m_collinearQuark.MellinConvolutionPDF(x, CDelta_iH_tilde, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate) + m_collinearQuark.MellinConvolutionPDF(x,CDelta_gH, &CollinearPDF<t_QuarkFlavour,t_HadronFlavour>::Evaluate));
}








#endif /* TMDPDF_h */
