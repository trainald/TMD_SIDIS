//
//  main.cpp
//  QCD_case
//
//  Created by Tommaso Rainaldi on 11/7/22.
//

#include <iostream>
#include <variant>
#include <fstream>
#include <string>
#include <vector>
#include <memory>


#include "TMDPDF.h"
#include "TMDFF.h"
#include "WtermFunctions.h"
#include "YtermFunctions.h"

#include "Evolution.h"
#include "Hardfactor.h"


//For timing the code
#include <chrono>

// Functions prototypes
void PrintHeader(std::ofstream *filename);
std::string WFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D);
std::string Filename_small(std::string f,double x, double z, double Q, double mu, double m_F);
std::string YFilename(std::string f,double x, double z, double Q, double mu);

std::string OPEWFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax);

// ##########################################################################################################

constexpr int Target  = UNDEFINED_HADRON;// unknown

// ==========================================================================================================

// BEGIN main

int main()
{
    std::cout << "Program BEGIN" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();// start counting
    
    
    // Initial parameters
    double x  = 0.1;// xBjorken for PDF
    double z  = 0.3;// zh for FF
    double Q  = 4.0;// Hard Scale
    double mu = 4.0;// Renormalization Scale
    
    double y  = 0.5; // for cross section
    
    //double m = 0.0001; // only for asymptotic
    double m_F = 0.4;// remember to change also ansatzs for gs.
    double m_D = 0.3;// remember to change also ansatzs for gs.
                                              //m_i   //m_g  //M_i
    InputTMDPDF_massParameters massParamsPDF = {m_F   , m_F ,  m_F};
    InputTMDFF_massParameters  massParamsFF  = {m_D*z , m_D*z ,m_D*z};
    
    
    
    bool NonPerturbative = true;// Flag for Input TMD PDF or FF to compute
    
    bool isLHAPDF = true;// if Collinear distributions are from LHAPDF or user defined
    
    bool PRINT_HEADER = true;// print header in txt file with info

    
    // NOTE: All PDFs and FFs must be initialized before any operation is executed
    
    // Initialize Collinear PDFs
    // Quarks
    CollinearPDF<UP,Target> colPDFUp(mu,Q,isLHAPDF);
    CollinearPDF<DOWN,Target> colPDFDown(mu,Q,isLHAPDF);
    CollinearPDF<STRANGE,Target> colPDFStrange(mu,Q,isLHAPDF);
    // Antiquarks
    CollinearPDF<UBAR,Target> colPDFUbar(mu,Q,isLHAPDF);
    CollinearPDF<DBAR,Target> colPDFDbar(mu,Q,isLHAPDF);
    CollinearPDF<SBAR,Target> colPDFSbar(mu,Q,isLHAPDF);
    // Boson
    CollinearPDF<GLUON,Target> colPDFGluon(mu,Q,isLHAPDF);// Collinear PDF for Gluon

    // Get input non-perturbative TMD PDF
    // Quarks
    TMDPDF<UP,GLUON,Target> TMDpdfUp(colPDFUp, colPDFGluon, x, massParamsPDF, NonPerturbative);
    TMDPDF<DOWN,GLUON,Target> TMDpdfDown(colPDFDown, colPDFGluon, x, massParamsPDF, NonPerturbative);
    TMDPDF<STRANGE,GLUON,Target> TMDpdfStrange(colPDFStrange, colPDFGluon, x, massParamsPDF, NonPerturbative);
    // Antiquarks
    TMDPDF<UBAR,GLUON,Target> TMDpdfUbar(colPDFUbar, colPDFGluon, x, massParamsPDF, NonPerturbative);
    TMDPDF<DBAR,GLUON,Target> TMDpdfDbar(colPDFDbar, colPDFGluon, x, massParamsPDF, NonPerturbative);
    TMDPDF<SBAR,GLUON,Target> TMDpdfSbar(colPDFSbar, colPDFGluon, x, massParamsPDF, NonPerturbative);
    // Get input perturbative TMD PDF
    // Quarks
    TMDPDF<UP,GLUON,Target> TMDpdfUp_pert(colPDFUp, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    TMDPDF<DOWN,GLUON,Target> TMDpdfDown_pert(colPDFDown, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    TMDPDF<STRANGE,GLUON,Target> TMDpdfStrange_pert(colPDFStrange, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    // Antiquarks
    TMDPDF<UBAR,GLUON,Target> TMDpdfUbar_pert(colPDFUbar, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    TMDPDF<DBAR,GLUON,Target> TMDpdfDbar_pert(colPDFDbar, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    TMDPDF<SBAR,GLUON,Target> TMDpdfSbar_pert(colPDFSbar, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    
    // Initialize Collinear FFs
    // Quarks
    CollinearFF<Target,UP> colFFUp(mu,Q, isLHAPDF);
    CollinearFF<Target,DOWN> colFFDown(mu,Q, isLHAPDF);
    CollinearFF<Target,STRANGE> colFFStrange(mu,Q, isLHAPDF);
    // Antiquarks
    CollinearFF<Target,UBAR> colFFUbar(mu,Q, isLHAPDF);
    CollinearFF<Target,DBAR> colFFDbar(mu,Q, isLHAPDF);
    CollinearFF<Target,SBAR> colFFSbar(mu,Q, isLHAPDF);
    // Boson
    CollinearFF<Target,GLUON> colFFGluon(mu,Q, isLHAPDF);
    
    // Get input non-perturbative TMD FF
    // Quarks
    TMDFF<Target,UP,GLUON> TMDffUp(colFFUp, colFFGluon, z, massParamsFF, NonPerturbative);
    TMDFF<Target,DOWN,GLUON> TMDffDown(colFFDown, colFFGluon, z, massParamsFF, NonPerturbative);
    TMDFF<Target,STRANGE,GLUON> TMDffStrange(colFFStrange, colFFGluon, z, massParamsFF, NonPerturbative);
    // Antiquarks
    TMDFF<Target,UBAR,GLUON> TMDffUbar(colFFUbar, colFFGluon, z, massParamsFF, NonPerturbative);
    TMDFF<Target,DBAR,GLUON> TMDffDbar(colFFDbar, colFFGluon, z, massParamsFF, NonPerturbative);
    TMDFF<Target,SBAR,GLUON> TMDffSbar(colFFSbar, colFFGluon, z, massParamsFF, NonPerturbative);
    // Get input perturbative TMD FF
    // Quarks
    TMDFF<Target,UP,GLUON> TMDffUp_pert(colFFUp, colFFGluon, z, massParamsFF, !NonPerturbative);
    TMDFF<Target,DOWN,GLUON> TMDffDown_pert(colFFDown, colFFGluon, z, massParamsFF, !NonPerturbative);
    TMDFF<Target,STRANGE,GLUON> TMDffStrange_pert(colFFStrange, colFFGluon, z, massParamsFF, !NonPerturbative);
    // Antiquarks
    TMDFF<Target,UBAR,GLUON> TMDffUbar_pert(colFFUbar, colFFGluon, z, massParamsFF, !NonPerturbative);
    TMDFF<Target,DBAR,GLUON> TMDffDbar_pert(colFFDbar, colFFGluon, z, massParamsFF, !NonPerturbative);
    TMDFF<Target,SBAR,GLUON> TMDffSbar_pert(colFFSbar, colFFGluon, z, massParamsFF, !NonPerturbative);
    


    
    // Set up the txt files to write on
    
    // W term --------------------------------------------------------------
    
    // Quarks
    // Up
    std::ofstream file_W_UP;
    file_W_UP.open (WFilename("UP", x, z, Q, mu, m_F, m_D));
//    file_W_UP.open (Filename_small("UP", x, z, Q, mu, m));
    file_W_UP.precision(8);
    file_W_UP.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_W_DOWN;
    file_W_DOWN.open (WFilename("DOWN", x, z, Q, mu, m_F, m_D));
//    file_W_DOWN.open (Filename_small("DOWN", x, z, Q, mu, m));
    file_W_DOWN.precision(8);
    file_W_DOWN.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_W_STRANGE;
    file_W_STRANGE.open (WFilename("STRANGE", x, z, Q, mu, m_F, m_D));
//    file_W_STRANGE.open (Filename_small("STRANGE", x, z, Q, mu, m));
    file_W_STRANGE.precision(8);
    file_W_STRANGE.setf(std::ios::fixed | std::ios::showpoint);

    // Antiquarks
    // UBAR
    std::ofstream file_W_UBAR;
    file_W_UBAR.open (WFilename("UBAR", x, z, Q, mu, m_F, m_D));
//    file_W_UBAR.open (Filename_small("UBAR", x, z, Q, mu, m));
    file_W_UBAR.precision(8);
    file_W_UBAR.setf(std::ios::fixed | std::ios::showpoint);
    // DBAR
    std::ofstream file_W_DBAR;
    file_W_DBAR.open (WFilename("DBAR", x, z, Q, mu, m_F, m_D));
//    file_W_DBAR.open (Filename_small("DBAR", x, z, Q, mu, m));
    file_W_DBAR.precision(8);
    file_W_DBAR.setf(std::ios::fixed | std::ios::showpoint);
    // UBAR
    std::ofstream file_W_SBAR;
    file_W_SBAR.open (WFilename("SBAR", x, z, Q, mu, m_F, m_D));
//    file_W_SBAR.open (Filename_small("SBAR", x, z, Q, mu, m));
    file_W_SBAR.precision(8);
    file_W_SBAR.setf(std::ios::fixed | std::ios::showpoint);


    // Sum of all the flavors
    std::ofstream file_W_ALL;
    file_W_ALL.open (WFilename("ALL", x, z, Q, mu, m_F, m_D));
//    file_W_ALL.open (Filename_small("ALL", x, z, Q, mu, m));
    file_W_ALL.precision(8);
    file_W_ALL.setf(std::ios::fixed | std::ios::showpoint);


    // Y term -----------------------------------------------------------------
    // UP
    std::ofstream file_Y_UP;
    file_Y_UP.open (YFilename("UP", x, z, Q, mu));
    file_Y_UP.precision(8);
    file_Y_UP.setf(std::ios::fixed | std::ios::showpoint);

    // DOWN
    std::ofstream file_Y_DOWN;
    file_Y_DOWN.open (YFilename("DOWN", x, z, Q, mu));
    file_Y_DOWN.precision(8);
    file_Y_DOWN.setf(std::ios::fixed | std::ios::showpoint);

    // STRANGE
    std::ofstream file_Y_STRANGE;
    file_Y_STRANGE.open (YFilename("STRANGE", x, z, Q, mu));
    file_Y_STRANGE.precision(8);
    file_Y_STRANGE.setf(std::ios::fixed | std::ios::showpoint);

    // UBAR
    std::ofstream file_Y_UBAR;
    file_Y_UBAR.open (YFilename("UBAR", x, z, Q, mu));
    file_Y_UBAR.precision(8);
    file_Y_UBAR.setf(std::ios::fixed | std::ios::showpoint);

    // DBAR
    std::ofstream file_Y_DBAR;
    file_Y_DBAR.open (YFilename("DBAR", x, z, Q, mu));
    file_Y_DBAR.precision(8);
    file_Y_DBAR.setf(std::ios::fixed | std::ios::showpoint);

    // SBAR
    std::ofstream file_Y_SBAR;
    file_Y_SBAR.open (YFilename("SBAR", x, z, Q, mu));
    file_Y_SBAR.precision(8);
    file_Y_SBAR.setf(std::ios::fixed | std::ios::showpoint);

    // Sum of all the flavors
    std::ofstream file_Y_ALL;
    file_Y_ALL.open (YFilename("ALL", x, z, Q, mu));
    file_Y_ALL.precision(8);
    file_Y_ALL.setf(std::ios::fixed | std::ios::showpoint);
    

    // Set initial parameters
    double qT_in    = 1e-3;// starting qT
    double qT       = qT_in;
    double qT_fin   = 4.0;// final qT
    int    qT_steps = 200;// number of points in qT space
    double qT_Delta = (qT_fin - qT_in)/qT_steps;// step size in qT space

    // Initialize workspace for gsl double integration
    gsl_integration_workspace * wsp1 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    gsl_integration_workspace * wsp2 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);

    
    // Initialize workspace for VEGAS
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    const gsl_rng_type *T;
    gsl_rng *r = 0;
    // Fill VEGAS parameters structure
    // Quarks
    my_Vegas_parameters<UP,GLUON,Target> VEGAS_params_UP = {TMDpdfUp_pert, TMDffUp_pert, qT};
    my_Vegas_parameters<DOWN,GLUON,Target> VEGAS_params_DOWN = {TMDpdfDown_pert, TMDffDown_pert, qT};
    my_Vegas_parameters<STRANGE,GLUON,Target> VEGAS_params_STRANGE = {TMDpdfStrange_pert, TMDffStrange_pert, qT};
    // Antiquarks
    my_Vegas_parameters<UBAR,GLUON,Target> VEGAS_params_UBAR = {TMDpdfUbar_pert, TMDffUbar_pert, qT};
    my_Vegas_parameters<DBAR,GLUON,Target> VEGAS_params_DBAR = {TMDpdfDbar_pert, TMDffDbar_pert, qT};
    my_Vegas_parameters<SBAR,GLUON,Target> VEGAS_params_SBAR = {TMDpdfSbar_pert, TMDffSbar_pert, qT};
    
    // Find charges
    double eUP      = colPDFUp.Get_PartonCharge();
    double eDOWN    = colPDFDown.Get_PartonCharge();
    double eSTRANGE = colPDFStrange.Get_PartonCharge();
    double eUBAR    = colPDFUbar.Get_PartonCharge();
    double eDBAR    = colPDFDbar.Get_PartonCharge();
    double eSBAR    = colPDFSbar.Get_PartonCharge();
    // make array of charges
    double charge[6] = {eUP, eDOWN, eSTRANGE, eUBAR, eDBAR, eSBAR};
    std::cout << charge[0] << charge[1] << charge[2] << charge[3] << charge[4] << charge[5] << std::endl;
    
    // Find the hard factor
    double Hard_Factor = Hardfactor(Q, mu, colFFUp.Get_AlphaStrong(mu));// chosen from FF !!!!!!!!!!!!!
    std::cout << "Hard factor = " << Hard_Factor << std::endl;

    // Declare output variables
    // W term
    double result_Wterm_exact[6];
    double result_Wterm_AsyLO[6];
    double result_Wterm_Correction[6];
    double result_Wterm_exact_bT[6];
    double result_Wterm_OPEbstar[6];
    // Y term
    double result_FixedOrder[6];
    
    

    // HEADER for txt files
    if(PRINT_HEADER)
    {
        // UP
    PrintHeader(&file_W_UP);
        // DOWN
    PrintHeader(&file_W_DOWN);
        // STRANGE
    PrintHeader(&file_W_STRANGE);
        // UBAR
    PrintHeader(&file_W_UBAR);
        // DBAR
    PrintHeader(&file_W_DBAR);
        // SBAR
    PrintHeader(&file_W_SBAR);
        // ALL
    PrintHeader(&file_W_ALL);
        
    }

    // *** Begin loop over qT ***
    for(int j = 0; j <= qT_steps; j++)
    {
        // For ALL flavours
        double Wterm_exact      = 0;
        double Wterm_AsyLO      = 0;
        double Wterm_Correction = 0;
        double Wterm_exact_bT   = 0;
        double Wterm_OPEbstar   = 0;
        
        double Yterm_FixedOrder = 0;
        
        qT = qT_in + j*qT_Delta;// update qT
        VEGAS_params_UP.qT = qT;// update qT in VEGAS structure
        VEGAS_params_DOWN.qT = qT;// update qT in VEGAS structure
        VEGAS_params_STRANGE.qT = qT;// update qT in VEGAS structure
        VEGAS_params_UBAR.qT = qT;// update qT in VEGAS structure
        VEGAS_params_DBAR.qT = qT;// update qT in VEGAS structure
        VEGAS_params_SBAR.qT = qT;// update qT in VEGAS structure
        
        

        // Compute the integrations
        // Up
        // W term
        result_Wterm_exact[0]      = Get_WtermSingleFlavour_exact(TMDpdfUp, TMDffUp, qT, wsp1, wsp2);
        result_Wterm_AsyLO[0]      = Get_WtermSingleFlavour_AsyLO(TMDpdfUp, TMDpdfUp_pert, TMDffUp, TMDffUp_pert, qT);
        result_Wterm_Correction[0] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_UP, s);
        result_Wterm_exact_bT[0]   = Get_WtermSingleFlavour_exact_bT(TMDpdfUp, TMDffUp, qT, wsp1);
        result_Wterm_OPEbstar[0]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFUp, TMDpdfUp_pert, colFFUp, TMDffUp_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[0]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFUp, TMDpdfUp_pert, colFFUp, TMDffUp_pert, wsp1);
        
        
        // Y term
        result_FixedOrder[0]       = colPDFUp.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFUp, colPDFGluon, colFFUp, colFFGluon, x, z, qT, y, wsp1);
        
        // UP
        file_W_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[0]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[0]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[0]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[0]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[0]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[0]
        << std::endl;
        
        file_Y_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[0]
        << std::endl;

        // Down
        result_Wterm_exact[1]      = Get_WtermSingleFlavour_exact(TMDpdfDown, TMDffDown, qT, wsp1, wsp2);
        result_Wterm_AsyLO[1]      = Get_WtermSingleFlavour_AsyLO(TMDpdfDown, TMDpdfDown_pert, TMDffDown, TMDffDown_pert, qT);
        result_Wterm_Correction[1] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_DOWN, s);
        result_Wterm_exact_bT[1]   = Get_WtermSingleFlavour_exact_bT(TMDpdfDown, TMDffDown, qT, wsp1);
        result_Wterm_OPEbstar[1]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFDown, TMDpdfDown_pert, colFFDown, TMDffDown_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[1]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFDown, TMDpdfDown_pert, colFFDown, TMDffDown_pert, wsp1);

        // Y term
        result_FixedOrder[1]       = colPDFDown.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFDown, colPDFGluon, colFFDown, colFFGluon, x, z, qT, y, wsp1);

        // DOWN
        file_W_DOWN << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[1]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[1]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[1]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[1]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[1]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[1]
        << std::endl;

        file_Y_DOWN << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[1]
        << std::endl;

        // Strange
        result_Wterm_exact[2]      = Get_WtermSingleFlavour_exact(TMDpdfStrange, TMDffStrange, qT, wsp1, wsp2);
        result_Wterm_AsyLO[2]      = Get_WtermSingleFlavour_AsyLO(TMDpdfStrange, TMDpdfStrange_pert, TMDffStrange, TMDffStrange_pert, qT);
        result_Wterm_Correction[2] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_STRANGE, s);
        result_Wterm_exact_bT[2]   = Get_WtermSingleFlavour_exact_bT(TMDpdfStrange, TMDffStrange, qT, wsp1);
        result_Wterm_OPEbstar[2]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFStrange, TMDpdfStrange_pert, colFFStrange, TMDffStrange_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[2]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFStrange, TMDpdfStrange_pert, colFFStrange, TMDffStrange_pert, wsp1);

        // Y term
        result_FixedOrder[2]       = colPDFStrange.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFStrange, colPDFGluon, colFFStrange, colFFGluon, x, z, qT, y, wsp1);

        // STRANGE
        file_W_STRANGE << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[2]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[2]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[2]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[2]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[2]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[2]
        << std::endl;

        file_Y_STRANGE << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[2]
        << std::endl;

        // Ubar
        result_Wterm_exact[3]      = Get_WtermSingleFlavour_exact(TMDpdfUbar, TMDffUbar, qT, wsp1, wsp2);
        result_Wterm_AsyLO[3]      = Get_WtermSingleFlavour_AsyLO(TMDpdfUbar, TMDpdfUbar_pert, TMDffUbar, TMDffUbar_pert, qT);
        result_Wterm_Correction[3] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_UBAR, s);
        result_Wterm_exact_bT[3]   = Get_WtermSingleFlavour_exact_bT(TMDpdfUbar, TMDffUbar, qT, wsp1);
        result_Wterm_OPEbstar[3]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFUbar, TMDpdfUbar_pert, colFFUbar, TMDffUbar_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[3]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFUbar, TMDpdfUbar_pert, colFFUbar, TMDffUbar_pert, wsp1);

        // Y term
        result_FixedOrder[3]       = colPDFUbar.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFUbar, colPDFGluon, colFFUbar, colFFGluon, x, z, qT, y, wsp1);

        // UBAR
        file_W_UBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[3]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[3]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[3]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[3]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[3]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[3]
        << std::endl;

        file_Y_UBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[3]
        << std::endl;

        // Dbar
        result_Wterm_exact[4]      = Get_WtermSingleFlavour_exact(TMDpdfDbar, TMDffDbar, qT, wsp1, wsp2);
        result_Wterm_AsyLO[4]      = Get_WtermSingleFlavour_AsyLO(TMDpdfDbar, TMDpdfDbar_pert, TMDffDbar, TMDffDbar_pert, qT);
        result_Wterm_Correction[4] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_DBAR, s);
        result_Wterm_exact_bT[4]   = Get_WtermSingleFlavour_exact_bT(TMDpdfDbar, TMDffDbar, qT, wsp1);
        result_Wterm_OPEbstar[4]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFDbar, TMDpdfDbar_pert, colFFDbar, TMDffDbar_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[4]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFDbar, TMDpdfDbar_pert, colFFDbar, TMDffDbar_pert, wsp1);

        // Y term
        result_FixedOrder[4]       = colPDFDbar.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFDbar, colPDFGluon, colFFDbar, colFFGluon, x, z, qT, y, wsp1);

        // DBAR
        file_W_DBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[4]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[4]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[4]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[4]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[4]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[4]
        << std::endl;

        file_Y_DBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[4]
        << std::endl;

        // Sbar
        result_Wterm_exact[5]      = Get_WtermSingleFlavour_exact(TMDpdfSbar, TMDffSbar, qT, wsp1, wsp2);
        result_Wterm_AsyLO[5]      = Get_WtermSingleFlavour_AsyLO(TMDpdfSbar, TMDpdfSbar_pert, TMDffSbar, TMDffSbar_pert, qT);
        result_Wterm_Correction[5] = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_SBAR, s);
        result_Wterm_exact_bT[5]   = Get_WtermSingleFlavour_exact_bT(TMDpdfSbar, TMDffSbar, qT, wsp1);
        result_Wterm_OPEbstar[5]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFSbar, TMDpdfSbar_pert, colFFSbar, TMDffSbar_pert, wsp1, wsp2);
//        result_Wterm_OPEbstar[5]   = Get_WtermSingleFlavour_OPEbstar_integrand(qT, mubstar(bstar(qT)),colPDFSbar, TMDpdfSbar_pert, colFFSbar, TMDffSbar_pert, wsp1);

        // Y term
        result_FixedOrder[5]       = colPDFSbar.Get_AlphaStrong(mu) * Get_TMDFixedOrder_unpolarized(colPDFSbar, colPDFGluon, colFFSbar, colFFGluon, x, z, qT, y, wsp1);

        // SBAR
        file_W_SBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact[5]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_AsyLO[5]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_Correction[5]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_exact_bT[5]
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[5]
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[5]
        << std::endl;

        file_Y_SBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_FixedOrder[5]
        << std::endl;


        // Sum contributions from all flavours
        for(int i = 0; i < sizeof(charge)/sizeof(charge[0]); i++)
        {
            Wterm_exact      += pow(charge[i],2) * result_Wterm_exact[i];
            Wterm_AsyLO      += pow(charge[i],2) * result_Wterm_AsyLO[i];
            Wterm_Correction += pow(charge[i],2) * result_Wterm_Correction[i];
            Wterm_exact_bT   += pow(charge[i],2) * result_Wterm_exact_bT[i];
            Wterm_OPEbstar   += pow(charge[i],2) * result_Wterm_OPEbstar[i];
            Yterm_FixedOrder += result_FixedOrder[i];
        }
        // ALL
        file_W_ALL << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_exact
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_AsyLO
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_Correction
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_exact_bT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_OPEbstar
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_OPEbstar
        << std::endl;

        file_Y_ALL << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Yterm_FixedOrder
        << std::endl;
    
        
        
        

        // Display progress
        std::cout << "#"+ std::string((j)*100.0/qT_steps, '*');// update loading bar
        std::cout << ">" ;
        std::cout << " Complete " << (double(j))*100.0/double(qT_steps) << " % " << std::endl;// update percentage
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(now - start);
        std::cout << elapsed_time.count() << " seconds" << std::endl;

    }
    
    std::cout << "Hard factor = " << Hard_Factor << std::endl;

    // Free all allocated space
    gsl_integration_workspace_free(wsp1);
    gsl_integration_workspace_free(wsp2);
    gsl_monte_vegas_free (s);


    // Duration of the program
    auto stop = std::chrono::high_resolution_clock::now();// start counting
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Program END" << std::endl;
    std::cout << "Executed in " << duration.count() << " seconds !!!" << std::endl;
    
    return 0;
}

// ##########################################################################################################

// END main


// W term
std::string WFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D)
{
    return "QCDWterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*10)) + "_mDoverz_0_" + std::to_string(int(m_D*10)) +".txt";
}

std::string Filename_small(std::string f,double x, double z, double Q, double mu, double m)
{
    return "QCDWterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_m_0_000" + std::to_string(int(m*10000)) +".txt";
}

// Y term
std::string YFilename(std::string f,double x, double z, double Q, double mu)
{
    return "QCDYterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) +".txt";
}

// OPE W term
std::string OPEWFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax)
{
    return "QCDOPEWterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*10)) + "_mDoverz_0_" + std::to_string(int(m_D*10)) + "_bmax_" + std::to_string(int(bmax)) +".txt";
}


// Print HEADER
void PrintHeader(std::ofstream *filename)
{
    *filename << "#"
    << std::right << std::setfill(' ') << std::setw(20) << "qT"
    << std::right << std::setfill(' ') << std::setw(20) << "W"
    << std::right << std::setfill(' ') << std::setw(20) << "AsyLO"
    << std::right << std::setfill(' ') << std::setw(20) << "Correction"
    << std::right << std::setfill(' ') << std::setw(20) << "W_bT"
    << std::right << std::setfill(' ') << std::setw(20) << "OPEbstar"
    << std::endl;
}





