//
//  main.cpp
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/21/22.
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


//For timing the code
#include <chrono>

// Functions prototypes
template<int t_QuarkFlavour, int t_BosonFlavour, int t_HadronFlavour>
void PrintInputTMD_PDF_FF_Coefficients(TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> tmdpdf,
                                       TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> tmdff);
void PrintHeader(std::ofstream *filename);

// ##########################################################################################################

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
                                              //m_i // m_g // M_i
    InputTMDPDF_massParameters massParamsPDF = {0.1,   0.1,   0.1};
    InputTMDFF_massParameters  massParamsFF  = {0.1*z,   0.1*z,   0.1*z};
    
    bool NonPerturbative = true;// Flag for Input TMD PDF or FF to compute
    
    bool isLHAPDF = true;// if Collinear distributions are from LHAPDF or user defined
    
    bool PRINT_HEADER = true;// print header in txt file with info
    
    // Constant expressions known at compiler time
    constexpr int Target  = UNDEFINED_HADRON;// unknown
    constexpr int Parton1 = UP;
    constexpr int Parton2 = DOWN;
    constexpr int Gluon   = GLUON;
    
    // NOTE: All PDFs and FFs must be initialized before any operation is executed
    
    // Store polymorphic types
//    struct CollinearPDF_list {CollinearPDF<Parton1,Target> colPDFUp(double,double,bool);// Collinear PDF for Parton1
//                              CollinearPDF<Parton2,Target> colPDFDown(double,double,bool);// Collinear PDF for Parton2
//                              CollinearPDF<Gluon,Target> colPDFGluon(double,double,bool);// Collinear PDF for Gluon};
//                              };
//    CollinearPDF_list myCollinearPDFlist = {colPDFUp(mu,Q,isLHAPDF), colPDFDown(mu,Q,isLHAPDF), colPDFGluon(mu,Q,isLHAPDF) };
    //std::variant< CollinearPDF<Gluon,Target>, CollinearPDF<Parton1,Target>,CollinearPDF<Parton2,Target> > colPDFGluon, colPDFUp, colPDFDown;
    
    // Initialize Collinear PDFs
    CollinearPDF<Parton1,Target> colPDFUp(mu,Q,isLHAPDF);// Collinear PDF for Parton1
    CollinearPDF<Parton2,Target> colPDFDown(mu,Q,isLHAPDF);// Collinear PDF for Parton2
    CollinearPDF<Gluon,Target> colPDFGluon(mu,Q,isLHAPDF);// Collinear PDF for Gluon
    // Get input non-perturbative TMD PDF
    TMDPDF<Parton1,Gluon,Target> TMDpdfUp(colPDFUp, colPDFGluon, x, massParamsPDF, NonPerturbative);
    TMDPDF<Parton2,Gluon,Target> TMDpdfDown(colPDFDown, colPDFGluon, x, massParamsPDF, NonPerturbative);
    // Get input perturbative TMD PDF
    TMDPDF<Parton1,Gluon,Target> TMDpdfUp_pert(colPDFUp, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    TMDPDF<Parton2,Gluon,Target> TMDpdfDown_pert(colPDFDown, colPDFGluon, x, massParamsPDF, !NonPerturbative);
    
    // Initialize Collinear FFs
    CollinearFF<Target,Parton1> colFFUp(mu,Q, isLHAPDF);// Collinear FF for Parton1
    CollinearFF<Target,Parton2> colFFDown(mu,Q, isLHAPDF);// Collinear FF for Parton2
    CollinearFF<Target,Gluon> colFFGluon(mu,Q, isLHAPDF);// Collinear FF for Gluon
    // Get input non-perturbative TMD FF
    TMDFF<Target,Parton1,Gluon> TMDffUp(colFFUp, colFFGluon, z, massParamsFF, NonPerturbative);
    TMDFF<Target,Parton2,Gluon> TMDffDown(colFFDown, colFFGluon, z, massParamsFF, NonPerturbative);
    // Get input perturbative TMD FF
    TMDFF<Target,Parton1,Gluon> TMDffUp_pert(colFFUp, colFFGluon, z, massParamsFF, !NonPerturbative);
    TMDFF<Target,Parton2,Gluon> TMDffDown_pert(colFFDown, colFFGluon, z, massParamsFF, !NonPerturbative);
    
    
//    // Set ansatzes for g functions
//    TMDpdfUp.Set_gAnsatz(false);
//    TMDffUp.Set_gAnsatz(false);


    
    
    // Set up the txt files to write on
    // Up
    std::ofstream file_W_UP;
    file_W_UP.open ("WtermUP.txt");
    file_W_UP.precision(8);
    file_W_UP.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_W_DOWN;
    file_W_DOWN.open ("WtermDOWN.txt");
    file_W_DOWN.precision(8);
    file_W_DOWN.setf(std::ios::fixed | std::ios::showpoint);
    
    // Set initial parameters
    double qT_in    = 0.0;// starting qT
    double qT       = qT_in;
    double qT_fin   = 10.0;// final qT
    int    qT_steps = 200;// number of points in qT space
    double qT_Delta = (qT_fin - qT_in)/qT_steps;// step size in qT space
    
    // Initialize workspace for gsl double integration
    gsl_integration_workspace * wsp1 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    gsl_integration_workspace * wsp2 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    
    std::cout << 4.0*M_PI*TMDffUp.EvaluateOPE_gq_bTstar_1(mubstar(bstar(0.5)), x, wsp1) << std::endl;
    
    // Initialize workspace for VEGAS
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    const gsl_rng_type *T;
    gsl_rng *r = 0;
    // Fill VEGAS parameters structure
    my_Vegas_parameters<Parton1,Gluon,Target> VEGAS_params_UP = {TMDpdfUp_pert, TMDffUp_pert, qT};
    my_Vegas_parameters<Parton2,Gluon,Target> VEGAS_params_DOWN = {TMDpdfDown_pert, TMDffDown_pert, qT};
    
    // Declare output variables
    // Up
    double result_Wterm_exact_UP;
    double result_Wterm_AsyLO_UP;
    double result_Wterm_Correction_UP;
    double result_Wterm_exact_bT_UP;
    double result_Wterm_OPEbstar_UP;
    // Down
    // Declare output variables
    double result_Wterm_exact_DOWN;
    double result_Wterm_AsyLO_DOWN;
    double result_Wterm_Correction_DOWN;
    double result_Wterm_exact_bT_DOWN;
    double result_Wterm_OPEbstar_DOWN;
    
    // HEADER for txt files
    if(PRINT_HEADER)
    {
        // UP
    PrintHeader(&file_W_UP);
        // DOWN
    PrintHeader(&file_W_DOWN);
    }
    
    // *** Begin loop over qT ***
    for(int j = 0; j <= qT_steps; j++)
    {
        qT = qT_in + j*qT_Delta;// update qT
        VEGAS_params_UP.qT = qT;// update qT in VEGAS structure
        VEGAS_params_DOWN.qT = qT;// update qT in VEGAS structure
        
        // Compute the integrations
        // Up
        result_Wterm_exact_UP      = Get_WtermSingleFlavour_exact(TMDpdfUp, TMDffUp, qT, wsp1, wsp2);
        result_Wterm_AsyLO_UP      = Get_WtermSingleFlavour_AsyLO(TMDpdfUp, TMDpdfUp_pert, TMDffUp, TMDffUp_pert, qT);
        result_Wterm_Correction_UP = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_UP, s);
        result_Wterm_exact_bT_UP   = Get_WtermSingleFlavour_exact_bT(TMDpdfUp, TMDffUp, qT, wsp1);
        result_Wterm_OPEbstar_UP   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFUp, TMDpdfUp_pert, colFFUp, TMDffUp_pert, wsp1, wsp2);
        
        // Down
        result_Wterm_exact_DOWN      = Get_WtermSingleFlavour_exact(TMDpdfDown, TMDffDown, qT, wsp1, wsp2);
        result_Wterm_AsyLO_DOWN      = Get_WtermSingleFlavour_AsyLO(TMDpdfDown, TMDpdfDown_pert, TMDffDown, TMDffDown_pert, qT);
        result_Wterm_Correction_DOWN = Get_WtermSingleFlavour_Correction(T, r, VEGAS_params_DOWN, s);
        result_Wterm_exact_bT_DOWN   = Get_WtermSingleFlavour_exact_bT(TMDpdfDown, TMDffDown, qT, wsp1);
        result_Wterm_OPEbstar_DOWN   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFDown, TMDpdfDown_pert, colFFDown, TMDffDown_pert, wsp1, wsp2);
        
        // Write results on external txt file
        // Up
        file_W_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_exact_UP
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_AsyLO_UP
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_Correction_UP
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_exact_bT_UP
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_OPEbstar_UP
        << std::endl;
        // Down
        file_W_DOWN << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_exact_DOWN
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_AsyLO_DOWN
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_Correction_DOWN
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_exact_bT_DOWN
        << std::right << std::setfill(' ') << std::setw(20) << result_Wterm_OPEbstar_DOWN
        << std::endl;
        
        // Display progress
        std::cout << "#"+ std::string((j)*100.0/qT_steps, '*');// update loading bar
        std::cout << ">" ;
        std::cout << " Complete " << (double(j))*100.0/double(qT_steps) << " % " << std::endl;// update percentage
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(now - start);
        std::cout << elapsed_time.count() << " seconds" << std::endl;
        
    }
    
    // Free all allocated space
    gsl_integration_workspace_free(wsp1);
    gsl_integration_workspace_free(wsp2);
    gsl_monte_vegas_free (s);
    
    
    std::cout << "alphaS = " << colPDFUp.Get_AlphaStrong(mu) << std::endl;
    // Print coefficients of TMD PDF and TMD FF
    PrintInputTMD_PDF_FF_Coefficients(TMDpdfUp, TMDffUp);// UP (NonPerturbative = true)
    PrintInputTMD_PDF_FF_Coefficients(TMDpdfDown, TMDffDown);// Down (NonPerturbative = true)
    
    // Duration of the program
    auto stop = std::chrono::high_resolution_clock::now();// start counting
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Program END" << std::endl;
    std::cout << "Executed in " << duration.count() << " seconds !!!" << std::endl;
    
    return 0;
}

// ##########################################################################################################

// END main





// Print coefficients
template<int t_QuarkFlavour, int t_BosonFlavour, int t_HadronFlavour>
void PrintInputTMD_PDF_FF_Coefficients(TMDPDF<t_QuarkFlavour,t_BosonFlavour,t_HadronFlavour> tmdpdf,
                                       TMDFF<t_HadronFlavour,t_QuarkFlavour,t_BosonFlavour> tmdff)
{
    double x = tmdpdf.Get_xBjorken();
    double z = tmdff.Get_z();
    double mupdf = tmdpdf.Get_RenormalizationScale();
    double muff  = tmdff.Get_RenormalizationScale();
    std::cout << "----------------------------" << std::endl;
    std::cout << " Input TMD PDF Coefficients " << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Af_" <<t_QuarkFlavour << "," << t_HadronFlavour << " (x = " << x << ", mu = " << mupdf << ") = " << tmdpdf.Get_A_PDFQuarkHadron() << std::endl;
    std::cout << "Bf_" <<t_QuarkFlavour << "," << t_HadronFlavour << " (x = " << x << ", mu = " << mupdf << ") = " << tmdpdf.Get_B_PDFQuarkHadron() << std::endl;
    std::cout << "Af_" <<t_BosonFlavour << "," << t_HadronFlavour << " (x = " << x << ", mu = " << mupdf << ") = " << tmdpdf.Get_A_PDFBosonHadron() << std::endl;
    std::cout << "Cf_" <<t_QuarkFlavour << "," << t_HadronFlavour << " (x = " << x << ", mu = " << mupdf << ") = " << tmdpdf.Get_C_PDFQuarkHadron() << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << " Input TMD FF Coefficients  " << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "AD_" <<t_HadronFlavour << "," << t_QuarkFlavour << " (z = " << z << ", mu = " << muff << ") = " << tmdff.Get_A_FFQuarkHadron() << std::endl;
    std::cout << "BD_" <<t_HadronFlavour << "," << t_QuarkFlavour << " (z = " << z << ", mu = " << muff << ") = " << tmdff.Get_B_FFQuarkHadron() << std::endl;
    std::cout << "AD_" <<t_HadronFlavour << "," << t_BosonFlavour << " (z = " << z << ", mu = " << muff << ") = " << tmdff.Get_A_FFBosonHadron() << std::endl;
    std::cout << "CD_" <<t_HadronFlavour << "," << t_QuarkFlavour << " (z = " << z << ", mu = " << muff << ") = " << tmdff.Get_C_FFQuarkHadron() << std::endl;
    std::cout << "----------------------------" << std::endl;
    
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




