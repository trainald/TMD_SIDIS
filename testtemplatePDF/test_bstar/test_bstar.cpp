//
//  main.cpp
//  test_bstar
//
//  Created by Tommaso Rainaldi on 2/24/23.
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
std::string W_bstarFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax);
std::string WTMD_bstarFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax);
void PrintHeader_TMDW(std::ofstream *filename);

// ##########################################################################################################

constexpr int Target  = UNDEFINED_HADRON;// unknown

// ==========================================================================================================

// BEGIN main

int main(int argc, const char * argv[]) {
    
    std::cout << "Program BEGIN" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();// start counting
    
    
    // Initial parameters
    double x  = 0.1;// xBjorken for PDF
    double z  = 0.3;// zh for FF
    double Q  = 4.0;// Hard Scale
    double mu = 4.0;// Renormalization Scale
    
    double bmax  = 0.1; // GeV-1. Change according to what is found in bstar.h
    
    //double m = 0.0001; // only for asymptotic
    double m_F = 0.25;// remember to change also ansatzs for gs.
    double m_D = 0.25;// remember to change also ansatzs for gs.
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
    
    // Quarks
    // Up
    std::ofstream file_W_bstar_UP;
    file_W_bstar_UP.open (W_bstarFilename("UP", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_UP.precision(8);
    file_W_bstar_UP.setf(std::ios::fixed | std::ios::showpoint);
    
    // Down
    std::ofstream file_W_bstar_DOWN;
    file_W_bstar_DOWN.open (W_bstarFilename("DOWN", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_DOWN.precision(8);
    file_W_bstar_DOWN.setf(std::ios::fixed | std::ios::showpoint);
    
    // Strange
    std::ofstream file_W_bstar_STRANGE;
    file_W_bstar_STRANGE.open (W_bstarFilename("STRANGE", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_STRANGE.precision(8);
    file_W_bstar_STRANGE.setf(std::ios::fixed | std::ios::showpoint);
    
    // Ubar
    std::ofstream file_W_bstar_UBAR;
    file_W_bstar_UBAR.open (W_bstarFilename("UBAR", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_UBAR.precision(8);
    file_W_bstar_UBAR.setf(std::ios::fixed | std::ios::showpoint);
    
    // Dbar
    std::ofstream file_W_bstar_DBAR;
    file_W_bstar_DBAR.open (W_bstarFilename("DBAR", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_DBAR.precision(8);
    file_W_bstar_DBAR.setf(std::ios::fixed | std::ios::showpoint);
    
    // Sbar
    std::ofstream file_W_bstar_SBAR;
    file_W_bstar_SBAR.open (W_bstarFilename("SBAR", x, z, Q, mu, m_F, m_D, bmax));
    file_W_bstar_SBAR.precision(8);
    file_W_bstar_SBAR.setf(std::ios::fixed | std::ios::showpoint);
    
    // Sum of all the flavors
    std::ofstream file_WTMD_bstar_ALL;
    file_WTMD_bstar_ALL.open (WTMD_bstarFilename("ALL", x, z, Q, mu, m_F, m_D, bmax));
    file_WTMD_bstar_ALL.precision(8);
    file_WTMD_bstar_ALL.setf(std::ios::fixed | std::ios::showpoint);
    
    PrintHeader_TMDW(&file_WTMD_bstar_ALL);
    
    
    // Set initial parameters
    double qT_in    = 1e-3;// starting qT
    double qT       = qT_in;
    double qT_fin   = 4.0;// final qT
    int    qT_steps = 500;// number of points in qT space
    double qT_Delta = (qT_fin - qT_in)/qT_steps;// step size in qT space
    
    double bT_in    = 1e-6;// starting bT
    double bT       = bT_in;
    double bT_fin   = 100.0;// final bT
    int    bT_steps = qT_steps;// number of points in bT space
    double bT_Delta = (bT_fin / bT_in);// step size in bT space
    
    // Initialize workspace for gsl double integration
    gsl_integration_workspace * wsp1 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    gsl_integration_workspace * wsp2 = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);

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
    double result_Wterm_OPEbstar[6];
    
    // *** Begin loop over qT ***
    for(int j = 0; j <= qT_steps; j++)
    {
        double Wterm_OPEbstar  = 0.0;
        
        // update qT
        qT = qT_in + j*qT_Delta;// update qT
        bT = bT * pow(bT_Delta,1.0/bT_steps);// update bT
        
        
        // Compute the W term with the OPE + bstar method
        
        // UP
        result_Wterm_OPEbstar[0]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFUp, TMDpdfUp_pert, colFFUp, TMDffUp_pert, wsp1, wsp2);
        
//        file_W_bstar_UP << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[0]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // DOWN
        result_Wterm_OPEbstar[1]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFDown, TMDpdfDown_pert, colFFDown, TMDffDown_pert, wsp1, wsp2);
        
//        file_W_bstar_DOWN << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[1]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // STRANGE
        result_Wterm_OPEbstar[2]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFStrange, TMDpdfStrange_pert, colFFStrange, TMDffStrange_pert, wsp1, wsp2);
        
//        file_W_bstar_STRANGE << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[2]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // UBAR
        result_Wterm_OPEbstar[3]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFUbar, TMDpdfUbar_pert, colFFUbar, TMDffUbar_pert, wsp1, wsp2);
        
//        file_W_bstar_UBAR << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[3]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // DBAR
        result_Wterm_OPEbstar[4]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFDbar, TMDpdfDbar_pert, colFFDbar, TMDffDbar_pert, wsp1, wsp2);
        
//        file_W_bstar_DBAR << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[4]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // SBAR
        result_Wterm_OPEbstar[5]   = Get_WtermSingleFlavour_OPEbstar(qT, colPDFSbar, TMDpdfSbar_pert, colFFSbar, TMDffSbar_pert, wsp1, wsp2);
        
//        file_W_bstar_SBAR << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << qT
//        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * result_Wterm_OPEbstar[5]
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp1)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp1)
//        << std::endl;
        
        // Sum contributions from all flavours
        for(int i = 0; i < sizeof(charge)/sizeof(charge[0]); i++)
        {
            Wterm_OPEbstar   += pow(charge[i],2) * result_Wterm_OPEbstar[i];
        }
        
        file_WTMD_bstar_ALL << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << Hard_Factor * Wterm_OPEbstar
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

    
    // Duration of the program
    auto stop = std::chrono::high_resolution_clock::now();// start counting
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Program END" << std::endl;
    std::cout << "Executed in " << duration.count() << " seconds !!!" << std::endl;
    
    return 0;
}


// OPE single flavor term bstar method
std::string W_bstarFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax)
{
    return "QCD_OPE_SingleFlavor_Wterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*100)) + "_mDoverz_0_" + std::to_string(int(m_D*100)) + "_10bmax_" + std::to_string(int(bmax*10)) +".txt";
}

std::string WTMD_bstarFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D, double bmax)
{
    return "QCD_OPE_WTMDterm" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*100)) + "_mDoverz_0_" + std::to_string(int(m_D*100)) + "_10bmax_" + std::to_string(int(bmax*10)) +".txt";
}

// Print HEADER
void PrintHeader_TMDW(std::ofstream *filename)
{
    *filename << "#"
    << std::right << std::setfill(' ') << std::setw(20) << "qT"
    << std::right << std::setfill(' ') << std::setw(20) << "W TMD"
    << std::endl;
}
