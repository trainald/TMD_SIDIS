//
//  main.cpp
//  TMDs
//
//  Created by Tommaso Rainaldi on 11/21/22.
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


std::string TMDFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D);
void PrintTMDHeader(std::ofstream *filename);


constexpr int Target  = UNDEFINED_HADRON;// unknown

int main(int argc, const char * argv[]) {
    
    std::cout << "Program BEGIN" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();// start counting
    
    
    // Initial parameters
    double x  = 0.1;// xBjorken for PDF
    double z  = 0.3;// zh for FF
    double Q  = 2.0;// Hard Scale
    double mu = 2.0;// Renormalization Scale
    
    double m_F = 0.9;// remember to change also ansatzs for gs.
    double m_D = 0.9;// remember to change also ansatzs for gs.
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
    
    
    std::cout << " AlphaS pdf = " << colPDFUp.Get_AlphaStrong(mu) << std::setw(20) << " AlphaS ff = " << colFFUp.Get_AlphaStrong(mu) << std::endl;
    
    // Quarks
    // Up
    std::ofstream file_TMD_UP;
    file_TMD_UP.open (TMDFilename("UP", x, z, Q, mu, m_F, m_D));
//    file_TMD_UP.open (Filename_small("UP", x, z, Q, mu, m));
    file_TMD_UP.precision(8);
    file_TMD_UP.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_TMD_DOWN;
    file_TMD_DOWN.open (TMDFilename("DOWN", x, z, Q, mu, m_F, m_D));
//    file_TMD_DOWN.open (Filename_small("DOWN", x, z, Q, mu, m));
    file_TMD_DOWN.precision(8);
    file_TMD_DOWN.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_TMD_STRANGE;
    file_TMD_STRANGE.open (TMDFilename("STRANGE", x, z, Q, mu, m_F, m_D));
//    file_TMD_STRANGE.open (Filename_small("STRANGE", x, z, Q, mu, m));
    file_TMD_STRANGE.precision(8);
    file_TMD_STRANGE.setf(std::ios::fixed | std::ios::showpoint);
    
    // Antiquarks
    // UBAR
    std::ofstream file_TMD_UBAR;
    file_TMD_UBAR.open (TMDFilename("UBAR", x, z, Q, mu, m_F, m_D));
//    file_TMD_UBAR.open (Filename_small("UBAR", x, z, Q, mu, m));
    file_TMD_UBAR.precision(8);
    file_TMD_UBAR.setf(std::ios::fixed | std::ios::showpoint);
    // DBAR
    std::ofstream file_TMD_DBAR;
    file_TMD_DBAR.open (TMDFilename("DBAR", x, z, Q, mu, m_F, m_D));
//    file_TMD_DBAR.open (Filename_small("DBAR", x, z, Q, mu, m));
    file_TMD_DBAR.precision(8);
    file_TMD_DBAR.setf(std::ios::fixed | std::ios::showpoint);
    // UBAR
    std::ofstream file_TMD_SBAR;
    file_TMD_SBAR.open (TMDFilename("SBAR", x, z, Q, mu, m_F, m_D));
//    file_TMD_SBAR.open (Filename_small("SBAR", x, z, Q, mu, m));
    file_TMD_SBAR.precision(8);
    file_TMD_SBAR.setf(std::ios::fixed | std::ios::showpoint);
    
    
    std::ofstream file_test;
    file_test.open ("test.txt");
    file_test.precision(8);
    file_test.setf(std::ios::fixed | std::ios::showpoint);
    
    // Set initial parameters
    double qT_in    = 1e-3;// starting qT
    double qT       = qT_in;
    double qT_fin   = 1.5*Q;// final qT
    int    qT_steps = 1000;// number of points in qT space
    double qT_Delta = (qT_fin - qT_in)/qT_steps;// step size in qT space
    
    double bT_in    = 1e-6;// starting bT
    double bT       = bT_in;
    double bT_fin   = 100.0;// final bT
    int    bT_steps = qT_steps;// number of points in bT space
//    double bT_Delta = (bT_fin - bT_in)/bT_steps;// step size in bT space
    double bT_Delta = (bT_fin / bT_in);// step size in bT space
    
    
    // HEADER for txt files
    if(PRINT_HEADER)
    {
        // UP
    PrintTMDHeader(&file_TMD_UP);
        // DOWN
    PrintTMDHeader(&file_TMD_DOWN);
        // STRANGE
    PrintTMDHeader(&file_TMD_STRANGE);
        // UBAR
    PrintTMDHeader(&file_TMD_UBAR);
        // DBAR
    PrintTMDHeader(&file_TMD_DBAR);
        // SBAR
    PrintTMDHeader(&file_TMD_SBAR);
        
    }
    
    // Allocate space for GSL numerical integration
    gsl_integration_workspace * wsp = gsl_integration_workspace_alloc(NUMERICAL_INTEGRATION_MAX_INTERVALS);
    
    for(int j = 0; j <= qT_steps; j++)
    {
        qT = qT_in + j*qT_Delta;// update qT
//        bT = bT_in * exp(j*bT_Delta);// update bT
        bT = bT * pow(bT_Delta,1.0/bT_steps);// update bT
        
//        file_test << std::scientific
//        << std::right << std::setfill(' ') << std::setw(20) << bT
//        << std::right << std::setfill(' ') << std::setw(20) << Q0bar(bT)
//        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_inbTspace_improved(x, bT,wsp)
//        << std::endl;
        
        std::cout << "bT = " << bT << std::endl;
        
//        std::cout << "mubstar = " << mubstar(bstar(bT)) << std::endl;
//        std::cout << "AlphaS(mustar) = " << colPDFUp.Get_AlphaStrong(mubstar(bstar(bT))) << std::endl;
        
        // UP
        file_TMD_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
        << std::endl;
        
        // DOWN
        file_TMD_DOWN << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
        << std::endl;
        
        // STRANGE
        file_TMD_STRANGE << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
        << std::endl;
        
        // UBAR
        file_TMD_UBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
        << std::endl;
        
        // DBAR
        file_TMD_DBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
        << std::endl;
        
        // SBAR
        file_TMD_SBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << qT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar_pert.Evaluate(x, qT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar_pert.Evaluate(z, qT)
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate_inbTspace(x, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate_inbTspace(z, bT)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate_inbTspace_improved(x, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate_inbTspace_improved(z, bT,wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate_Quark_inbTspace_bstarMethod(x, bT, wsp)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate_Quark_inbTspace_bstarMethod(z, bT, wsp)
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
    gsl_integration_workspace_free(wsp);
    
    // Duration of the program
    auto stop = std::chrono::high_resolution_clock::now();// start counting
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Program END" << std::endl;
    std::cout << "Executed in " << duration.count() << " seconds !!!" << std::endl;
    
    return 0;
}




std::string TMDFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D)
{
    return "QCD_TMDs" + f + "_x_0_" + std::to_string(int(x*10)) + "_z_0_" + std::to_string(int(z*10)) +  "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*10)) + "_mDoverz_0_" + std::to_string(int(m_D*10)) +".txt";
}


// Print HEADER
void PrintTMDHeader(std::ofstream *filename)
{
    *filename << "#"
    << std::right << std::setfill(' ') << std::setw(20) << "qT"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDPDF"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDFF"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDPDF_pert"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDFF_pert"
    << std::right << std::setfill(' ') << std::setw(20) << "bT"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDPDF bT"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDFF bT"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDPDF bT impr"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDFF bT impr"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDPDF bT OPE"
    << std::right << std::setfill(' ') << std::setw(20) << "TMDFF bT OPE"
    << std::endl;
}
