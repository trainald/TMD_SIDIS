//
//  main.cpp
//  CutoffCollinears
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



//For timing the code
#include <chrono>


std::string CutoffCollinearFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D);
void PrintCutoffCollinearHeader(std::ofstream *filename);


constexpr int Target  = UNDEFINED_HADRON;// unknown

int main(int argc, const char * argv[])
{
    
    std::cout << "Program BEGIN" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();// start counting
    
    
    // Initial parameters
    double x  = 1E-2;// xBjorken for PDF
    double z  = 1E-2;// zh for FF
    double Q  = 4.0;// Hard Scale
    double mu = 4.0;// Renormalization Scale
    
    double kc = 1.0*Q;
    
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
    
    std::cout << "AlphaS PDF U  = " << colPDFUp.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF U  = " << colFFUp.Get_AlphaStrong(mu) << std::endl;
    std::cout << "AlphaS PDF D  = " << colPDFDown.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF D  = " << colFFDown.Get_AlphaStrong(mu) << std::endl;
    std::cout << "AlphaS PDF S  = " << colPDFStrange.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF S  = " << colFFStrange.Get_AlphaStrong(mu) << std::endl;
    std::cout << "AlphaS PDF Ub = " << colPDFUbar.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF Ub = " << colFFUbar.Get_AlphaStrong(mu) << std::endl;
    std::cout << "AlphaS PDF Db = " << colPDFDbar.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF Db = " << colFFDbar.Get_AlphaStrong(mu) << std::endl;
    std::cout << "AlphaS PDF Sb = " << colPDFSbar.Get_AlphaStrong(mu) << std::setw(10)<< "  AlphaS FF Sb = " << colFFSbar.Get_AlphaStrong(mu) << std::endl;
    
    // Quarks
    // Up
    std::ofstream file_CutoffCollinear_UP;
    file_CutoffCollinear_UP.open (CutoffCollinearFilename("UP", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_UP.open (Filename_small("UP", x, z, Q, mu, m));
    file_CutoffCollinear_UP.precision(8);
    file_CutoffCollinear_UP.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_CutoffCollinear_DOWN;
    file_CutoffCollinear_DOWN.open (CutoffCollinearFilename("DOWN", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_DOWN.open (Filename_small("DOWN", x, z, Q, mu, m));
    file_CutoffCollinear_DOWN.precision(8);
    file_CutoffCollinear_DOWN.setf(std::ios::fixed | std::ios::showpoint);
    // Down
    std::ofstream file_CutoffCollinear_STRANGE;
    file_CutoffCollinear_STRANGE.open (CutoffCollinearFilename("STRANGE", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_STRANGE.open (Filename_small("STRANGE", x, z, Q, mu, m));
    file_CutoffCollinear_STRANGE.precision(8);
    file_CutoffCollinear_STRANGE.setf(std::ios::fixed | std::ios::showpoint);
    
    // Antiquarks
    // UBAR
    std::ofstream file_CutoffCollinear_UBAR;
    file_CutoffCollinear_UBAR.open (CutoffCollinearFilename("UBAR", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_UBAR.open (Filename_small("UBAR", x, z, Q, mu, m));
    file_CutoffCollinear_UBAR.precision(8);
    file_CutoffCollinear_UBAR.setf(std::ios::fixed | std::ios::showpoint);
    // DBAR
    std::ofstream file_CutoffCollinear_DBAR;
    file_CutoffCollinear_DBAR.open (CutoffCollinearFilename("DBAR", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_DBAR.open (Filename_small("DBAR", x, z, Q, mu, m));
    file_CutoffCollinear_DBAR.precision(8);
    file_CutoffCollinear_DBAR.setf(std::ios::fixed | std::ios::showpoint);
    // UBAR
    std::ofstream file_CutoffCollinear_SBAR;
    file_CutoffCollinear_SBAR.open (CutoffCollinearFilename("SBAR", x, z, Q, mu, m_F, m_D));
//    file_CutoffCollinear_SBAR.open (Filename_small("SBAR", x, z, Q, mu, m));
    file_CutoffCollinear_SBAR.precision(8);
    file_CutoffCollinear_SBAR.setf(std::ios::fixed | std::ios::showpoint);
    
    // Set initial parameters
    double x_in    = x;// starting qT
    x       = x_in;
    double x_fin   = 1.0;// final qT
    int    x_steps = 200;// number of points in qT space
    double x_Delta = (x_fin - x_in)/x_steps;// step size in qT space
    
    
    // HEADER for txt files
    if(PRINT_HEADER)
    {
        // UP
    PrintCutoffCollinearHeader(&file_CutoffCollinear_UP);
        // DOWN
    PrintCutoffCollinearHeader(&file_CutoffCollinear_DOWN);
        // STRANGE
    PrintCutoffCollinearHeader(&file_CutoffCollinear_STRANGE);
        // UBAR
    PrintCutoffCollinearHeader(&file_CutoffCollinear_UBAR);
        // DBAR
    PrintCutoffCollinearHeader(&file_CutoffCollinear_DBAR);
        // SBAR
    PrintCutoffCollinearHeader(&file_CutoffCollinear_SBAR);
        
    }
    
    
    for(int j = 0; j <= x_steps; j++)
    {
        
        x = x_in + j*x_Delta;// update qT
        z = x;
        
        // UPDATES
        // UP
        TMDpdfUp.Set_xBjorken(x);
        TMDpdfUp.Set_A_PDFQuarkHadron();
        TMDpdfUp.Set_B_PDFQuarkHadron();
        TMDpdfUp.Set_A_PDFBosonHadron();
        TMDpdfUp.Set_C_PDFQuarkHadron();
        
        TMDffUp.Set_z(z);
        TMDffUp.Set_A_FFQuarkHadron();
        TMDffUp.Set_B_FFQuarkHadron();
        TMDffUp.Set_A_FFBosonHadron();
        TMDffUp.Set_C_FFQuarkHadron();
        
        // DOWN
        TMDpdfDown.Set_xBjorken(x);
        TMDpdfDown.Set_A_PDFQuarkHadron();
        TMDpdfDown.Set_B_PDFQuarkHadron();
        TMDpdfDown.Set_A_PDFBosonHadron();
        TMDpdfDown.Set_C_PDFQuarkHadron();
        
        TMDffDown.Set_z(z);
        TMDffDown.Set_A_FFQuarkHadron();
        TMDffDown.Set_B_FFQuarkHadron();
        TMDffDown.Set_A_FFBosonHadron();
        TMDffDown.Set_C_FFQuarkHadron();
        
        // STRANGE
        TMDpdfStrange.Set_xBjorken(x);
        TMDpdfStrange.Set_A_PDFQuarkHadron();
        TMDpdfStrange.Set_B_PDFQuarkHadron();
        TMDpdfStrange.Set_A_PDFBosonHadron();
        TMDpdfStrange.Set_C_PDFQuarkHadron();
        
        TMDffStrange.Set_z(z);
        TMDffStrange.Set_A_FFQuarkHadron();
        TMDffStrange.Set_B_FFQuarkHadron();
        TMDffStrange.Set_A_FFBosonHadron();
        TMDffStrange.Set_C_FFQuarkHadron();
        
        // UBAR
        TMDpdfUbar.Set_xBjorken(x);
        TMDpdfUbar.Set_A_PDFQuarkHadron();
        TMDpdfUbar.Set_B_PDFQuarkHadron();
        TMDpdfUbar.Set_A_PDFBosonHadron();
        TMDpdfUbar.Set_C_PDFQuarkHadron();
        
        TMDffUbar.Set_z(z);
        TMDffUbar.Set_A_FFQuarkHadron();
        TMDffUbar.Set_B_FFQuarkHadron();
        TMDffUbar.Set_A_FFBosonHadron();
        TMDffUbar.Set_C_FFQuarkHadron();
        
        // DBAR
        TMDpdfDbar.Set_xBjorken(x);
        TMDpdfDbar.Set_A_PDFQuarkHadron();
        TMDpdfDbar.Set_B_PDFQuarkHadron();
        TMDpdfDbar.Set_A_PDFBosonHadron();
        TMDpdfDbar.Set_C_PDFQuarkHadron();
        
        TMDffDbar.Set_z(z);
        TMDffDbar.Set_A_FFQuarkHadron();
        TMDffDbar.Set_B_FFQuarkHadron();
        TMDffDbar.Set_A_FFBosonHadron();
        TMDffDbar.Set_C_FFQuarkHadron();
        
        // SBAR
        TMDpdfSbar.Set_xBjorken(x);
        TMDpdfSbar.Set_A_PDFQuarkHadron();
        TMDpdfSbar.Set_B_PDFQuarkHadron();
        TMDpdfSbar.Set_A_PDFBosonHadron();
        TMDpdfSbar.Set_C_PDFQuarkHadron();
        
        TMDffSbar.Set_z(z);
        TMDffSbar.Set_A_FFQuarkHadron();
        TMDffSbar.Set_B_FFQuarkHadron();
        TMDffSbar.Set_A_FFBosonHadron();
        TMDffSbar.Set_C_FFQuarkHadron();
        
        // UP
        file_CutoffCollinear_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFUp.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFUp.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUp.Delta_MSbar_iH_1(x, mu, colPDFUp.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUp.Delta_MSbar_Hi_1(z,mu,colFFUp.Get_AlphaStrong(mu))
        << std::endl;
        
        // DOWN
        file_CutoffCollinear_DOWN << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFDown.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFDown.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDown.Delta_MSbar_iH_1(x, mu, colPDFDown.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDown.Delta_MSbar_Hi_1(z,mu,colFFDown.Get_AlphaStrong(mu))
        << std::endl;
        
        // STRANGE
        file_CutoffCollinear_STRANGE << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFStrange.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFStrange.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfStrange.Delta_MSbar_iH_1(x, mu, colPDFStrange.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffStrange.Delta_MSbar_Hi_1(z,mu,colFFStrange.Get_AlphaStrong(mu))
        << std::endl;
        
        // UBAR
        file_CutoffCollinear_UBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFUbar.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFUbar.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfUbar.Delta_MSbar_iH_1(x, mu, colPDFUbar.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffUbar.Delta_MSbar_Hi_1(z,mu,colFFUbar.Get_AlphaStrong(mu))
        << std::endl;
        
        // DBAR
        file_CutoffCollinear_DBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFDbar.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFDbar.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfDbar.Delta_MSbar_iH_1(x, mu, colPDFDbar.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffDbar.Delta_MSbar_Hi_1(z,mu,colFFDbar.Get_AlphaStrong(mu))
        << std::endl;
        
        // SBAR
        file_CutoffCollinear_SBAR << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Evaluate_CutoffIntegral(x,kc)
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Evaluate_CutoffIntegral(z,kc)
        << std::right << std::setfill(' ') << std::setw(20) << colPDFSbar.Evaluate(x)
        << std::right << std::setfill(' ') << std::setw(20) << colFFSbar.Evaluate(z)
        << std::right << std::setfill(' ') << std::setw(20) << TMDpdfSbar.Delta_MSbar_iH_1(x, mu, colPDFSbar.Get_AlphaStrong(mu))
        << std::right << std::setfill(' ') << std::setw(20) << TMDffSbar.Delta_MSbar_Hi_1(z,mu,colFFSbar.Get_AlphaStrong(mu))
        << std::endl;
        
        
        
        // Display progress
        std::cout << "#"+ std::string((j)*100.0/x_steps, '*');// update loading bar
        std::cout << ">" ;
        std::cout << " Complete " << (double(j))*100.0/double(x_steps) << " % " << std::endl;// update percentage
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(now - start);
        std::cout << elapsed_time.count() << " seconds" << std::endl;
    }
    
    // Duration of the program
    auto stop = std::chrono::high_resolution_clock::now();// start counting
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Program END" << std::endl;
    std::cout << "Executed in " << duration.count() << " seconds !!!" << std::endl;
    
    return 0;
}

// main ends




std::string CutoffCollinearFilename(std::string f,double x, double z, double Q, double mu, double m_F, double m_D)
{
    return "QCD_CutoffCollinear" + f + "_Q_" + std::to_string(int(Q)) + "_mu_" + std::to_string(int(mu)) + "_mF_0_" + std::to_string(int(m_F*10)) + "_mDoverz_0_" + std::to_string(int(m_D*10)) +".txt";
}


// Print HEADER
void PrintCutoffCollinearHeader(std::ofstream *filename)
{
    *filename << "#"
    << std::right << std::setfill(' ') << std::setw(20) << "qT"
    << std::right << std::setfill(' ') << std::setw(20) << "CutoffCollinearPDF"
    << std::right << std::setfill(' ') << std::setw(20) << "CutoffCollinearFF"
    << std::right << std::setfill(' ') << std::setw(20) << "MsbarCollinearPDF"
    << std::right << std::setfill(' ') << std::setw(20) << "MsbarCollinearFF"
    << std::right << std::setfill(' ') << std::setw(20) << "Delta_CutoffCollinearPDF"
    << std::right << std::setfill(' ') << std::setw(20) << "Delta_CutoffCollinearFF"
    << std::endl;
}
