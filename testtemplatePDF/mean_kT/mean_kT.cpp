//
//  main.cpp
//  mean_kT
//
//  Created by Tommaso Rainaldi on 11/4/22.
//

#include <iostream>
#include <fstream>
#include <string>

#include "TMDPDF.h"
#include "TMDFF.h"

int main(int argc, const char * argv[]) {
    
    
    double Q  = 2.0;// Hard Scale
    double mu = 2.0;// Renormalization Scale
    
    // Set initial parameters
    double x_in    = 0.35;// starting qT
    double x       = x_in;
    double z       = x_in;
    double x_fin   = 0.99;// final qT
    int    x_steps = 200;// number of points in qT space
    double x_Delta = (x_fin - x_in)/x_steps;// step size in qT space
                                              //m_i // m_g // M_i
    InputTMDPDF_massParameters massParamsPDF = {0.2,   0.2,   0.2};
    InputTMDFF_massParameters  massParamsFF  = {0.2,   0.2,   0.2};
    
    bool NonPerturbative = true;// Flag for Input TMD PDF or FF to compute
    
    bool isLHAPDF = true;// if Collinear distributions are from LHAPDF or user defined
    
    
    // Constant expressions known at compiler time
    constexpr int Target  = UNDEFINED_HADRON;// unknown
    constexpr int Parton1 = UP;
    constexpr int Parton2 = DOWN;
    constexpr int Gluon   = GLUON;
    
    
    
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
//    CollinearFF<Target,Parton1> colFFUp(mu,Q, isLHAPDF);// Collinear FF for Parton1
//    CollinearFF<Target,Parton2> colFFDown(mu,Q, isLHAPDF);// Collinear FF for Parton2
//    CollinearFF<Target,Gluon> colFFGluon(mu,Q, isLHAPDF);// Collinear FF for Gluon
//    // Get input non-perturbative TMD FF
//    TMDFF<Target,Parton1,Gluon> TMDffUp(colFFUp, colFFGluon, z, massParamsFF, NonPerturbative);
//    TMDFF<Target,Parton2,Gluon> TMDffDown(colFFDown, colFFGluon, z, massParamsFF, NonPerturbative);
//    // Get input perturbative TMD FF
//    TMDFF<Target,Parton1,Gluon> TMDffUp_pert(colFFUp, colFFGluon, z, massParamsFF, !NonPerturbative);
//    TMDFF<Target,Parton2,Gluon> TMDffDown_pert(colFFDown, colFFGluon, z, massParamsFF, !NonPerturbative);
    
    
    // Set up the txt files to write on
    // Up
    std::ofstream file_meankT_UP;
    file_meankT_UP.open ("meankT_UP.txt");
    file_meankT_UP.precision(8);
    file_meankT_UP.setf(std::ios::fixed | std::ios::showpoint);
    
    double mean_kT_PDF;
    double mean_kT_PDF_pert;
    double mean_kT_FF;
    double mean_kT_FF_pert;
    
    double kT2cutoff = Q*Q/16.0;
    
    
    for(int j = 0; j <= x_steps; j++)
    {
        x = x_in + j*x_Delta;// update qT
        z = x;
        // PDFs
        TMDpdfUp.Set_xBjorken(x);
        TMDpdfUp.Set_A_PDFQuarkHadron();
        TMDpdfUp.Set_B_PDFQuarkHadron();
        TMDpdfUp.Set_A_PDFBosonHadron();
        TMDpdfUp.Set_C_PDFQuarkHadron();
        
        TMDpdfUp_pert.Set_xBjorken(x);
        TMDpdfUp_pert.Set_A_PDFQuarkHadron();
        TMDpdfUp_pert.Set_B_PDFQuarkHadron();
        TMDpdfUp_pert.Set_A_PDFBosonHadron();
        TMDpdfUp_pert.Set_C_PDFQuarkHadron();
        
        
        // FFs
//        TMDffUp.Set_z(z);
//        TMDffUp.Set_A_FFQuarkHadron();
//        TMDffUp.Set_B_FFQuarkHadron();
//        TMDffUp.Set_A_FFBosonHadron();
//        TMDffUp.Set_C_FFQuarkHadron();
//
//        TMDffUp_pert.Set_z(z);
//        TMDffUp_pert.Set_A_FFQuarkHadron();
//        TMDffUp_pert.Set_B_FFQuarkHadron();
//        TMDffUp_pert.Set_A_FFBosonHadron();
//        TMDffUp_pert.Set_C_FFQuarkHadron();
        
        
        
        
        mean_kT_PDF      = TMDpdfUp.Get_mean_squared_transverse_momentum(kT2cutoff);
        std::cout << TMDpdfUp.Get_A_PDFQuarkHadron() << TMDpdfUp.Af_iH(TMDpdfUp.Get_xBjorken(), mu) << std::endl;
        mean_kT_PDF_pert = TMDpdfUp_pert.Get_mean_squared_transverse_momentum(kT2cutoff);
//        mean_kT_FF       = TMDffUp.Get_mean_squared_transverse_momentum(kT2cutoff);
//        mean_kT_FF_pert  = TMDffUp_pert.Get_mean_squared_transverse_momentum(kT2cutoff);
        
        // Write results on external txt file
        // Up
        file_meankT_UP << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << x
        << std::right << std::setfill(' ') << std::setw(20) << mean_kT_PDF
        << std::right << std::setfill(' ') << std::setw(20) << mean_kT_PDF_pert
        << std::right << std::setfill(' ') << std::setw(20) << mean_kT_FF
        << std::right << std::setfill(' ') << std::setw(20) << mean_kT_FF_pert
        << std::endl;
        
    }
    
    
    return 0;
}
