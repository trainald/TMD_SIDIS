//
//  main.cpp
//  Q0bar_Evolution
//
//  Created by Tommaso Rainaldi on 2/3/23.
//


#include <iostream>
#include <fstream>
#include <string>

#include "Q0barEvolution.h"


int main(int argc, const char * argv[]) {
    
    
//    double Q  = 4.0;// Hard Scale
//    double mu = 4.0;// Renormalization Scale

    double bT;
    double bT_in  = 1E-3;// initial b_T
    double bT_fin = 1.5;// final b_T
    int    bT_steps = 200;// number of points in qT space
    double bT_Delta = (bT_fin - bT_in)/bT_steps;// step size in bT space
    
    
    // Set up the txt files to write on
    std::ofstream file_Q0barEvo;
    file_Q0barEvo.open ("Q0barEvo.txt");
    file_Q0barEvo.precision(8);
    file_Q0barEvo.setf(std::ios::fixed | std::ios::showpoint);
    

    double EQ0bar2Q0;
    
    std::cout << "bT" << std::setw(20)<< "E_Q0bar2Q0(bT)" << std::setw(20) << "Q0bar" << std::endl;
    for(int j = 0; j <= bT_steps; j++)
    {
        bT = bT_in + j*bT_Delta;// update bT

        std::cout << bT << std::setw(20)<< E_Q0bar2Q0(bT) << std::setw(20)<< Q0bar(bT) <<std::endl;
    
        
        // Write results on external txt file
        file_Q0barEvo << std::scientific
        << std::right << std::setfill(' ') << std::setw(20) << bT
        << std::right << std::setfill(' ') << std::setw(20) << E_Q0bar2Q0(bT)
        << std::right << std::setfill(' ') << std::setw(20) << Q0bar(bT)
        << std::endl;
        
    }
    
    
    return 0;
}




