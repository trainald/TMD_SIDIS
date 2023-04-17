//
//  main.cpp
//  test_LHAPDF
//
//  Created by Tommaso Rainaldi on 10/31/22.
//

#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "YtermFunctions.h"

using namespace std;
int main(int argc, const char * argv[]) {
    // See what sets are available
    for (const string& p : LHAPDF::paths())
        cout << p << endl;
     
    cout << "@" << LHAPDF::findFile("lhapdf.conf") << "@" << endl;
     
    cout << "-------------------------------" << endl;
    cout << "List of available PDFs:" << endl;
    cout << "-------------------------------" << endl;
      for (const string& s : LHAPDF::availablePDFSets())
        cout << "# " << s << endl;
    cout << "-------------------------------" << endl;
    
    LHAPDF::PDF* colpdf;
    colpdf = LHAPDF::mkPDF("cteq66", 0);
    //const LHAPDF::PDF* pdf = LHAPDF::mkPDF("cteq66", 0);
    double Q = 2;
    double alpha = colpdf->alphasQ2(Q*Q);
    std::cout << alpha << std::endl;
    
    
    
//    double x = LHAPDF::sin(1.0);
//    std::cout<< x;
    
    
    
    return 0;
}
