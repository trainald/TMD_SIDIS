//
//  Flavour.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/21/22.
//

#ifndef Flavour_h
#define Flavour_h

enum Flavour {
    // Quarks
    DOWN    = 1, DBAR = -1,
    UP      = 2, UBAR = -2,
    STRANGE = 3, SBAR = -3,
    CHARM   = 4, CBAR = -4,
    BOTTOM  = 5, BBAR = -5,
    TOP     = 6, TBAR = -6,
    // Bosons
    GLUON   = 0, PHOTON = 7,
    // Hadrons (not standard)
    PROTON  = (2*UP + DOWN)*100, NEUTRON = (2*DOWN + UP)*100,
    //undefined (not standard)
    UNDEFINED_HADRON = 100, UNDEFINED_PARTON = 50
};



#endif /* Flavour_h */



