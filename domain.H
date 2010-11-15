/*
 *
 * TERRA grid Definition - the Icosahedron Domain
 *
 * 
 */



#ifndef _GEOMORPH_ICOSAHEDRON_H
#define _GEOMORPH_ICOSAHEDRON_H

#include <math.h>
#include <stdio.h>

class Domain {
public:

    enum Type { NEAREST, LINEAR };
    
    int id,nr,mt;
    bool northern;
    
    double *xn;
    
    Domain();
    Domain(int, bool);
    Domain(int, bool, int, int);
    
    int    defineDomain       (int, int, int);
    int    idx                (int, int, int, int);
    bool   grdgen             ();
    int    getValue           (int, int );
    bool   getNearestDataValue(int);
    bool   importData();

};

#endif
