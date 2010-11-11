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

    int id,mt,nr;
    bool northern;

    double *xn;

    Domain(int, bool);
    Domain(int, bool, int, int);

    int  idx(int, int, int, int);
    bool grdgen();


};

#endif
