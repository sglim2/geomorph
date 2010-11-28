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

#include "data.H"

class Domain {
public:

    enum Type { NEAREST, LINEAR };
    
    int id,nr,mt;
    bool northern;
    
    double *xn,*yn,*zn;
    double *V;
    
    Domain();
    Domain(int, bool);
    Domain(int, bool, int, int);
    
    bool   defineDomain        (int, int, int);
    int    idx                 (int, int, int);
    bool   grdgen              (double);
    int    getValue            (Data*, int, int );
    double getNearestDataValue (Data*, int);
    double getNearestDataValue2(Data*, int);
    bool   importData          (Data*);
    bool   exportMVIS          (FILE *, int, int, int);
    int    sqrti               (int);
    bool   rotate3d            (double*, double*, double*, double, double, double, double);

};

#endif
