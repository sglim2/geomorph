/*
 *
 * TERRA grid Definition
 *
 * 
 */


#ifndef _GEOMORPH_GRID_H
#define _GEOMORPH_GRID_H

#include <math.h>
#include <stdio.h>

#include "domain.h"

class Grid {
public:

    int       mt,nt,nd;
    int       nproc;
    int       nr;
    int       idmax;
    int       npts;
    double    rmax,rmin;
    Domain    *domains;
    
    Grid();
    Grid(int, int, int);
    
    bool   genGrid    (double);
    int    suggestGrid(int);
    int    xnProc     (int);
    int    xnProc     (int, int, int, int, int);
    int    idx        (int, int, int, int, int);
    bool   importData (Data*);
    bool   exportGrid (Data*);
    bool   exportMVIS (Data*);

};

extern"C" {
    void grdgen_(double *, int *);
}

#endif
