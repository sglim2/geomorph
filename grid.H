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

#include "domain.H"

class Grid {
public:

    int       mt,nt,nd;
    int       nproc;
    int       nr;
    int       idmax;
    int       npts;
    double    rmax,rmin;
    double    *x,*y,*z,*V;
    Domain    *domains;
    
    Grid();
    Grid(int, int, int);
    
    bool genGrid();
    int  findGrid(int);
    int  xnProc(int idx);
    int  xnProc(int r, int id, int i1, int j1, int xyz);
    int  idx(int r, int id, int i1, int j1, int xyz);
    
};

extern"C" {
    void grdgen_(double *, int *);
}

#endif
