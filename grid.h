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
#include <time.h>

#include "domain.h"

class Grid {
public:

    int       suffix;
    int       mt,nt,nd;
    bool      mtSet,ntSet,ndSet; // GUI only
    int       nproc;
    int       nr;
    int       idmax;
    int       npts;
    double    rmax,rmin;
    Domain    *domains;
    
    Grid();
    Grid(int, int, int);
    
    bool   genGrid     (double);
    int    suggestGrid (int);
    int    idx         (int, int, int, int, int);
    bool   importData  (Data*);
    bool   importMVIS  (char*, double);
    bool   importTERRA (char*, int);
    bool   exportGrid  (Data*);
    bool   exportMVIS  (Data*);
    bool   exportTERRA (Data*, int);

};

extern"C" {
    void grdgen_(double *, int *);
}

#endif
