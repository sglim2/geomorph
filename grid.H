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

class Grid {
public:

  int       mt,nt,nd;
  int       nproc;
  int       nr;
  int       npts;
  double    rmax,rmin;
  double    *x,*y,*z,*V; 


  Grid();
  Grid(int, int, int);

  bool genGrid();
  int findGrid(int);
};

#endif
