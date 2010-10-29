/*
 *
 * TERRA grid Definition
 *
 * 
 */


#ifndef _GEOMORPH_GRID_H
#define _GEOMORPH_GRID_H


class Grid {
public:

  int       mt,nt,nd;
  double    *x,*y,*z,*V;
   
    

  Grid();
  Grid(int, int, int);

  bool genGrid();
};

#endif
