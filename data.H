/*
 *
 * Data Manipulation
 *
 * 
 */



#ifndef _GEOMORPH_DATA_H
#define _GEOMORPH_DATA_H

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "grid.H"

class Data {
public:

  enum Type { UNDEF, MVIS, TERRA, MITP };

  char *    infile;
  char *    outfile;
  int       intype;
  int       outtype;
  long int  nval;
  double    *x,*y,*z,*V;
  double    Vmean,Vmax,Vmin;
  Grid      grid;



  Data ();
  Data (char *, int);

  bool  Read ();
  bool  mitpRead ();
  char* intypeConverter();
  char* outtypeConverter();
  bool  getStats();
    
};

#endif
