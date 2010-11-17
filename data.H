/*
 *
 * Data Manipulation
 *
 * 
 */



#ifndef _GEOMORPH_DATA_H
#define _GEOMORPH_DATA_H

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
    int     nlat,nlng,ndpth;
  double    minR,maxR;
  double    * lyrs;
  long int  nval;
  int       nlayr;
  double    *x,*y,*z,*V;
  double    Vmean,Vmax,Vmin;

  Data ();
  Data (char *, int);

  bool   Read ();
  bool   mitpRead ();
  double mitpDepth2Radius(double);

  char* intypeConverter();
  char* outtypeConverter();

  bool  getStats();
    
  bool  findBoundary();
  bool  findLayers();
};

#endif
