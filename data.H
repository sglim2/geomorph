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
  Grid      grid;



  Data ();
  Data (char *, int);

  bool Read ();

};

#endif
