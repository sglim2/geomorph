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
    
    enum Type { UNDEF, MVIS, TERRA, MITP, FILT };
    
    char *    infile;
    char *    outfile;
    int       intype;
    int       outtype;
    int       filtinstart,filtinend,filtinnumfiles;
    int       filtoutstart,filtoutend,filtoutnumfiles;
    int       nlat,nlng,ndpth;
    double    minR,maxR;
    double    * lyrs;
    long int  nval;
    int       nlayr;
    double    *x,*y,*z,*V;
    double    Vmean,Vmax,Vmin;
    double    a,cmb;
    
    Data ();
    Data (char *, int);
    
    bool   Read ();
  
    // MITP
    bool   mitpRead ();
    double mitpDepth2Radius(double);

    // FILT
    bool   filtRead ();
    double filtDepth2Radius(double);
    
    char*  intypeConverter();
    char*  outtypeConverter();
    
    bool   getStats();
    
    bool   findBoundary();
    bool   findLayers();
};

#endif
