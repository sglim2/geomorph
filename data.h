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

class Data {
public:
    
    enum { UNDEF, MVIS, TERRA, MITP, FILT };
    enum { NEAREST, NEAREST2, LINEAR, CUBIC };

    char *    infile;
    char *    indir;  // GUI only
    char *    outfile;
    char *    outdir;  // GUI only
    int       intype;
    int       outtype;
    int       interp;  // interpolation routine to be used
    int       filtinstart,filtinend,filtinnumfiles;
    int       filtoutstart,filtoutend,filtoutnumfiles;
    int       nlat,nlng,ndpth,nvalpershell;
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
    char*  interpConverter();
    
    bool   getStats();
    
    bool   findBoundary();
    bool   findLayers();
};

#endif
