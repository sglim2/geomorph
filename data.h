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

// we need a forward declaration of the Grid class
// due to a cycle of dependency.
class Grid;

#include "grid.h"

class Data {
public:
    
    enum { MVIS, TERRA_CC, TERRA_CV , MITP, FILT, UNDEF };
    enum { NEAREST, NEAREST2, LINEAR, CUBIC };

    char *    infile;
    bool      infileSet; // GUI only
    char *    indir;  // GUI only
    bool      indirSet;  // GUI only
    char *    outfile;
    bool      outfileSet; // GUI only
    char *    outdir;  // GUI only
    bool      outdirSet;  // GUI only
    int       intype;
    bool      intypeSet;  // GUI only
    int       outtype;
    bool      outtypeSet; // GUI only
    int       interp;  // interpolation routine to be used
    bool      interpSet; // GUI only
    int       filtinstart,filtinend,filtinnumfiles;
    bool      filtinstartSet,filtinendSet,filtinnumfilesSet;
    int       filtoutstart,filtoutend,filtoutnumfiles;
    int       nlat,nlng,ndpth,nvalpershell;
    double    minR,maxR;
    double    * lyrs;
    long int  nval;
    int       nlayr;
    double    *x,*y,*z,*V;
    double    Vmean,Vmax,Vmin;
    double    a,cmb;
    Grid      *mvis;
    
    Data ();
    Data (char *, int);
    
    bool   Read ();
  
    // MITP
    bool   mitpRead ();
    double mitpDepth2Radius(double);

    // FILT
    bool   filtRead ();
    double filtDepth2Radius(double);

    // MVIS
    bool   mvisRead();

    // TERRA
    bool   terraRead();

    char*  intypeConverter();
    char*  outtypeConverter();
    char*  interpConverter();
    
    bool   getStats();
    bool   getStatsData();
    bool   getStatsGrid();
    
    bool   findBoundary();
    bool   findLayers();
};

#endif
