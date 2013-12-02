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
    
    enum { MVIS, TERRA_CC, TERRA_CV , MITP, FILT, GYPSUMP, GYPSUMS, UNDEF };
    enum { NEAREST, NEAREST2, LINEAR, CUBIC };

    char *    infile;      /**< Base-name of input file */
    bool      infileSet;   /**< GUI only. Set to true when infile has been set via the GUI */
    char *    indir;       /**< GUI only. input base-name is split between path and file within the GUI */
    bool      indirSet;    /**< GUI only. Set to true when indir has been set via the GUI */
    char *    outfile;     /**< Base-name of output file */ 
    bool      outfileSet;  /**< GUI only. Set to true when outfile has been set via the GUI */
    char *    outdir;      /**< GUI only. output base-name is split between path and file within the GUI */
    bool      outdirSet;   /**< GUI only. Set to true when outdir has been set via the GUI */
    int       intype;      /**< input type - MVIS,TERRA_CV,FILT,etc... */
    bool      intypeSet;   /**< GUI only. Set to true when intype has been set via the GUI */
    int       outtype;     /**< output type - MVIS,TERRA_CV,etc... */
    bool      outtypeSet;  /**< GUI only. Set to true when outtype has been set via the GUI */
    int       interp;      /**< interpolation routine to be used */
    bool      interpSet;   /**< GUI only. Set to true when interp has been set via the GUI */
    int       filtinstart;
    int       filtinend;
    int       filtinnumfiles; 
    bool      filtinstartSet;    /**< GUI only. Set to true when filtinstart has been set via the GUI */
    bool      filtinendSet;      /**< GUI only. Set to true when filtinend has been set via the GUI */
    bool      filtinnumfilesSet; /**< GUI only. Set to true when filtinnumfiles has been set via the GUI */
    int       filtoutstart;
    int       filtoutend;
    int       filtoutnumfiles;
    int       gypsuminnumfiles;
    char *    gypsumlatloninfile;
    char *    gypsumdepthinfile;
    bool      gypsumlatloninfileSet,gypsumdepthinfileSet;
    int       nlat,nlng,ndpth,nvalpershell;
    double    minR,maxR;
    double    * lyrs;
    long int  nval;
    int       nlayr;
    double    *x,*y,*z,*V;
    double    Vmean,Vmax,Vmin;
    double    a,cmb;
    bool      cmbinSet;
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

    // GYPSUM
    bool   gypsumRead ();
    double gypsumDepth2Radius(double);

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
