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
    int       filtinstart;       /**< The start depth of the FILT input data */
    int       filtinend;         /**< The end depth of the FILT input data */
    int       filtinnumfiles;    /**< The number of datafiles for the FILT input data */
    bool      filtinstartSet;    /**< GUI only. Set to true when filtinstart has been set via the GUI */
    bool      filtinendSet;      /**< GUI only. Set to true when filtinend has been set via the GUI */
    bool      filtinnumfilesSet; /**< GUI only. Set to true when filtinnumfiles has been set via the GUI */
    int       filtoutstart;      /**< The start depth of the FILT output - not implemented */
    int       filtoutend;        /**< The end depth of the FILT output - not implemented */
    int       filtoutnumfiles;   /**< The number of outfiles used for FILT output - not implemeted */
    int       gypsuminnumfiles;   /**< The number of files containing the gypsum data-set */
    char *    gypsumlatloninfile; /**< The path of the gyspum lat-lon file */
    char *    gypsumdepthinfile;  /**< The path of the gypsum depth file */
    bool      gypsumlatloninfileSet; /**< GUI only. Set to true when gypsumlatloninfile has been set via the GUI - not implemented*/
    bool      gypsumdepthinfileSet;  /**< GUI only. Set to true when gypsumdepthinfile has been set via the GUI - not implemented*/
    int       nlat;             /**< The number of latitude points in the input data */
    int       nlng;             /**< The number of longitude points in the input data */
    int       ndpth;            /**< The number of depth points in the input data */
    int       nvalpershell;     /**< The number of values per shell in the input data - usually nlat x nlng */
    double    minR;             /**< The minimum radius found in the input data*/
    double    maxR;             /**< The maximum radius found in the input data */
    double    * lyrs;           /**< The pointer to the input data lyrs */
    long int  nval;             /**< The number of grid-points in the input data */
    int       nlayr;            /**<  The number of radial layers in the input data */
    double    *x;               /**< The x-coordinates of the input data */
    double    *y;               /**< The y-coordinates of the input data */
    double    *z;               /**< The z-corrdinates of the input data */
    double    *V;               /**< The values of the input data */
    double    Vmean;            /**< The statistical mean on the input values */
    double    Vmax;             /**< The maximum value of the input data */
    double    Vmin;             /**< The minimum value of the input data */
    double    a;                /**< The surface Radius */
    double    cmb;              /**< The Core-Mantle boundary in Earth-Radius units */
    bool      cmbinSet;         /**< GUI only. Set to true when cmb has been set via the GUI */
    Grid      *mvis;            /**< The geomorph grid definition*/
    
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
