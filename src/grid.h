/**
 \file 

 \brief
  Grid class definition 
 */

#ifndef _GEOMORPH_GRID_H
#define _GEOMORPH_GRID_H

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "domain.h"

class Grid {
public:

    int       suffix;    /**< The value used for the filename suffices. Converted to a 
                              2-digit number with leading-zeroes */
    bool      suffixSet; /**< GUI only. set to true once the suffix value has been set via the GUI */
    int       mt;        /**< The mt value */
    int       nt;        /**< The nt value */
    int       nd;        /**< The nd value */
    bool      mtSet;     /**< GUI only. set to true once the mt value has been set via the GUI */
    bool      ntSet;     /**< GUI only. set to true once the nt value has been set via the GUI */
    bool      ndSet;     /**< GUI only. set to true once the nd value has been set via the GUI */
    int       nproc;     /**< In MVIS/TERRA_CC/TERRA_CV format a nproc value is defined 
                              which is the number of processors used during the calculation.
                              This is calculated from the relationship nproc=pow((mt/nt),2)*10/nd */
    int       nr;        /**< Number of radial divisions, giving the number of radial layers as nr+1 */
    int       idmax;     /**< The number of domains in the grid. Always 10.*/
    int       npts;      /**< The total number of grid points, across all domains*/
    Domain    *domains;  /**< Each grid is split in to 10 domains */
    
    Grid();
    Grid(int, int, int);
    
    bool   genGrid     (double); 
    int    suggestGrid (int);
    int    idx         (int, int, int, int, int);
    bool   importData  (Data*);
    bool   importMVIS  (char*, double);
    bool   importTERRA (char*, int);
    bool   exportGrid  (Data*);
    bool   exportMVIS  (Data*);
    bool   exportTERRA (Data*, int);

};

extern"C" {
    void grdgen_(double *, int *);
}

#endif
