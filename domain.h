/**
 \file 

 \brief
  Domain class definition 
 */

#ifndef _GEOMORPH_ICOSAHEDRON_H
#define _GEOMORPH_ICOSAHEDRON_H

#include <math.h>
#include <stdio.h>

//#include "main.h"
#include "data.h"

class Domain {
public:
    
    int  id;       /**< The domain id 0..9 */ 
    int  nr;       /**< Number of radial divisions, number of layers = nr+1 */
    int  mt;       /**< The mt value of the TERRA grid */
    bool northern; /**< true for northern hemisphere, false for southern */
    
    double *xn;  /**< The x-position of all grid points in the domain */
    double *yn;  /**< The y-position of all grid points in the domain */
    double *zn;  /**< The z-position of all grid points in the domain */
    double *V;   /**< The value of the grid points in the domain */
    double *P;   /**< Not yet implemented */
    double *vel; /**< Not yet implemented */

    double Vmax; /**< The max value of V in the domain */
    double Vmin; /**< The minimum value of V in the domain */

    Domain();
    Domain(int, bool);
    Domain(int, bool, int, int);
    
    bool   defineDomain(int, int, int);

    int    idx(int, int, int);
    bool   midpt(double*, double*, double*, 
		 double, double, double, 
		 double, double, double);
    bool   grdgen(double);
    bool   grdgen2(double);

    int    getValue(Data*, int, int );
    double getNearestDataValue(Data*, int);
    double getNearestDataValue2_mitp(Data*, int);
    double getNearestDataValue2_filt(Data*, int);

    int    getLinearDataValue(Data*);

    bool   importData(Data*);
    bool   importData_gpu(Data*);

    bool   importMVIS(FILE *, int, int, double);
    bool   importTERRA(FILE *, int, int, int, int);
    bool   exportMVIS(FILE *, int, int);
    bool   exportTERRA(FILE *, int, int, int, int, long int &);

    int    sqrti(int);
    bool   rotate3d(double*, double*, double*, 
		    double, double, double, double);
    bool   crossProduct(double, double, double, 
		        double*, double*, double*);
    double dotProduct(double, double, double, 
		      double, double, double);
    double getAngle3d(double, double, double, 
		      double, double, double);
    bool   normalise(double *, double *, double *);

};

extern"C" {
  //    double getNearestDataValue2_gpu(int, double, double, int, int, 
  //				  double *, double *, double *,
  //				  double *, double *, double *, double *,
  //				  int);
  
    bool importData_c_gpu(int, int, 
			  int, double, double, int, int, 
			  double *, double *, double *, double *,
			  double *, double *, double *, double *);
}

#endif
