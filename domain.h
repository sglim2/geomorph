/*
 *
 * TERRA grid Definition - the Icosahedron Domain
 *
 * 
 */



#ifndef _GEOMORPH_ICOSAHEDRON_H
#define _GEOMORPH_ICOSAHEDRON_H

#include <math.h>
#include <stdio.h>

//#include "main.h"
#include "data.h"

class Domain {
public:
    
    int id,nr,mt;
    bool northern;
    
    double *xn,*yn,*zn;
    double *V;

    double Vmax,Vmin;

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

    bool   importData(Data*);
    bool   importData_gpu(Data*);

    bool   importMVIS(FILE *, int, int, int, double);
    bool   exportMVIS(FILE *, int, int, int);
    bool   exportTERRA(FILE *, int, int, int, int, int, long int &);

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
