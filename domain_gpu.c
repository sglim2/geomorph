/*
 *  This routine is a copy of the Domain::getNearestDataValue2() algorithm,
 *  but written in C (rather than C++). This allows the option of compiling
 *  with the PGI Accelerator, which is not available in C++ (at the time of
 *  writing).
 */

#include <math.h>

double getNearestDataValue2_gpu(int _ndpth, double _minR, double _maxR, int _nlat, int _nlng, 
				double *_xn, double *_yn, double *_zn,
				double *_x,  double *_y , double *_z , double *_V,
				int index)
{
    double localVeryLarge=1E+99;
    int nr=0, di, ir; 
    
    double gx,gy,gz;
    gx = _xn[index];
    gy = _yn[index];
    gz = _zn[index];
    
    // what's our current radius
    double rad=sqrt(gx*gx + gy*gy + gz*gz);

    // find closest radial layer in Data::dptr nr
    double dR=localVeryLarge;
    double dataR=0.;
    for ( ir=0 ; ir<_ndpth ; ir++ ){
      dataR = _minR + ir*(_maxR - _minR)/_ndpth;
      if (fabs(dataR - rad) < dR) {
	dR = fabs(dataR - rad) ;
	nr = _ndpth - ir;  // reverse ordering.
      }
    }

    double dataV;
    dataV=0.;

    double d2 = 1.E+99;
    double xd,yd,zd;
    double tmpd2;
    
    // loop over all Data values within layer nr
    for ( di=nr*_nlat*_nlng; di<nr*_nlat*_nlng + _nlat*_nlng; di++ ){
      xd=_x[di] - _xn[index];
      yd=_y[di] - _yn[index];
      zd=_z[di] - _zn[index];
      tmpd2=xd*xd + yd*yd + zd*zd;
      if ( tmpd2 < d2 ) {
	d2 = tmpd2;
	dataV = _V[di];
      }
    }

    return dataV;
}
