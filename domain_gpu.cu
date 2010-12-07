
#include <math.h>

////////////////////////////////////////
// idx_gpu()
//
int idx_gpu(int r, int i2, int i1, int nr, int mt)
{
  int   idx;
  int   rbase;
  int   i2base;
  int   i1base;
  idx=0;

  rbase = r  * (mt+1)*(mt+1);
  i2base = i2 * (mt+1);
  i1base = i1;
  
  idx = rbase + i2base + i1base;

   return idx;
}

////////////////////////////////////////
// getNearestDataValue2_gpu()
//
double getNearestDataValue2_gpu(int _ndpth, double _minR, double _maxR, int _nlat, int _nlng, 
				double *_xn, double *_yn, double *_zn,
				double *_x,  double *_y , double *_z , double *_V,
				int index)
{
    double rad,dR,dataR;
    double d2;
    double xd,yd,zd;
    double tmpd2;
    double dataV;
    double localVeryLarge;
    int nr, di, ir; 
    
    double gx,gy,gz;
    gx = _xn[index];
    gy = _yn[index];
    gz = _zn[index];

    localVeryLarge=1E+99;

    rad=sqrt(gx*gx + gy*gy + gz*gz);

    nr = 0;
    dR = localVeryLarge;
    dataR = 0.;
    for ( ir=0 ; ir<_ndpth ; ir++ ){
      dataR = _minR + ir*(_maxR - _minR)/_ndpth;
      if (fabs(dataR - rad) < dR) {
	dR = fabs(dataR - rad) ;
	nr = _ndpth - ir;
      }
    }

    dataV = 0.;

    d2 = 1.E+99;
    
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


////////////////////////////////////////
// importData_c_gpu()
//
extern "C" bool importData_c_gpu(int nr, int mt, 
		      int _ndpth, double _minR, double _maxR, int _nlat, int _nlng, 
		      double *_xn, double *_yn, double *_zn, double *_Vn,
		      double *_x,  double *_y , double *_z , double *_V)
{
  int ri,i2,i1,index;
  
  for ( ri=0 ; ri < nr ; ri++ ){
	for ( i2 = 0 ; i2 < mt+1 ; i2++ ){
	    for ( i1 = 0 ; i1 < mt+1 ; i1++ ) {
	      index=idx_gpu(ri,i2,i1,nr,mt);
	      _Vn[index] = getNearestDataValue2_gpu(_ndpth, _minR, _maxR, 
						    _nlat, _nlng,
						    _xn, _yn, _zn,
						    _x,  _y,  _z, _V,
						    index);
	    }
	}
    }

  return 0;
}



