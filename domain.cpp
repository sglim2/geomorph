/*
 *
 * Data Manipulation within Domains
 *
 * 
 */

#include "main.h"
#include "domain.h"

/*
 * utilities
 */


////////////////////////////////////////
// construction/destruction
////////////////////////////////////////

////////////////////////////////////////
// Domain::Domain
// 
//
Domain::Domain()
    : id(), nr(), mt(), northern(), xn(), yn(), zn(), V(), Vmax(), Vmin()
{
    Vmin = veryLarge;
    Vmax = verySmall;
}

////////////////////////////////////////
// Domain::Domain
// 
//
Domain::Domain(int _id, bool _northern)
    : id(), nr(), mt(), northern(), xn(), yn(), zn(), V(), Vmax(), Vmin()
{
    if ( _id >= 10 ) {
	printf("Domain::Domain fail");
    }
    
    id=_id;
    northern=_northern;

    Vmin = veryLarge;
    Vmax = verySmall;
}

////////////////////////////////////////
// Domain::Domain
// 
//
Domain::Domain(int _id, bool _northern, int _mt, int _nr)
  : id(), nr(), mt(), northern(), xn(), yn(), zn(), V(), Vmax(), Vmin()
{
    if (_id >= 10) {
	printf("Domain::Domain fail");
    }
    
    id=_id;
    northern=_northern;
    mt = _mt;
    nr = _nr;

    Vmin = veryLarge;
    Vmax = verySmall;

    xn = new double[nr*(mt+1)*(mt+1)];
    yn = new double[nr*(mt+1)*(mt+1)];
    zn = new double[nr*(mt+1)*(mt+1)];
    V  = new double[nr*(mt+1)*(mt+1)];
    vel= new double[nr*(mt+1)*(mt+1)*3];
    P  = new double[nr*(mt+1)*(mt+1)];
}


////////////////////////////////////////
// operators
////////////////////////////////////////



////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Domain::defineDomain
//
//
bool Domain::defineDomain(int _id, int _nr, int _mt)
{
    
    id = _id;
    nr = _nr;
    mt = _mt;
    
    // not needed??
    //    id < 5 ? northern = true : northern = false;

    // not needed??
    //    id < 5 ? northern = true : northern = false;

    xn = new double[nr*(mt+1)*(mt+1)];
    yn = new double[nr*(mt+1)*(mt+1)];
    zn = new double[nr*(mt+1)*(mt+1)];
    V  = new double[nr*(mt+1)*(mt+1)];

    return 0;
}

////////////////////////////////////////
// Domain::idx
//
// given:
//   xyz....x=0,y=1,z=2
//   i1.....domain indices
//   i2.....domain indices
//   r......radial layer  
// returns the index value for array xn.
// Array xn has structure:
//    xn[nr][mt+1][mt+1][xyz]
//
int Domain::idx(int r, int i2, int i1)
{
  int nerror = 0;

  if (  r  > nr    || r     < 0 ) nerror=1; 
  if ( i2  > mt    || i2    < 0 ) nerror=1;
  if ( i1  > mt    || i1    < 0 ) nerror=1;


  int idx=0;
  
  int   rbase = r  * (mt+1)*(mt+1);
  int  i2base = i2 * (mt+1);
  int  i1base = i1;
  
  idx = rbase + i2base + i1base;

  if (nerror!=0){
      printf("Domain::idx Error: ir = %d; i2 = %d; i1 = %d\n",r,i2,i1);
  }

  return idx;
}

////////////////////////////////////////
// Domain::getNearestDataValue
//
// Given the index of the geomorph grid, searches for the 'V' value of the
// nearest spatial point of the 'Data' data.
double Domain::getNearestDataValue(Data *dptr, int index)
{
  
    // search Data::x,y,z for nearest spatial point
    double dataV;
    dataV=0.;

    double d2 = 1.E+99;
    double xd,yd,zd;

    // loop over all Data values
    for ( int di = 0 ; di < dptr->nval ; di++ ){
      xd=dptr->x[di] - xn[index];
      yd=dptr->y[di] - yn[index];
      zd=dptr->z[di] - zn[index];
      if ( xd*xd + yd*yd + zd*zd < d2 ) {
	d2 = xd*xd + yd*yd + zd*zd;
	dataV = dptr->V[di];
      }
    }

    return dataV;
}

////////////////////////////////////////
// Domain::getNearestDataValue2
//
// Given the index of the geomorph grid, searches for the 'V' value of the
// nearest spatial point in the Data class. A more optimised algorithm than
// getNearestDataValue.
double Domain::getNearestDataValue2_mitp(Data *dptr, int index)
{
    int nr=0; // our layer of interest in Data:dptr

    double gx,gy,gz;
    gx = xn[index];
    gy = yn[index];
    gz = zn[index];
    
    // what's our current radius
    double rad=sqrt(gx*gx + gy*gy + gz*gz);

    // find closest radial layer in Data::dptr nr
    double dR=veryLarge;
    double dataR=0.;
    for ( int ir=0 ; ir<dptr->ndpth ; ir++ ){
      dataR = dptr->minR + ir*(dptr->maxR - dptr->minR)/dptr->ndpth;
      if (fabs(dataR - rad) < dR) {
	dR = fabs(dataR - rad) ;
	nr = dptr->ndpth - ir;  // reverse ordering.
      }
    }

    double dataV;
    dataV=0.;

    double d2 = 1.E+99;
    double xd,yd,zd;
    double tmpd2;
    // loop over all Data values within layer nr
    for (int di=nr*dptr->nlat*dptr->nlng; di<nr*dptr->nlat*dptr->nlng + dptr->nlat*dptr->nlng; di++ ){
      xd=dptr->x[di] - xn[index];
      yd=dptr->y[di] - yn[index];
      zd=dptr->z[di] - zn[index];
      tmpd2=xd*xd + yd*yd + zd*zd;
      if ( tmpd2 < d2 ) {
	d2 = tmpd2;
	dataV = dptr->V[di];
        if ( dataV < Vmin ) Vmin = dataV;
        if ( dataV > Vmax ) Vmax = dataV;
      }
    }

    return dataV;
}

////////////////////////////////////////
// Domain::getNearestDataValue2_filt
//
// Given the index of the geomorph grid, searches for the 'V' value of the
// nearest spatial point in the Data class for FILT data. As with
// getNearestDataValue2-mitp, a more optimised algorithm than
// getNearestDataValue.
double Domain::getNearestDataValue2_filt(Data *dptr, int index)
{
    int nr=0; 

    double gx,gy,gz;
    gx = xn[index];
    gy = yn[index];
    gz = zn[index];
    
    // what's our current radius
    double rad=sqrt(gx*gx + gy*gy + gz*gz);

    // find closest radial layer in Data::dptr nr
    double dR=veryLarge;
    double dataR=0.;
    for ( int ir=0 ; ir<dptr->ndpth ; ir++ ){
      dataR = dptr->minR + ir*(dptr->maxR - dptr->minR)/dptr->ndpth;
      if (fabs(dataR - rad) < dR) {
	dR = fabs(dataR - rad) ;
	nr = dptr->ndpth - ir;  // reverse ordering.
      }
    }

    double dataV;
    dataV=0.;

    double d2 = 1.E+99;
    double xd,yd,zd;
    double tmpd2;
    // loop over all Data values within layer nr
    for (int di=nr*dptr->nvalpershell; di<nr*dptr->nvalpershell + dptr->nvalpershell; di++ ){
      xd=dptr->x[di] - xn[index];
      yd=dptr->y[di] - yn[index];
      zd=dptr->z[di] - zn[index];
      tmpd2=xd*xd + yd*yd + zd*zd;
      if ( tmpd2 < d2 ) {
	d2 = tmpd2;
	dataV = dptr->V[di];
      }
    }

    return dataV;
}

////////////////////////////////////////
// Domain::getLinearDataValue
//
// Converts one grid to another, with up-scaling or down-scaling, if necessary.
//
// mtin = mtout -> 1. simple copy of grid
//
// mtin < mtout -> 1. straight copy of grid points into new grid (relevant
//                    points only) 
//                 2. cycle through each layer (having previously defined
//                    points) on output grid, performing a bilinear
//                    interpoloation of all undefined points on that layer.
//                 3. Cycle through all adjacent layers with defined points,
//                    performing a linearinterpolation between layers for all
//                    layer points.
//
// mtin > mtout -> 1. copy a sub-set of input grid points to the new grid.
//
int Domain::getLinearDataValue(Data *dptr)
{
  int index1=0;
  int index2=0;
  int i22 = 0;
  int i12 = 0;
  int ri2 = 0;
  double * Vptr=V;
  int ScaleFactor=mt/dptr->mvis->domains[id].mt;
  
  if ( dptr->mvis->mt == mt ) {
      // straight copy
      for ( int ri=0 ; ri < nr ; ri++) {
	  for ( int i2 = 0 ; i2 < mt + 1 ; i2++) {
	      for ( int i1 = 0 ; i1 < mt + 1 ; i1++) {
		  index1 = dptr->mvis->domains[id].idx(ri,i2,i1);
		  index2 = index1;
		  Vptr[index2] = dptr->mvis->domains[id].V[index1];
	      }
	  }
      } 
  } // if mtin==mtout
  
  if ( dptr->mvis->mt <  mt ) {
      // Up-scaling

      // 1. Copy existing data...
      for ( int ri1=0 ; ri1 < dptr->mvis->domains[id].nr ; ri1++) {
	  for ( int i21 = 0 ; i21 < dptr->mvis->domains[id].mt + 1 ; i21++) {
	      for ( int i11 = 0 ; i11 < dptr->mvis->domains[id].mt + 1 ; i11++) {
		  index1 = dptr->mvis->domains[id].idx(ri1,i21,i11);
		  ri2 = ri1 * ScaleFactor;
		  i22 = i21 * ScaleFactor;
		  i12 = i11 * ScaleFactor;
		  index2 = idx(ri2,i22,i12);
		  Vptr[index2] = dptr->mvis->domains[id].V[index1];
	      }
	  }
      }

      // 2. Cycle through each defined layer, bilinear-interpolating between
      //    all existing points.
      for ( int ri1=0 ; ri1 < dptr->mvis->domains[id].nr ; ri1++) {
	  ri2 = ri1 * ScaleFactor;

	  for ( int i21 = 0 ; i21 < dptr->mvis->domains[id].mt +1 ; i21++) {
	      i22 = i21*ScaleFactor;

	      for ( int i11 = 0 ; i11 < dptr->mvis->domains[id].mt ; i11++) {

		  int A = i11 * ScaleFactor; // starting i1 point
		  int B = (i11+1) * ScaleFactor; // ending i1 point 
		  for ( i12 = A+1 ; i12 < B ; i12++) {
		      index2 = idx(ri2,i22,i12);
		      V[index2] = V[idx(ri2,i22,A)] + 
          			  (i12-A) * 
	        		  ( V[idx(ri2,i22,B)] - V[idx(ri2,i22,A)]) / ScaleFactor;
		  } 
	      }
	  }

	  for ( int i21 = 0 ; i21 < dptr->mvis->domains[id].mt ; i21++) {
	      int A = i21 * ScaleFactor;
	      int B = (i21+1) * ScaleFactor;
	      
	      for ( int i12 = 0 ; i12 < mt+1 ; i12++) { 
		  
		  for ( i22 = A+1 ; i22 < B ; i22++) {
		      index2 = idx(ri2,i22,i12);
		      V[index2] = V[idx(ri2,A,i12)] + 
          			  (i22-A) * 
	        		  ( V[idx(ri2,B,i12)] - V[idx(ri2,A,i12)]) / ScaleFactor;
		      
		  }
	      }
	  }
      } // radial layers

      // 3. Cycle through all adjacent layers with defined data-values. For
      //    each point withon theses layers, interpolate lineraly for all
      //    sandwiched layers.
      for ( int ri1=0 ; ri1 < dptr->mvis->domains[id].nr - 1 ; ri1++) {
	  int A = ri1 * ScaleFactor;
	  int B = (ri1+1) * ScaleFactor;

	  for ( int i22 = 0 ; i22 < mt+1 ; i22++) { 
	      for ( int i12 = 0 ; i12 < mt+1 ; i12++) { 
		  
		  for ( ri2 = A+1 ; ri2 < B ; ri2++) {
		      index2 = idx(ri2,i22,i12);
		      V[index2] = V[idx(A,i22,i12)] + 
          			  (ri2-A) * 
			          ( V[idx(B,i22,i12)] - V[idx(A,i22,i12)]) / ScaleFactor;
		      
		  }
	      }
	  }

      }
  } // if mtin < mtout

  if ( dptr->mvis->mt > mt ) {
      ScaleFactor = dptr->mvis->mt / mt;

      for ( int ri2=0 ; ri2 < nr ; ri2++) {
	  int ri1 = ri2 * ScaleFactor;
	  for ( int i22 = 0 ; i22 < mt+1 ; i22++) { 
	      int i21 = i22 * ScaleFactor;
	      for ( int i12 = 0 ; i12 < mt+1 ; i12++) { 
		  int i11 = i12 * ScaleFactor;

		  index1 = dptr->mvis->domains[id].idx(ri1,i21,i11);
		  index2 = idx(ri2,i22,i12);
		  V[index2] = dptr->mvis->domains[id].V[index1];
	      }
	  }
      }
  }// if mtin > mtout

  return 0;
}

////////////////////////////////////////
// Domain::getValue
//
// Given the index of the geomorph grid, uses the 'interp' method to calculate
// its value from that contained in Data::x,y,z, and V.
int Domain::getValue(Data * dptr, int index, int interp)
{

  if ( interp == dptr->NEAREST ){
    // brute-force...
    if (dptr->intype == dptr->FILT) {
      V[index] = getNearestDataValue(dptr,index);
    }
    if (dptr->intype == dptr->MITP) {
      V[index] = getNearestDataValue(dptr,index);
    }
  }else if ( interp == dptr->NEAREST2 ){
    // a slightly more intelligent routine...
    if (dptr->intype == dptr->FILT) {
      V[index] = getNearestDataValue2_filt(dptr,index);
    }
    if (dptr->intype == dptr->MITP) {
      V[index] = getNearestDataValue2_mitp(dptr,index);
      // the slightly more intelligent routine, plus GPU computation
      //V[index] = getNearestDataValue2_gpu(dptr->ndpth, dptr->minR, dptr->maxR, 
      //                                    dptr->nlat, dptr->nlng,
      //				    xn, yn, zn,
      //				    dptr->x,  dptr->y,  dptr->z, dptr->V,
      //
    }
  }else if ( interp == dptr->LINEAR ){          
  }else if ( interp == dptr->CUBIC ){
  }
  
  return 0;
}

////////////////////////////////////////
// Domain::importData
//
// For each xn element, find the best match data value V using the defined
// algorithm
//
bool Domain::importData(Data *dptr)
{
    int index=0;

    // Here we decide which algorithm to use depending on input type. It moves
    // data found in dptr to the current domain.

    // ...Non-terra based grids
    if ( dptr->intype == dptr->MITP || dptr->intype == dptr->FILT ) {

      for ( int ri=0 ; ri < nr ; ri++){
	printf("......Layer %d\n",ri);
	for ( int i2 = 0 ; i2 < mt+1 ; i2++) {
	  for ( int i1 = 0 ; i1 < mt+1 ; i1++) {
	    index = idx(ri,i2,i1);
	    getValue(dptr, index, dptr->interp);
	  }
	}
      }
    }else{
      // ... must be intype == MVIS or TERRA.
      getLinearDataValue(dptr);
    }

    return 0; // success
}

////////////////////////////////////////
// Domain::importData_gpu
//
// A simple interfaceto the gpu-enabled code
//
bool Domain::importData_gpu(Data *dptr)
{

  importData_c_gpu(nr,mt,
		   dptr->ndpth, dptr->minR, dptr->maxR, 
		   dptr->nlat, dptr->nlng,
		   xn, yn, zn, V,
		   dptr->x,  dptr->y,  dptr->z, dptr->V);
  
    return 0; // success
}

////////////////////////////////////////
// Domain::sqrti
//
// Returns the (floored int) integer square root of 'a'.
//   2 <= sqrt(a) <= 2^15
//
int Domain::sqrti(int a)
{
    int i;

    if (a<1)  return -1;

    for ( i=1; i<=32768; i++){
	if (i*i > a ){
	    i--;
	    return i;
	}
    }

    return -1; // fail
}

////////////////////////////////////////
// Domain::importMVIS
//
// import data from MVIS files
// 
bool Domain::importMVIS(FILE * fptr, int proc, int nt, double cmb)
{
    
    float buf;
    int   index=0;
    int   index1=0;
    int   i1start,i1end,i2start,i2end;
    div_t divresult;

    divresult = div(proc, mt/nt);

    i1start = nt * divresult.rem;
    i1end   = i1start + nt;
    
    i2start = nt * divresult.quot ;
    i2end   = i2start + nt;

    // cycle through our 'process' points
    for ( int i2=i2start ; i2<=i2end ; i2++ ){
	for ( int i1=i1start ; i1<=i1end ; i1++ ){

	    index1 = idx(0,i2,i1);
	    fscanf(fptr,"%f",&buf);
	    xn[index1] = buf;
	    fscanf(fptr,"%f",&buf);
	    yn[index1] = buf; 
	    fscanf(fptr,"%f",&buf); 
	    zn[index1] = buf;

	    for ( int ir=nr-1 ; ir>=0 ; ir-- ){
		index = idx(ir,i2,i1);
		fscanf(fptr,"%f",&buf);
		V[index] = buf;
		
		xn[index] = xn[index1] - xn[index1]*ir*(1-cmb)/(nr-1);
		yn[index] = yn[index1] - yn[index1]*ir*(1-cmb)/(nr-1);
		zn[index] = zn[index1] - zn[index1]*ir*(1-cmb)/(nr-1);

	    }
	}
    }
    
    return 0; // success
}

////////////////////////////////////////
// Domain::importTERRA
//
// import Data::dptr in the TERRA format (convection or circulation  model)
// 
bool Domain::importTERRA(FILE * fptr, int proc, int nt, int ir, int tvp)
{
    float buf;
    int   index=0;
    int   i1start,i1end,i2start,i2end;
    div_t divresult;

    divresult = div(proc, mt/nt);

    i1start = nt * divresult.rem;
    i1end   = i1start + nt;
    
    i2start = nt * divresult.quot ;
    i2end   = i2start + nt;

    for ( int i2=i2start ; i2<=i2end ; i2++ ){
      for ( int i1=i1start ; i1<=i1end ; i1++ ){
	
	if (tvp == 0) {
	  index = idx(ir,i2,i1);
	  fscanf(fptr,"%f",&buf);
	  V[index] = buf;
	}else if (tvp == 1) {
	  for ( int xyz=0 ; xyz<3 ; xyz++) {
	    fscanf(fptr,"%f",&buf);
	    vel[3*index+xyz] = buf;
	  }
	}else if (tvp == 2) {
	  fscanf(fptr,"%f",&buf);
	  P[index] = buf;
	}else {
	  return 1; // failure
	}
      } // i1
    } // i2
    
    return 0; // success
}

////////////////////////////////////////
// Domain::exportMVIS
//
// export Data::dptr in the MVIS format.
// 
bool Domain::exportMVIS(FILE * fptr, int proc, int nt)
{
    
    int   index=0;
    int   i1start,i1end,i2start,i2end;
    div_t divresult;

    divresult = div(proc, mt/nt);

    i1start = nt * divresult.rem;
    i1end   = i1start + nt;
    
    i2start = nt * divresult.quot ;
    i2end   = i2start + nt;

    // cycle through our 'process' points
    for ( int i2=i2start ; i2<=i2end ; i2++ ){
	for ( int i1=i1start ; i1<=i1end ; i1++ ){
	    index = idx(0,i2,i1);
	    fprintf(fptr,"%16.8E\t%16.8E\t%16.8E\n",xn[index],yn[index],zn[index]);
	    for ( int ir=nr-1 ; ir>=0 ; ir-- ){
		index = idx(ir,i2,i1);
		fprintf(fptr,"%16.8E\n",V[index]);
	    }
	}
    }
    
    return 0; // success
}

////////////////////////////////////////
// Domain::exportTERRA
//
// export Data::dptr in the TERRA format (convection or circulation  model)
// 
bool Domain::exportTERRA(FILE * fptr, int proc, int nt, int ir, int tvp, long int &colcntr)
{
   
    int   index=0;
    int   i1start,i1end,i2start,i2end;
    div_t divresult;

    divresult = div(proc, mt/nt);

    i1start = nt * divresult.rem;
    i1end   = i1start + nt;
    
    i2start = nt * divresult.quot ;
    i2end   = i2start + nt;

    for ( int i2=i2start ; i2<=i2end ; i2++ ){
	for ( int i1=i1start ; i1<=i1end ; i1++ ){
	    
	    if (tvp == 0) {
		index = idx(ir,i2,i1);
		fprintf(fptr,"%10.3f",V[index]);
		if ( colcntr%15 == 0 )  fprintf(fptr,"\n"); // print in columns of 15
		colcntr++;
	    }else if (tvp == 1) {
		for ( int xyz=0 ; xyz<3 ; xyz++) {
		    fprintf(fptr,"%10.3E",0.);
		    if ( colcntr%15 == 0 )  fprintf(fptr,"\n"); // print in columns of 15
		    colcntr++;
		}
	    }else if (tvp == 2) {
		fprintf(fptr,"%10.3E",0.);
		if ( colcntr%15 == 0 )  fprintf(fptr,"\n"); // print in columns of 15
		colcntr++;
	    }else {
		return 1; // failure
	    }
	} // i1
    } // i2

    return 0; // success
}

////////////////////////////////////////
// Domain::getAngle3d
//
// Given two vectors Aand B, getAngle3d returns the angle A0B.
//
double Domain::getAngle3d(double ax, double ay, double az, double bx, double by, double bz)
{
    double maga = sqrt( ax*ax + ay*ay + az*az);
    double magb = sqrt( bx*bx + by*by + bz*bz);

    return acos( dotProduct(ax,ay,az,bx,by,bz)/maga/magb );
}

////////////////////////////////////////
// Domain::rotate3d
//
// Given the point to be rotated, the rotation axis _x,_y_z, and the rotation
// angle, rotates the point about the axis by phi radians. The resulting point
// is returned in px,py,pz.
//  Assumes the input point (px^2 + py^2 +  pz^2) = 1.
bool Domain::rotate3d(double *px, double *py, double *pz, double rx, double ry, double rz, double phi)
{
  double R[3][3], A[3], B[3];

  B[0] = *px;
  B[1] = *py;
  B[2] = *pz;

//printf("px; py; pz = %12.8E;  %12.8E;  %12.8E \n",px, py, pz);
//printf("B0; B1; B2 = %12.8E;  %12.8E;  %12.8E \n",B[0], B[1], B[2]);

  A[0] = 0.;
  A[1] = 0.;
  A[2] = 0.;

  R[0][0] = cos(phi)    + rx*rx*(1-cos(phi));
  R[0][1] = rx*ry*(1-cos(phi)) - rz*sin(phi);
  R[0][2] = rx*rz*(1-cos(phi)) + ry*sin(phi);
  
  R[1][0] = ry*rx*(1-cos(phi)) + rz*sin(phi);
  R[1][1] = cos(phi)    + ry*ry*(1-cos(phi));
  R[1][2] = ry*rz*(1-cos(phi)) - rx*sin(phi);

  R[2][0] = rz*rx*(1-cos(phi)) - ry*sin(phi);
  R[2][1] = rz*ry*(1-cos(phi)) + rx*sin(phi);
  R[2][2] = cos(phi)    + rz*rz*(1-cos(phi));

  // Rotate the vertices by Ry
  for ( int ix = 0 ; ix < 3 ; ix++ ){
    for ( int iy = 0 ; iy < 3 ; iy++ ){
      A[ix] += R[ix][iy]*B[iy];
    }
  }	
  
  px[0] = A[0];
  py[0] = A[1];
  pz[0] = A[2];
  
  return 0; // success

}

////////////////////////////////////////
// Domain::normalise
//
// Returns the dot product of the vectors p1 x p2
//
bool Domain::normalise(double *ax, double *ay, double *az)
{
    double maga = sqrt(ax[0]*ax[0] + ay[0]*ay[0] + az[0]*az[0]);
    ax[0]=ax[0]/maga;
    ay[0]=ay[0]/maga;
    az[0]=az[0]/maga;
    
    return 0;
}

////////////////////////////////////////
// Domain::dotProduct
//
// Returns the dot product of the vectors p1 x p2
//
double Domain::dotProduct(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z)
{
    return p1x*p2x + p1y*p2y + p1z*p2z; 
}

////////////////////////////////////////
// Domain::crossProduct
//
// Overwrites p2* with the cross-product of  p1 x p2
//
bool Domain::crossProduct(double p1x, double p1y, double p1z, double *p2x, double *p2y, double *p2z)
{
    double ax,ay,az;
    
    ax = p1y*p2z[0] - p1z*p2y[0];
    ay = p1z*p2x[0] - p1x*p2z[0];
    az = p1x*p2y[0] - p1y*p2x[0];

    p2x[0] = ax;
    p2y[0] = ay;
    p2z[0] = az;

    return 0; //success
}

////////////////////////////////////////
// Domain::grdgen
//
//
bool Domain::grdgen(double cmb)
{
    
    int    index;
    double x0,y0,z0;

    double Tbx,Tby,Tbz,Tcx,Tcy,Tcz;

    double a,tau,rho,u,v,Beta,phi;
    double Ry[3][3], A[12][3], Ad[12][3];

    tau = (sqrt(5) + 1)/2;

    a=1.;
    rho=tau-1;
    u=a/(sqrt(1+pow(rho,2)));
    v=rho*u;
    phi = 2 * asin( 1/sqrt(tau*sqrt(5)) );
    
    Beta= atan(v/u);
    Ry[0][0] = cos(Beta);
    Ry[0][1] = 0;
    Ry[0][2] = -sin(Beta);

    Ry[1][0] = 0;
    Ry[1][1] = 1;
    Ry[1][2] = 0;

    Ry[2][0] = sin(Beta);
    Ry[2][1] = 0;
    Ry[2][2] = cos(Beta);

    for ( int i = 0 ; i < 12 ; i++ ){
      Ad[i][0]=0;
      Ad[i][1]=0;
      Ad[i][2]=0;
    }
  
    // Set up the intial regular icosahedral grid
    // Note: this grid needs to be rotated about y-axis in order to match the TERRA grid.
    A[0][0] = v;  A[0][1] = 0; A[0][2] = u; 
    A[1][0] = u;  A[1][1] = v; A[1][2] = 0; 
    A[2][0] = 0;  A[2][1] = u; A[2][2] = v;   
    A[3][0] =-v;  A[3][1] = 0; A[3][2] = u;  
    A[4][0] = 0;  A[4][1] =-u; A[4][2] = v;  
    A[5][0] = u;  A[5][1] =-v; A[5][2] = 0;
    A[6][0] = v;  A[6][1] = 0; A[6][2] =-u; 
    A[7][0] = 0;  A[7][1] = u; A[7][2] =-v;
    A[8][0] =-u;  A[8][1] = v; A[8][2] = 0;
    A[9][0] =-u;  A[9][1] =-v; A[9][2] = 0;
    A[10][0]= 0;  A[10][1]=-u; A[10][2]=-v;  
    A[11][0]=-v;  A[11][1]= 0; A[11][2]=-u;   

    // Rotate the vertices by Ry
    for ( int i = 0 ; i < 12 ; i++ ){
	for ( int ix = 0 ; ix < 3 ; ix++ ){
	    for ( int iy = 0 ; iy < 3 ; iy++ ){
		Ad[i][ix] += Ry[ix][iy]*A[i][iy];
	    }
	}	
    }
   
    // Vertex Points........
    // Northern Hemisphere
    if (id<5) {
	// North Pole - 0,0
	index = idx(0, 0, 0);
	xn[index] = Ad[0][0];
	yn[index] = Ad[0][1];
	zn[index] = Ad[0][2];
	// mt,0
	index = idx(0, mt, 0);
	xn[index] = Ad[id+1][0];
	yn[index] = Ad[id+1][1];
	zn[index] = Ad[id+1][2]; 
	// 0,mt
	index = idx(0, 0, mt);
	if (id == 0) {
	    xn[index] = Ad[id+5][0];
	    yn[index] = Ad[id+5][1];
	    zn[index] = Ad[id+5][2]; 
	}else{
	    xn[index] = Ad[id][0];
	    yn[index] = Ad[id][1];
	    zn[index] = Ad[id][2]; 
	}
	// mt,mt
	index = idx(0, mt, mt);
	xn[index] = Ad[id+6][0];
	yn[index] = Ad[id+6][1];
	zn[index] = Ad[id+6][2]; 
	
	// Southern Hemisphere
    }else{
	// South Pole - 0,0
	index = idx(0, 0, 0);
	xn[index] = Ad[11][0];
	yn[index] = Ad[11][1];
	zn[index] = Ad[11][2];
	// mt,0
	index = idx(0, mt, 0);
	if (id == 9) {
	    xn[index] = Ad[id-3][0];
	    yn[index] = Ad[id-3][1];
	    zn[index] = Ad[id-3][2]; 
	}else{
	    xn[index] = Ad[id+2][0];
	    yn[index] = Ad[id+2][1];
	    zn[index] = Ad[id+2][2]; 
	}
	// 0,mt
	index = idx(0, 0, mt);
	xn[index] = Ad[id+1][0];
	yn[index] = Ad[id+1][1];
	zn[index] = Ad[id+1][2]; 
	// mt,mt
	index = idx(0, mt, mt);
	xn[index] = Ad[id-4][0];
	yn[index] = Ad[id-4][1];
	zn[index] = Ad[id-4][2]; 
    }
    
    // Edge Points......... (assumes unit radius)
    // i2,i1=0
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, 0);
	
	Tbx = xn[idx(0,0,0)];
	Tby = yn[idx(0,0,0)];
	Tbz = zn[idx(0,0,0)];

	Tcx = xn[idx(0, mt, 0)];
	Tcy = yn[idx(0, mt, 0)];
	Tcz = zn[idx(0, mt, 0)];
	
	phi = getAngle3d(Tbx,Tby,Tbz,Tcx,Tcy,Tcz);
	crossProduct(Tbx,Tby,Tbz,&Tcx,&Tcy,&Tcz);
	normalise(&Tcx,&Tcy,&Tcz);
	rotate3d(&Tbx,&Tby,&Tbz,Tcx,Tcy,Tcz, i2*phi/mt );
	normalise(&Tbx,&Tby,&Tbz); // just in case
	
	xn[index] = Tbx;
	yn[index] = Tby;
	zn[index] = Tbz;
    }
    // i2=0,i1
    for ( int i1 = 1 ; i1 < mt ; i1++ ){

	index=idx(0, 0, i1);
	
	Tbx = xn[idx(0, 0, 0)];
	Tby = yn[idx(0, 0, 0)];
	Tbz = zn[idx(0, 0, 0)];

	Tcx = xn[idx(0, 0, mt)];
	Tcy = yn[idx(0, 0, mt)];
	Tcz = zn[idx(0, 0, mt)];
	
	phi = getAngle3d(Tbx,Tby,Tbz,Tcx,Tcy,Tcz);
	crossProduct(Tbx,Tby,Tbz,&Tcx,&Tcy,&Tcz);
	normalise(&Tcx,&Tcy,&Tcz);
	rotate3d(&Tbx,&Tby,&Tbz,Tcx,Tcy,Tcz, i1*phi/mt );
	normalise(&Tbx,&Tby,&Tbz); // just in case

	xn[index] = Tbx;
	yn[index] = Tby;
	zn[index] = Tbz;
    }
    // i2,i1=mt
    for ( int i2 = 1 ; i2 < mt ; i2++ ){

	index=idx(0, i2, mt);
	
	Tbx = xn[idx(0, 0, mt)];
	Tby = yn[idx(0, 0, mt)];
	Tbz = zn[idx(0, 0, mt)];

	Tcx = xn[idx(0, mt, mt)];
	Tcy = yn[idx(0, mt, mt)];
	Tcz = zn[idx(0, mt, mt)];
	
	phi = getAngle3d(Tbx,Tby,Tbz,Tcx,Tcy,Tcz);
	crossProduct(Tbx,Tby,Tbz,&Tcx,&Tcy,&Tcz);
	normalise(&Tcx,&Tcy,&Tcz);
	rotate3d(&Tbx,&Tby,&Tbz,Tcx,Tcy,Tcz, i2*phi/mt );
	normalise(&Tbx,&Tby,&Tbz); // just in case
	
	xn[index] = Tbx;
	yn[index] = Tby;
	zn[index] = Tbz;
    }
    // i1,i2=mt
    for ( int i1 = 1 ; i1 < mt ; i1++ ){

	index=idx(0, mt, i1);
	
	Tbx = xn[idx(0, mt, 0)];
	Tby = yn[idx(0, mt, 0)];
	Tbz = zn[idx(0, mt, 0)];

	Tcx = xn[idx(0, mt, mt)];
	Tcy = yn[idx(0, mt, mt)];
	Tcz = zn[idx(0, mt, mt)];
	
	phi = getAngle3d(Tbx,Tby,Tbz,Tcx,Tcy,Tcz);
	crossProduct(Tbx,Tby,Tbz,&Tcx,&Tcy,&Tcz);
	normalise(&Tcx,&Tcy,&Tcz);
	rotate3d(&Tbx,&Tby,&Tbz,Tcx,Tcy,Tcz, i1*phi/mt );
	normalise(&Tbx,&Tby,&Tbz); // just in case
	
	xn[index] = Tbx;
	yn[index] = Tby;
	zn[index] = Tbz;
    }
   
    // Everywhere inbetween...
    //
    // There's a bug here somwhere. It seems this algoritm (edge point
    // rotation) doesn't match the TERRA-Grid.
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
        for ( int i1 = 1 ; i1 < mt ; i1++ ){


	    index=idx(0, i2, i1);

	    Tbx = xn[idx(0, 0, i1)];
	    Tby = yn[idx(0, 0, i1)];
	    Tbz = zn[idx(0, 0, i1)];
	    
	    Tcx = xn[idx(0, mt, i1)];
	    Tcy = yn[idx(0, mt, i1)];
	    Tcz = zn[idx(0, mt, i1)];
	    
	    phi = getAngle3d(Tbx,Tby,Tbz,Tcx,Tcy,Tcz);
	    crossProduct(Tbx,Tby,Tbz,&Tcx,&Tcy,&Tcz);
	    normalise(&Tcx,&Tcy,&Tcz);
	    rotate3d(&Tbx,&Tby,&Tbz,Tcx,Tcy,Tcz, i2*phi/mt );
	    normalise(&Tbx,&Tby,&Tbz); // just in case
	    
	    xn[index] = Tbx;
	    yn[index] = Tby;
	    zn[index] = Tbz;
	}
    }
    
    // generate radial points (we already know ir=0)
    for ( int ir = 1 ; ir < nr ; ir++){
      index=idx(ir, 0, 0);
      for ( int i2 = 0 ; i2 < mt+1 ; i2++ ){
	for ( int i1 = 0 ; i1 < mt+1 ; i1++ ){
	  x0 = xn[idx(0,i2,i1)];
	  y0 = yn[idx(0,i2,i1)];
	  z0 = zn[idx(0,i2,i1)];
	  
	  xn[index] = x0 - x0*(a-cmb)*ir/nr;
	  yn[index] = y0 - y0*(a-cmb)*ir/nr;
	  zn[index] = z0 - z0*(a-cmb)*ir/nr;
	  index++;
	}
      }
    }

    return 0;
}

////////////////////////////////////////
// Domain::midpt
//
//
bool Domain::midpt(double *x, double *y, double *z,
		   double x1, double y1, double z1,
		   double x2, double y2, double z2)
{
  x[0] = x1 + x2;
  y[0] = y1 + y2;;
  z[0] = z1 + z2;

  double norm = 1./sqrt(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);

  x[0] = norm*x[0];
  y[0] = norm*y[0];
  z[0] = norm*z[0];

  return 0;
}


////////////////////////////////////////
// Domain::grdgen2
//
// Should be used over grdgen().
//
bool Domain::grdgen2(double cmb)
{
    
    int    index,index1,index2;
    double x0,y0,z0;

    double a,tau,rho,u,v,Beta;
    double Ry[3][3], A[12][3], Ad[12][3];

    int    lvt = int (1.45*log(mt));
    int    m,l,l2,i1,i2;

    tau = (1+sqrt(5))/2;

    a=1.;
    tau=0.5*(sqrt(5)+1);
    rho=tau-1;
    u=a/(sqrt(1+pow(rho,2)));
    v=rho*u;
    
    Beta= atan(v/u);
    Ry[0][0] = cos(Beta);
    Ry[0][1] = 0;
    Ry[0][2] = -sin(Beta);

    Ry[1][0] = 0;
    Ry[1][1] = 1;
    Ry[1][2] = 0;

    Ry[2][0] = sin(Beta);
    Ry[2][1] = 0;
    Ry[2][2] = cos(Beta);

    for ( int i = 0 ; i < 12 ; i++ ){
      Ad[i][0]=0;
      Ad[i][1]=0;
      Ad[i][2]=0;
    }
  
    // Set up the intial regular icosahedral grid
    // Note: this grid needs to be rotated about y-axis in order to match the TERRA grid.
    A[0][0] = v;  A[0][1] = 0; A[0][2] = u; 
    A[1][0] = u;  A[1][1] = v; A[1][2] = 0; 
    A[2][0] = 0;  A[2][1] = u; A[2][2] = v;   
    A[3][0] =-v;  A[3][1] = 0; A[3][2] = u;  
    A[4][0] = 0;  A[4][1] =-u; A[4][2] = v;  
    A[5][0] = u;  A[5][1] =-v; A[5][2] = 0;
    A[6][0] = v;  A[6][1] = 0; A[6][2] =-u; 
    A[7][0] = 0;  A[7][1] = u; A[7][2] =-v;
    A[8][0] =-u;  A[8][1] = v; A[8][2] = 0;
    A[9][0] =-u;  A[9][1] =-v; A[9][2] = 0;
    A[10][0]= 0;  A[10][1]=-u; A[10][2]=-v;  
    A[11][0]=-v;  A[11][1]= 0; A[11][2]=-u;   

    // Rotate the vertices by Ry
    for ( int i = 0 ; i < 12 ; i++ ){
	for ( int ix = 0 ; ix < 3 ; ix++ ){
	    for ( int iy = 0 ; iy < 3 ; iy++ ){
		Ad[i][ix] += Ry[ix][iy]*A[i][iy];
	    }
	}	
    }
   
    // Vertex Points........
    // Northern Hemisphere
    if (id<5) {
	// North Pole
	index = idx(0, 0, 0);
	xn[index] = Ad[0][0];
	yn[index] = Ad[0][1];
	zn[index] = Ad[0][2];
	// mt,0
	index = idx(0, mt, 0);
	xn[index] = Ad[id+1][0];
	yn[index] = Ad[id+1][1];
	zn[index] = Ad[id+1][2]; 
	// 0,mt
	index = idx(0, 0, mt);
	if (id == 0) {
	    xn[index] = Ad[id+5][0];
	    yn[index] = Ad[id+5][1];
	    zn[index] = Ad[id+5][2]; 
	}else{
	    xn[index] = Ad[id][0];
	    yn[index] = Ad[id][1];
	    zn[index] = Ad[id][2]; 
	}
	// mt,mt
	index = idx(0, mt, mt);
	xn[index] = Ad[id+6][0];
	yn[index] = Ad[id+6][1];
	zn[index] = Ad[id+6][2]; 
	
	// Southern Hemisphere
    }else{
	// South Pole
	index = idx(0, 0, 0);
	xn[index] = Ad[11][0];
	yn[index] = Ad[11][1];
	zn[index] = Ad[11][2];
	// mt,0
	index = idx(0, mt, 0);
	if (id == 9) {
	    xn[index] = Ad[id-3][0];
	    yn[index] = Ad[id-3][1];
	    zn[index] = Ad[id-3][2]; 
	}else{
	    xn[index] = Ad[id+2][0];
	    yn[index] = Ad[id+2][1];
	    zn[index] = Ad[id+2][2]; 
	}
	// 0,mt
	index = idx(0, 0, mt);
	xn[index] = Ad[id+1][0];
	yn[index] = Ad[id+1][1];
	zn[index] = Ad[id+1][2]; 
	// mt,mt
	index = idx(0, mt, mt);
	xn[index] = Ad[id-4][0];
	yn[index] = Ad[id-4][1];
	zn[index] = Ad[id-4][2]; 
    }

    
    for ( int k = 0; k < lvt ; k++){
      
      m  = int (pow(2,k)+0.1);
      l  = mt/m;
      l2 = l/2;


      // rows of diamond--
      for (int j1=0 ; j1 < m + 1 ; j1++ ){
	for (int j2=0 ; j2 < m  ; j2++ ){

	  i1 = (j1)*l;
	  i2 = (j2)*l + l2;
	  index = idx(0,i2,i1);
	  index1= idx(0,i2-l2,i1);
	  index2= idx(0,i2+l2,i1);
	  midpt(
		&xn[index], &yn[index], &zn[index],
                xn[index1],yn[index1],zn[index1],
                xn[index2],yn[index2],zn[index2]);
	  
	} // j2
      } // j1


      // columns of diamond--
      for (int j1=0 ; j1 < m + 1 ; j1++ ){
	for (int j2=0 ; j2 < m  ; j2++ ){

	  i1 = (j2)*l + l2;
	  i2 = (j1)*l;
	  index = idx(0,i2,i1);
	  index1= idx(0,i2,i1-l2);
	  index2= idx(0,i2,i1+l2);
	  midpt(
		&xn[index], &yn[index], &zn[index],
                xn[index1],yn[index1],zn[index1],
                xn[index2],yn[index2],zn[index2]);
	  
	} // j2
      } // j1


      // diagonals of diamond--
      for (int j1=0 ; j1 < m  ; j1++ ){
	for (int j2=0 ; j2 < m  ; j2++ ){

	  i1 = (j1)*l + l2;
	  i2 = (j2)*l + l2;
	  index = idx(0,i2,i1);
	  index1= idx(0,i2+l2,i1-l2);
	  index2= idx(0,i2-l2,i1+l2);
	  midpt(
		&xn[index], &yn[index], &zn[index],
                xn[index1],yn[index1],zn[index1],
                xn[index2],yn[index2],zn[index2]);
	  
	} // j2
      } // j1

    } // k


    // generate radial points (we already know ir=0)
    for ( int ir = 1 ; ir < nr ; ir++){
      index=idx(ir, 0, 0);
      for ( int i2 = 0 ; i2 < mt+1 ; i2++ ){
	for ( int i1 = 0 ; i1 < mt+1 ; i1++ ){
	  x0 = xn[idx(0,i2,i1)];
	  y0 = yn[idx(0,i2,i1)];
	  z0 = zn[idx(0,i2,i1)];
	  
	  xn[index] = x0 - x0*(a-cmb)*ir/nr;
	  yn[index] = y0 - y0*(a-cmb)*ir/nr;
	  zn[index] = z0 - z0*(a-cmb)*ir/nr;
	  index++;
	}
      }
    }
    

    return 0;
}
