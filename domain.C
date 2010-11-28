/*
 *
 * Data Manipulation within Domains
 *
 * 
 */


#include "main.H"
#include "domain.H"

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
    : id(), northern(), xn(), mt(), nr()
{
}

////////////////////////////////////////
// Domain::Domain
// 
//
Domain::Domain(int _id, bool _northern)
    : id(), northern(), xn(), mt(), nr()
{
    if ( _id >= 10 ) {
	printf("Domain::Domain fail");
    }
    
    id=_id;
    northern=_northern;
}

////////////////////////////////////////
// Domain::Domain
// 
//
Domain::Domain(int _id, bool _northern, int _mt, int _nr)
  : id(), northern(), xn(), mt(), nr()
{
    if (_id >= 10) {
	printf("Domain::Domain fail");
    }
    
    id=_id;
    northern=_northern;
    mt = _mt;
    nr = _nr;

    xn = new double[nr*(mt+1)*(mt+1)];
    yn = new double[nr*(mt+1)*(mt+1)];
    zn = new double[nr*(mt+1)*(mt+1)];
    V  = new double[nr*(mt+1)*(mt+1)];
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
    
    id < 5 ? northern = true : northern = false;

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

  if (nerror!=0) printf("Domain::idx Error: ir = %d; i2 = %d; i1 = %d\n",r,i2,i1);

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
double Domain::getNearestDataValue2(Data *dptr, int index)
{
    int nr; // our layer of interest in Data:dptr (or in fact either side of this!)

    double gx,gy,gz;
    gx = xn[index];
    gy = yn[index];
    gz = zn[index];
    
    // what's our current radius
    double rad=sqrt(gx*gx + gy*gy + gz*gz);

    // find closest radial layer in Data::dptr nr
    double dR=veryLarge;
    double dataR=0.;
//    double tmpR=0;   // only needed for temporary printf statement
    for ( int ir=0 ; ir<dptr->ndpth ; ir++ ){
      dataR = dptr->minR + ir*(dptr->maxR - dptr->minR)/dptr->ndpth;
      if (fabs(dataR - rad) < dR) {
	dR = fabs(dataR - rad) ;
	nr = dptr->ndpth - ir;  // reverse ordering.
//	tmpR = dataR;
      }
    }

//    printf("grad=%12.8g\tdrad=%12.8g\tnr=%d\n",rad,tmpR,nr);

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
      }
    }

    return dataV;
}

////////////////////////////////////////
// Domain::getValue
//
// Given the index of the geomorph grid, uses the 'interp' method to calculate
// its value from that contained in Data::x,y,z, and V.
int Domain::getValue(Data * dptr, int index, int interp)
{
  switch (interp) {
  case NEAREST:
    //   V[index] = getNearestDataValue(dptr,index);
    V[index] = getNearestDataValue2(dptr,index);
    break;
  case LINEAR:
    //do something else;
    break;
  default:
    break;
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
    int interp=NEAREST;
    
    for ( int ri=0 ; ri < nr ; ri++){
    // layer 12 only...
    // for ( int ri=12 ; ri < 13 ; ri++){
	printf("......Layer %d\n",ri);
	for ( int i2 = 0 ; i2 < mt+1 ; i2++) {
	    for ( int i1=0 ; i1 < mt+1 ; i1++) {
		index=idx(ri,i2,i1);
		getValue(dptr, index, interp)	;
	    }
	}
    }

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
// Domain::exportMVIS
//
// export Data::dptr in the MVIS format.
// 
bool Domain::exportMVIS(FILE * fptr, int nproc, int proc, int nt)
{
    
    int   index=0;
    int   i1start,i1end,i2start,i2end;
    div_t divresult;

    divresult = div(proc, mt/nt);

//    i1start = nt * ( proc % (mt/nt) );
    i1start = nt * divresult.rem;
    i1end   = i1start + nt;
    
    i2start = nt * divresult.quot ;
    i2end   = i2start + nt;

    // cycle through our 'process' points
    
    for ( int i2=i2start ; i2<=i2end ; i2++ ){
	for ( int i1=i1start ; i1<=i1end ; i1++ ){
	    index = idx(0,i2,i1);
	    fprintf(fptr,"%16.8E\t%16.8E\t%16.8E\n",xn[index],yn[index],zn[index]);
	    for ( int ir=0 ; ir<nr ; ir++ ){
		index = idx(ir,i2,i1);
		fprintf(fptr,"%16.8E\n",V[index]);
	    }
	}
    }
    
    return 0; // success
}

////////////////////////////////////////
// Domain::rotate3d
//
// Given the point to be rotated, the rotation axis _x,_y_z, and the rotation
// angle, rotates the point about the axis by phi radians. The resulting point
// is returned in px,py,pz.
//  Assumes the input point (px^2 + py^2 +  pz^2) = 1.
bool Domain::rotate3d(double *px, double *py, double *pz, double rx,double ry, double rz, double phi)
{
  double R[3][3], A[3], B[3];

  B[0] = *px;
  B[1] = *py;
  B[2] = *pz;

  printf("px; py; pz = %12.8E;  %12.8E;  %12.8E \n",px, py, pz);
  printf("B0; B1; B2 = %12.8E;  %12.8E;  %12.8E \n",B[0], B[1], B[2]);

  A[0] = 0.;
  A[1] = 0.;
  A[2] = 0.;

  R[0][0] = cos(phi) + rx*rx*(1-cos(phi));
  R[0][1] = rx*ry*(1-cos(phi)) - rz*sin(phi);
  R[0][2] = rx*rz*(1-cos(phi)) + ry*sin(phi);
  
  R[1][0] = ry*rx*(1-cos(phi)) + rz*sin(phi);
  R[1][1] = cos(phi) + ry*ry*(1-cos(phi));
  R[1][2] = ry*rz*(1-cos(phi)) - rx*sin(phi);

  R[2][0] = rz*rx*(1-cos(phi)) - ry*sin(phi);
  R[2][1] = rz*ry*(1-cos(phi)) + rx*sin(phi);
  R[2][2] = cos(phi) + rz*rz*(1-cos(phi));

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
// Domain::grdgen
//
//
bool Domain::grdgen(double cmb)
{
    
    int    index;
    double R,x0,y0,z0;

    double TA,Tax,Tay,Taz,Tbx,Tby,Tbz,Tcx,Tcy,Tcz,F,E;

    double a,tau,rho,u,v,Beta;
    double Ry[3][3], A[12][3], Ad[12][3];

    tau = (1+sqrt(5))/2;
    TA  = asin( 1/(sqrt(tau*sqrt(5))) );

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
    
    // Edge Points......... (assumes unit radius)
    // i2,i1=0
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, 0);
	
	Tbx = xn[idx(0,i2-1,0)];
	Tby = yn[idx(0,i2-1,0)];
	Tbz = zn[idx(0,i2-1,0)];

	Tcx = xn[index];
	Tcy = yn[index];
	Tcz = zn[index];


	/*
	xn[index] =   
	    xn[idx(0, 0, 0)] + i2*(xn[idx(0, mt, 0)] - xn[idx(0, 0, 0)])/mt;
	yn[index] =   
	    yn[idx(0, 0, 0)] + i2*(yn[idx(0, mt, 0)] - yn[idx(0, 0, 0)])/mt;
	zn[index] =   
	    zn[idx(0, 0, 0)] + i2*(zn[idx(0, mt, 0)] - zn[idx(0, 0, 0)])/mt;
	*/
	R = sqrt(pow(xn[index],2) +  
	         pow(yn[index],2) +
	         pow(zn[index],2) );
	xn[index] = 1/R * xn[index];
	yn[index] = 1/R * yn[index];
	zn[index] = 1/R * zn[index];
    }
    // i1,i2=0
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	index=idx(0, 0, i1);
	xn[index] =   
	    xn[idx(0, 0, 0)] + i1*(xn[idx(0, 0, mt)] - xn[idx(0, 0, 0)])/mt;
	yn[index] =   
	    yn[idx(0, 0, 0)] + i1*(yn[idx(0, 0, mt)] - yn[idx(0, 0, 0)])/mt;
	zn[index] =   
	    zn[idx(0, 0, 0)] + i1*(zn[idx(0, 0, mt)] - zn[idx(0, 0, 0)])/mt;
	
	R = sqrt(pow(xn[index],2) +  
	         pow(yn[index],2) +
	         pow(zn[index],2) );
	xn[index] = 1/R * xn[index];
	yn[index] = 1/R * yn[index];
	zn[index] = 1/R * zn[index];
    }
    // i2,i1=mt
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, mt);
	xn[index] =   
	    xn[idx(0, 0, mt)] + i2*(xn[idx(0, mt, mt)] - xn[idx(0, 0, mt)])/mt;
	yn[index] =   
	    yn[idx(0, 0, mt)] + i2*(yn[idx(0, mt, mt)] - yn[idx(0, 0, mt)])/mt;
	zn[index] =   
	    zn[idx(0, 0, mt)] + i2*(zn[idx(0, mt, mt)] - zn[idx(0, 0, mt)])/mt;
	
	R = sqrt(pow(xn[index],2) +  
	         pow(yn[index],2) +
	         pow(zn[index],2) );
	xn[index] = 1/R * xn[index];
	yn[index] = 1/R * yn[index];
	zn[index] = 1/R * zn[index];
    }
    // i1,i2=mt
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	index=idx(0, mt, i1);
	xn[index] =   
	    xn[idx(0, mt, 0)] + i1*(xn[idx(0, mt, mt)] - xn[idx(0, mt, 0)])/mt;
	yn[index] =   
	    yn[idx(0, mt, 0)] + i1*(yn[idx(0, mt, mt)] - yn[idx(0, mt, 0)])/mt;
	zn[index] =   
	    zn[idx(0, mt, 0)] + i1*(zn[idx(0, mt, mt)] - zn[idx(0, mt, 0)])/mt;
	
	R = sqrt(pow(xn[index],2) +  
	         pow(yn[index],2) +
	         pow(zn[index],2) );
	xn[index] = 1/R * xn[index];
	yn[index] = 1/R * yn[index];
	zn[index] = 1/R * zn[index];
    }
   
    // Everywhere inbetween
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, i1);
	xn[index] =   
	    xn[idx(0, i1, 0)] + i2*(xn[idx(0, i1, mt)] - xn[idx(0, i1, 0)])/mt;  
	yn[index] =   
	    yn[idx(0, i1, 0)] + i2*(yn[idx(0, i1, mt)] - yn[idx(0, i1, 0)])/mt;
	zn[index] =   
	    zn[idx(0, i1, 0)] + i2*(zn[idx(0, i1, mt)] - zn[idx(0, i1, 0)])/mt;
	
	R = sqrt(pow(xn[index],2) +  
	         pow(yn[index],2) +
	         pow(zn[index],2) );
	xn[index] = 1/R * xn[index];
	yn[index] = 1/R * yn[index];
	zn[index] = 1/R * zn[index];
	}
    }
    
//    printf("cmb = %12.8g\n",cmb);
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
//	  xn[index] = x0 - (a-cmb)*ir/nr;
//	  yn[index] = y0 - (a-cmb)*ir/nr;
//	  zn[index] = z0 - (a-cmb)*ir/nr;
	  index++;
	}
      }
    }

    return 0;
}
