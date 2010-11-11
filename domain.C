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

    xn = new double[nr*(mt+1)*(mt+1)*3];
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
int Domain::defineDomain(int _id, int _nr, int _mt)
{
    
    id = _id;
    nr = _nr;
    mt = _mt;
    
    id < 5 ? northern = true : northern = false;

    xn = new double[nr*(mt+1)*(mt+1)*3];
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
int Domain::idx(int r, int i2, int i1, int xyz)
{
  int nerror = 0;

  if (  r  > nr    || r     < 0 ) nerror=1; 
  if ( i1  > mt    || i1    < 0 ) nerror=1;
  if ( i2  > mt    || i2    < 0 ) nerror=1;
  if ( xyz > 3     || xyz   < 0 ) nerror=1;

  int idx=0;
  
  int   rbase = r  * (mt+1)*(mt+1)*3;
  int  i2base = i2 * (mt+1)*3;
  int  i1base = i1 * 3;
  int xyzbase = xyz;
  
  idx = rbase + i2base + i1base + xyzbase;

  return idx;
}

////////////////////////////////////////
// Domain::grdgen
//
//
bool Domain::grdgen()
{
    
    int    nn,index;
    double R,x0,y0,z0;

    double a,tau,rho,u,v,Beta;
    double Ry[3][3], A[12][3], Ad[12][3];

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
	index = idx(0, 0, 0, 0);
	xn[index] = Ad[0][0]; index++;
	xn[index] = Ad[0][1]; index++;
	xn[index] = Ad[0][2];
	// mt,0
	index = idx(0, mt, 0, 0);
	xn[index] = Ad[id+1][0]; index++;
	xn[index] = Ad[id+1][1]; index++;
	xn[index] = Ad[id+1][2]; 
	// 0,mt
	index = idx(0, 0, mt, 0);
	if (id == 0) {
	    xn[index] = Ad[id+5][0]; index++;
	    xn[index] = Ad[id+5][1]; index++;
	    xn[index] = Ad[id+5][2]; 
	}else{
	    xn[index] = Ad[id][0]; index++;
	    xn[index] = Ad[id][1]; index++;
	    xn[index] = Ad[id][2]; 
	}
	// mt,mt
	index = idx(0, mt, mt, 0);
	xn[index] = Ad[id+6][0]; index++;
	xn[index] = Ad[id+6][1]; index++;
	xn[index] = Ad[id+6][2]; 
	
	// Southern Hemisphere
    }else{
	// South Pole
	index = idx(0, 0, 0, 0);
	xn[index] = Ad[11][0]; index++;
	xn[index] = Ad[11][1]; index++;
	xn[index] = Ad[11][2];
	// mt,0
	index = idx(0, mt, 0, 0);
	if (id == 9) {
	    xn[index] = Ad[id-3][0]; index++;
	    xn[index] = Ad[id-3][1]; index++;
	    xn[index] = Ad[id-3][2]; 
	}else{
	    xn[index] = Ad[id+2][0]; index++;
	    xn[index] = Ad[id+2][1]; index++;
	    xn[index] = Ad[id+2][2]; 
	}
	// 0,mt
	index = idx(0, 0, mt, 0);
	xn[index] = Ad[id+1][0]; index++;
	xn[index] = Ad[id+1][1]; index++;
	xn[index] = Ad[id+1][2]; 
	// mt,mt
	index = idx(0, mt, mt, 0);
	xn[index] = Ad[id-4][0]; index++;
	xn[index] = Ad[id-4][1]; index++;
	xn[index] = Ad[id-4][2]; 
    }
    
    // Edge Points......... (assumes unit radius)
    // i2,i1=0
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, 0, 0);
	xn[index+0] =   
	    xn[idx(0, 0, 0, 0)] + i2*(xn[idx(0, mt, 0, 0)] - xn[idx(0, 0, 0, 0)])/mt;
	xn[index+1] =   
	    xn[idx(0, 0, 0, 1)] + i2*(xn[idx(0, mt, 0, 1)] - xn[idx(0, 0, 0, 1)])/mt;
	xn[index+2] =   
	    xn[idx(0, 0, 0, 2)] + i2*(xn[idx(0, mt, 0, 2)] - xn[idx(0, 0, 0, 2)])/mt;
	
	R = sqrt(pow(xn[index+0],2) +  
	         pow(xn[index+1],2) +
	         pow(xn[index+2],2) );
	xn[index+0] = 1/R * xn[index+0];
	xn[index+1] = 1/R * xn[index+1];
	xn[index+2] = 1/R * xn[index+2];
    }
    // i1,i2=0
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	index=idx(0, 0, i1, 0);
	xn[index+0] =   
	    xn[idx(0, 0, 0, 0)] + i1*(xn[idx(0, 0, mt, 0)] - xn[idx(0, 0, 0, 0)])/mt;
	xn[index+1] =   
	    xn[idx(0, 0, 0, 1)] + i1*(xn[idx(0, 0, mt, 1)] - xn[idx(0, 0, 0, 1)])/mt;
	xn[index+2] =   
	    xn[idx(0, 0, 0, 2)] + i1*(xn[idx(0, 0, mt, 2)] - xn[idx(0, 0, 0, 2)])/mt;
	
	R = sqrt(pow(xn[index+0],2) +  
	         pow(xn[index+1],2) +
	         pow(xn[index+2],2) );
	xn[index+0] = 1/R * xn[index+0];
	xn[index+1] = 1/R * xn[index+1];
	xn[index+2] = 1/R * xn[index+2];
    }
    // i2,i1=mt
    for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, mt, 0);
	xn[index+0] =   
	    xn[idx(0, 0, mt, 0)] + i2*(xn[idx(0, mt, mt, 0)] - xn[idx(0, 0, mt, 0)])/mt;
	xn[index+1] =   
	    xn[idx(0, 0, mt, 1)] + i2*(xn[idx(0, mt, mt, 1)] - xn[idx(0, 0, mt, 1)])/mt;
	xn[index+2] =   
	    xn[idx(0, 0, mt, 2)] + i2*(xn[idx(0, mt, mt, 2)] - xn[idx(0, 0, mt, 2)])/mt;
	
	R = sqrt(pow(xn[index+0],2) +  
	         pow(xn[index+1],2) +
	         pow(xn[index+2],2) );
	xn[index+0] = 1/R * xn[index+0];
	xn[index+1] = 1/R * xn[index+1];
	xn[index+2] = 1/R * xn[index+2];
    }
    // i1,i2=mt
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	index=idx(0, mt, i1, 0);
	xn[index+0] =   
	    xn[idx(0, mt, 0, 0)] + i1*(xn[idx(0, mt, mt, 0)] - xn[idx(0, mt, 0, 0)])/mt;
	xn[index+1] =   
	    xn[idx(0, mt, 0, 1)] + i1*(xn[idx(0, mt, mt, 1)] - xn[idx(0, mt, 0, 1)])/mt;
	xn[index+2] =   
	    xn[idx(0, mt, 0, 2)] + i1*(xn[idx(0, mt, mt, 2)] - xn[idx(0, mt, 0, 2)])/mt;
	
	R = sqrt(pow(xn[index+0],2) +  
	         pow(xn[index+1],2) +
	         pow(xn[index+2],2) );
	xn[index+0] = 1/R * xn[index+0];
	xn[index+1] = 1/R * xn[index+1];
	xn[index+2] = 1/R * xn[index+2];
    }
   
    // Everywhere inbetween
    for ( int i1 = 1 ; i1 < mt ; i1++ ){
	for ( int i2 = 1 ; i2 < mt ; i2++ ){
	index=idx(0, i2, i1, 0);
	xn[index+0] =   
	    xn[idx(0, i1, 0, 0)] + i2*(xn[idx(0, i1, mt, 0)] - xn[idx(0, i1, 0, 0)])/mt;  
	xn[index+1] =   
	    xn[idx(0, i1, 0, 1)] + i2*(xn[idx(0, i1, mt, 1)] - xn[idx(0, i1, 0, 1)])/mt;
	xn[index+2] =   
	    xn[idx(0, i1, 0, 2)] + i2*(xn[idx(0, i1, mt, 2)] - xn[idx(0, i1, 0, 2)])/mt;
	
	R = sqrt(pow(xn[index+0],2) +  
	         pow(xn[index+1],2) +
	         pow(xn[index+2],2) );
	xn[index+0] = 1/R * xn[index+0];
	xn[index+1] = 1/R * xn[index+1];
	xn[index+2] = 1/R * xn[index+2];
	}
    }
    
    // genrerate radial points
    
    for ( int ir = 1 ; ir < nr ; ir++){
      index=idx(ir, 0, 0, 0);
      for ( int i2 = 0 ; i2 < mt+1 ; i2++ ){
	for ( int i1 = 0 ; i1 < mt+1 ; i1++ ){
	  x0 = xn[idx(0,i2,i1,0)];
	  y0 = xn[idx(0,i2,i1,1)];
	  z0 = xn[idx(0,i2,i1,2)];
	  
	  xn[index] = x0 - x0*(a-cmb)*ir/nr;
	  index++;
	  xn[index] = y0 - y0*(a-cmb)*ir/nr;
	  index++;
	  xn[index] = z0 - z0*(a-cmb)*ir/nr;
	  index++;

	}
      }
    }

}
