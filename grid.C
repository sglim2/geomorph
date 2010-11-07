/*
 *
 * Data Manipulation
 *
 * 
 */


#include "grid.H"

/*
 * utilities
 */

////////////////////////////////////////
// construction/destruction
////////////////////////////////////////


Grid::Grid()
    : mt(), nt(), nd(), nr(), idmax(), nproc(), rmax(), rmin()
{
  idmax = 10;
}


////////////////////////////////////////
// Grid::Grid
// 
// Assume nr=mt/2
//  i.e # radial layers = nr+1
//
Grid::Grid(int _mt, int _nt, int _nd)
  : mt(), nt(), nd(), nr(), idmax(), nproc(), rmax(), rmin()
{
  idmax = 10;
  mt=_mt;
  nt=_nt;
  nd=_nd;
  
  nr=mt/2;
  nproc=pow((mt/nt),2)*10/nd;
  
}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Grid:xnProc
//
// given idx, returns TERRA processor number
//   
int Grid::xnProc(int idx)
{
  int xnproc=0;

  return idx;
}

////////////////////////////////////////
// Grid:xnProc
//
// given:
//   id.....diamond
//   xyz....x=0,y=1,z=2
//   i1.....i
//   j1.....j
//   r......radial layer  
// returns the TERRA processor number.
// Array xn has structure:
//    xn[nr+1][id][mt+1][mt+1][xyz]
//
int Grid::xnProc(int r, int id, int i1, int j1, int xyz)
{
  
  int xnproc=0;
  
  


  return xnproc;

}

////////////////////////////////////////
// Grid::idx
//
// given:
//   id.....diamond
//   xyz....x=0,y=1,z=2
//   i1.....i
//   j1.....j
//   r......radial layer  
// returns the index value for array xn.
// Array xn has structure:
//    xn[nr+1][id][mt+1][mt+1][xyz]
//
int Grid::idx(int r, int id, int i1, int j1, int xyz)
{
  int idmax  = 10;
  int nerror = 0;

  if (  r  > nr    || r     < 0 ) nerror=1; 
  if ( id  > idmax || idmax < 0 ) nerror=1;
  if ( i1  > mt+1  || i1    < 0 ) nerror=1;
  if ( j1  > mt+1  || j1    < 0 ) nerror=1;
  if ( xyz > 3     || xyz   < 0 ) nerror=1;

  int idx=0;
  
  int   rbase = r  * idmax*(mt+1)*(mt+1)*3;
  int  idbase = id * (mt+1)*(mt+1)*3;
  int   ibase = i1 * (mt+1)*3;
  int   jbase = j1 * 3;
  int xyzbase = xyz;
  
  idx = rbase + idbase + ibase + jbase + xyzbase;


  return idx;

}

////////////////////////////////////////
// Grid::gridfind
//
// A basic comparison of any given grid structire to the TERRA
// grid. Returns the mt value of the closest matching TERRA grid
// comparing the total number of grid points across the two grid
// structures.
//
int Grid::findGrid(int npts)
{
    int _npts=0;
    int _mt=1;      // eventually provides the closest matching mt value

    int _mt_old=0,_npts_old=0;

    while (_npts <= npts){
	_mt_old=_mt;
	_mt *= 2;
	_npts_old=_npts;
	_npts = (10*_mt*_mt + 2)*(_mt/2 + 1);
    }
    

    printf("_mt=%d\n",_mt);
    printf("_npts=%d\n",_npts);
    if (  npts - (_npts + _npts_old)/2 < 0  )  {
	_mt=_mt_old;
    }

    return _mt;
}

bool Grid::grdgen(double * xn, int mt)
{

  double fifthpi,w,phi,sinw,cosw,sgn;
  int    id,nn,index;

  fifthpi = 0.4*asin(1.);
  w       = 2.0*acos(1./(2.*sin(fifthpi)));
  cosw    = cos(w);
  sinw    = sin(w);
  nn      = pow(mt+1,2)*10;

  for ( id=0 ; id<10 ; id++ ){

    id < 5 ? sgn = 1. : sgn = -1.;
    
    phi = (2*((id+1)%5) + (id)/5)*fifthpi;

    index = idx(0, id, 0, 0, 0);
    
    xn[index] =  0.;
    index++;
    xn[index] =  0.;
    index++;
    xn[index] =  sgn;

    index = idx(0, id, 0, mt, 0);
    xn[index] =  sinw*cos(phi);
    index++;
    xn[index] =  sinw*sin(phi);
    index++;
    xn[index] =  cosw*sgn;

    index = idx(0, id, mt, 0, 0);
    xn[index] =  sinw*cos(phi + fifthpi + fifthpi);
    index++;
    xn[index] =  sinw*sin(phi + fifthpi + fifthpi);
    index++;
    xn[index] =  cosw*sgn;

    index = idx(0, id, mt, mt, 0);
    xn[index] =  sinw*cos(phi + fifthpi);
    index++;
    xn[index] =  sinw*sin(phi + fifthpi);
    index++;
    xn[index] = -cosw*sgn;
 
  }


}

bool Grid::genGrid()
{
  

    npts  = 10*(mt+1)*(mt+1)*(mt/2 + 1);
    nr    = mt/2+1;
    nproc = (mt/nt)*(mt/nt)*10/nd; 
    
    printf("TERRA-grid Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  mt        =  %12d                  |\n"        , mt);
    printf("|  nt        =  %12d                  |\n"        , nt);;
    printf("|  nd        =  %12d                  |\n"        , nd);
    printf("|  npts      =  %12d                  |\n"        , npts);
    printf("|  nr        =  %12d                  |\n"        , nr);
    printf("|  nproc     =  %12d                  |\n"        , nproc);
    printf("+----------------------------------------------+\n");
    
//    x = new double[npts];
//    y = new double[npts];
//    z = new double[npts];
//    V = new double[npts];

    // mt (and subsequently npts) for test purposes only
    mt = 1;
    npts  = (10*mt*mt + 2)*(mt/2 + 1);


    //    xn = new double[(mt+1)*(mt+1)*10*3];
    xn = new double[nr*10*(mt+1)*(mt+1)*3];

    //    printf("\nnproc_total = %d\n",nr*10*(mt+1)*(mt+1)*3);

    //    grdgen_(xn,&mt);    // fortran
    grdgen(&xn[0],mt);            // c++

    /*
    for ( int i=0 ; i<10 ; i++){
      printf("\nid %d\n",i);
      for ( int yi=0 ; yi<mt+1 ; yi++){
	for ( int xi=0 ; xi<mt+1 ; xi++){
	  printf("%12.6g\t",xn[0*10*(mt+1)*(mt+1) + i*(mt+1)*(mt+1) + yi*(mt+1) + xi]);
	}
      }
    }            
    printf("\n");

    printf("%d\n",2*10*(mt+1)*(mt+1) + 9*(mt+1)*(mt+1) + mt*(mt+1) + mt + 1);

    printf("%d\n",(mt+1)*(mt+1)*10*3);
    */
    
    int index=0;
    for ( int i=0 ; i<10 ; i++){
      for ( int xi=0 ; xi<mt+1 ; xi++){
	for ( int yi=0 ; yi<mt+1 ; yi++){
	  for ( int k=0 ; k<3 ; k++){
	    index = idx(0,i,xi,yi,k);
	    printf("%12.6g ",xn[index]);
	  }
	  printf("\n");
	}
      }      
      printf("\n");
    }

    
    for ( int i=0 ; i < ( 10*(mt+1)*(mt+1)*3 ) ; i++){
      printf("%d\t%12.8g\n",i,xn[i]);
    }
    

    return 0;

}
