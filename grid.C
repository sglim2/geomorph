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
    : domains(), mt(), nt(), nd(), nr(), idmax(), nproc(), rmax(), rmin()
{
  idmax = 10;
  domains = new Domain[idmax];
}


////////////////////////////////////////
// Grid::Grid
// 
// Assume nr=mt/2
//  i.e # radial layers = nr+1
//
Grid::Grid(int _mt, int _nt, int _nd)
  : domains(), mt(), nt(), nd(), nr(), idmax(), nproc(), rmax(), rmin()
{
  idmax = 10;
  mt=_mt;
  nt=_nt;
  nd=_nd;
  
  nr=mt/2;
  nproc=pow((mt/nt),2)*10/nd;
  
  domains = new Domain[idmax];
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
int Grid::xnProc(int r, int id, int j1, int i1, int xyz)
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
int Grid::idx(int r, int id, int i2, int i1, int xyz)
{
  int idmax  = 10;
  int nerror = 0;

  if (  r  > nr    || r     < 0 ) nerror=1; 
  if ( id  > idmax || idmax < 0 ) nerror=1;
  if ( i1  > mt    || i1    < 0 ) nerror=1;
  if ( i2  > mt    || i2    < 0 ) nerror=1;
  if ( xyz > 3     || xyz   < 0 ) nerror=1;

  int idx=0;
  
  int   rbase = r  * idmax*(mt+1)*(mt+1)*3;
  int  idbase = id * (mt+1)*(mt+1)*3;
  int  i2base = i2 * (mt+1)*3;
  int  i1base = i1 * 3;
  int xyzbase = xyz;
  
  idx = rbase + idbase + i2base + i1base + xyzbase;


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

////////////////////////////////////////
// Grid::genfind
//
bool Grid::genGrid()
{
  
    Domain * dptr;
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
    
    // mt (and subsequently npts) for test purposes only
//    mt = 8;
//    npts  = (10*mt*mt + 2)*(mt/2 + 1);
//    nr = 1;

    // make sure all domains are set up correctly
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	dptr->defineDomain(id,nr,mt);
	dptr++;
    }
    
    // generate grid points within each domain
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	dptr->grdgen(); 
	dptr++;
    }
      
    /*
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	printf("id = %d\n",id);
	for ( int i=0 ; i <  nr*(mt+1)*(mt+1)*3  ; i++){
	    printf("%d\t%12.8g\n",i,dptr->xn[i]);
	}
	dptr++;
    }
    */

    
    // print outer shell (i.e nr=0)
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	for ( int i=0 ; i < (mt+1)*(mt+1)*3 ; i+=3){
	    printf("%12.8g\t%12.8g\t%12.8g\n",dptr->xn[i+0],dptr->xn[i+1],dptr->xn[i+2]);
	}
	dptr++;
    }
    

    return 0;
}
