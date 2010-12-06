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
    mt = 0;
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
  // not yet implemented (or needed!)
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

  if (nerror!=0) printf("Grid::idx Error.");

  return idx;

}

////////////////////////////////////////
// Grid::importData
//
// imports Data::data into all domains.
//
bool Grid::importData(Data * dptr)
{
    int _idmax = 10;

#ifndef GEOMORPH_GPU
#pragma omp parallel for
    for (int i=0 ; i<_idmax ; i++){
	printf("...Domain %d\n",i);
	if (domains[i].importData(dptr)){
	    printf("Error in Domain::importData()");
	    //	    return 1; // fail
	}
    }
#pragma end parallel for
#else
    // if we're using the GPU we need to avoid moving in and out of the GPU
    // unit. For each domain, we move onto the GPU and stay there - for all
    // points in the domain.
    
    for (int i=0 ; i<_idmax ; i++){
	printf("...Domain %d\n",i);
	if (domains[i].importData_gpu(dptr)){
	    printf("Error in Domain::importData_gpu()");
	}
    }
#endif
    
    return 0; // success
}

////////////////////////////////////////
// Grid::exportData
//
// exports data in desired format.
//
// Due to the necessary sequential output, this must be done in serial for
// each process.
//
bool Grid::exportMVIS(Data * dptr)
{
    // assuming nd=10...for now

    FILE * fptr;
    char outfile[256];
    
    int  fail=0;

    // we now cycle over 'processors'
    for ( int proc=0 ; proc < nproc ; proc++ ){
		
	// open a file per 'process'
	sprintf(outfile,"%s.%04d.01", dptr->outfile, proc);
	printf("outfile = %s\n",outfile);
        
	// open file
	fptr=fopen(outfile,"a");
	if (fptr==NULL){
	}

	// call each domain which is part of our 'process'  (if nd=10, this means all of them)
	for (int i=0 ; i < idmax*nd/10 ; i++){
	    
	    // call our domain export routine
	    if ( domains[i].exportMVIS(fptr, nproc, proc, nt) ){
		printf("Error in Domain::exportMVIS()");
	    }
	}

	// close file, ready for re-assigning to a new 'process'
	fclose(fptr);	    
    }

    return fail; 
}

////////////////////////////////////////
// Grid::exportData
//
// exports data in desired format.
//
// Note, we cannot use a switch statement here since it defies the standard
// (and throws an error on pgi/gnu): 
//    the -> operator is not allowed in an integral constant expression.
//
bool Grid::exportGrid(Data * dptr)
{
    
    if ( dptr->outtype == dptr->MVIS ) {
	if (exportMVIS(dptr) ){
	    printf("Error in exportMVIS\n");
	    return 1; // fail
	}
    }else
	if ( dptr->outtype == dptr->TERRA ) {
//	    if (exportTERRA(dptr) ){
//		printf("Error in exportTERRA\n");
//		return 1; // fail
//	    }
	}else 
	    if ( dptr->outtype == dptr->MITP ) {
//		if (exportMITP(dptr) ){
//		    printf("Error in exportMITP\n");
//		    return 1; // fail
//		}
	    }else 
		if ( dptr->outtype == dptr->UNDEF ) {
		    return 1; // fail
		}else 
		    return 1; // fail

    return 0; // success
}

////////////////////////////////////////
// Grid::suggestGrid
//
// A basic comparison of any given grid structire to the TERRA
// grid. Returns the mt value of the closest matching TERRA grid
// comparing the total number of grid points across the two grid
// structures.
//
int Grid::suggestGrid(int npts)
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
// Grid::genGrid
//
bool Grid::genGrid(double _cmb)
{
    int fail = 0;
    
    Domain * dptr;
    npts  = 10*(mt+1)*(mt+1)*(mt/2 + 1);
    nr    = mt/2+1;
    nproc = (mt/nt)*(mt/nt)*10/nd; 
    
    printf("TERRA-grid Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  mt        =  %12d                   |\n"        , mt);
    printf("|  nt        =  %12d                   |\n"        , nt);;
    printf("|  nd        =  %12d                   |\n"        , nd);
    printf("|  npts      =  %12d                   |\n"        , npts);
    printf("|  nr        =  %12d                   |\n"        , nr);
    printf("|  nproc     =  %12d                   |\n"        , nproc);
    printf("+----------------------------------------------+\n");
    
    // make sure all domains are set up correctly
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	if (dptr->defineDomain(id,nr,mt)){
	    printf("Error in Domain::defineDomain()");
//	    fail = 1; // commented to enable auto-parallelization on PGI compiler
	}
	dptr++;
    }
    
    // generate grid points within each domain
    dptr = domains;
    for ( int id = 0 ; id < idmax ; id++ ){
	if (dptr->grdgen2(_cmb)){
	    printf("Error in Domain::grdgen");
//	    fail = 1; // commented to enable auto-parallelization on PGI compiler
	}
	dptr++;
    }
      
    return fail; 
}
