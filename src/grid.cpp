/**
 \file
 
  \brief
  Grid class constructors/destructors, methods and operators. The geomorph native grid structure. 
 */

#include "grid.h"
#include "main.h"


////////////////////////////////////////
// construction/destruction
////////////////////////////////////////

Grid::Grid()
/** 
   Construction of class Grid, with no mt defined.
 */
  : suffix(), suffixSet(), mt(), nt(), nd(), mtSet(), ntSet(), ndSet(), nproc(),  nr(), idmax(), npts(), domains()
{
  mt = 0;
  mtSet=false;
  ntSet=false;
  ndSet=false;
  suffixSet=false;
  idmax = 10;
  domains = new Domain[idmax];
  for (int i=0 ; i<idmax ; i++) domains[i].id=i;
}


////////////////////////////////////////
// Grid::Grid
Grid::Grid(int _mt, int _nt, int _nd)
/** 
   Construction of class Grid, with mt,nt,nd defined.
   We assume nr=mt/2, giving the number of radial layers as nr+1. This is standard in TERRA models.
*/
  : suffix(), suffixSet(), mt(), nt(), nd(), mtSet(), ntSet(), ndSet(), nproc(),  nr(), idmax(), npts(), domains()
{
  idmax = 10;
  suffixSet=false;
  mt=_mt;
  nt=_nt;
  nd=_nd;
  mtSet=true;
  ntSet=true;
  ndSet=true;
  
  nr=mt/2;
  nproc=pow((mt/nt),2)*10/nd;
  
  domains = new Domain[idmax];
  for (int i=0 ; i<idmax ; i++) domains[i].id=i;
}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Grid::idx
//
/**
  Each domain has one-dimensiional arrays which represent points in
  3-dimensional space. This routine takes the 3-dimensional space parameters
  and converts them into an index of the one-dimensional array.

  \param r the radial layer index
  \param id the domain index
  \param i2 the first shell coordinate
  \param i1 the second shell coordinate
  \param xyz we are asking for either the x,y, or z index (1,2, or 3 respectively)
  \return idx the one-demensional index representing the 3-d space.
 */
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
/**
 Imports Data::data into all domains.
 \param dptr A pointer to the input data

 return 0 for success, 1 for failure.
*/
bool Grid::importData(Data * dptr)
{
    int _idmax = 10;
    bool fail=false;

#ifndef GEOMORPH_GPU
#pragma omp parallel for
    for (int i=0 ; i<_idmax ; i++){
	printf("...Domain %d\n",i);
	if (domains[i].importData(dptr)){
	    printf("Error in Domain::importData()");
            //      fail=true;
            //      return fail;
	}
    }

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
    
    return fail; // success if false
}

////////////////////////////////////////
// Grid::importMVIS
/**
 Imports data from MVIS files.
 \param dinfile Basename of the mvis input file. Further parameters will be suffixed on to this file to make up the fill filename.
 \param cmb The Core-Mantle boundary in Earth-Radius units.
 \return fail Boolean return value, false for success, true for failure.
*/
bool Grid::importMVIS(char * dinfile, double cmb)
{
    FILE * fptr;
    char infile[256];
    
    bool  fail=false;

    // we now cycle over 'processors' (or files).
    for ( int proc=0 ; proc < nproc ; proc++ ){

	// open a file per 'process'
	strcpy(infile,dinfile);
        sprintf(infile,"%s.%04d.%02d", infile, proc, suffix);
	printf("infile = %s\n",infile);
        
	// open file
	fptr=fopen(infile,"r");
	if (fptr==NULL){
	}

	// call each domain which is part of our 'process'  
        // (if nd=10, this means all of them; nd=5, 1 hemisphere only)
	for (int i=0 ; i < nd ; i++){
	    
	  // call our domain export routine
	  if ( nd == 10 ){
	      if ( domains[i].importMVIS(fptr, proc, nt, cmb) ){
		  printf("Error in Domain::importMVIS()");
	      }
	  } // if nd == 10
	  
	  if ( nd == 5 ) {
	      if ( proc < nproc/2 ){
		  if ( domains[i].importMVIS(fptr, proc, nt, cmb) ){
		      printf("Error in Domain::importMVIS()");
		  }
	      }else{
		  if ( domains[i+nd].importMVIS(fptr, proc-nproc/2, nt, cmb) ){
		      printf("Error in Domain::importMVIS()");
		  }
	      }
	  } // if nd == 5
	  
	} // for i
	// close file, ready for re-assigning to a new 'process'
	fclose(fptr);	    
    }
    
    return fail; 
}

////////////////////////////////////////
// Grid::exportMVIS
/**
 Exports data in the desired format. Due to the necessary sequential output, this must be done in serial for
 each process.
 \param dptr A pointer to the input data.
 \return fail Boolean return value, false for success, true for failure.
*/
bool Grid::exportMVIS(Data * dptr)
{
    FILE * fptr;
    char outfile[256];
    
    int  fail=0;

    // we now cycle over 'processors'
    for ( int proc=0 ; proc < nproc ; proc++ ){
		
	// open a file per 'process'
        sprintf(outfile,"%s.%04d.%02d", dptr->outfile, proc,suffix);
	printf("outfile = %s\n",outfile);
        
	// open file
	fptr=fopen(outfile,"a");
	if (fptr==NULL){
	}

	// Call each domain which is part of our 'process'.
	// (if nd=10, this means all of them; nd=5, 1 hemisphere only)
	for (int i=0 ; i < nd ; i++){
	    
	  // call our domain export routine
	  if ( nd == 10 ){
	    if ( domains[i].exportMVIS(fptr, proc, nt) ){
		printf("Error in Domain::exportMVIS()");
	    }
	  } // if nd == 10

	  if ( nd == 5 ) {
	    if ( proc < nproc/2 ){
	      if ( domains[i].exportMVIS(fptr, proc, nt) ){
		printf("Error in Domain::exportMVIS()");
	      }
	    }else{
	      if ( domains[i+nd].exportMVIS(fptr, proc-nproc/2, nt) ){
		printf("Error in Domain::exportMVIS()");
	      }
	    }
	  } // if nd == 5

	} // for i
	// close file, ready for re-assigning to a new 'process'
	fclose(fptr);	    
    }

    return fail; 
}

////////////////////////////////////////
// Grid::importTERRA
/**
 Imports data in desired format (TERRA Convection model, or Circulation
 model).  Due to the necessary sequential input, this must be done in serial
 for each process.
 \param dinfile Basename of the TERRA input file. Further parameters will be suffixed on to this file to make up the fill filename.
 \param terratype Distinction between TERRA types - Convection (CV, terratype=0, default) or Circulation (CC, terratype=1) model. 
 \return fail Boolean return value, false for success, true for failure.
*/
bool Grid::importTERRA(char * dinfile, int terratype)
{
    FILE * fptr;
    char infile[256];
    int  tvppmax;   // temp/vel/pressure/plate-history
    int  fail=0;

    // terratype=0(default) - > convection model
    // terratype=1 - > circulation model (i.e. with plate-histories)
    terratype == 1 ? tvppmax=4 : tvppmax=3;

    // we now cycle over 'processors'
    for ( int proc=0 ; proc < nproc ; proc++ ){
      
      // open a file per 'process'
      strcpy(infile,dinfile);
      sprintf(infile,"%s.%04d.%02d", infile, proc,suffix);
      printf("infile = %s\n",infile);
      
      // open file
      fptr=fopen(infile,"r");
      if (fptr==NULL){
	  printf("Error opening file %s\n",infile);
      }
      
      for ( int tvpp=0 ; tvpp<tvppmax ; tvpp++ ){ // temp/vel/pressure/plate-history output
	
	// read initial blurb and throw away
	char * tmpbuf = new char[256];
	float buf;
	for ( int i=0 ; i<5 ; i++){
	  fgets(tmpbuf , 256 , fptr);
	}
	
	// read radii of layers and throw away
	for (int i=1 ; i <= nr ; i++){
	  fscanf(fptr,"%f", &buf );
	}
	    
	// read propr array and throw away
	for (int i=1 ; i <= 20 ; i++){
	  fscanf(fptr,"%f", &buf );
	}
	
	if (tvpp != 3) { // i.e. not plate-histories
	  
	  // cycle over layers
//        for ( int ir=nr-1 ; ir>=0 ; ir-- ){
	  for ( int ir=0 ; ir<nr ; ir++ ){
	    
	    // call each domain which is part of our 'process'  (if nd=10, this means all of them; nd=5, 1 hemisphere only)
	    for (int i=0 ; i < nd ; i++){
	      
	      // call our domain export routine
	      if ( nd == 10 ){
		if ( domains[i].importTERRA(fptr, proc, nt, ir, tvpp) ){
		  printf("Error in Domain::importTERRA()");
		}
	      } // if nd == 10
	      
	      if ( nd == 5 ) {
		if ( proc < nproc/2 ){
		  if ( domains[i].importTERRA(fptr, proc, nt, ir, tvpp) ){
		    printf("Error in Domain::importTERRA()");
		  }
		}else{
		  if ( domains[i+nd].importTERRA(fptr, proc-nproc/2, nt, ir, tvpp) ){
		    printf("Error in Domain::importTERRA()");
		  }
		}
	      } // if nd == 5
	      
	    } // for i (domain)
	  } // ir
	}else{ 
	  
	  // read plate histories and throw away
	  for (int i=1 ; i <= (nt+1)*(nt+1)*nd*1*2 ; i++){
	    fscanf(fptr,"%f", &buf );
	  } // for i
	}
      } // tvpp
      
      // close file, ready for re-assigning to a new 'process'
      fclose(fptr);	    
      
    } // proc
    
    return fail; 
} // importTERRA


////////////////////////////////////////
// Grid::exportTERRA
/**
 Exports data in desired format (TERRA Convection model, or Circulation model).  Due to the necessary sequential output, this must be done in serial for each process.
 \param dptr A pointer to the input data.
 \param terratype Distinction between TERRA types - Convection (CV, terratype=0, default) or Circulation (CC, terratype=1) model. 
 \return fail Boolean return value, false for success, true for failure.
*/
bool Grid::exportTERRA(Data * dptr, int terratype)
{
    FILE * fptr;
    char outfile[256];
    int tvppmax;   // temp/vel/pressure/plate-history
    int  fail=0;

    // terratype=0(default) - > convection model
    // terratype=1 - > circulation model (i.e. with plate-histories)
    terratype == 1 ? tvppmax=4 : tvppmax=3;

    // we now cycle over 'processors'
    for ( int proc=0 ; proc < nproc ; proc++ ){
		
	// open a file per 'process'
        sprintf(outfile,"%s.%04d.%02d", dptr->outfile, proc,suffix);
	printf("outfile = %s\n",outfile);
        
	// open file
	fptr=fopen(outfile,"a");
	if (fptr==NULL){
	}
	
	// get time
	time_t rawtime;
	time ( &rawtime );
	
	for ( int tvpp=0 ; tvpp<tvppmax ; tvpp++ ){ // temp/vel/pressure/plate-history output
	    
	    // write initial blurb
	    fprintf(fptr,"%5d%5d\n",nr-1,nt);
	    fprintf(fptr,"CASE 001, in-type = %s\n",dptr->intypeConverter());
	    fprintf(fptr,"GEOMORPH-GENERATED CONVERSION\n");
	    fprintf(fptr,"-\n");
	    fprintf(fptr,"%s",ctime(&rawtime));  // new-line taken care of with ctime()
	    
	    // write radii of layers
	    for (int i=1 ; i <= nr ; i++){
		fprintf(fptr,"%15.8E", EarthRadKM*1000*(dptr->maxR-i*(dptr->maxR-dptr->minR)/(nr-1)) );
		if ( i%10 == 0 || i==nr ) fprintf(fptr,"\n"); // print in columns of 10
	    }
	    
	    // write propr array
	    for (int i=1 ; i <= 20 ; i++){
		if ( i==5 ){
                  fprintf(fptr,"%15.8E", 0.001);
                }else{
                  fprintf(fptr,"%15.8E", 0. );
                }
		if ( i%10 == 0 || i==20 ) fprintf(fptr,"\n"); // print in columns of 10
	    }


	    if (tvpp != 3) { // i.e. no plate-histories

		// cycle over layers
		long int colcntr=1;
//	        for ( int ir=nr-1 ; ir>=0 ; ir-- ){
		for ( int ir=0 ; ir<nr ; ir++ ){
		    
		    // call each domain which is part of our 'process'  (if nd=10, this means all of them; nd=5, 1 hemisphere only)
		    for (int i=0 ; i < nd ; i++){
			
			// call our domain export routine
			if ( nd == 10 ){
			    if ( domains[i].exportTERRA(fptr, proc, nt, ir, tvpp, colcntr) ){
				printf("Error in Domain::exportTERRA()");
			    }
			} // if nd == 10
			
			if ( nd == 5 ) {
			    if ( proc < nproc/2 ){
				if ( domains[i].exportTERRA(fptr, proc, nt, ir, tvpp, colcntr) ){
				    printf("Error in Domain::exportTERRA()");
				}
			    }else{
				if ( domains[i+nd].exportTERRA(fptr, proc-nproc/2, nt, ir, tvpp, colcntr) ){
				    printf("Error in Domain::exportTERRA()");
				}
			    }
			} // if nd == 5
			
		    } // for i (domain)
		} // ir
	    }else{ // tvpp != 3
		long int colcntr2=1;
		for (int i=1 ; i <= (nt+1)*(nt+1)*nd*1*2 ; i++){
		    fprintf(fptr,"%10.3E", 0. );
		    if ( colcntr2%15 == 0 ) fprintf(fptr,"\n"); // print in columns of 15
		    colcntr2++;
		} // for i
	    }
	    if (tvpp!=1) fprintf(fptr,"\n");
	} // tvpp
	
	// close file, ready for re-assigning to a new 'process'
	fclose(fptr);	    

    } // proc

    return fail; 
} // exportTERRA

////////////////////////////////////////
// Grid::exportGrid
/**
 Exports data in the desired format. 
 \param dptr A pointer to the input data.
 \return fail Boolean return value, false for success, true for failure.
*/
bool Grid::exportGrid(Data * dptr)
{

  bool fail=false;
  /* Note, we cannot use a switch statement here since it defies the standard
    (and throws an error on pgi/gnu):
    the -> operator is not allowed in an integral constant expression.
  */
    if ( dptr->outtype == dptr->MVIS ) {
	if (exportMVIS(dptr) ){
	    printf("Error in exportMVIS\n");
	    fail=true; // fail
	    return fail;
	}
    }else
	if ( dptr->outtype == dptr->TERRA_CV ) {
	    if (exportTERRA(dptr,0) ){
		printf("Error in exportTERRA()\n");
		fail=true; // fail
		return fail;
	    }
	}else
	    if ( dptr->outtype == dptr->TERRA_CC ) {
		if (exportTERRA(dptr,1) ){
		    printf("Error in exportTERRA()\n");
		    fail=true; // fail
		    return fail;
		}
	    }else 
		if ( dptr->outtype == dptr->MITP ) {
//		if (exportMITP(dptr) ){
//		    printf("Error in exportMITP\n");
//                  fail=true; // fail
//		    return fail;
//		}
		}else 
		    if ( dptr->outtype == dptr->UNDEF ) {
		      fail=true; // fail
		      return fail;
		    }else 
		fail=true; // fail
		return fail;
    
    return fail; // success it false
}

////////////////////////////////////////
// Grid::suggestGrid
/**
 A basic comparison of any given grid structure to the TERRA grid. Returns the
 mt value of the closest matching TERRA grid comparing the total number of
 grid points across the two grid structures.
 \param npts Total number of points of the comparison grid-structure.
 \return _mt The closest matching TERRA-grid in terms of mt.
*/
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
    
//    printf("_mt=%d\n",_mt);
//    printf("_npts=%d\n",_npts);
    if (  npts - (_npts + _npts_old)/2 < 0  )  {
	_mt=_mt_old;
    }

    return _mt;
}

////////////////////////////////////////
// Grid::genGrid
/**
 Defines each of the 10 domains of the TERRA grid-format based on pre-loaded mt,nt, and nd values.
 \param _cmb Core-Mantle boundary in Earth-Radius units.
 \return fail Boolean return value, false for success, true for failure.
 */
bool Grid::genGrid(double _cmb)
{
    int fail = 0;
    
    Domain * dptr;
    npts  = 10*(mt+1)*(mt+1)*(mt/2 + 1);
    nr    = mt/2+1;
    nproc = (mt/nt)*(mt/nt)*10/nd; 
    
    printf("Grid Statistics....\n");
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
