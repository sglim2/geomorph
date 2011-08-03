/*
 *
 * Data Manipulation
 *
 * 
 */

#include "main.h"
#include "data.h"

/*
 * utilities
 */

////////////////////////////////////////
// construction/destruction
////////////////////////////////////////


Data::Data()
  : infile(), outfile(), intype(), outtype(), nval(), x(), y(), z(), V()
{

  intype=UNDEF;
  outtype=UNDEF;
  
  indirSet = false;
  outdirSet = false;

}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Data::findBoundary
// An alternative (and faster) method would be to calculate these values as
// the data is read in. That way we share a loop over nval, and (for mitp data
// at least) is in the native format (i.e. no sqrt necessary).
bool Data::findBoundary()
{
  double tmpcmb=0.;
  double tmpa=0.;
  cmb = veryLarge;
  a   = verySmall;

  for (int in=0 ; in<nval ; in++){
    tmpcmb = sqrt(x[in]*x[in] + y[in]*y[in] + z[in]*z[in]);
    tmpa   = sqrt(x[in]*x[in] + y[in]*y[in] + z[in]*z[in]);
    tmpcmb < cmb ? cmb = tmpcmb : cmb=cmb;
    tmpa   > a   ? a   = tmpa   : a=a;
  }

  return 0;// success
}

////////////////////////////////////////
// Data::findLayers 
// As with Data::findBoundary(), an alernative method would be to calculate
// these values as the data is read in. Also, this method isn't very robust,
// and accurracy depends upon nbins.
bool Data::findLayers()
{
    nlayr = 0;
    int nbins=10000;
    int * bins = new int[nbins];
    //  double * rads = new double[nval];
    double rad;
    
    for ( int n=0 ; n<nbins ; n++) bins[n]=0;
    
    for (int in=0 ; in<nval ; in++){
	//    rads[in] = sqrt (x[in]*x[in] + y[in]*y[in] + z[in]*z[in] );
	rad      = sqrt (x[in]*x[in] + y[in]*y[in] + z[in]*z[in] );
	bins[ int(floor(nbins*(rad-cmb)/(a-cmb))) ]++;
    }
    
    
    for ( int n=0 ; n<nbins ; n++) {
//      printf("%d\t%d\n",n,bins[n]);
	if (bins[n] > 10){
	    nlayr++;
	}
    }
    
    return 0;
}


////////////////////////////////////////
// Read
bool Data::Read()
{
    switch (intype) {
	case UNDEF :
	    printf("Error.Intype undefined");
	    return 1; //fail
	case MVIS :
	    return mvisRead();
	case TERRA_CC :
	    // do something
	    break;
	case TERRA_CV :
	    // do something
	    break;
	case MITP :
	    return mitpRead();
	case FILT :
	    return filtRead();
	default :
	    printf("Error: intype undefined.");
	    return 1; //fail
    }

    return 0;    
}


////////////////////////////////////////
// mitpDepth2Radius
// --------
//
//
double Data::mitpDepth2Radius(double dpth)
{

 return (EarthRadKM - dpth) / EarthRadKM;
 
}


////////////////////////////////////////
// filtDepth2Radius
// --------
//
//
double Data::filtDepth2Radius(double dpth)
{

    return mitpDepth2Radius(dpth);
 
}

////////////////////////////////////////
// mitpRead
// --------
//   input data:
//     lat  ( degrees: -90  :  90 )
//     lng  ( degrees:   0  : 360 )
//     dpth ( km  (0=crust) )
//
//   Convert on input to:
//     lat  ( degrees:  -90 :  90 )
//     lng  ( degrees: -180 : 180 )
//     dpth ( EarthRadii  (0=centre) )
//
//   Convert to x,y,z (Earth Radius units) using formulae:
//    x = - cos(lat) * cos(lng)
//    y =   sin(lat) 
//    z =   cos(lat) * sin(lng)
//
//
// Input data is expected to be in the format:
//  line1: 
//    headers
//  subsequent lines: 
//    data in the order: 
//      Lat Long Depth Value
//  Sort order is:
//    increasing  Depth; 
//    increasing  Long;
//    then increasing  Lat.
bool Data::mitpRead()
{ 
  int ferror;

  double lat ;
  double lng ;
  double dpth;

  nval = 0;
  
  printf("Input File/base: %s\n",infile);
  FILE * fptr=fopen(infile,"r");
  if (fptr==NULL){
    return 1; //fail
  }
  
  char buf[255]; 
  while( fgets(buf,sizeof(buf),fptr) != NULL) {
      nval++;
  } 
  
  // discard the header record
  nval--;
  // define arrays
  x = new double[nval];
  y = new double[nval];
  z = new double[nval];
  V = new double[nval];
 
  // We have nval, now let's find nlat, nlng and ndpth..........
  // move to the first data-line
  fclose(fptr);
  fptr=fopen(infile,"r");
  if (fptr==NULL){
      return 1; //fail
  }
  // discard headers
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);

  lat = -veryLarge;
  lng = -veryLarge;
  dpth= -veryLarge;
  minR= +veryLarge;
  maxR= -veryLarge;

  // initialize limits
  nlat = 0;
  nlng = 0;
  ndpth = 0;

  for ( int i = 0 ; i < nval ; i++ ){
      ferror = fscanf(fptr,"%s", buf);
      if (atof(buf) > lat + quiteSmall ) {
	  nlat++;
	  lat = atof(buf);
      }
      
      ferror = fscanf(fptr,"%s", buf);
      if (atof(buf) > lng + quiteSmall ) {
	  nlng++;
	  lng = atof(buf);
      }

      ferror = fscanf(fptr,"%s", buf);
      if (atof(buf) > dpth + quiteSmall ) {
	  ndpth++;
	  dpth = atof(buf);
	  if ( mitpDepth2Radius(dpth) < minR ) minR = mitpDepth2Radius(dpth) ;
	  if ( mitpDepth2Radius(dpth) > maxR ) maxR = mitpDepth2Radius(dpth) ;
      }

      //discard last column
      ferror = fscanf(fptr,"%s", buf);
  }
  
  
  // define more arrays
  //lyrs = new double[ndpth];
  // collect layers radii
  //  .... to be finished

  //Collect Data..........

  // move to the first data-line
  fclose(fptr);
  fptr=fopen(infile,"r");
  if (fptr==NULL){
      return 1; //fail
  }
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);
  ferror = fscanf(fptr,"%s", buf);

  // Collect the rest of the data
  for ( int i = 0 ; i<nval ; i++ ){
      // Collect data
      ferror = fscanf(fptr,"%s", buf);  lat = atof (buf) * pi/180. ;
      ferror = fscanf(fptr,"%s", buf);  lng = atof (buf) * pi/180. ;
      ferror = fscanf(fptr,"%s", buf);  dpth = mitpDepth2Radius(atof(buf));
      ferror = fscanf(fptr,"%s", buf);  V[i] = atof (buf) ;

      // Convert to xyz
      // We may have issues here....
      // It's possible x,y,z values are swapped somewhere along the line.
      //attempt 1
//      x[i] = - 1. * cos(lat) * cos(lng) ;
//      y[i] =   1. * sin(lat) ;
//      z[i] =   1. * cos(lat) * sin(lng) ;
      // attempt 2
//      x[i] = - dpth * cos(lat) * cos(lng) ;
//      y[i] =   dpth * sin(lat) ;
//      z[i] =   dpth * cos(lat) * sin(lng) ;
      // attempt 3
      x[i] =   dpth * cos(lat) * cos(lng) ;
      z[i] =   dpth * sin(lat) ;
      y[i] =   dpth * cos(lat) * sin(lng) ;
  }

  fclose(fptr);
  return 0; //success
}


////////////////////////////////////////
// Dat::mvisRead()
// --------
//
bool Data::mvisRead()
{ 
  // this->mvis is already created, but we need to generate the domains
  mvis->genGrid(cmb);
  
  // for mvis input, our data isn't kept in Data->x,y,z,V, but
  // this->mvis->domains[].
  mvis->importMVIS(infile, cmb);
  
  int ferror;


  char buf[255]; 
  char tmpinfile[64];
  nval = 0;


  

  return 0; //success
}


////////////////////////////////////////
// filtRead
// --------
// Filt data consists of multiple files, each continaing a 
// single layer of data.
bool Data::filtRead()
{ 

  int ferror;

  double lat ;
  double lng ;
  double dpth;

  int file_nval=0;
  char buf[255]; 
  char tmpinfile[64];
  nval = 0;

  // find nval....
  for (int dpth=filtinstart; dpth<=filtinend ; dpth+=(filtinend-filtinstart)/(filtinnumfiles-1) ) {

      sprintf(tmpinfile,"%s.%dkm.xyz",infile,dpth);
//      printf("%s\n",tmpinfile);
 
      FILE * fptr=fopen(tmpinfile,"r");
      if (fptr==NULL){
	  printf("Cannot open file %s for reading\n",tmpinfile);
	  return 1; //fail
      }

      file_nval=0;
      while( fgets(buf,sizeof(buf),fptr) != NULL) {
	  file_nval++;
      } 

      // no header to discard
      // file_nval--;
      nval += file_nval;

      fclose(fptr);
  } // for dpth

  
  // define arrays based on nval
  x = new double[nval];
  y = new double[nval];
  z = new double[nval];
  V = new double[nval];

  // find minR, maxR, ndpth, nvalpershell
  minR = filtDepth2Radius(filtinend);
  maxR = filtDepth2Radius(filtinstart);
  ndpth = filtinnumfiles;
  nvalpershell = file_nval;  // should be the same across files.

  int nval_counter = 0;
  
  // Collect Data....
  for (int dpthkm=filtinstart; dpthkm<=filtinend ; dpthkm+=(filtinend-filtinstart)/(filtinnumfiles-1) ) {
      sprintf(tmpinfile,"%s.%dkm.xyz",infile,dpthkm);
 
      FILE * fptr=fopen(tmpinfile,"r");
      if (fptr==NULL){
	  return 1; //fail
      }
      
      // find file_nval again...
      file_nval=0;
      while( fgets(buf,sizeof(buf),fptr) != NULL) {
	  file_nval++;
      } 
      
      // move to beginning of file
      fclose(fptr);
      fptr=fopen(tmpinfile,"r");
      if (fptr==NULL){
	  return 1; //fail
      }
      
      lat = -veryLarge;
      lng = -veryLarge;
      dpth= -veryLarge;
      
      for ( int i = nval_counter ; i<nval_counter+file_nval ; i++ ){
	  // Collect data
          ferror = fscanf(fptr,"%s", buf);  lng = atof (buf) * pi/180. ;
          ferror = fscanf(fptr,"%s", buf);  lat = atof (buf) * pi/180. ;
	                           dpth = filtDepth2Radius(dpthkm);
          ferror = fscanf(fptr,"%s", buf);  V[i] = atof (buf) ;
	  
	  // Convert to xyz
	  x[i] =   dpth * cos(lat) * cos(lng) ;
	  z[i] =   dpth * sin(lat) ;
	  y[i] =   dpth * cos(lat) * sin(lng) ;
      }
      nval_counter+=file_nval;
      fclose(fptr);
  } // for dpth

  return 0; //success
}


////////////////////////////////////////
// intypeConverter
char* Data::intypeConverter()
{
    char * buf = new char[16];
    switch (intype) {
	    case UNDEF: strcpy(buf, "UNDEF");
		break;
	    case MVIS : strcpy(buf, "MVIS");
		break;
	    case TERRA_CC: strcpy(buf, "TERRA_CC");
		break;
	    case TERRA_CV: strcpy(buf, "TERRA_CV");
		break;
	    case MITP : strcpy(buf, "MITP");
		break;
	    case FILT : strcpy(buf, "FILT");
		break;
	    default: strcpy(buf, "UNDEF");
	}
    return buf;
}

////////////////////////////////////////
// outtypeConverter
char* Data::outtypeConverter()
{

    char * buf = new char[16];
    switch (this->outtype) {
	    case UNDEF: strcpy(buf, "UNDEF");
		break;
	    case MVIS : strcpy(buf, "MVIS");
		break;
	    case TERRA_CC: strcpy(buf, "TERRA_CC");
		break;
	    case TERRA_CV: strcpy(buf, "TERRA_CV");
		break;
	    case MITP : strcpy(buf, "MITP");
		break;
	    case FILT : strcpy(buf, "FILT");
		break;
	    default: strcpy(buf, "UNDEF");
	}
    return buf;
}

////////////////////////////////////////
// outtypeConverter
char* Data::interpConverter()
{

    char * buf = new char[16];
    switch (this->interp) {
	    case NEAREST : strcpy(buf, "NEAREST");
		break;
	    case NEAREST2 : strcpy(buf, "NEAREST2");
		break;
	    case LINEAR : strcpy(buf, "LINEAR");
		break;
	    case CUBIC : strcpy(buf, "CUBIC");
		break;
	    default: strcpy(buf, "Oops!");
	}
    return buf;
}

////////////////////////////////////////
// getStats
bool Data::getStats()
{
    double max = -veryLarge,
           min =  veryLarge;
	    
    Vmean = 0;

    for ( int i=0; i<nval ; i++){
	if ( V[i] > max )  max=V[i] ;
	if ( V[i] < min )  min=V[i] ;
	Vmean += V[i];
    }
    
    Vmax   = max;
    Vmin   = min;
    Vmean /= nval;

  // find a and cmb
  if (findBoundary()){
      printf("Error in Data::findBoudary().");
      return 1;
  }
  // find the number of layers
  if (findLayers()){
      printf("Error in Data::findLayers().");
      return 1;
  }
  

    printf("Input Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  nvals        =  %12ld                  |\n"       , nval);
    printf("|  nlayr(~)     =  %12d                  |\n"        , nlayr);
    printf("|  nlat         =  %12d                  |\n"        , nlat);
    printf("|  nlng         =  %12d                  |\n"        , nlng);
    printf("|  ndpth        =  %12d                  |\n"        , ndpth);
    printf("|  nvalpershell =  %12d                  |\n"        , nvalpershell);
    printf("|  minR         =  %12.8g                  |\n"      , minR);
    printf("|  maxR         =  %12.8g                  |\n"      , maxR);
    printf("|  V(mean)      =  %12.8g                  |\n"      , Vmean);
    printf("|  Vmin         =  %12.8g                  |\n"      , Vmin);
    printf("|  Vmax         =  %12.8g                  |\n"      , Vmax);
    printf("|  a            =  %12.8g                  |\n"      , a);
    printf("|  cmb          =  %12.8g                  |\n"      , cmb);
    printf("+----------------------------------------------+\n");
    
    
    return 0; //success
}
