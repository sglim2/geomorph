/**
 \file
 
  \brief
  Data class constructors/destructors, methods and operators. A flexible input grid structure. 
 */

#include "main.h"
#include "data.h"


////////////////////////////////////////
// construction/destruction
////////////////////////////////////////

/**
 Construction of class Data
 
 The Data class will hold the input data of any input type and move it 
 into a geomorph grid within the class. 
 */
Data::Data()
  : infile(), outfile(), intype(), outtype(), nval(), x(), y(), z(), V(), intypeSet(), outtypeSet(), interpSet()
{

  intype=UNDEF;
  outtype=UNDEF;
  
  indirSet = false;
  outdirSet = false;

  // FILT data...
  filtinstart    = 0;
  filtinend      = 0;
  filtinnumfiles = 0;
  nvalpershell   = 0;
  filtinstartSet = false;
  filtinendSet = false;
  filtinnumfilesSet = false;

  // MITP data...
  nlat = 0;
  nlng = 0;

  // GUI 
  infileSet = false;
  outfileSet = false;
  intypeSet = false;
  outtypeSet = false;
  interpSet = false;

  cmbinSet = false;
}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Data::findBoundary
/** 
 Finds the rdial boundaries, cmb and a.
 An alternative (and faster) method would be to calculate these values as
 the data is read in. That way we share a loop over nval, and (for mitp data
 at least) is in the native format (i.e. no sqrt necessary).
*/
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
/**
 As with Data::findBoundary(), an alernative method would be to calculate
 these values as the data is read in. Also, this method isn't very robust,
 and accurracy depends upon nbins.

 THIS ROUTINE IS INACCURATE AND SHOULD BE AVOIDED!
*/
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
// Data::Read
/**
  Function controlling the read functions, depending on input type
 */
bool Data::Read()
{
    switch (intype) {
	case UNDEF :
	    printf("Error.Intype undefined");
	    return 1; //fail
	case MVIS :
	    return mvisRead();
            break;
	case TERRA_CC :
	    return terraRead();
	    break;
	case TERRA_CV :
	    return terraRead();
	    break;
	case MITP :
	    return mitpRead();
	    break;
	case FILT :
	    return filtRead();
	    break;
	case GYPSUMP :
	    return gypsumRead();
	    break;
	case GYPSUMS :
	    return gypsumRead();
	    break;
	default :
	    printf("Error: intype undefined.");
	    return 1; //fail
    }
    return 0;    
}


////////////////////////////////////////
// Data::mitpDepth2Radius
/**
 Given a depth, returns a radius for MITP data
 \param dpth The mitp depth
 */
double Data::mitpDepth2Radius(double dpth)
{
 return (EarthRadKM - dpth) / EarthRadKM;
}

////////////////////////////////////////
// Data::filtDepth2Radius
/**
 Given a depth, returns a radius for FILT data
 \param dpth The mitp depth
 */
double Data::filtDepth2Radius(double dpth)
{
    return mitpDepth2Radius(dpth);
}

////////////////////////////////////////
// gypsumDepth2Radius
/**
 Given a depth, returns a radius for GYPSUM data
 \param dpth The mitp depth
 */
double Data::gypsumDepth2Radius(double dpth)
{
    return mitpDepth2Radius(dpth);
}

////////////////////////////////////////
// mitpRead
/**
  Reads MITP data.
  
  The MITP data has lat values in the range ( degrees: -90  :  90 ),
  lng values in the range ( degrees:   0  : 360 )
  and dpth as km  (0=crust)

  Convert on input to lat  ( degrees:  -90 :  90 ),
  lng  ( degrees: -180 : 180 ), and 
  dpth ( EarthRadii  (0=centre) )

  Convert to x,y,z (Earth Radius units) using formulae:
    x = - cos(lat) * cos(lng)
    y =   sin(lat) 
    z =   cos(lat) * sin(lng)

 Input data is expected to be in the format:
  line1: 
    headers
  subsequent lines: 
    data in the order: 
      Lat Long Depth Value
  Sort order is:
    increasing  Depth; 
    increasing  Long;
    then increasing  Lat.
*/
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
/**
  Reads in MVIS data
 */
bool Data::mvisRead()
{ 
  // this->mvis is already created, but we need to generate the domains
  mvis->genGrid(cmb);
  
  // for mvis input, our data isn't kept in Data->x,y,z,V, but
  // this->mvis->domains[].
  mvis->importMVIS(infile, cmb);
  
  return 0; //success
}

////////////////////////////////////////
// Data::terraRead()
/**
 Reads in TERRA data
 */
bool Data::terraRead()
{ 
  // this->mvis is already created, but we need to generate the domains
  mvis->genGrid(cmb);
  
  // for mvis input, our data isn't kept in Data->x,y,z,V, but
  // this->mvis->domains[].
  mvis->importTERRA(infile, 0); // always terra_cc for now.
  
  return 0; //success
}

////////////////////////////////////////
// Data::filtRead
/**
 Reads in FLIT data.
 
 Filt data consists of multiple files, each containing a single layer of data.
*/
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
// Data::gypsumRead
/**
 Reagds in GYPSUM data.

 GYPSUM data consists of multiple files, each contianing a 
 single layer of data.
 You can have P- or S- data. A secondary file defines the latlon 
 values. A third file defines the depth of the layers in km.
*/
bool Data::gypsumRead(){

  int ferror;

  double * lat ;
  double * lon ;
  double * dpth;

  int file_nval=0; // number of lines in file
  char buf[255]; 
  char tmpinfile[64];
  nval = 0; 

  // We know for gypsum how many layers there are:
  dpth = new double[gypsuminnumfiles];

  // get layer depths #######################
  FILE * fptr=fopen(gypsumdepthinfile,"r");
  if (fptr==NULL){
    printf("Cannot open file %s for reading\n",gypsumdepthinfile);
    return 1; //fail
  }  
  
  for ( int i = 0 ; i<gypsuminnumfiles ; i++ ){
    // Collect latlon data
    ferror = fscanf(fptr,"%s", buf); dpth[i] = atof (buf) ;
  }

  //printf("Depths: "); for (int i=0;i<gypsuminnumfiles;i++) printf("%6.2f ",dpth[i]);printf("\n");

  fclose(fptr);
  // End of - get layer depths #######################
 
  // Obtain latlon values  #######################
  fptr=fopen(gypsumlatloninfile,"r");
  if (fptr==NULL){
    printf("Cannot open file %s for reading\n",gypsumlatloninfile);
    return 1; //fail
  }

  while( fgets(buf,sizeof(buf),fptr) != NULL) {
    file_nval++;
  }

  lat = new double[file_nval];
  lon = new double[file_nval];

  // move back to the beginning of the file
  fclose(fptr);
  fptr=fopen(gypsumlatloninfile,"r");
  if (fptr==NULL){
    printf("Cannot open file %s for reading\n",gypsumlatloninfile);
    return 1; //fail
  }

  for ( int i = 0 ; i<file_nval ; i++ ){
    // Collect latlon data
    ferror = fscanf(fptr,"%s", buf); lat[i] = atof (buf) * pi/180. ;
    ferror = fscanf(fptr,"%s", buf); lon[i] = atof (buf) * pi/180. ;
  }
  fclose(fptr);
  // End of -  Obtain latlon values  #######################

  
  // define arrays based on nval
  nval = file_nval*gypsuminnumfiles;
  x = new double[nval];
  y = new double[nval];
  z = new double[nval];
  V = new double[nval];

  // find minR, maxR, ndpth, nvalpershell
  minR = gypsumDepth2Radius(dpth[gypsuminnumfiles-1]);
  maxR = gypsumDepth2Radius(dpth[0]);
  ndpth = gypsuminnumfiles;
  nvalpershell = file_nval;  // should be equal across files.

  //int nval_counter = 0;
  
  // Collect Data....
  for (int f=0 ; f<ndpth ; f++ ) {
      sprintf(tmpinfile,"%s.%02d.txt",infile,f+1);
 
      fptr=fopen(tmpinfile,"r");
      if (fptr==NULL){
          printf("Cannot open file %s for reading\n",tmpinfile);
	  return 1; //fail
      }
     
      // test
      //printf("Radii: f=%02d %6.2f\n",f,gypsumDepth2Radius(dpth[f]));

      for ( int i = 0 ; i<nvalpershell ; i++ ){
	  // Collect data
          ferror = fscanf(fptr,"%s", buf);  V[f*nvalpershell+i] = atof (buf) ;
	  
	  //printf("%7.5f\n",V[f*nvalpershell+i]);
          //if (i==0){ // test to make sure our data is shifting nicely along the layers
	    //printf("f=%02d i=%05d V=%8.6f",f,i,V[f*nvalpershell+i]);
	  //}

	  // Convert to xyz
	  x[f*nvalpershell+i] =   gypsumDepth2Radius(dpth[f]) * cos(lat[i]) * cos(lon[i]) ;
	  z[f*nvalpershell+i] =   gypsumDepth2Radius(dpth[f]) * sin(lat[i]) ;
	  y[f*nvalpershell+i] =   gypsumDepth2Radius(dpth[f]) * cos(lat[i]) * sin(lon[i]) ;

	  //printf("data xyz test  %7.5g  %7.5g  %7.5g  %7.5g  %7.5g  %7.5g  %7.5g  %7.5g\n",x[i],y[i],z[i],x[i]*x[i]+y[i]*y[i]+z[i]*z[i],lat[i],lon[i],sin(lat[i]),cos(lat[i]),sin(lon[i]),cos(lon[i]));
      }
      fclose(fptr);
  } // for dpth


  printf("data xyz test x y z radius\n");
  for (int ii=0;ii<nval;ii++){
    printf("data xyz test %7.5g %7.5g %7.5g %7.5g\n",x[ii],y[ii],z[ii],x[ii]*x[ii]+y[ii]*y[ii]+z[ii]*z[ii]);
  }

  return 0; //success
}


////////////////////////////////////////
// Data::intypeConverter
/**
 Defines the input data type. 
 \return buf character string describing the input type
 */
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
	    case GYPSUMP : strcpy(buf, "GYPSUMP");
		break;
	    case GYPSUMS : strcpy(buf, "GYPSUMS");
     	        break;
	    default: strcpy(buf, "UNDEF");
	}
    return buf;
}

////////////////////////////////////////
// Data::outtypeConverter
/**
 Defines the output data type. 
 \return buf character string describing the output type
 */
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
		break;
	}
    return buf;
}

////////////////////////////////////////
// Data::interpConverter
/**
 Defines the interpolation routine to be used during conversion. 
 \return buf character string describing the interp selection
 */
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
// Data::getStats
/**
 Redirects to various functions based on input type 
 */
bool Data::getStats()
{
  if ( intype == MITP || intype == FILT || intype == GYPSUMP || intype == GYPSUMS ){
	return getStatsData();
    }else{
	return getStatsGrid();
    }
    
    return 0;
}

////////////////////////////////////////
// Data::getStatsGrid
/**
  Outputs input statistics to stdout 
 */
bool Data::getStatsGrid()
{
    double max = -veryLarge,
           min =  veryLarge;

    nval = mvis->npts;
    ndpth = mvis->nr;
	
    for ( int d=0 ; d<10 ; d++){
	for ( int i=0; i<nval/10 ; i++){
	    if ( mvis->domains[d].V[i] > max )  max=mvis->domains[d].V[i] ;
	    if ( mvis->domains[d].V[i] < min )  min=mvis->domains[d].V[i] ;
	    Vmean += mvis->domains[d].V[i];
	}
    }
    Vmax   = max;
    Vmin   = min;
    Vmean /= nval;
    for ( int d=0 ; d<10 ; d++){
	mvis->domains[d].Vmax = max;
	mvis->domains[d].Vmin = min;
    }

    printf("Input Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  nvals        =  %12ld                  |\n"       , nval);
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

    return 0;
}

////////////////////////////////////////
//Data:: getStatsData
/**
 Prints out statistical data
 */
bool Data::getStatsData()
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
  //  if (findLayers()){
  //      printf("Error in Data::findLayers().");
  //      return 1;
  //  }
  

    printf("Input Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  nvals        =  %12ld                  |\n"       , nval);
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
