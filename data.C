/*
 *
 * Data Manipulation
 *
 * 
 */

#include "main.H"
#include "data.H"

/*
 * utilities
 */

////////////////////////////////////////
// construction/destruction
////////////////////////////////////////


Data::Data()
  : infile(), outfile(), intype(), outtype(), nval(), x(), y(), z(), V(), grid()
{

  // Create new (empty) Grid instance
  Grid();

  intype=UNDEF;
  outtype=UNDEF;

}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////

////////////////////////////////////////
// Read
bool Data::Read()
{
    switch (intype) {
	case UNDEF :
	    printf("Error.Intyoue undefined");
	    return 1; //fail
//	    break;  // unreachable!!
	case MVIS :
	    // do something
	    break;
	case TERRA :
	    // do something
	    break;
	case MITP :
	    return mitpRead();
//	    break; // unreachable!!
	default :
	    printf("Error.Intyoue undefined");
	    return 1; //fail
//	    break;  // unreachable!!
    }

    return 0;    
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
bool Data::mitpRead()
{
  nval = 0;
  
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
  double lat ;
  double lng ;
  double dpth;

  x = new double[nval];
  y = new double[nval];
  z = new double[nval];
  V = new double[nval];

  // move tothe first data-line
  fclose(fptr);
  fptr=fopen(infile,"r");
  if (fptr==NULL){
      return 1; //fail
  }
  fscanf(fptr,"%s", &buf);
  fscanf(fptr,"%s", &buf);
  fscanf(fptr,"%s", &buf);
  fscanf(fptr,"%s", &buf);

  for ( int i = 0 ; i<nval ; i++ ){
      // Collect data
      fscanf(fptr,"%s", &buf);  lat = atof (buf) * pi/180. ;
      fscanf(fptr,"%s", &buf);  lng = (atof (buf) - 180.) * pi/180.  ;
      fscanf(fptr,"%s", &buf);  dpth = (EarthRadKM - atof (buf)) / EarthRadKM;
      fscanf(fptr,"%s", &buf);  V[i] = atof (buf) ;

      // Convert to xyz
      x[i] = - 1. * cos(lat) * cos(lng) ;
      y[i] =   1. * sin(lat) ;
      z[i] =   1. * cos(lat) * sin(lng) ;

  }

  
  for ( int i = 0 ; i<20 ; i++ ){
      printf ("%12.8g\t%12.8g\t%12.8g\t%12.8g\n", x[i],y[i],z[i],V[i]);
  }

  fclose(fptr);

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
	    case TERRA: strcpy(buf, "TERRA");
		break;
	    case MITP : strcpy(buf, "MITP");
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
	    case TERRA: strcpy(buf, "TERRA");
		break;
	    case MITP : strcpy(buf, "MITP");
		break;
	    default: strcpy(buf, "UNDEF");
	}
    return buf;
}

////////////////////////////////////////
// getStats
bool Data::getStats()
{
    double veryLarge=1E+99;

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

    printf("Some Stats....\n");
    printf("+----------------------------------------------+\n");
    printf("|  nvals      =  %12d                  |\n"        , nval);
    printf("|  V(mean)    =  %12.8g                  |\n"      , Vmean);
    printf("|  Vmin       =  %12.8g                  |\n"      , Vmin);
    printf("|  Vmax       =  %12.8g                  |\n"      , Vmax);
    printf("+----------------------------------------------+\n");
    
    
    return 0; //success
}

////////////////////////////////////////
// defineGrid
bool Data::defineGrid()
{
   
    grid.x = new double[ grid.nt+1 * grid.nt+1 * 10 ];
    grid.y = new double[ grid.nt+1 * grid.nt+1 * 10 ];
    grid.z = new double[ grid.nt+1 * grid.nt+1 * 10 ];
    grid.V = new double[ grid.nt+1 * grid.nt+1 * 10 ];
  
    return 0; // success
}

////////////////////////////////////////
// destroyGrid
bool Data::destroyGrid()
{
   
    delete [] grid.x;
    delete [] grid.y;
    delete [] grid.z;
    delete [] grid.V;

    return 0; // success
}
