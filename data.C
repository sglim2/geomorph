/*
 *
 * Data Manipulation
 *
 * 
 */


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

bool Data::Read()
{
  int linecount=0;
  
  FILE * fptr=fopen(infile,"r");
  if (fptr==NULL){
    return 1; //fail
  }
  
  char buf[255]; 
  while( fgets(buf,sizeof(buf),fptr) != NULL) {
    linecount++;
  } 
  

  printf("linecount = %d\n", linecount);

  // get rid of the header
  linecount--;

  // define arrays
  x = new double[linecount];
  y = new double[linecount];
  z = new double[linecount];
  V = new double[linecount];

  fseek ( fptr , 1L , SEEK_SET );
  fin.read(seat[0],3*linecount);

  //  for ( int i = 0 ; i<linecount ; i++ ){
    

  } 

  fclose(fptr);
  return 0; //success
}

////////////////////////////////////////
// methods
////////////////////////////////////////
