/*
 *
 * GeoMorph
 *
 * A grid conversion tool for popular geodynamics and seismic tomography data
 * 
 * Exmaple:
 * $ > ./geomorph --mt 256 --nt 16 --nd 10 --infile ../data/MITP08.txt --outfile mvis001 --intype mitp --outtype mvis
 */

#include "main.H"


////////////////////////////////////////
// return the basename of a path
static char* gm_basename(char* path)
{
  char* name = &path[0];
  char* ptr  = &path[strlen(path) - 1];

  while (ptr > path) {
    if (*(ptr - 1) == '\\') {
      name = ptr;
      break;
    }
    ptr--;
  }

  return name;
}

////////////////////////////////////////
// print usage instructions
bool gm_usage()
{
  
  printf("usage: %s [options] [args]\n\n", programName);
  
  printf("Available Commands:\n");
  printf(" --help . . . . . . print help-page and exit\n");
  printf(" --intype . . . . . input type [ MVIS | TERRA | MITP ]\n");
  printf(" --outtype. . . . . output type [ MVIS | TERRA | MITP ]\n");
  printf(" --infile,  -i. . . input filename\n");
  printf(" --outfile, -o. . . output filename base\n");
  printf(" --mt . . . . . . . \n");
  printf(" --nt . . . . . . . \n");
  printf(" --nd . . . . . . . \n");


  return 0;
}

////////////////////////////////////////
// process the commandline arguments
bool gm_processCommandLine(int argc, char* argv[])
{
  /*
    --help . . . . . . print help-page and exit
    --intype . . . . . input type [ MVIS | TERRA | MITP ]
    --outtype. . . . . output type [ MVIS | TERRA | MITP ]
    --infile,  -i. . . input filename
    --outfile, -o. . . output filename base
    --mt . . . . . . . 
    --nt . . . . . . . 
    --nd . . . . . . . 
   */

  programName = gm_basename(argv[0]);

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      gm_usage();
      exit(0);
    }

    if (strcmp(argv[i], "--intype") == 0) {
      i++;
      if (strcmp(argv[i], "mvis") == 0) {
	data->intype = data->MVIS;
      }else if (strcmp(argv[i], "terra") == 0) {
	data->intype = data->TERRA;
      }else if (strcmp(argv[i], "mitp") == 0) {
	data->intype = data->MITP;
      }else {
	gm_usage();
	exit(0);
      }
    }
    if (strcmp(argv[i], "--outtype") == 0) {
      i++;
      if (strcmp(argv[i], "mvis") == 0) {
	data->outtype = data->MVIS;
      }else if (strcmp(argv[i], "terra") == 0) {
	data->outtype = data->TERRA;
      }else if (strcmp(argv[i], "mitp") == 0) {
	data->outtype = data->MITP;
      }else {
	gm_usage();
	exit(0);
      }
    }

    if (strcmp(argv[i], "--infile") == 0) {
      i++;
      data->infile = new char[strlen(argv[i])];
      strcpy(data->infile,argv[i]);
    }
    if (strcmp(argv[i], "--outfile") == 0) {
      i++;
      data->outfile = new char[strlen(argv[i])];
      strcpy(data->outfile,argv[i]);
    }
  
    if (strcmp(argv[i], "--mt") == 0) {
      i++;
      data->grid.mt=atoi(argv[i]);
    }
    if (strcmp(argv[i], "--nt") == 0) {
      i++;
      data->grid.nt=atoi(argv[i]);
    }
    if (strcmp(argv[i], "--nd") == 0) {
      i++;
      data->grid.nd=atoi(argv[i]);
    }

  }

  printf("\nInput Arguments:\n");
  printf(" intype   %s(%d) \n", data->intypeConverter(), data->intype);
  printf(" outtype  %s(%d) \n", data->outtypeConverter(), data->outtype);
  printf(" infile   %s \n", data->infile);
  printf(" outfile  %s \n", data->outfile);
  printf(" mt       %d \n", data->grid.mt); 
  printf(" nt       %d \n", data->grid.nt); 
  printf(" nd       %d \n", data->grid.nd);
  printf("\n");

  return 0;
}

////////////////////////////////////////
// main
int main(int argc, char* argv[])
{

  // Create a new (empty) Data instance
  data = new Data;
  
  // process command-line and fill some data objects
  if (gm_processCommandLine(argc, argv)) {
    gm_usage();
    exit(-1);
  }

  // read input and convert to geomorph x,y,z,V structure
  //  if (data->Read()){
  //  printf("Error reading input file.");
  //};

  // get input stats
  if (data->getStats()){
    printf("Error computing Stats.");
  }

  // find best grid to match input data
  int my_mt= data->grid.findGrid(data->nval);
  printf("best matching mt value = %d\n", my_mt);
  if ( data->grid.mt != my_mt ) {
      printf("Overriding any user-specified mt value.\nUsing mt=%d\n", my_mt);
      data->grid.mt = my_mt; 
  }

  if (data->grid.genGrid()){
    printf("Error computing TERRA stats.\n");
  }

  // generate grid 
//  data->defineGrid();
  
  //convert to terra-grid

  // convert to output format



  data->destroyGrid();
  delete [] data;

 return 0;

}
