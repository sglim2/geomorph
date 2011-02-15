/*
 * +==========================================================+
 * |                                                  __      |
 * |   ____   ____  ____   _____   ____  _____  ____ |  |__   |
 * |  / __ \ / __ \/  _ \ /     \ /  _ \/  __ \/___ \|  |  \  |
 * | / /_/  >  ___(  <_> )  Y Y  (  <_> )  | \/  |_> >   Y  \ |
 * | \___  / \____ \____/|__|_|__/\____/|__|  |   __/|___|__/ |
 * |/_____/                                   |__|            |
 * +==========================================================+
 *
 * A grid conversion tool for popular geodynamics and seismic tomography data
 *
 * Exmaple 1:
 * =========
 * $ > ./geomorph --mt 256 --nt 16 --nd 10 --infile ../data/MITP08.txt --outfile mvis001 --intype mitp --outtype mvis --interp nearest2
 *
 * Exmaple 2:
 * =========
 * $ > ./geomorph --mt 32 --nt 16 --nd 10 --infile ../data/Filt --outfile mvis001 --intype filt --filtinstart 50 --filtinend 2850 --filtinnum 57 --outtype mvis --interp nearest2 
 *
 */

#ifndef GEO_TUI_
#include <QtGui/QApplication>
#include "geomorphmainwindow.h"
#endif

#include "main.h"

// globals
Grid* grid=0;
Data* data=0;

// GUI already has an object called 'data'
Data * gdata=0;

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
  printf(" --interp . . . . . interpolation routine [ nearest | nearest2 ]\n");


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

  if (!data){
    printf("Error: Data::data no initialized.\n");
    exit(0);
  }

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
      }else if (strcmp(argv[i], "filt") == 0) {
        data->intype = data->FILT;
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
      }else if (strcmp(argv[i], "filt") == 0) {
        data->outtype = data->FILT;
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

    if (strcmp(argv[i], "--filtinstart") == 0) {
      i++;
      data->filtinstart = atoi(argv[i]);
    }
    if (strcmp(argv[i], "--filtinend") == 0) {
      i++;
      data->filtinend = atoi(argv[i]);
    }
    if (strcmp(argv[i], "--filtinnum") == 0) {
      i++;
      data->filtinnumfiles = atoi(argv[i]);
    }

    if (strcmp(argv[i], "--filtoutstart") == 0) {
      i++;
      data->filtoutstart = atoi(argv[i]);
    }
    if (strcmp(argv[i], "--filtoutend") == 0) {
      i++;
      data->filtoutend = atoi(argv[i]);
    }
    if (strcmp(argv[i], "--filtoutnum") == 0) {
      i++;
      data->filtoutnumfiles = atoi(argv[i]);
    }

    if (strcmp(argv[i], "--mt") == 0) {
      i++;
      grid->mt=atoi(argv[i]);
    }
    if (strcmp(argv[i], "--nt") == 0) {
      i++;
      grid->nt=atoi(argv[i]);
    }
    if (strcmp(argv[i], "--nd") == 0) {
      i++;
      grid->nd=atoi(argv[i]);
    }

    if (strcmp(argv[i], "--interp") == 0) {
      i++;
      if (strcmp(argv[i], "nearest") == 0) {
        data->interp = data->NEAREST;
      }else if (strcmp(argv[i], "nearest2") == 0) {
        data->interp = data->NEAREST2;
/*      }else if (strcmp(argv[i], "linear") == 0) {
        data->interp = data->LINEAR;
      }else if (strcmp(argv[i], "cubic") == 0) {
        data->interp = data->CUBIC;   */
      }else {
        gm_usage();
        exit(0);
      }
    }

  }

  printf("\nInput Arguments:\n");
  printf(" intype      %s(%d) \n", data->intypeConverter(), data->intype);
  printf(" outtype     %s(%d) \n", data->outtypeConverter(), data->outtype);
  printf(" infile      %s \n", data->infile);
  printf(" outfile     %s \n", data->outfile);
  printf(" filtinstart %d \n", data->filtinstart);
  printf(" filtinend   %d \n", data->filtinend);
  printf(" filtinnum   %d \n", data->filtinnumfiles);
  printf(" mt          %d \n", grid->mt);
  printf(" nt          %d \n", grid->nt);
  printf(" nd          %d \n", grid->nd);
  printf(" interp      %s(%d) \n", data->interpConverter(), data->interp);
  printf("\n");

  return 0;
}

////////////////////////////////////////
// writeGrid
//
// writes the 'rad' layer of the geomorph grid only.
//
bool writeGrid()
{
  int _idmax=10;
  int rad=0;

    // test some data
    Domain * dptr = grid->domains;
    FILE * Xptr=fopen("outX","w");
    if (Xptr==NULL){
        return 1; //fail
    }
    FILE * Yptr=fopen("outY","w");
    if (Yptr==NULL){
        return 1; //fail
    }
    FILE * Zptr=fopen("outZ","w");
    if (Zptr==NULL){
        return 1; //fail
    }
    FILE * Cptr=fopen("outC","w");
    if (Cptr==NULL){
        return 1; //fail
    }
    dptr = grid->domains;
    for ( int id = 0 ; id < _idmax ; id++ ){
        for ( int i2=0 ; i2<(dptr[id].mt+1) ; i2++ ){
            for ( int i1=0 ; i1<(dptr[id].mt+1) ; i1++ ){
                fprintf(Xptr,"%12.8g\t",dptr[id].xn[dptr[id].idx(rad,i2,i1)]);
            }
            fprintf(Xptr,"\n");
        }
    }
    dptr = grid->domains;
    for ( int id = 0 ; id < _idmax ; id++ ){
        for ( int i2=0 ; i2<(dptr[id].mt+1) ; i2++ ){
            for ( int i1=0 ; i1<(dptr[id].mt+1) ; i1++ ){
                fprintf(Yptr,"%12.8g\t",dptr[id].yn[dptr[id].idx(rad,i2,i1)]);
            }
            fprintf(Yptr,"\n");
        }
    }
    dptr = grid->domains;
    for ( int id = 0 ; id < _idmax ; id++ ){
        for ( int i2=0 ; i2<(dptr[id].mt+1) ; i2++ ){
            for ( int i1=0 ; i1<(dptr[id].mt+1) ; i1++ ){
                fprintf(Zptr,"%12.8g\t",dptr[id].zn[dptr[id].idx(rad,i2,i1)]);
            }
            fprintf(Zptr,"\n");
        }
    }
    dptr = grid->domains;
    for ( int id = 0 ; id < _idmax ; id++ ){
        for ( int i2=0 ; i2<(dptr[id].mt+1) ; i2++ ){
            for ( int i1=0 ; i1<(dptr[id].mt+1) ; i1++ ){
                fprintf(Cptr,"%12.8g\t",dptr[id].V[dptr[id].idx(rad,i2,i1)]);
            }
            fprintf(Cptr,"\n");
        }
    }

    return 0; // success
}

////////////////////////////////////////
// process
//
void geo_process()
{

  // read input and convert to geomorph x,y,z,V structure
  if (data->Read()){
      printf("Error reading input file.\n");
  }

  // get input stats
  if (data->getStats()){
    printf("Error computing Stats.\n");
  }

  // find best grid to match input data
  int my_mt = grid->suggestGrid(data->nval);
  printf("best matching mt value = %d\n", my_mt);
//  if ( grid->mt != my_mt ) {
//  printf("Overriding any user-specified mt value.\nUsing mt=%d\n", my_mt);
//      grid->mt = my_mt;
//  }
  if ( grid->mt == 0 ){
      grid->mt = my_mt;
  }else if ( grid->mt != my_mt ){
      printf("Warning: User-specified value on mt is not optimal. Try using mt=%d\n", my_mt);
  }

  printf("Generating grid....\n");
  if (grid->genGrid(data->cmb)){
    printf("Error computing TERRA stats.\n");
  }

  // import Data::data to geomorph grid
  printf("Converting Data....\n");
  if (grid->importData(data)){
    printf("Error importing Data into geomorph grid.\n");
  }

  // write a single shell (for use in matlab, or alternative)
//  if (writeGrid()){
//    printf("Error writing Grid data.\n");
//  }

  printf("Exporting Data....\n");
  if (grid->exportGrid(data)){
    printf("Error exporting Data failed.\n");
    //    return 1; // fail
  }

  // delete [] data;
  // delete [] grid;
}



int main(int argc, char *argv[])
{

    grid = new Grid;
    data = new Data;

#ifndef GEO_TUI_
    gdata= data;
#endif

    // process command-line and fill some data objects
    if (gm_processCommandLine(argc, argv)) {
      gm_usage();
      exit(-1);
    }


#ifdef GEO_TUI_
    geo_process();
#else
    QApplication geo(argc, argv);
    GeomorphMainWindow geoWin;
    geoWin.show();
    return geo.exec();
#endif

}
