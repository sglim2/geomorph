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
 * $ > ./geomorph --mt 256 --nt 16 --nd 10 --suffix 01 --infile ../data/MITP08.txt --outfile mvis001 --intype mitp --outtype mvis --interp nearest2
 *
 * Exmaple 2:
 * =========
 * $ > ./geomorph --mt 32 --nt 16 --nd 10 --suffix 01 --infile ../data/Filt --outfile mvis001 --intype filt --filtinstart 50 --filtinend 2850 --filtinnum 57 --outtype mvis --interp nearest2
 *
 * Exmaple 3:
 * =========
 * $ > ./geomorph --mt 32 --nt 8 --nd 5 --suffix 01 --infile ../data/MITP08.txt --outfile c001 --intype mitp --outtype terra_cv --interp nearest2
 *
 * Exmaple 4:
 * =========
 * $ > ./geomorph --mt 32 --nt 8 --nd 5 --suffix 01 --infile ../data/MITP08.txt --outfile c002 --intype mitp --outtype terra_cc --interp nearest2
 *
 * Exmaple 5:
 * =========
 * $ > ./geomorph --mt 16 --nt 8 --nd 10 --suffix 01 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile ../data/mvis002 --outfile mvis002 --intype mvis --outtype mvis --interp linear
 *
 * Exmaple 6:
 * =========
 * $ > ./geomorph --mt 32 --nt 8 --nd 10 --suffix 00 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile ../data/c002 --outfile mvis001 --intype terra_cv --outtype mvis --interp linear
 *
 * Exmaple 7:
 * =========
 * $ > ./geomorph --mt 64 --nt 16 --nd 10 --suffix 00 --mtin 32 --ntin 8 --ndin 10 --suffixin 01 --cmbin 0.54940 --infile ../data/c002 --outfile c003 --intype terra_cv --outtype terra_cv --interp linear
 */

#ifndef GEO_TUI_
#include <QtGui/QApplication>
#include "geomorphmainwindow.h"
#endif

#include "main.h"

// globals
Grid* grid=0;
Data* data=0;

bool previewAvailable = false;  // GUI only
bool gridAvailable = false;     // GUI only
bool triangleAvailable = false; // GUI only
int  previewLayer = 0;          // GUI only
bool autoRotate = false;        // GUI only
bool radialAverage = false;     // GUI only

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
  printf(" --intype . . . . . input type [ MVIS | TERRA_CC | TERRA_CV | MITP ]\n");
  printf(" --outtype. . . . . output type [ MVIS | TERRA_CC | TERRA_CV | MITP ]\n");
  printf(" --infile . . . . . input filename\n");
  printf(" --outfile. . . . . output filename base\n");
  printf(" --mt . . . . . . . Desired MT value for output\n");
  printf(" --nt . . . . . . . Desired NT value for output\n");
  printf(" --nd . . . . . . . Desired ND value for output\n");
  printf(" --mtin . . . . . . Input MT value for MVIS/TERRA input-types\n");
  printf(" --ntin . . . . . . Input NT value for MVIS/TERRA input-types\n");
  printf(" --ndin . . . . . . Input ND value for MVIS/TERRA input-types\n");
  printf(" --cmbin. . . . . . Core/Mantle boudary value for MVIS/TERRA input-types\n");
  printf(" --suffix . . . . . Output file-name suffix for MVIS/TERRA output-types\n");
  printf(" --suffixin . . . . Input file-name suffix for MVIS/TERRA input-types\n");
  printf(" --filtinstart. . . FILT-type start depth (in km)\n");
  printf(" --filtinend. . . . FILT-type end depth (in km)\n");
  printf(" --filtinnum. . . . FILT-type number of data files\n");
  printf(" --interp . . . . . interpolation routine [ nearest | nearest2 | linear ]\n");
  printf("\n");
  printf("Inerpolation routines...\n");
  printf("  nearest   - used with intype = MITP or FILT \n");
  printf("  nearest2  - used with intype = MITP or FILT \n");
  printf("  linear    - used with intype = MVIS, TERRA_CV, or TERRA_CC \n");
  printf("\n");
  printf("\n");
  printf("Example Usage...\n");
  printf("\n");
  printf("Converts MIT-P08 data to mvis format:\n");
  printf("./geomorph --mt 256 --nt 16 --nd 10 --suffix 01 \\\n");
  printf("           --infile ../data/MITP08.txt --intype mitp\\\n");
  printf("           --outfile mvis001 --outtype mvis \\\n");
  printf("           --interp nearest2\n");
  return 0;
}

////////////////////////////////////////
// process the commandline arguments
bool gm_processCommandLine(int argc, char* argv[])
{
  /*
    --help . . . . . . print help-page and exit
    --intype . . . . . input type [ MVIS | TERRA_CC | TERRA_CV | MITP ]
    --outtype. . . . . output type [ MVIS | TERRA_CC | TERRA_CV | MITP ]
    --infile . . . . . input filename
    --outfile. . . . . output filename base
    --mt . . . . . . . output mt value
    --nt . . . . . . . output nt value
    --nd . . . . . . . output nd value
    --mtin . . . . . . input mt value
    --ntin . . . . . . input nt value
    --ndin . . . . . . input nd value
    --cmbin. . . . . . Core/Mantle boudary value for MVIS/TERRA input-types
    --suffix . . . . . output suffix (for MVIS and TERRA) 
    --suffixin . . . . input suffix (for MVIS and TERRA) 
    --filtinstart. . . FILT-type start depth (in km)
    --filtinend. . . . FILT-type end depth (in km)
    --filtinnum. . . . FILT-type number of data files
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
      data->intypeSet=true;
      if (strcmp(argv[i], "mvis") == 0) {
        data->intype = data->MVIS;
      }else if (strcmp(argv[i], "terra_cc") == 0) {
        data->intype = data->TERRA_CC;
      }else if (strcmp(argv[i], "terra_cv") == 0) {
        data->intype = data->TERRA_CV;
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
      data->outtypeSet=true;
      if (strcmp(argv[i], "mvis") == 0) {
        data->outtype = data->MVIS;
      }else if (strcmp(argv[i], "terra_cc") == 0) {
        data->outtype = data->TERRA_CC;
      }else if (strcmp(argv[i], "terra_cv") == 0) {
        data->outtype = data->TERRA_CV;
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
      data->infile = new char[strlen(argv[i])+1];
      strcpy(data->infile,argv[i]);
      data->infileSet=true;
    }
    if (strcmp(argv[i], "--outfile") == 0) {
      i++;
      data->outfile = new char[strlen(argv[i])+1];
      strcpy(data->outfile,argv[i]);
      data->outfileSet=true;
    }

    if (strcmp(argv[i], "--filtinstart") == 0) {
      i++;
      data->filtinstart = atoi(argv[i]);
      data->filtinstartSet = true;
    }
    if (strcmp(argv[i], "--filtinend") == 0) {
      i++;
      data->filtinend = atoi(argv[i]);
      data->filtinendSet = true;
    }
    if (strcmp(argv[i], "--filtinnum") == 0) {
      i++;
      data->filtinnumfiles = atoi(argv[i]);
      data->filtinnumfilesSet = true;
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

    if (strcmp(argv[i], "--mtin") == 0) {
      i++;
      data->mvis->mt=atoi(argv[i]);
      data->mvis->mtSet=true;
    }
    if (strcmp(argv[i], "--ntin") == 0) {
      i++;
      data->mvis->nt=atoi(argv[i]);
      data->mvis->ntSet=true;
    }
    if (strcmp(argv[i], "--ndin") == 0) {
      i++;
      data->mvis->nd=atoi(argv[i]);
      data->mvis->ndSet=true;
    }

    if (strcmp(argv[i], "--suffixin") == 0) {
      i++;
      data->mvis->suffix=atoi(argv[i]);
      data->mvis->suffixSet=true;
    }
 

    if (strcmp(argv[i], "--mt") == 0) {
      i++;
      grid->mt=atoi(argv[i]);
      grid->mtSet=true;
    }
    if (strcmp(argv[i], "--nt") == 0) {
      i++;
      grid->nt=atoi(argv[i]);
      grid->ntSet=true;
    }
    if (strcmp(argv[i], "--nd") == 0) {
      i++;
      grid->nd=atoi(argv[i]);
      grid->ndSet=true;
    }

    if (strcmp(argv[i], "--suffix") == 0) {
      i++;
      grid->suffix=atoi(argv[i]);
      grid->suffixSet=true;
    }
 

    if (strcmp(argv[i], "--cmbin") == 0) {
      i++;
      data->cmb=atof(argv[i]);
      data->cmbinSet=true;
    }

    if (strcmp(argv[i], "--interp") == 0) {
      i++;
      data->interpSet = true;
      if (strcmp(argv[i], "nearest") == 0) {
        data->interp = data->NEAREST;
      }else if (strcmp(argv[i], "nearest2") == 0) {
        data->interp = data->NEAREST2;
      }else if (strcmp(argv[i], "linear") == 0) {
        data->interp = data->LINEAR;
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
  printf(" mtin        %d \n", data->mvis->mt);
  printf(" ntin        %d \n", data->mvis->nt);
  printf(" ndin        %d \n", data->mvis->nd);
  printf(" cmbin       %6.5f \n", data->cmb);
  printf(" suffix      %d \n", grid->suffix);
  printf(" suffixin    %d \n", data->mvis->suffix);
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

  if ( grid->mt == 0 ){
      grid->mt = my_mt;
  }else if ( grid->mt != my_mt ){
      printf("Warning: User-specified value on mt is not optimal. Try using mt=%d\n", my_mt);
  }

  if (grid->genGrid(data->cmb)){
    printf("Error computing grid statistics.\n");
  }

  // import Data::data to geomorph grid
  if (grid->importData(data)){
    printf("Error importing Data into geomorph grid.\n");
  }

  // write a single shell (for use in matlab, or alternative)
//  if (writeGrid()){
//    printf("Error writing Grid data.\n");
//  }

  if (grid->exportGrid(data)){
    printf("Error exporting Data failed.\n");
  }

  // delete [] data;
  // delete [] grid;
}



int main(int argc, char *argv[])
{

    grid = new Grid;
    data = new Data;
    data->mvis = new Grid;

#ifndef GEO_TUI_
    gdata = data;
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
#ifdef _WIN32
    geoWin.setWindowIcon(QIcon("globe1.png"));
#endif
    geoWin.show();
    return geo.exec();
#endif

}
