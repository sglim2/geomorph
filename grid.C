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
    : mt(), nt(), nd(), nr(), nproc(), rmax(), rmin()
{

}


////////////////////////////////////////
// Grid::Grid
// 
// Assume nr=mt/2
//  i.e # radial layers = nr+1
//
Grid::Grid(int _mt, int _nt, int _nd)
    : mt(), nt(), nd(), nr(), nproc(), rmax(), rmin()
{
  mt=_mt;
  nt=_nt;
  nd=_nd;
  
  nr=mt/2;
  nproc=pow((mt/nt),2)*10/nd;
  
}

////////////////////////////////////////
// operators
////////////////////////////////////////

////////////////////////////////////////
// methods
////////////////////////////////////////


////////////////////////////////////////
// Grid::gridfind
//
// A basic comparison of any given gris structire to the TERRA
// grid. Returns the mt value of the closest matching TERRA grid
// comparing the total number of grid points acroos the two grid
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

bool Grid::genGrid()
{

    npts  = (10*mt*mt + 2)*(mt/2 + 1);
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
    
//    x = new double[npts];
//    y = new double[npts];
//    z = new double[npts];
//    V = new double[npts];

    // mt =1 for test purposes only
    mt = 2;
    xn = new double[(mt+1)*(mt+1)*10*3];

    grdgen_(xn,&mt);

    for ( int i=0 ; i<10 ; i++){
      printf("\nid %d\n",i);
      for ( int yi=0 ; yi<mt+1 ; yi++){
	for ( int xi=0 ; xi<mt+1 ; xi++){
	  printf("%12.6g\t",xn[0*10*(mt+1)*(mt+1) + i*(mt+1)*(mt+1) + yi*(mt+1) + xi]);
	}
      }
    }            
    printf("\n");

    printf("%d\n",2*10*(mt+1)*(mt+1) + 9*(mt+1)*(mt+1) + mt*(mt+1) + mt + 1);

    printf("%d\n",(mt+1)*(mt+1)*10*3);

    for ( int k=0 ; k<3 ; k++){
      for ( int i=0 ; i<10 ; i++){
	for ( int yi=0 ; yi<mt+1 ; yi++){
	  for ( int xi=0 ; xi<mt+1 ; xi++){
	    printf("%12.6g ",xn[k*10*(mt+1)*(mt+1) + i*(mt+1)*(mt+1) + yi*(mt+1) + xi]);
	  }
	}
      }      
      printf("\n");
    }

    return 0;

}
