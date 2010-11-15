/*
 *
 * Global Definitions
 *
 * 
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "grid.H"
#include "data.H"

// constants
const double pi=3.141592653589793;

// Earth Data
const double EarthRadKM = 6371;
const double CoreRadiusKM = 3400;

// globals
static char* programName = 0;
static Data* data=0;
static Grid* grid=0;
static double a;
static double cmb;
