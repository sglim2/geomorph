/*
 *
 * Global Definitions
 *
 * 
 */


#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "grid.h"
#include "data.h"

// constants
const double pi=3.141592653589793;

// Earth Data
const double EarthRadKM = 6371;
const double CoreRadiusKM = 3400;

// globals
static char* programName = 0;
//static Data* data=0;

//static double a;
//static double cmb;
static double veryLarge=1E+99;
static double verySmall=1E-99;
static double quiteLarge=1E+5;
static double quiteSmall=1E-5;


#endif // MAIN_H
