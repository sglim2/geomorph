/**
  \file 

  \brief 
  Global Definitions
*/


#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "grid.h"
#include "data.h"

// constants
const double pi=3.141592653589793; /**< Pi */

// Earth Data
const double EarthRadKM = 6371;    /**< A defined Earth Radius */
const double CoreRadiusKM = 3400;  /**< A defined Core-Mantle boundary radius */

// globals
static char* programName = 0;   /**< programName will be set to the program name */
//static Data* data=0;

//static double a;
//static double cmb;
static double veryLarge=1E+99;  /**< A very large number */
static double verySmall=1E-99;  /**< A very small number */
static double quiteLarge=1E+5;  /**< A quite large number */
static double quiteSmall=1E-5;  /**< A quite small number */


#endif // MAIN_H
