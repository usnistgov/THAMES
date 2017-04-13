/**
@file global.h
@brief Declare and assign global variables for THAMES.

*/
#ifndef GLOBALH
#define GLOBALH

/*
#define DEBUG
#define VERBOSE
*/

// Flag to indicate exiting the program
const int QUIT_PROGRAM = 1;

// Flag to indicate simulation of hydration only
const int HYDRATION = 2;

// Flag to indicate simulation of leaching only
const int LEACHING = 3;

// Flag to indicate simulation of external sulfate only
const int SULFATE_ATTACK = 4;


// Flag to indicate kinetic equations handled within code (NOT USED)
const int INTERNAL_KINETICS = 0;

// Flag to indicate kinetics is done a priori outside the code (NOT USED)
const int EXTERNAL_KINETICS = 1;


// Normal return condition flag (NOT USED)
const int RETURN_NORMAL = 0;


// Flag to indicate an element is out of bounds (NOT USED)
const int DB_EOB = 1;

// Flag to indicate call to GEM failed (NOT USED)
const int GEMRUN_ERROR = 1;

// Flag to indicate that a lattice element is out of bounds (NOT USED)
const int MESH_EOB = 1;

// Flag to indicate that an array element is out of bounds (NOT USED)
const int ARRAY_EOB = 1;

// Flag to indicate that a file could not be opened (NOT USED)
const int FILE_OPEN_ERROR = 2;

// Flag to indicate unexpected end of file (NOT USED)
const int PREMATURE_EOF = 3;

// Flag to indicate datum is bad or of wrong kind
const int INVALID_INPUT = 4;


// Special phase ids that are important and must always be the same value
const int VOIDID = 0;
const int WATERID = 1;
const int C3SID = 2;
const int C2SID = 3;
const int C3AID = 4;
const int C4AFID = 5;
const int NUMCLINKERPHASES = 4;

// The number of face, edge, and corner neighbors to a cubic lattice site
const unsigned int NUM_NEAREST_NEIGHBORS = 6;
const unsigned int NUM_SECONDNEAREST_NEIGHBORS = 12;
const unsigned int NUM_THIRDNEAREST_NEIGHBORS = 8;

// Maximum allowed string length (NOT USED)
const int MAXSTRING = 128;

// Reference temperature for kinetic calculations [K]
const double REFTEMP = 298.15;

// Reference lattice resolution [micrometers]
const double REFRES = 4.0;

// Ideal gas constant [J/mol/K]
const double GASCONSTANT = 8.314;

// Flags for different kinds of probability distributions (NOT USED)
const int DELTA = +1;
const int UNIFORM = +2;
const int GAUSSIAN = +3;
const int LOGNORMAL = +4;

//
// Next constant specifies boundary condition configuration
// 0 = periodic boundaries everywhere
// 1 = periodic in y and z, greased in x
// 2 = periodic in x and z, greased in y
// 3 = periodic in x and y, greased in z
//

const int BC = 0;

// Growth mode constants
const unsigned int DLA = 0;
const unsigned int WMC = 1;
const unsigned int IWMC = 2;

// Saturation value for colors
const double COLORSATVAL = 255.0;

// Flags for data formats (NOT USED)
const int GEMSFORMAT = 0;
const int EXTERNALFORMAT = 0;


#include "Exceptions.h"
#include <stdexcept>
#endif
