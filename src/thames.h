/**
@file thames.h
@brief Header for the main THAMES program.

*/

/**
@mainpage THAMES Overview

@author Jeffrey W. Bullard, NIST
@author Pan Feng, Southeast University
@author Zachary C. Grasley, Texas A&M University
@author Xiaodan Li, Oklahoma State University
@date 2009--2017

@section intro_sec Introduction

THAMES is a model of 3-D microstructure development in cement paste.
The overall philosophy of the model is that, if a kinetic model of
clinker phase consumption can be developed, then that model can be
used to establish the change in concentration of dissolved species
in solution.  This latter data are then used in a thermodynamic
speciation model to predict the solution speciation <b>and</b>
the equilibrium mass fractions of hydration products.  That is,
the kinetic model establishes initial conditions on the solution
at each time step, and these initial conditions are used to solve
for chemical equilibrium of hydration products and pore solution
chemistry at that time step.

This overall approach has been used with considerable success by
several research groups.  Lothenbach and
coworkers [1--4] have used the kinetic model of Parrot and Killoh [5]
for clinker phase consumption to establish the initial conditions on
the pore solution at each time step, and then used the GEMS
geochemical model [6,7] to predict speciation and hydration product
mass fractions.  

THAMES uses this same modeling philosophy, but adds the extra
dimension that the microstructure itself is simulated as a function of 
time, using simple nucleation and dissolution/growth rules for each phase.
The 3-D microstructure is in the form of a digital image (called a
lattice), with each lattice site assuming a unique unsigned integer
value corresponding to a valid phase.  Data about the phases are stored
in a database which must also be input at the beginning of a
simulation.  The microstructure representation used by THAMES is similar to
the CEMHYD3D cellular automaton model, developed primarily by
Bentz [8], between 1994 and 2002 at the National
Institute of Standards and Technology (NIST).  In turn, CEMHYD3D
is based on other numerical modeling work on cement microstructures
undertaken at NIST starting in 1988.

The overall flow of THAMES is shown in the figure.
To achieve this flow, THAMES needs to link to the GEM3K kernel library.

@image html THAMES-Concrete-Modules-small.png
@image latex THAMES-Concrete-Modules-scaled.png

Data about dissolution, nucleation, and growth behavior of each phase is
stored in a database which is read at the beginning of the simulation,
along with the initial 3-D paste microstructure and the kinetic model
used. The kinetic model takes information about the phase surface
areas, degree of reaction, etc. to predict the consumption of clinker
phases during each time step.  The corresponding solution chemistry
is used as an initial condition for a solution speciation thermodynamic
model (GEMS) which predicts the equilibrium solution composition and
mass fractions of each hydrated phase.  The microstructure evolver
updates the microstructure, and the process is repeated until the desired
final time is achieved.

THAMES is written in C++ and is heavily object-oriented.  We believe
that this feature makes it easier to encapsulate and operate on data.
Fewer specialized functions are required as a result, and much more of
the tedious details can be conveniently hidden from the casual user if
desired.  At the same time, modification to the codes becomes easier because
each class can be examined and changed rather independently from the
rest of the programs.

This document contains just the main executable program and its associated
header file.  The details of the various classes that are defined for
THAMES are contained in library files.

@section ref_sec References

-# Lothenbach, B., Winnefeld, F., Thermodynamic modelling of the hydration of Portland cement, Cement and Concrete Research 36 (2006) 209--226.
-# Lothenbach, B., Wieland, E., A thermodynamic approach to the hydration of sulphate-resisting Portland cement, Waste Management 26 (2006) 706--719.
-# Lothenbach, B., Matschei, T., Moeschner, G., Glasser, F.P., Thermodynamic modelling of the effect of temperature on the hydration and porosity of Portland cement, Cement and Concrete Resaerch 38 (2008) 1--18.
-# Lothenbach, B., Le Saout, G., Gallucci, E., Scrivener, K., Influence of limestone on the hydration of Portland cements, Cement and Concrete Research 38 (2008) 848--860.
-# Parrot, L.J., Killoh, D.C., Prediction of cement hydration, British Ceramics Proceedings 35 (1984) 41--53.
-# Kulik, D.A., Wagner, T., Dmytrieva, S.V., Kosakowski, G., Hingerl, F.F., Chudnenko, K.V., Berner, U., GEM-Selektor geochemical modeling package: revised algorithm and GEMS3K numerical kernel for coupled simulation codes, Computational Geosciences 17 (2013) 1--24.
-# Wagner, T., Kulik, D.A., Hingerl, F.F., Dmytrieva, S.V., GEM-Selektor geochemical modeling package: TSolMod library and data interface for multicomponent phase models, Canadian Mineralogist 50 (2012) 701--723.
-# Bentz, D.P., CEMHYD3D: A three-dimensional cement hydration and microstructural development modelling package, version 2.0, NISTIR 6485, U.S. Department of Commerce, April, 2000.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include "thameslib/Controller.h"
#include "thameslib/ChemicalSystem.h"
#include "thameslib/Solution.h"
#include "thameslib/Interface.h"
#include "thameslib/Site.h"
#include "thameslib/Lattice.h"
#include "thameslib/KineticModel.h"
#include "thameslib/AppliedStrain.h"
#include "thameslib/ThermalStrain.h"
#include "thameslib/StrainEnergy.h"

/**
@brief The vector of component elastic energies.

The strainenergy vector is passed to the GEM3K library to modify the
Gibbs energy of formation of the varioud dependent components (DCs) as
a result of elastic deformation, either by an applied load or as a result
of phase transformation misfit strain.
*/
vector<double> strainenergy;

/**
@brief Write the formatted report file listing job properties and input.

Almost all the actual formatted output is done by the ChemicalSystem object
through its `writeChemSys` method.

@param jobRoot is the root name of the THAMES simulation
@param itmie is the start time of the job
@param mfname is name of the initial microstructure image file
@param parfilename is the name of the input parameter file
@param csname is the name of the GEM chemical system definition (CSD) file
@param csys is a pointer to the ChemicalSystem object for the simulation
@param ctr is a pointer to the Controller object for the simulation

If a file is not present, the file name should be given as an empty string.
*/
void writeReport (const string &jobroot,
                  struct tm *itime,
                  const string &mfname,
                  const string &parfilename,
                  const string &csname,
                  ChemicalSystem *csys,
                  Controller *ctr);

using namespace std;
