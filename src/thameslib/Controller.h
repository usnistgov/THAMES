/**
@file Controller.h
@brief Declaration of the Controller class.
*/

#ifndef CONTROLLERH
#define CONTROLLERH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>
#include "Lattice.h"
#include "KineticModel.h"
#include "ThermalStrain.h"
#include "Site.h"
#include "global.h"

using namespace std;

/**
@class Controller
@brief Controls the running of simulation iterations.

The `Controller` class is the hub for THAMES simulations.  It has pointers
to all the other major objects that are instantiated by THAMES, including
the KineticModel, the ThermalStrain, the ChemicalSystem, and the Lattice
that are associated with the system.

The `Controller` object is also responsible for running each iteration
of the simulation and deciding which modules to run based on whether
hydration, leaching, or sulfate attack is desired.

In THAMES, the `Lattice` class can be thought of as the system.
By itself, it can identify what materials it contains, the properties of
the materials, the temperature, current age, degree of hydration, etc.
However, the `Lattice` class does not possess a <i>driver</i> to
control the microstructure development.  This latter
functionality is contained in the `Controller` class, which operates
directly on the `Lattice` to determine how to modify the lattice at
each time step specified in the phase input file.

THAMES keeps track of physical and chemical data
about individual material phases, including
specific gravity, internal porosity, composition, molar volume, etc.  All the
phases are stored in a phase database.

The ultimate objective of the `Controller` class is to cycle through
the phase input file data, and to modify the lattice accordingly at each
time step.
*/

class Controller {

protected:
    
string jobroot_;                    /**< Root name for all output files */
Lattice *lattice_;                  /**< Pointer to microstructure lattice object */
KineticModel *kineticmodel_;        /**< Pointer to kinetic model object */
ThermalStrain *thermalstr_;         /**< Pointer to the finite element model object */

/**
@brief Stores the moles of independent components dissolved.

@warning This member appears to not be used and may be obsolete
@todo Investigate whether to delete this member
*/
vector<double> molesdissolved_;

double imgfreq_;                    /**< Frequency to output microstructure image */
ChemicalSystem *chemsys_;           /**< Pointer to `ChemicalSystem` object */
Solution *solut_;                   /**< Pointer to the `Solution` object */
vector<double> time_;               /**< List of simulation times for each iteration */
double statfreq_;                   /**< Frequency to output statistics */
	
private:
    
double sattack_time_;               /**< Simulation time at which to begin sulfate attack,
                                            in days */
double leach_time_;                 /**< Simulation time at which to begin leaching,
                                            in days */
long int damagecount_;              /**< Number of pixels in the lattice that are damaged */

public:
    
/**
@brief The constructor.

This is the only Controller constructor provided.  It requires that all the
auxiliary objects be defined already, including

    - The lattice object
    - The kinetic model object
    - The chemical system object (interface between GEM and THAMES
    - The solution object
    - The finite element model for tracking strain and stress

@param msh is a pointer to the already-instantiated `Lattice` object
@param km is a pointer to the already-instantiated `KineticModel` object
@param cs is a pointer to the already-instantiated `ChemicalSystem` object
@param solut is a pointer to the already-instantiated `Solution` object
@param thmstr is a pointer to the already-instantiated `ThermalStrain` object
@param parfilename is the name of the input parameter file
@param jobname is the root name to give to all output files
*/
Controller (Lattice *msh,
            KineticModel *km,
            ChemicalSystem *cs,
            Solution *solut,
            ThermalStrain *thmstr,
            const string &parfilename,
            const string &jobname);
    
/**
@brief Run a computational iteration.

This method launches one computational iteration, which includes

    - Consulting the kinetic model to determine the number of moles of
        each independent component (IC) to add or subtract from the system
    - Running the GEM thermodynamic calculation
    - Running the finite element code (optionally) to update stress and strain states
    - Updating the lattice to reflect the new microstructure

@param statfilename is the name of the file to store phase statistics
@param choice is an int flag to specify whether simulating hydration, leaching, or
sulfate attack
*/
void doCycle (const string &statfilename,
              int choice);

/**
@brief Calculate the state of the system (called by doCycle).

This method calculates the change in state of the system during a cycle, including

    - Consulting the kinetic model to get IC moles dissolved or added
    - (Optionally) determining IC moles to add from an external sulfate solution
    - Launching a thermodynamic calculation
    - (Optionally) determining the AFt saturation index for crystallization pressure
            calculations
    - Updating the microstructure
    - Outputting the microstructure phase volume fractions and other data

@param time is the simulation time [days]
@param dt is the change in simulation time used by the kinetic model [days]
@param isfirst is true iff this is the first state calculation (initialization)
*/
void calculateState (double time,
                     double dt,
                     bool isfirst);
 
/**
@brief Parse the input XML file specifying Controller parameters to use.

Controller parameters that need to be input are

    - Length of time to calculate [days]
    - Frequency to output microstructure images

@param docname is the name of the XML input file containing the Controller parameters
*/
void parseDoc (const string &docname);
    
/**
@brief Set the simulation time at which to begin sulfate attack simulation.

@param sattacktime is the simulation time to begin sulfate attack [days]
*/
void setSattack_time (const double sattacktime)
{
    sattack_time_ = sattacktime;
}

/**
@brief Get the simulation time at which to begin sulfate attack simulation.

@return the simulation time to begin sulfate attack [days]
*/
double getSattack_time (void) const
{
    return sattack_time_;
}
                                                  
};      // End of Controller class
#endif
