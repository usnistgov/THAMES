/**
@file Solution.h
@brief Declare the Solution class.

@note Perhaps this is really a derived class from the ChemicalSystem.
@todo Investigate whether to make this a derived class.

*/
#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#include "GEM3K/node.h"
#include "GEM3K/io_arrays.h"

#include "global.h"

using namespace std;

/**
@class Solution
@brief Declares the Solution class.

This class keeps track of the composition of the aqueous solution, and
has a lot of member variables that enable direct communication with the
GEM3K library.

The solution is kept separate from the ChemicalSystem because we
often want to just know the driving force (saturation index) for
dissolution or precipitation of solid phases.  Therefore we may
want to set the upper limit of solid phases to a very low number
approximating zero.

@note This class's members consist almost entirely of a subset of the
members of the ChemicalSystem class.

@todo Find out why we should not make the ChemicalSystem object a member
of this class and then access the members that way.
*/
class Solution {

private:

TNode *node_;                       /**< Pointer to the GEM3K Tnode object */

double T_;                          /**< System-wide temperature [K] */
double P_;                          /**< System-wide pressure [Pa] */
double Vs_;                         /**< System total volume [m<sup>3</sup>] */
double Ms_;                         /**< System total mass [kg] */
double Gs_;                         /**< System total Gibbs energy [J] */
double Hs_;                         /**< System total enthalpy [J] */
double ionicstrength_;              /**< Solution ionic strength [mol/kgw] */
double pH_;                         /**< Solution pH */
double pe_;                         /**< Solution pe */
double Eh_;                         /**< Solution Eh [volts] */

long int nodehandle_;               /**< integer flag used to identify a node */
long int nodestatus_;               /**< integer flag used to identify node's status */
long int iterdone_;                 /**< number of iterations performed in the most
                                                recent GEM calculation on the node */

unsigned int ICnum_;                /**< Number of independent components (IC) */
unsigned int DCnum_;                /**< Number of dependent components (DC) */
unsigned int phasenum_;             /**< Number of GEM phases in the CSD */
unsigned int solutionphasenum_;     /**< Number of GEM solution phases in the CSD;
                                              solution phases are non-stoichiometric */

double *ICmoles_;                   /**< List of number of moles of each IC in system */
double *ICresiduals_;               /**< List of errors in IC moles for mass balance */
double *ICchempot_;                 /**< List of chemical potentials of each IC, in
                                            the GEM dual solution */
double *DCmoles_;                   /**< List of moles of each DC */
double *DCactivitycoeff_;           /**< List of activity coefficients for each DC */

double *phasemoles_;                /**< List of moles of each phase in the system */
double *phasemass_;                 /**< List of mass of each phase in the system [kg] */
double *phasevolume_;               /**< List of volume of each phase in the system
                                            [m<sup>3</sup> */

double *phasestoich_;               /**< List of amount of moles of each IC in a
                                                GEM CSD phase (pointer form) */
/**
@brief Solid stoichiometry list for communicating with GEM-IPM.

@warning Not sure how this variable is used
*/
double *solidstoich_;  

double *carrier_;                   /**< List of moles of carrier (solvent) in
                                            multicomponent asymmetric phases */
double *surfacearea_;               /**< List of specific surface area of each
                                            phase, in m<sup>2</sup>/kg */
double *DCupperlimit_;              /**< List of upper bound on moles of each DC */
double *DClowerlimit_;              /**< List of lower bound on moles of each DC,
                                            generally non-zero for numerical
                                            stability of the thermodynamic
                                            calculations */
vector<string> ICname_;             /**< Names of ICs in the GEM CSD */
vector<string> DCname_;             /**< Names of DCs in the GEM CSD */
vector<string> phasename_;          /**< Names of phases in the GEM CSD */

/**
@brief Saturation index of each phase in the GEM CSD.

The departure from equilibrium between a given solid phase and an aqueous solution
is characterized by the saturation index, `SI_`, which is defined as the activity
product for the dissolution reaction divided by the equilibrium value of that activity
product.  This variable stores the current SI for each solid phase in the GEM CSD.
*/
vector<double> SI_;

double crystrain_;                  /**< Assigned strain from FE model */

public:

/**
@brief Constructor.

The class members can all be assigned by reading GEM3K input files,
which are passed to the constructor.

@param GEMfilename is the name of the GEM DCH file
@param GEMdbrname is the name fo the GEM DBR (data bridge) file
*/
Solution (const string &GEMfilename,
          const string &GEMdbrname);

/**
@brief Destructor.

A destructor is needed to free the memory allocated by the constructor.
*/
~Solution ();

/**
@brief Calculate the thermodynanmic equilibrium state of the solution.

This method calculates the thermodynamic equilibrium state of the solution
under the assumption that no solids are allowed to precipitate or dissolve.
As a result, GEM3K calculates and stores the saturation indices with respect
to each dependent component in the GEM chemical system definition (CSD).

@todo Find out how GEM knows not to precipitate or dissolve in this version.

@param isfirst is true if this is the first equilibrium calculation, false otherwise
*/
void calculateState (bool isfirst);

/**
@brief Get the list of multicomponent phase volumes.

@return a pointer to the array of phase volumes [m<sup>3</sup>]
*/
double *getPhasevolume () const
{
  return phasevolume_; 
}

/**
@brief Get the volume of one multicomponent phase.

@param idx is the index of the phase in the array of phase volumes
@return the volume [m<sup>3</sup>]
*/
double getPhasevolume (const unsigned int idx)
{
  if (idx >= phasenum_) {
      cout << "idx beyonds the range of phasenum_" << endl;
      exit(1);
  }
  return phasevolume_[idx];
}

/**
@brief Get the total number of phases in the thermodynamic system.

@return the total number of defined phases
*/
int getPhasenum ()
{
    return phasenum_;
}

/**
@brief Get the name of a phase specified by its index in the phase array.

@param idx is the index of the phase in the array of phase names
@return the phase name
*/
string &getPhasename (const unsigned int idx)
{
    if (idx >= phasename_.size()) {
        cout << "idx beyonds the range of phasename_" << endl;
        exit(1);
    }
    return (string &)phasename_[idx];
}

/**
@brief Get the list of multicomponent phase moles.

@return a pointer to the array of phase moles
*/
double *getPhasemoles () const
{
    return phasemoles_;
}

/**
@brief Get the list of all multicomponent phase names.

@return a pointer to the array of phase names
*/
vector<string> getPhasename () const
{
    return phasename_;
}

/**
@brief Get the number of iterations needed for equilibration by GEM.

@return the number of iterations executed by GEM
*/
long int getIterdone ();

/**
@brief Set the number of moles of a given independent component (IC) in the solution.

@param idx is the index of the IC in the array of ICs
@param val is the number of moles to set for this component
*/
void setICmoles (const unsigned int idx,
                 const double val)
{
    if (idx >= ICnum_) {
        cout << "index beyond the range of ICnum, so exit the program." << endl;
        exit(1);
    } else {
        ICmoles_[idx] = val;
    }
    return;
}

/**
@brief Set the number of moles of all independent components (ICs) in the solution.

@param val is an array of mole values, one for each IC
*/
void setICmoles (vector<double> val)
{
    for (int i = 0; i < ICnum_; i++) {
      ICmoles_[i] = val[i];
    }

    return;
}

/**
@brief Get the list of all independent component (IC) moles in the solution.

@return a pointer to the array of IC moles
*/
double *getICmoles ()
{
    return ICmoles_;
}

/**
@brief Get the number of moles of an independent component (IC) specified by its index.

@todo Perform bounds checking and throw an error if out of bounds.

@param idx is the index of the IC to query
@return the number of moles of that IC in the system
*/
double getICmoles (const unsigned int idx)
{
    return ICmoles_[idx];
}

/**
@brief Set the saturation index of each dependent component.
*/
void setSI ()
{
    SI_.clear();
    double *Falp;
    Falp = (node_->ppmm())->Falp;
    
    for (int i = 0; i < phasenum_; i++) {
        cout << "log10(SI) for " << phasename_[i] << " is: "
             << Falp[i] << endl;
        double si = pow(10,Falp[i]);
        SI_.push_back(si);
    }
    return;
}

/**
@brief Get the list of all saturation indices.

@return the vector of saturation indices, one for each solid phase
*/
vector<double> getSI ()
{
    return SI_;
}

/**
@brief Get the saturation index of a phase based on its index in the array.

@todo Perform bounds checking and throw an error if out of bounds.

@param phaseid is the position of the queried phase in the array of SIs
@return the saturation index of that phase
*/
double getSI (int phaseid)
{
    return SI_[phaseid];
}
   
/**
@brief Get a pointer to the GEM node doing the equilibrium calculations.

@return the GEM node doing the calculations
*/
TNode *getNode ()
{
    return node_;
}  
 
/**
@brief Calculate the site strain based on local crystallization pressure.

The calculation of local crystallization pressure happens in this method,
based on the saturation index of ettringite.  Therefore, the method is
currently restricted to problems of sulfate attack, and it probably can
be generalized somewhat.

Assuming that the saturation index \f$\beta\f$ is known for ettringite,
then the crystallization pressure should be the difference in the Laplace pressure
between the large pore entrance, \f$r_{cp}\f$, and the size of the average gel
pore, \f$r_{gp}\f$.  This pressure difference is

\f[
    p_c = 2 \gamma \left[ \frac{1}{r_{gp} - \delta} - \frac{1}{r_{cp} - \delta} \right]
\f]

where we subtract \f$\delta\f$, the liquid film thickness, so that the terms in
the denominator are the radii of curvature of the actual crystal in these two
locations.

In THAMES, we <i>assume</i> that the largest pore entrance to the gel porosity
is about half the size of a lattice site.  The usual dimension of a lattice site
in THAMES is 1 micrometer--- <b>although this is not necessary</b>--- so we assume
that \f$r_{cp} = 500\f$ nm.

The Thompson-Freundlich relation tells us the size of a crystal that is in equilibrium
with a supersaturated solution with saturation index, relative to that crystal,
of \f$\beta\f$.  The condition of equilibrium is that the chemical potential of
formula units in the crystal be equal to the chemical potential of dissolved formula
units in the solution:

\f[
    V_c \gamma \kappa_{cr} = R T \ln \frac{Q}{K} = R T \ln \beta
\f]

where \f$V_c\f$ is the stress-free molar volume of the crystal, \f$\gamma\f$ is the
surface energy of the crystal-liquid interface, \f$\kappa_{cr}\f$ is the mean
curvature of the crystal in equilibrium with the solution, <i>R</i> and <i>T</i>
are the gas constant and absolute temperature, respectively, and <i>Q</i> and
<i>K</i> are the activity product and equilibrium constant of the crystal dissociation
reaction.  Therefore, the mean curvature is

\f[
    \kappa_{cr} = \frac{R T}{V_c \gamma} \ln \beta
\f]

And since \f$\kappa_{cr} \equiv 2/r_{cr}\f$, we have \f$r_{cr} = 2/\kappa_{cr}\f$
But \f$r_{gp}\f$ is the radius of curvature of the gel pore, not the crystal, so
we must add the liquid film thickness: \f$r_{gp} = r_{cr} + \delta\f$.

With the crystallization pressure, <i>p</i><sub>c</sub>, calculated, we must
now estimate the local strain.  Ordinarily, this could be done only with detailed
knowledge of the local pore structure.  However, to make the calculation tractable,
we adopt the assumptions of poromechanics, in which case the strain is estimated
by a volume-averaged equation,

\f[
    \epsilon_x \approx \left( \frac{1}{3 K_p} - \frac{1}{3 K_s} \right)
    \left[ \phi p_c + pl - p_{atm} \right]
\f]

where \f$\epsilon_x\f$ is the linear strain (strain is assumed to be isotropic),
<i>K</i><sub>p</sub> and <i>K</i><sub>s</sub> are the volume averaged bulk and
shear moduli of the porous body, \f$\phi\f$ is the local porosity, <i>p</i><sub>l</sub>
is the hydrostatic pressure in the liquid film, and <i>p</i><sub>atm</sub> is
atmospheric pressure.

@todo Consider renaming this function to calcCrystallizationStrain.

@param ettrSI is the saturation index of the growing crystal
@param porevolfrac is the local porosity
@param Kp is the effective bulk modulus of the porous body
@param Ks is the effective shear modulus of the porous body
@return the calculated crystallization strain
*/
double calculateCrystrain (double ettrSI,
                           double porevolfrac,
                           double Kp,
                           double Ks);

};      // End of the Solution class

#endif
