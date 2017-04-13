/**
@file Interface.h
@brief Declaration of the Interface class.

@section Introduction
THAMES reads an input table of the mass fraction of
each phase at each calculated time.  Regardless of how this
table is generated, it must be stored by THAMES to use
to modify the lattice at each time increment.

This document describes a class called `Interface`, which is
primarily used to store the list of sites of each phase that are
at an interface with one or more other phases.

@section Description 
The `Interface` class is basically an STL vector data structure with
some functions to access or modify that data.

@todo Add more exception handling besides what is in the constructor,
especially for sorting operations and removing or adding elements to vectors.
*/

#ifndef INTERFACEH
#define INTERFACEH

#include <vector>
#include <algorithm>
#include <typeinfo>

#include "Site.h"
#include "Isite.h"
#include "global.h"

using namespace std;

///
/// The ChemicalSystem class needs to be declared explicitly because
/// it is used as a member in the Interface class.
///

class ChemicalSystem;

/**
@class Declare the Interface class for handling and sorting interface voxels.

The Interface object is a collection of voxels of all the same phase that share
at least one face with a different type of phase voxel.  We have containers for
storing the list of voxels and for sorting the list in descending order of potential
for dissolution.
*/
class Interface {

private:
    
    unsigned int phaseid_;              /**< The phase id of the voxels at this interface */
    ChemicalSystem *chemsys_;           /**< The `ChemicalSystem` object for the simulation */
    vector<Isite> growth_sites_;        /**< The list of all sites eligible for
                                                adjacent growth */
    vector<Isite> dissolution_sites_;   /**< The list of sites eligible for self-dissolution */
    RanGen *rg_;                        /**< The random number generator object */

public:
    
/**
@brief The default constructor, initializing members to empty or zero values.

@note NOT USED.
*/
Interface ();

/**
@brief Overloaded constructor initializing the random number generator.

@param rg is a pointer to the random number generator object to assign
*/
Interface (RanGen *rg);

/**
@brief Overloaded constructor, initializing all members to prescribed values.

@param csys is a pointer to the ChemicalSystem object being used in the simulation
@param rg is a pointer to the random number generator
@param gv is the list of pointers to growth sites adjacent to the interface for this phase
@param dv is the list of pointers to dissolution sites of this interface
@param pid is the integer id of the phase associated with this interface
*/
Interface (ChemicalSystem *csys,
           RanGen *rg,
           vector<Site *> gv,
           vector<Site *> dv,
           unsigned int pid);
    
/**
@brief Destructor for the Interface class.

*/
~Interface ();
    

/**
@brief Gets the integer phase id associated with this interface.

@note NOT USED?

@return the integer id for the phase associated with this interface
*/
unsigned int getPhaseId () const
{
    return (unsigned int)(phaseid_);
}
    
/**
@brief Gets the list of sites where growth of this phase can occur adjacent to the interface.

@return the vector of Isite objects where growth can occur
*/
vector<Isite> getGrowthSites ()
{
    return growth_sites_;
}

/**
@brief Gets the list of sites where dissolution of this phase can occur at the interface.

@return the vector of Isite objects where dissolution can occur
*/
vector<Isite> getDissolutionSites ()
{
    return dissolution_sites_;
}

/**
@brief Gets the number of sites where growth can occur adjacent to the interface.

@note NOT USED.

@return the number of potential growth sites adjacent to the interface
*/
unsigned long int getGrowthNumSites ()
{
    return growth_sites_.size();
}

/**
@brief Gets the number of sites where dissolution can occur at the interface.

@note NOT USED.

@return the number of potential dissolution sites adjacent to the interface
*/
unsigned long int getDissolutionNumSites ()
{
    return dissolution_sites_.size();
}
    
/**
@brief Add a site to the list of sites where growth can occur adjacent to the interface.

@param loc is a pointer to the site to add to the list of growth sites
@return true if the site was added successfully, false otherwise
*/
bool addGrowthSite (Site *loc);

/**
@brief Add a site to the list of sites where dissolution can occur at the interface.

@param loc is a pointer to the site to add to the list of dissolution sites
@return true if the site was added successfully, false otherwise
*/
bool addDissolutionSite (Site *loc);
    
/**
@brief Sort the list of growth sites in descending order of potential for growth event

@param ste is the list of sites to sort for growth potential
@param pid is the phase that could grow at these sites
@return true if the list was sorted successfully, false otherwise
*/
bool sortGrowthSites (vector<Site> &ste,
                      unsigned int pid);

/**
@brief Sort the list of dissolution sites in descending order of potential for dissolution event

@param ste is the list of sites to sort for dissolution potential
@param pid is the phase that could dissolve at these sites
@return true if the list was sorted successfully, false otherwise
*/
bool sortDissolutionSites (vector<Site> &ste,
                           unsigned int pid);
    
/**
@brief Remove a site from the list of sites where growth can occur adjacent to the interface.

@todo Add possibility of verbose output.

@param loc is a pointer to the site to remove from the list of growth sites
@return true if the site was removed successfully, false otherwise
*/
bool removeGrowthSite(Site *loc);

/**
@brief Remove a site from the list of sites where dissolution can occur adjacent
to the interface.

@param loc is a pointer to the site to remove from the list of growth sites
@param verbose will generate more output to standard out if set to true
@return true if the site was removed successfully, false otherwise
*/
bool removeDissolutionSite(Site *loc, bool verbose);

};      // End of Interface class

#endif

///
/// The functions below are used to aid in comparison of one site to another, by
/// which means lists of the sites can be sorted.
///

#ifndef CMPFUNCS
#define CMPFUNCS

/**
@brief Compare two sites, returning true is the first site is "less than" the second.

The comparison is made on the basis of what THAMES loosely calls the
<i>weighted mean curvature</i>, (wmc).  A site with high wmc is a site where dissolution of
a phase is likely to occur, and growth of another phase is unlikely to occur.
Conversely, a site with a low wmc is a site where growth of a phase is likely to
occur but dissolution of a phase is unlikely to occur.

@param s1 is a pointer to the first site in the comparison
@param s2 is a pointer to the second site in the comparison
@return true if the first site has lower wmc than the second, false otherwise
*/
bool cmp (const Site *s1,
          const Site *s2);


/**
@brief Sort two sites based on their affinity for a given phase.

The comparison is made on the basis of what THAMES loosely calls the
<i>affinity</i>.  A site with high affinity is a site where growth of a phase
is more likely to occur because of an affinity between it and the interface.

@param s1 is the first site in the comparison
@param s2 is the second site in the comparison
@return true if the first site has <i>greater</i> affinity than the second, false otherwise
*/
bool affinitySort(const Isite s1, const Isite s2);

#endif
