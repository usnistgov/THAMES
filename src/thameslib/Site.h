/**
@file Site.h
@brief Declaration of the Site class.

*/

#ifndef SITEH
#define SITEH

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "global.h"
#include "ChemicalSystem.h"
#include "RanGen.h"

using namespace std;

/**
@class Site
@brief Handle behavior of individual lattice sites.

The Site class manages changes in the phase id of a lattice site,
and whether or not it is a dissolution site or a growth site.  The site
can store a flag that determines whether or not it is damaged by
deterioration, and stores a list of its neighbor sites by their index
number.

*/
class Site {

protected:
    
unsigned int x_;                    /**< x-coordinate in mesh coordinate frame */
unsigned int y_;                    /**< y-coordinate in mesh coordinate frame */
unsigned int z_;                    /**< y-coordinate in mesh coordinate frame */
unsigned long int id_;              /**< Unique id in the 1D array of all sites */
unsigned int phaseid_;              /**< The microstructure phase assignment */
ChemicalSystem *chemsys_;           /**< Pointer to simulation's ChemicalSystem object */
vector<unsigned int> dissolution_;  /**< Vector of phases that can dissolve at this site */
vector<unsigned int> growth_;       /**< Vector of phases that can grow at this site */
double strfreevolume_;              /**< Stress-free volume of the site */
double truevolume_;                 /**< Actual volume of site, accounting for stress */
bool damage_;                       /**< True if site is damaged, false otherwise */
vector<Site *> nb_;                 /**< List of site ids that are neighbors to this site */

/**
@brief Ranking of potential for dissolution if the site is an interface site.

The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
potential for a site's phase to dissolve based on its local environment.  It is
kind of like mean curvature, but potentially can be weighted to account for
crystalline anisotropy.

THAMES is not currently able to quantitatively calculate wmc as defined in
the metallurgy literature:

    - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
        growth, Acta metall. mater. 40 [7] (1992) 1443--1474. 
    - Taylor, J.E., II--Mean curvature and weighted mean curvature,
        Acta metall. mater. 40 [7] (1992) 1475--1485.

Instead, THAMES uses the digital template method to calculate a quantity that
is roughly linearly proportional to mean curvature:

    - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
        methods for computing interfacial mean curvature, Comput. Mater. Sci.
        4 (1995) 103--116.

This provides a ranking of dissolution potential only.
*/
double wmc_; 
double expstrain_;                  /**< Assigned expansion strain by phase
                                            constrained transformation or an
                                            applied load */

public:
    
/**
@brief Default constructor.

@note NOT USED.
*/
Site ();

/**
@brief Overloaded constructor.

This constructor takes arguments for all the member variables so that they can
be assigned during construction.

The (x,y,z) coordinates and the (x,y,z) dimensions of the lattice are used
to calculate and assign a unique index number for the site in a 1D array
for the lattice.

@param xp is the x-coordinate of the site in the mesh coordinate frame
@param yp is the y-coordinate of the site in the mesh coordinate frame
@param zp is the z-coordinate of the site in the mesh coordinate frame
@param xs is the number of sites in the x dimension of the lattice
@param ys is the number of sites in the y dimension of the lattice
@param zs is the number of sites in the z dimension of the lattice
@param neigh is the number of adjacent sites to be considered as neighbors
@param csys is a pointer to the simulation's ChemicalSystem object
*/
Site (unsigned int xp,
      unsigned int yp,
      unsigned int zp,
      unsigned int xs,
      unsigned int ys,
      unsigned int zs,
      unsigned int neigh,
      ChemicalSystem *csys);
    
/**
@brief Get a pointer to a given site in the site's neighborhood.

The neighbors are stored in a 1D vector of site index values.  The construction
of the neighbor table happens in the Lattice class.

@param pos is the index of the neighbor in the neighbor table
@return a pointer to the neighboring site
*/
Site *nb (const unsigned int pos) const 
{
    if (pos >= nb_.size()) {
        throw EOBException("Site","nb","nb_",nb_.size(),pos);
    }
    return nb_[pos];
}

/**
@brief Get the size of the neighbor table (number of neighbors).

@param dist is the distance (site dimensions) away from the site to consider
@return the number of neighbors within this distance
*/
unsigned int nbSize (int dist=3) const 
{
    switch (dist) {
        case 0:
            return 1;
            break;
        case 1:
            return NUM_NEAREST_NEIGHBORS;
            break;
        case 2:
            return (NUM_NEAREST_NEIGHBORS + NUM_SECONDNEAREST_NEIGHBORS);
            break;
        default:
            return nb_.size();
            break;
    }
    return nb_.size();
}

/**
@brief Set a neighbor position with a particular site pointer.

This method is not used to push neighbors onto the neighbor vector.  The creation
of the neighbor table happens elsewhere in the Lattice class.

@param i is the index in the allocated neighbor table.
@param neigh is a pointer to the site that is to be assigned to index i of the neighbor table
*/
void setNb (unsigned int i,
            Site *neigh) 
{
    if (i >= nb_.size()) throw EOBException("Site","setNb","nb_",nb_.size(),i);
    nb_[i] = neigh;
    return;
}
    
/**
@brief Get the index number of the site (position in the 1D Lattice vector).

@return the index number of the site
*/
unsigned long int getId () const
{
    return id_;
}

/**
@brief Get the microstructure phase id number assigned to the site.

@return the microstructure phase id number
*/
unsigned int getPhaseId () const
{
    return phaseid_;
}

/**
@brief Set the microstructure phase id number assigned to the site.

@param pid is the microstructure phase id number to assign
*/
void setPhaseId (const unsigned int pid)
{
    phaseid_ = pid;
    strfreevolume_ = 1.0;
    return;
}
    
/**
@brief Get the x-coordinate of the site in the mesh coordinate frame.

@return the x-coordinate of the site in the mesh coordinate frame
*/
unsigned int getX () const
{
    return ((unsigned int)x_);
}
    
/**
@brief Get the y-coordinate of the site in the mesh coordinate frame.

@return the y-coordinate of the site in the mesh coordinate frame
*/
unsigned int getY () const
{
    return ((unsigned int)y_);
}
    
/**
@brief Get the z-coordinate of the site in the mesh coordinate frame.

@return the z-coordinate of the site in the mesh coordinate frame
*/
unsigned int getZ () const
{
    return ((unsigned int)z_);
}
    
/**
@brief Get the "weighted mean curvature" of the site.

The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
potential for a site's phase to dissolve based on its local environment.  It is
kind of like mean curvature, but potentially can be weighted to account for
crystalline anisotropy.

THAMES is not currently able to quantitatively calculate wmc as defined in
the metallurgy literature:

    - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
        growth, Acta metall. mater. 40 [7] (1992) 1443--1474. 
    - Taylor, J.E., II--Mean curvature and weighted mean curvature,
        Acta metall. mater. 40 [7] (1992) 1475--1485.

Instead, THAMES uses the digital template method to calculate a quantity that
is roughly linearly proportional to mean curvature:

    - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
        methods for computing interfacial mean curvature, Comput. Mater. Sci.
        4 (1995) 103--116.

This provides a ranking of dissolution potential only.

@return the relative potential for dissolution at this site
*/
double getWmc (void) const
{
    return wmc_;
}

/**
@brief Set the "weighted mean curvature" of the site.

@param wmcval is the value of wmc_ to assign to the site
*/
void setWmc (double wmcval)
{
    wmc_ = wmcval;
}

/**
@brief Increment the "weighted mean curvature" of the site.

@param dwmcval is the increment to make to the existing value of wmc
*/
void dWmc (double dwmcval)
{
    wmc_ += dwmcval;
}

/**
@brief Calculate the "weighted mean curvature" of the site.

The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
potential for a site's phase to dissolve based on its local environment.  It is
kind of like mean curvature, but potentially can be weighted to account for
crystalline anisotropy.

THAMES is not currently able to quantitatively calculate wmc as defined in
the metallurgy literature:

    - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
        growth, Acta metall. mater. 40 [7] (1992) 1443--1474. 
    - Taylor, J.E., II--Mean curvature and weighted mean curvature,
        Acta metall. mater. 40 [7] (1992) 1475--1485.

Instead, THAMES uses the digital template method to calculate a quantity that
is roughly linearly proportional to mean curvature:

    - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
        methods for computing interfacial mean curvature, Comput. Mater. Sci.
        4 (1995) 103--116.

This provides a ranking of dissolution potential only.
*/
void calcWmc(void);
    
/**
@brief Designate the site as a potential dissolution site for a particular phase.

@param pid is the microstructure phase id of the phase that can dissolve at the site
*/
void setDissolutionSite (unsigned int pid)
{
    dissolution_.clear();
    dissolution_.push_back(pid);
    growth_.clear();
}

/**
@brief Designate the site as a potential growth site for a particular phase.

@param pid is the microstructure phase id of the phase that can grow at the site
*/
void setGrowthSite (unsigned int pid)
{
    vector<unsigned int>::iterator start = growth_.begin();
    vector<unsigned int>::iterator end = growth_.end();
    vector<unsigned int>::iterator p = find(start,end,pid);
    if (p == growth_.end()) growth_.push_back(pid);
    dissolution_.clear();
}

/**
@brief Remove a phase from the list of phases that can dissolve from the site.

Only one phase can be at a site at a time, so the list of dissolution sites
is guaranteed to have only one member (if there is a solid phase there) or
zero members (if it is a water of void site).

@param pid is the id of the microstructure phase to remove
*/
void removeDissolutionSite (unsigned int pid)
{
    if (dissolution_.size() > 0 && dissolution_[0] == pid) dissolution_.clear();
}

/**
@brief Remove a phase from the list of phases that can grow at the site.

If a site is currently occupied by aqueous solution, then any of several solid
phases may be able to grow there.  Therefore, the list of phases in the `growth_`
vector can have multiple members.

@param pid is the id of the microstructure phase to remove
*/
void removeGrowthSite (unsigned int pid)
{
    vector<unsigned int>::iterator start = growth_.begin();
    vector<unsigned int>::iterator end = growth_.end();
    vector<unsigned int>::iterator p = find(start,end,pid);
    if (p != growth_.end()) growth_.erase(p);
}

/**
@brief Get the entire list of all phases that can grow at the site.

@return the list of ids of all microstructure phases that can grow at the site
*/
vector<unsigned int> getGrowthPhases () const
{
    return growth_;
}

/**
@brief Get the entire list of all phases that can grow at the site.

Only one phase can be at a site at a time, so the list of dissolution sites
is guaranteed to have only one member (if there is a solid phase there) or
zero members (if it is a water of void site).

@note NOT USED.

@return the list of ids of all microstructure phases that can grow at the site
*/
vector<unsigned int> getDissolutionPhases () const
{
    return dissolution_;
}
    
/**
@brief Find out if the site is designated as damaged by some kind of deterioration.

@todo Change the method name to isDamaged.

@return true if the site is damaged (damage_ flag set), or false otherwise
*/
bool IsDamage ()
{
    return damage_;
}

/**
@brief Designate the site as damaged by some kind of deterioration.

*/
void setDamage ()
{
    damage_ = true;
}
    
/**
@brief Set the stress-free volume of the site.

@note NOT USED.
@todo Change the method name to something like setStressfreevolume.

@param vol is the stress-free volume of the site to assign, normalized by reference site volume
*/
void setStrfreevolume (double vol)
{
    if (vol < 0) {
        cout << "in the setStrfreevolume function...volume should not be negative."
             << endl;
        exit(1);
    } else {
        strfreevolume_ = vol;
    }
}

/**
@brief Get the stress-free volume of the site.

@note NOT USED.
@todo Change the method name to something like getStressfreevolume.

@return the stress-free volume of the site, normalized by reference site volume
*/
double getStrfreevolume ()
{
    return strfreevolume_;
}

/**
@brief Get the true volume of the site.

@return the actual volume of the site, normalized by the strain-free site volume
*/
double getTruevolume ()
{
    return truevolume_;
}

/**
@brief Set the true volume of the site.

@param vol is the actual volume of the site, normalized by the strain-free site volume
*/
void setTruevolume (double vol)
{
    if (vol < 0) {
        cout << "in the setTruevolume function...volume should not be negative." 
             << endl;
        exit(1);
    } else {
        truevolume_ = vol;
    }
    return;
}

/**
@brief Set the expansion strain at the site.

@param val is the isotropic (<i>i.e.</i>, volumetric) strain at the site
*/
void setExpansionStrain (double val) 
{
    if (val > expstrain_) {
        expstrain_ = val; 
    } 
    return;
}

/**
@brief Get the expansion strain at the site.

@return the isotropic (<i>i.e.</i>, volumetric) strain at the site
*/
double getExpansionStrain ()
{
    return expstrain_;
}
    
/**
@brief One site is equal to another iff their wmc values are equal.

This kind of function is used to provide a comparison for sorting a list of Site objects
*/
friend bool operator==(const Site &s1, const Site &s2)
{
    return (s1.getWmc() == s2.getWmc());
}
    
/**
@brief One site is less than another iff its wmc value is less.

This kind of function is used to provide a comparison for sorting a list of Site objects
*/
friend bool operator<(const Site &s1, const Site &s2)
{
    return (s1.getWmc() < s2.getWmc());
}

};      // End of the Site class

#endif
