/**
@file Lattice.h
@brief Declaration of the Lattice class for storing the 3D microstructure

THAMES defines a Lattice class that is instantiated to a Lattice
object at the beginning of the program's execution.  The lattice defines the
three-dimensional environment within which a cement paste microstructure
exists, hydrates, and possibly deteriorates.
*/

#ifndef LATTICEH
#define LATTICEH

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <typeinfo>
#include <algorithm>
#include <cmath>
#include <map>
#include <ctime>

#include "Site.h"
#include "Isite.h"
#include "ChemicalSystem.h"
#include "Interface.h"
#include "global.h"
#include "RanGen.h"
#include "AppliedStrain.h"

const string VERSIONSTRING("Version:");
const string IMGSIZESTRING("Image_Size:");
const string IMGRESSTRING("Image_Resolution:");
const string XSIZESTRING("X_Size:");
const string YSIZESTRING("Y_Size:");
const string ZSIZESTRING("Z_Size:");

using namespace std;

/**
@class Lattice
@brief Defines and stores the 3D microstructure as a discrete lattice of voxel sites.

*/
class Lattice {

private:
    
string version_;                            /**< THAMES version for header information */
string jobroot_;                            /**< The root name for output files */
unsigned int xdim_;                         /**< Number of sites in the x dimension */
unsigned int ydim_;                         /**< Number of sites in the y dimension */
unsigned int zdim_;                         /**< Number of sites in the z dimension */
double resolution_;                         /**< Voxel edge length [micrometers] */
RanGen *rg_;                                /**< Pointer to random number generator object */
vector<Site> site_;                         /**< 1D list of Site objects (site = voxel) */
unsigned long int numsites_;                /**< Total number of sites */
const unsigned int siteneighbors_;          /**< Number of neighbor sites to a given site */
ChemicalSystem *chemsys_;                   /**< Pointer to simulation's ChemicalSystem */
Solution *solut_;                           /**< Pointer to the simulation's Solution */
AppliedStrain *FEsolver_;                   /**< Pointer to simulation's FE elastic solver */
vector<Interface> interface_;               /**< List of the different interface objects
                                                    in the microstructure */
vector<double> volumefraction_;             /**< Array of volume fractions of each 
                                                    microstructure phase */
vector<double> count_;                      /**< Number of sites of each different type */
double ettrSI_;                             /**< Current saturation index of AFt phase,
                                                    used only for sulfate attack simulation */
map<int,vector<double> > expansion_;        /**< Map of expansion strain of each voxel */
map<int,vector<int> > expansion_coordin_;   /**< Map of coordinates of sites with 
                                                    local expansion strain */
double waterchange_;                        /**< How much water must be added or subtracted
                                                    due to hydration or deterioration */
double time_;                               /**< The current simulation time [days] */
double temperature_;                        /**< The current simulation temperature [K] */
double oldtemp_;                            /**< The temperature in the previous
                                                    time step [K] */
double sattack_time_;                       /**< Simulation time at which to begin
                                                    simulation of sulfate attack [days] */
double leach_time_;                         /**< Simulation time at which to begin
                                                    simulation of leaching [days] */
double surfacearea_;                        /**< Total surface area [m<sup>2</sup>] */

bool deptheffect_;                          /**< Whether or not PNG images should have depth effect */

public:
    
/**
@brief Constructor without input microstructure file name.

This constructor simply initializes the dimensions and time to zero, sets
the temperature to the globally defined reference temperature, and
sets the lattice resolution to the globally defined reference value.

@note Not currently used in THAMES

@param cs is a pointer to the ChemicalSystem object for the simulation
@param solut is a pointer to the Solution object for the simulation
*/
Lattice (ChemicalSystem *cs,
         Solution *solut);

/**
@brief Overloaded constructor with input microstructure file name.

This constructor initializes the dimensions and time to zero, sets
the temperature to the globally defined reference temperature, and
sets the lattice resolution to the globally defined reference value.
Afterward, the input microstructure file is opened and read, so that
the voxel phase assignments can be made at each site.

@param cs is a pointer to the ChemicalSystem object for the simulation
@param solut is a pointer to the Solution object for the simulation
@param fname is the name of the file containing the microstructure data
*/
Lattice (ChemicalSystem *cs,
        Solution *solut,
        const string &fname);
    
/**
@brief Destructor.

This destructor clears out the `interface_` and `site_` vectors, and
also deletes the allocated memory for the random number generator object,
since this is the class that allocated the memory for that object.
*/
~Lattice ();
    
/**
@brief Set the number of sites in the x dimension.

@param x is the number of sites in the x dimension
*/
void setXDim (const unsigned int x)
{
    xdim_ = x;
    numsites_ = (long int)(xdim_ * ydim_ * zdim_);
}

/**
@brief Get the number of sites in the x dimension.

@return the number of sites in the x dimension
*/
unsigned int getXDim () const
{
    return xdim_;
}

/**
@brief Set the number of sites in the y dimension.

@param y is the number of sites in the y dimension
*/
void setYDim (const unsigned int y)
{
    ydim_ = y;
    numsites_ = (long int)(xdim_ * ydim_ * zdim_);
}

/**
@brief Get the number of sites in the y dimension.

@return the number of sites in the y dimension
*/
unsigned int getYDim () const
{
    return ydim_;
}

/**
@brief Set the number of sites in the z dimension.

@param z is the number of sites in the z dimension
*/
void setZDim (const unsigned int z)
{
    zdim_ = z;
    numsites_ = (long int)(xdim_ * ydim_ * zdim_);
}

/**
@brief Get the number of sites in the z dimension.

@return the number of sites in the z dimension
*/
unsigned int getZDim () const
{
    return zdim_;
}
    
/**
@brief Get the total number of lattice sites.

The lattice is rectangular, so the total number of sites is
`xdim_ * ydim_ * zdim_`, but we store this value as a class member to
save having to compute it multiple times.

@return the total number of lattice sites
*/
unsigned long int getNumsites () const
{
    return numsites_;
}
    
/**
@brief Get the volume fraction of a given microstructure phase.

This is simply the number of sites with a given phase divided by the
total number of sites.

@note NOT USED.

@param i is the index of the microstructure phase
@return the total number of lattice sites
*/
double getVFrac (unsigned int i)
{
    try {
        if (numsites_ == 0) {
            throw FloatException("Lattice","getVFrac",
                                 "Divide by zero (numsites_)");
        }
        return ((double)(count_.at(i))/(double)(numsites_));
    }
    catch (FloatException flex) {
        flex.printException();
        exit(1);
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","getVFrac","count_",
                        count_.size(),i);
        ex.printException();
        exit(1);
    }
}
    
/**
@brief Get the number of neighbor sites each site has.

This is simply the number of sites with a given phase divided by the
total number of sites.

@note NOT USED.

@return the number of neighbor sites each site has
*/
unsigned int getSiteneighbors () const
{
    return siteneighbors_;
}
    
/**
@brief Set the lattice resolution [micrometers].

The lattice resolution is the physical length associated with the edge
length of a site.

@param res is the lattice resolution [micrometers]
*/
void setResolution (const double res);
    
/**
@brief Get the lattice resolution [micrometers].

The lattice resolution is the physical length associated with the edge
length of a site.

@note NOT USED.

@return the lattice resolution [micrometers]
*/
double getResolution () const
{
    return resolution_;
}
    
/**
@brief Set the simulation time [days].

@note NOT USED.

@param tval is the simulation time [days]
*/
void setTime (const double tval)
{
    time_ = tval;
}
    
/**
@brief Get the simulation time [days].

@note NOT USED.

@return the simulation time [days]
*/
double getTime () const
{
    return time_;
}

/**
@brief Get the simulation time at which to start sulfate attack simulation [days].

@note NOT USED.

@return the simulation time at which to start sulfate attack [days]
*/
double getSattack_time () const
{
    return sattack_time_;
}

/**
@brief Set the simulation time at which to start sulfate attack simulation [days].

@param sattacktime is the simulation time at which to start sulfate attack [days]
*/
void setSattack_time (const double sattacktime)
{
    sattack_time_ = sattacktime;
}

/**
@brief Get the simulation time at which to start leaching simulation [days].

@note NOT USED.

@return the simulation time at which to start leaching [days]
*/
double getLeach_time () const
{
    return leach_time_;
}

/**
@brief Set the simulation time at which to start leaching simulation [days].

@param leachtime is the simulation time at which to start leaching [days]
*/
void setLeach_time (const double leachtime)
{
    leach_time_ = leachtime;
} 
    
/**
@brief Set the lattice temperature [K].

@param tmp is the temperature [K]
*/
void setTemperature (const double tmp)
{
    temperature_ = tmp;
}
    
/**
@brief Get the lattice temperature [K].

@return the temperature [K]
*/
double getTemperature () const
{
    return temperature_;
}
    
/**
@brief Get the version of THAMES

@note NOT USED.

@return the version number as a string
*/
const string &getVersion () const
{
    return version_;
}
   
/**
@brief Set the root name for simulation output files.

@param jobname is the root name for simulation output files
*/
void setJobroot (string jobname)
{
    jobroot_ = jobname;
}
 
/**
@brief Add a site at location (xp,yp,zp) to the lattice.

The site is checked for valid coordinates.  If valid a new Site object
is created and pushed back onto the class's `site_` vector.

@param xp is the x coordinate of the site to add
@param yp is the y coordinate of the site to add
@param zp is the z coordinate of the site to add
*/
void addSite (const unsigned int xp,
              const unsigned int yp,
              const unsigned int zp);
    
/**
@brief Get the x coordinate of a site with a given index in the 1D `site_` array.

@param i is the index of the site in the class's `site_` array
@return the x coordinate
*/
unsigned int getX (const unsigned long int i) const
{
    try {
        if (i >= site_.size()) {
            throw EOBException("Lattice","getX","site_",
                               site_.size(),(unsigned int)i);
        }
        return (site_[i].getX());
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the y coordinate of a site with a given index in the 1D `site_` array.

@param i is the index of the site in the class's `site_` array
@return the y coordinate
*/
unsigned int getY (const unsigned long int i) const
{
    try {
        if (i >= site_.size()) {
            throw EOBException("Lattice","getY","site_",
                               site_.size(),(unsigned int)i);
        }
        return (site_[i].getY());
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the z coordinate of a site with a given index in the 1D `site_` array.

@param i is the index of the site in the class's `site_` array
@return the x coordinate
*/
unsigned int getZ (const unsigned long int i) const
{
    try {
        if (i >= site_.size()) {
            throw EOBException("Lattice","getY","site_",
                               site_.size(),(unsigned int)i);
        }
        return (site_[i].getZ());
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
}
    
/**
@brief Get a site's index in the 1D `site_` array, given its (x,y,z) coordinates.

@param ix is the x coordinate of the site
@param iy is the x coordinate of the site
@param iz is the x coordinate of the site
@return the index of the site in the `site_` array
*/
unsigned long int getIndex (int ix,
                            int iy,
                            int iz) const;
    
/**
@brief Get the collection of site indices neighboring a given site.

@param sitenum is the index of the site in question
@param size is the maximum distance defining the neighborhood [sites]
@return a list of site indices for all neighbors within the maximum distance
*/
vector<unsigned long int> getNeighborhood (const unsigned long int sitenum,
                                           const int size);
    

/**
@brief Get a pointer to a Site object at a given index in the `site_` array.

@param index is the index of the Site object in the `site_` array
@return a pointer to the Site object in question
*/
Site *getSite (int index)
{
    try {
        if (index >= site_.size()) {
            throw EOBException("Lattice","getSite","site_",
                           site_.size(),(unsigned int)index);
        }
        return &site_[index];
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
}

/**
@brief Designate a site as damaged.

The site to be damaged is specified by its index in the `site_` array.

@note NOT USED.

@param index is the index of the Site object in the `site_` array
*/
void setDamage (int index)
{
    try {
        if (index >= site_.size()) {
            throw EOBException("Lattice","setDamage","site_",
                           site_.size(),(unsigned int)index);
        }
        site_[index].setDamage();
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }

}


/**
@brief Change the wmc (weighted mean curvature) of a site by a prescribed amount.

@param index is the index of the Site object in the `site_` array
@param dwmcval is the increment to add to the wmc
*/
void dWmc(int index,
          double dwmcval)
{
    try {
        if (index >= site_.size()) {
            throw EOBException("Lattice","dWmc","site_",
                           site_.size(),(unsigned int)index);
        }
        site_[index].setDamage();
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }

}

    
/**
@brief Master method to locate the interfaces for each phase in the microstructure.

*/
void findInterfaces ();
    
/**
@brief Add a prescribed number of sites of a given phase to the microstructure.

This method gets a list of all the <i>potential</i> growth sites, which
already is sorted in descending order of growth affinity.  The list is
then visited one site at a time, switching the id of the phase at the site.
Once this is done, the lists of growth sites and dissolution sites are
updated to account for the new local geometry.

@param phaseid is the id of the microstructure phase to add
@param numtoadd is the number of sites to switch to this phase
@return the actual number of sites that were changed
*/
long int growPhase (unsigned int phaseid,
                    long int numtoadd);
    
/**
@brief Remove a prescribed number of sites of a given phase from the microstructure.

This method gets a list of all the <i>potential</i> dissolution sites, which
already is sorted in descending order of <i>growth</i> affinity.  The list is
then visited in reverse order one site at a time, switching the id of the phase
at the site.  Once this is done, the lists of growth sites and dissolution
sites are updated to account for the new local geometry.

@param phaseid is the id of the microstructure phase to remove
@param numtoadd is the number of sites to switch from this phase
@return the actual number of sites that were changed
*/
long int dissolvePhase (unsigned int phaseid,
                        long int numtoadd);
    
/**
@brief Remove the water from a prescribed number of solution-filled sites.

This method constructs a list of all the <i>potential</i> void sites, based
on whether there are multiple connected solution-filled sites in a cluster.
The list is then sorted essentially by the effective pore size.  Only then
is the list visited and the prescribed number of sites switched to void.

@param numsites is the number of sites to switch to void
@return the actual number of sites that were changed
*/
long int emptyPorosity (long int numsites);
    
/**
@brief Add water to a prescribed number of empty pore sites.

This method constructs a list of all the void sites, based
on whether there are multiple connected void sites in a cluster.
The list is then sorted essentially by the effective pore size.  Only then
is the list visited and the prescribed number of sites switched to void.

@param numsites is the number of sites to switch from void to water
@return the actual number of sites that were changed
*/
long int fillPorosity(long int numsites);
    
/**
@brief Count the number of solution sites within a box centered on a given site.

@param boxsize is the linear dimension of the cubic box neighborhood
@param siteid is the index of the queried site in the `site_` array
@return the number of solution-filled sites in the box neighborhood
*/
int countBox(int boxsize,
             unsigned long int siteid);
    
/**
@brief Check whether a linear coordinate is outside the lattice boundaries.

If a given coordinate is outside the lattice boundaries, then the additive
adjustment is returned that will locate the equivalent site within the lattice,
assuming periodic boundary conditions.

@param pos is the linear coordinate to check
@param size is the dimension of the lattice in that dimension (number of sites)
@return the additive adjustment to locate the equivalent coordinate within
the lattice
*/
int checkBC (int pos,
             int size)
{
    if (pos >= size) return(-size);
    if (pos < 0) return(size);
    return(0);
}
    
/**
@brief Get a pointer to the ChemicalSystem object for the simulation.

@return a pointer to the ChemicalSystem object for the simulation
*/
ChemicalSystem* getChemsys () const
{
    return chemsys_;
}
    
/**
@brief Set the phase id of a given site, specified by a pointer to the Site object.

@param s is a pointer to the Site object
@param i is the phase index to set at that site
*/
void setPhaseId (Site *s,
                 const unsigned int i)
{
    try {
          count_.at(s->getPhaseId())--;        
          s->setPhaseId(i);
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","setPhaseId","count_",
                        count_.size(),i);
        ex.printException();
        exit(1);
    }
    return;
}

/**
@brief Set the phase id of a given site, specified by the site's index number.

@param sitenum is the index of the site in the `site_` array
@param i is the phase index to set at that site
*/
void setPhaseId (const long int sitenum,
                 const unsigned int i)
{
    try {
        if (i == chemsys_->getMicid("DAMAGE")) {
          site_.at(sitenum).setPhaseId(i);
          count_.at(i) += 1;
        } else {
          count_.at(site_.at(sitenum).getPhaseId())--;
          site_.at(sitenum).setPhaseId(i);
        }
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","setPhaseId","count_",
                        count_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the phase id of a given site, specified by the site's index number.

@param sitenum is the index of the site in the `site_` array
@return the microstructure phase id at the site
*/
int getPhaseId (const long int sitenum)
{
    try {
        Site *ste;
        ste = &site_[sitenum];
        return (ste->getPhaseId());
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","getPhaseId","sitenum",
                        numsites_,sitenum);
        ex.printException();
        exit(1);
    }
}
    
/**
@brief Add a site to the list of sites where dissolution of a given phase can occur.

@param loc is a pointer to the Site object to add to the list of potential dissolution sites
@param pid is the microstructure phase id
*/
void addDissolutionSite (Site *loc,
                         unsigned int pid);
    
/**
@brief Add a site to the list of sites where growth of a given phase can occur.

@param loc is a pointer to the Site object to add to the list of potential growth sites
@param pid is the microstructure phase id
*/
void addGrowthSite (Site *loc,
                    unsigned int pid);
    
/**
@brief Remove a site from the list of sites where dissolution of a given phase can occur.

@param loc is a pointer to the Site object to remove from the list of potential dissolution sites
@param pid is the microstructure phase id
*/
void removeDissolutionSite (Site *loc,
                            unsigned int pid);
    
/**
@brief Remove a site from the list of sites where growth of a given phase can occur.

@param loc is a pointer to the Site object to remove from the list of potential growth sites
@param pid is the microstructure phase id
*/
void removeGrowthSite (Site *loc,
                       unsigned int pid);
    
/**
@brief Master method to update a microstructure during after a given time interval.

Updating the microstructure includes determining how many sites of each phase to
add and subtract from the lattice, determining which site locations will be used
to do that, and then actually causing the switches in phase id to happen at those sites.
The interfaces and lists of dissolution and growth sites are updated accordingly, too.

@param time is is the simulation time [days]
@param isfirst is true if this is the first microstructure update, false otherwise
*/
void changeMicrostructure (double time,
                           bool isfirst);
    
/**
@brief Write the 3D microstructure to a file.

The microstructure output file will indicate the phase id at each site.

@param curtime is the current time in days
@param root is the root name of the output file to create
*/
void writeLattice (double curtime, const string &root);
    
/**
@brief Write the 3D microstructure to a file.

The damage output file is binary, each site either being damaged or not.

@param curtime is the current time in days
@param root is the root name of the output file to create
*/
void writeDamageLattice (double curtime, const string &root);
    

/**
@brief Write the 3D microstructure to a png file that can be immediately rendered.

@param curtime is the current time in days
@param root is the root name of the png output file to create
*/
void writeLatticePNG (double curtime, const string &root);
    
/**
@brief Write the 3D microstructure to a png file that can be immediately rendered.

The damage output file is binary, each site either being damaged or not.

@param curtime is the current time in days
@param root is the root name of the png output file to create
*/
void writeDamageLatticePNG (double curtime, const string &root);
    
/**
@brief Create files of sequential slices of the microstructure in the x direction.

The slices are individual PPM files of 2D (y,z) microstructure slices,
written back to back, in the same file.  Once created, the files are each
converted to GIFs using a system call to Imagemagick, and then the GIFs are
converted to an animated GIF file using a system call to gifsicle.

@note NOT USED.

@warning This method currently depends on system calls
@warning This method currently depends on having Imagemagick installed
@warning This method currently depends on having gifsicle installed

@todo Remove the dependence on system calls, Imagemagick, and gifsicle

@param root is the root name of the png output file to create
*/
void makeMovie (const string &root);
    
/**
@brief Set the expansion strain components of a site specified by its index.

This function changes the strain components of a site already in the
list of expansion sites.  If the prescribed site is not already in the 
list of expansion sites, then the site will be added to that list.

@param index is the index of the site in the `site_` array
@param val is the vector of expansion strain components to set
*/
void setExpansion (int index,
                   vector<double> val)
{
    map<int,vector<double> >::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
        p->second = val;
    } else {
        expansion_.insert(make_pair(index,val));
    }
}

/**
@brief Get the expansion strain components of a site specified by its index.

@param index is the index of the site in the `site_` array
@return the vector of expansion strain components to set
*/
vector<double> getExpansion (int index)
{
    string msg;
    map<int,vector<double> >::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
        return p->second;
    } else {
        msg = "Could not find expansion_ match to index provided";
        throw EOBException("Lattice","getExpansion",msg,expansion_.size(),0);
    }
}

/**
@brief Get the saturation index of ettringite.

The saturation index is the ratio of the activity product for the assumed
dissolution reaction in the GEM database to the equilibrium value of that
activity product (<i>i.e.</i>, the equilibrium constant).

@todo Change the function name to something like getEttringiteSI.
@note NOT USED.

@return the saturation index of ettringite.
*/
double getEttrSI ()
{
    return ettrSI_;
}

/**
@brief Set the saturation index of ettringite.

The saturation index is the ratio of the activity product for the assumed
dissolution reaction in the GEM database to the equilibrium value of that
activity product (<i>i.e.</i>, the equilibrium constant).

@todo Change the function name to something like setEttringiteSI.

@param val is the saturation index of ettringite.
*/
void setEttrSI (double val)
{
    ettrSI_ = val;
    return;
}

/**
@brief Get the expansion strain components for all strained sites in the lattice.

@return the map of the strain components, keyed to the site index numbers
*/
map<int, vector<double> > getExpansion ()
{
    return expansion_;
}

/**
@brief Get the coordinates of local region for calculating expansion stress.

This gets the coordinates of the center site of a box in the lattice within which
the expansion strain is calculated in the ThermalStrain model due to local
crystallization pressure.

@todo Change the function name to something like getExpansionSiteCoordinates.

@param index is the index of a site that has crystallization pressure
@return the (x,y,z) coordinates of the site
*/
vector<int> getExpansionCoordin (int index) 
{
    string msg;
    map<int,vector<int> >::iterator p = expansion_coordin_.find(index);
    if (p != expansion_coordin_.end()) {
        return p->second;
    } else {
        msg = "Could not find expansion_coordin_ match to index provided";
        throw EOBException("Lattice","getExpansionCoordin",msg,expansion_coordin_.size(),0);
    }
}

/**
@brief Set the coordinates of local site for calculating expansion stress.

This gets the coordinates of the center site of a box in the lattice within which
the expansion strain is calculated in the ThermalStrain model due to local
crystallization pressure.

@note NOT USED (commented in Controller)

@todo Change the function name to something like setExpansionSiteCoordinates

@param index is the index of a site that has crystallization pressure
@param coordin is the (x,y,z) triple of the site's coordinates
@return the (x,y,z) coordinates of the site
*/
void setExpansionCoordin (int index,
                          vector<int> coordin)
{
    string msg;
    map<int,vector<int> >::iterator p = expansion_coordin_.find(index);
    if (p == expansion_coordin_.end()) {
        expansion_coordin_.insert(make_pair(index,coordin));
    }
}
    
/**
@brief Get the number of sites of water that must be added after a time step.

@note Currently only used in sulfate attack simulations.

@return the amount of water that must be added [site units]
*/
double getWaterchange (void) const
{
    return waterchange_;
}

/**
@brief Set the number of sites of water that must be added after a time step.

@note NOT USED.

@param the number of sites of water that must be added [site units]
*/
void setWaterchange (double waterchangeval)
{
    waterchange_ = waterchangeval;
}

/**
@brief Increment the number of sites of water that must be added after a time step.

@note NOT USED.

@param the extra number of sites of water that must be added [site units]
*/
void dWaterchange (double dwaterchangeval)
{
    waterchange_ += dwaterchangeval;
}
    
/**
@brief Implement conversion of Al-bearing phases to ettringite.

This method is intended only for simulating sulfate attack.  The method locates
all the Al-bearing phases that are driven to transform. It also calculates the volume
of free space adjacent to this site to determine whether crystallization
pressure will arise.  If so, the method calculates the crystallization
pressure and crystallization stress-free strain.  It then applies the expansion strain
so that the new stress field can be calculated by the ThermalStrain FE model object.

@todo Consider breaking this method into smaller pieces for maintenance and
readability.

@param alphaseid is the microstructure phase id of the Al-bearing phase to dissolve
@param netsitesAlphaseid is the number of Al-bearing sites to dissolve
@param ettrid is the microstructure phase id of ettringite
@param netsitesEttrid is the number of ettringite sites to grow
@return vector (na,ne) where na is the number of Al-bearing sites actually changed,
and ne is the number of ettringite sites actually grown
*/
vector<int> transform (int alphaseid,
                       int netsitesAlphaseid,
                       int ettrid,
                       int netsitesEttrid,
                       double volumeratio);

/**
@brief Set a pointer to the AppliedStrain object for the simulation.

@param elas is a pointer to the AppliedStrain object for the simulation
*/
void setFEsolver (AppliedStrain *AppliedStrainSolver)
{
    FEsolver_ = AppliedStrainSolver;
}

/**
@brief Write a contiguous subvolume of the lattice to a file.

@param fname is the file name to which the subvolume will be written
@param centerste is a pointer to the central site within the subvolume
@param size is the extent of the subvolume in each direction away from the center site
@return a list of the site indices belonging to the subvolume that was written
*/
vector<unsigned long int> writeSubVolume (string fname,
                                          Site *centerste,
                                          int size);

/**
@brief Assign isotropic expansion strain at a set of prescribed sites.

This function changes the strain components of a site already in the
list of expansion sites.  If the prescribed site is not already in the 
list of expansion sites, then the site will be added to that list.

@todo Consider changing the name of this method to applyExpansion

@param alnb is the collection of site indices to which strain will be assigned
@param exp is the isotropic expansion strain to set
*/
void applyExp (vector<unsigned long int> alnb,
               double exp);

/**
@brief Estimate the <i>internal</i> surface area of a phase with the aqueous solution.

@param phaseid is the id of the microstructure phase
@return the estimated surface area [site face units]
*/
double getSurfaceArea (int phaseid);

};      // End of Lattice class
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
