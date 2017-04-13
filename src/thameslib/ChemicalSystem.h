/**
@file ChemicalSystem.h
@brief Declaration of the ChemicalSystem base class.

A `ChemicalSystem` object contains the data
relevant to the properties of every phase known to THAMES,
such as composition, growth characteristics, molar volume, specific
gravity, heat capacity, etc.  
This class divides matter into three major types:
   -# <b>Independent Components</b> (ICs), which are basically chemical elements.
   -# <b>Dependent Components</b> (DCs), which are composed of one or more ICs.
       DCs include ion complexes, pure condensed phases, and pure vapors.
   -# <b>Phases</b>, which are composed of one or more DCs.

The operation of THAMES depends heavily upon the design of the
`ChemicalSystem` class because the phase in each cell governs
the types of chemical reactions that can occur.  Therefore, careful attention must be
paid to making this class as useful and, at the same time, easy to use as
possible.
*/

#ifndef CHEMSYSH
#define CHEMSYSH

#include "global.h"
#include "utils.h"
#include "valid.h"
#include "Solution.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include "GEM3K/node.h"
#include "GEM3K/io_arrays.h"
#include "myconfig.h"
#include <iomanip>
#include <algorithm>

using namespace std;

#ifndef CHEMSYSDATASTRUCT
#define CHEMSYSDATASTRUCT

/**
@struct PhaseData
@brief Stores data about each phase possible in the system for ease of parsing the input files.

In THAMES, phases are identified either thermodynamically--- in the 
GEM data repository--- or microstructurally.  A microstructural phase can
be one, or a combination of more than on, thermodynamically defined phase.

The structure is defined to make it easier to parse the input file for the GEM
chemical system definitions (CSD). It is not used elsewhere in the code.  In fact,
the same members are identified as class variables in the `ChemicalSystem` class.

Most of the members have self-evident meanings:
    - `id` is the unique integer id of the microstructure phase
    - `randomgrowth` determines how much randomness is associated with the
       growth of the phase, rather than being determined only by mean curvature
    - `k2o`, `na2o`, `mgo`, and `so3` are the mass fractions of potassium,
       sodium, magnesium, and sulfur oxides dissolved within the phase.
    - `porosity` is the volume fraction of internal porosity in the phase,
       (e.g., C-S-H)
    - `red`, `green`, and `blue` are the rgb values for the color that the
       phase will have in simulated color micrographs
    - `gray` is the grayscale index the phase will have in simulated backscattered
        electron micrographs
    - `thamesname` is the name the phaes will have in the microstructure
    - `gemphasename` is a vector of the GEM (thermodynamic) phases making up the
       THAMES phase, given by name
    - `gemdcname` is a vector of the GEM DCs (dependent components) making up the
       THAMES phase, given by name
    - `gemphasedcmembers` is a vector of the GEM (thermodynamic) phase DCs making
       up the THAMES phase, given by their GEM DC id number
    - `gemphaseid` is a vector of the GEM (thermodynamic) phases making up
       the THAMES phase, given by their GEM phase id number
    - `gemdcid` is ???
    - `gtmplt` is a vector of the growth templates for the phase
    - `atmpvec` the vector for the effect of pressure on the thermodynamic parameters ???
    - `colors` is ???
*/

struct PhaseData {
    int id;
    double randomgrowth;
    double k2o;
    double na2o;
    double mgo;
    double so3;
    double porosity;
    double red;
    double green;
    double blue;
    double gray;
    string thamesname;
    vector<char * > gemphasename,gemdcname;
    vector<int> gemphasedcmembers;
    vector<int> gemphaseid;
    vector<int> gemdcid;
    vector<int> gtmplt;
    vector<int> atmpvec;
    vector<double> colors;
};
#endif

/**
@class ChemicalSystem
@brief Handles the tracking of phases and communications between GEM and THAMES.

A `ChemicalSystem` object contains the data relevant to the properties of
every phase known to a given run of THAMES, such as composition, growth
characteristics, molar volume, specific gravity, heat capacity, etc.

@section intro_sec Materials
The class divides matter into three major types:
    -# <b>Independent Components</b> (ICs), which are basically chemical
        elements or a unit of charge
    -# <b>Dependent Components</b> (DCs), which are composed of one or
        more ICs.  DCs include ion complexes, pure condensed phases, and
        pure vapors
    -# <b>Phases</b>, which are composed of one or more DCs.

@section Phases
One of the biggest challenges with ChemicalSystem is that it is an interface between
the thermodynamic representation of phases, stored in the associated GEM Chemical
System Definition (CSD), and the phases as they are represented in the 3D microstructure.
To accommodate this, ChemicalSystem keeps track of two entirely different lists of
phases, and keeps track of how the members of the list are associated with each other.
For example, the THAMES microstructure includes a phase called C-S-H, which represents
the calcium silicate hydrate product.  However, the GEM CSD defines several different
phases, which are different compositional and structural end members of a non-ideal
solid solution that makes up C-S-H.  Therefore, THAMES must keep a list of all the
GEM phases that are collectively called C-S-H.  The same thing is true for certain of the
AFm and AFt phases.

@section methods Methods
THAMES communicates frequently with the GEM-IPM library, and so a lot of its
data and methods are devoted to storing and communicating GEM data, and translating
that GEM data into THAMES data that are used to create and change microstructure.

*/

class ChemicalSystem {

Solution *solut_;                       /**< Pointer to a Solution object for the system */
unsigned int micphasenum_;              /**< Total number of material components that the
                                                microstructure can contain.  The GEM
                                                chemical system definition (CSD) will
                                                list all of the possible phases in the
                                                system, but we may want to combine two or
                                                more phases into a single microstructure
                                                component */
unsigned int micimpuritynum_;           /**< Number of impurities, as oxides, dissolved
                                                in the clinker phases that will be tracked
                                                Typically the value will be 4 (potassium,
                                                sodium, magnesium, and sulfur oxides) */
unsigned int ICnum_;                    /**< Number of independent components (IC) */
unsigned int DCnum_;                    /**< Number of dependent components (DC) */
unsigned int phasenum_;                 /**< Number of GEM phases in the CSD */
unsigned int solutionphasenum_;         /**< Number of GEM solution phases in the CSD;
                                              solution phases are non-stoichiometric */
vector<string> micphasename_;           /**< Names of phases identified in a THAMES
                                                microstructure */
vector<string> ICname_;                 /**< Names of ICs in the GEM CSD */
vector<string> DCname_;                 /**< Names of DCs in the GEM CSD */
vector<string> phasename_;              /**< Names of phases in the GEM CSD */
vector<int> micid_;                     /**< Unique ids of THAMES microstructure phases */
int c3sid_;                             /**< Specific id number for C3S phase, treated
                                                specially as one of the clinker phases as
                                                a matter of convenience in kinetic model */
int c2sid_;                             /**< Specific id number for C2S phase, treated
                                                specially as one of the clinker phases as
                                                a matter of convenience in kinetic model */
int c3aid_;                             /**< Specific id number for C3A phase, treated
                                                specially as one of the clinker phases as
                                                a matter of convenience in kinetic model */
int c4afid_;                            /**< Specific id number for C4AF phase, treated
                                                specially as one of the clinker phases as
                                                a matter of convenience in kinetic model */
int gypsumid_;                          /**< Specific id number for gypsum, treated
                                                specially as one of the soluble portland cement
                                                phases as a matter of convenience
                                                in kinetic model */
vector<double> randomgrowth_;             /**< One real number for each microstructure phase,
                                          that indicates the tendency for growth in random
                                          directions (ballistic or diffusion-limited
                                          aggregation) as opposed to compact growth */
vector<double> ICmolarmass_;              /**< One molar mass for each IC [g/mol] */
vector<double> DCmolarmass_;              /**< One molar mass for each DC [g/mol] */
vector<double> phasemolarmass_;           /**< One molar mass for each GEM phase [g/mol] */
vector<vector<int> > growthtemplate_;   /**< A list of the phases on which a given phase
                                                is allowed to grow; one list for each phase */
vector<vector<int> > affinity_;         /**< A list of the microstructure phases with which
                                                a given microstructure phase has an affinity to
                                                associate when growing */
map<int,vector<int> > micphasemembers_; /**< A list of all the CSD phase ids that are
                                                associated with a given microstructure phase */
map<int,vector<int> > micDCmembers_;     /**< A list of all the CSD DC ids that are
                                                associated with a given microstructure phase */
map<int,vector<int> > phaseDCmembers_;   /**< A list of all the CSD DC ids that are
                                                associated with a given CSD phase */

/**
@brief Volume fraction of each GEM CSD phase associated with a THAMES phase.

@warning This variable may not be used
*/
map<int,vector<double> > micphasemembersvolfrac_;

vector<double> porosity_;                 /**< The internal porosity of a given phase,
                                                such as C-S-H (dimensionless) */
vector<double> k2o_;                      /**< Mass fraction of K<sub>2</sub>O dissolved in
                                                each clinker phase, in units of
                                                g per 100 g of the phase */
vector<double> na2o_;                     /**< Mass fraction of Na<sub>2</sub>O dissolved in
                                                each clinker phase, in units of
                                                g per 100 g of the phase */
vector<double> mgo_;                      /**< Mass fraction of MgO dissolved in
                                                each clinker phase, in units of
                                                g per 100 g of the phase */
vector<double> so3_;                      /**< Mass fraction of SO<sub>3</sub> dissolved in
                                                each clinker phase, in units of
                                                g per 100 g of the phase */
vector<double> grayscale_;                /**< A number on [0,255] giving the relative
                                                grayscale brightness of the THAMES
                                                phases in a backscattered electron image */
vector<vector<double> > color_;           /**< A list of <r,g,b> values specifying the
                                                color of the THAMES phases in a false
                                                color micrograph */
map<string,int> micidlookup_;           /**< Map that returns the vector index of the
                                                microstructure phase name */
map<string,int> ICidlookup_;            /**< Map that returns the vector index of the
                                                IC name */
map<string,int> DCidlookup_;            /**< Map that returns the vector index of the
                                                DC name */
map<string,int> phaseidlookup_;         /**< Map that returns the vector index of the
                                                GEM CSD phase name */

map<int,int> mic2kinetic_;              /**< Map that returns the kinetic model id of
                                                a given microstructure phase id */
map<int,int> kinetic2mic_;              /**< Map that returns the microstructure phase id
                                                of a given kinetic model phase id */
map<int,int> mic2thermo_;               /**< Map that returns the GEM CSD phase id
                                                of a given microstructure phase id */
map<int,int> thermo2mic_;               /**< Map that returns the microstructure phase id
                                                of a given GEM CSD phase id */
map<int,vector<int> > mic2phase_;       /**< Map that returns the GEM CSD phase for
                                                a given microstructure phase */
map<int,vector<int> > mic2DC_;          /**< Map that returns the GEM DC for
                                                a given microstructure phase */
vector<vector<double> > DCstoich_;        /**< List of amount of moles of each IC in a DC */
vector<int> kineticphase_;              /**< List of GEM CSD phases that are NOT under
                                                thermodynamic control (clinker, gypsum, etc) */
vector<int> thermophase_;               /**< List of GEM CSD phases that ARE under
                                                thermodynamic control */
double *phasestoich_;                   /**< List of amount of moles of each IC in a
                                                GEM CSD phase (pointer form) */
/**
@brief Solid stoichiometry list for communicating with GEM-IPM.

@warning Not sure how this variable is used
*/
double *solidstoich_;           

/**
@brief Solution stoichiometry list for communicating with GEM-IPM.

@warning Not sure how this variable is used
*/
double *solutphasestoich_;

/**
@brief Solution solid (?) stoichiometry list for communicating with GEM-IPM.

@warning Not sure how this variable is used
*/
double *solutsolidstoich_;

vector<vector<double> > vphasestoich_;    /**< List of amount of moles of each IC in
                                                a given GEM CSD phase (vector form) */

double *ICmoles_;                       /**< List of number of moles of each IC in system */
double *ICresiduals_;                   /**< List of errors in IC moles for mass balance */
double *ICchempot_;                     /**< List of chemical potentials of each IC, in
                                                the GEM dual solution */
double *DCmoles_;                       /**< List of moles of each DC */
double *DCactivitycoeff_;               /**< List of activity coefficients for each DC */

double *phasemoles_;                    /**< List of moles of each phase in the system */
double *phasemass_;                     /**< List of mass of each phase in the system */
double *phasevolume_;                   /**< List of volume of each phase in the system */
double *solutphasemoles_;               /**< List of moles of each solution phase in
                                                the system */
double *solutphasemass_;                /**< List of mass of each solution phase in
                                                the system */
double *solutphasevolume_;              /**< List of volume of each phase in the system */
double *ophasemoles_;                   /**< List of moles of each phase in the system
                                                in the previous time step */
double *ophasemass_;                    /**< List of mass of each phase in the system
                                                in the previous time step */
double *ophasevolume_;                  /**< List of volume of each phase in the system
                                                in the previous time step */
double *carrier_;                       /**< List of moles of carrier (solvent) in
                                                multicomponent asymmetric phases */
double *surfacearea_;                   /**< List of specific surface area of each
                                                phase, in m<sup>2</sup>/kg */
double *DCupperlimit_;                  /**< List of upper bound on moles of each DC */
double *DClowerlimit_;                  /**< List of lower bound on moles of each DC,
                                                generally non-zero for numerical
                                                stability of the thermodynamic
                                                calculations */

/**
@brief List of one-letter class codes for each IC, specifying its kind.

Each independent component is exactly one of six kinds in GEM database:

    - e = any chemical element except hydrogen
    - h = hydrogen
    - a = component with unspecified stoichiometry, such as Nit
    - i = isotope of a chemical element
    - f = formula unit
    - z = electric charge
*/
vector<char> ICclasscode_;
                         
/**
@brief List of one-letter class codes for each DC, specifying its kind.

Each dependent component is exactly one of 15 kinds in GEM database:

    - O = a single-component phase
    - I = ideal end member of a solid solution phase (Raoult)
    - J = minor component in a solid solution (Henry)
    - M = majority component in a solid solution (Raoult)
    - S = any aqueous species except for H<sup>+</sup> or electrons
    - T = H<sup>+</sup>
    - E = electron
    - K = surface complex represented as an aqueous species
    - W = water, H<sub>2</sub>O
    - L = other components of a solvent, such as alcohol
    - G = gas
    - V = steam
    - C = carbon dioxide, CO<sub>2</sub>
    - H = hydrogen gas, H<sub>2</sub>
    - N = nitrogen gas, N<sub>2</sub>
*/
vector<char> DCclasscode_;

/**
@brief List of one-letter class codes for each phase defined in GEM CSD, specifying its kind.

Each phase in the GEM CSD is exactly one of seven kinds:

    - a = aqueous electrolyte
    - g = mixture of gases
    - f = supercritical fluid
    - l = non-electrolyte liquid (melt)
    - x = dispersed solid with ion exchange in an aqueous system
    - s = condensed solid solution phase
    - d = dispersed multicomponent solid phase
*/
vector<char> phaseclasscode_;

double T_;                              /**< System-wide temperature [K] */
double P_;                              /**< System-wide pressure [Pa] */
double Vs_;                             /**< System total volume [m<sup>3</sup>] */
double Ms_;                             /**< System total mass [kg] */
double Gs_;                             /**< System total Gibbs energy [J] */
double Hs_;                             /**< System total enthalpy [J] */
double ionicstrength_;                  /**< Solution ionic strength [mol/kgw] */
double pH_;                             /**< Solution pH */
double pe_;                             /**< Solution pe */
double Eh_;                             /**< Solution Eh [volts] */

TNode *node_;                           /**< Pointer to a GEM3K TNode object */

long int nodehandle_;                   /**< integer flag used to identify a node */
long int nodestatus_;                   /**< integer flag used to identify node's status */
long int iterdone_;                     /**< number of iterations performed in the most
                                                recent GEM calculation on the node */

/**
@brief Time to begin sulfate solution exposure, in days.

THAMES enables one to turn of the hydration simulation (i.e., kinetic dissolution of
clinker phases according to empirical rate laws in an otherwise closed system) 
and begin a simulation of sulfate attack on that hydrated system.  This variable
designates the time (in days) at which this switch should happen.
*/
double sattack_time_;

/**
@brief Time to begin leaching simulation, in days.

THAMES enables one to turn of the hydration simulation (i.e., kinetic dissolution of
clinker phases according to empirical rate laws in an otherwise closed system) 
and begin a simulation of leaching by a low-pH solution.  This variable
designates the time (in days) at which this switch should happen.
*/
double leach_time_;

vector<double> micphasevolfrac_;        /**< Change in volume fraction of each phase
                                                in the microstructure */
vector<double> micphasemass_;           /**< Absolute mass of each microstructure phase */
vector<double> micphasemassdissolved_;  /**< Absolute mass dissolved of each microstructure
                                                phase */
vector<double> micphasevolume_;         /**< Absolute volume of each microstructure phase */

double mictotvolume_;                   /**< Absolute volume of the microstructure */
double mictotinitvolume_;               /**< Initial absolute volume of the microstructure */
double micvoidvolume_;                  /**< Absolute volume of void space in microstrucxture */
double micvoidvolfrac_;                 /**< Volune fraction of void space in microstructure */

/**
@brief Saturation index of each phase in the GEM CSD.

The departure from equilibrium between a given solid phase and an aqueous solution
is characterized by the saturation index, `SI_`, which is defined as the activity
product for the dissolution reaction divided by the equilibrium value of that activity
product.  This variable stores the current SI for each solid phase in the GEM CSD.
*/
vector<double> SI_;

public:

/**
@brief Constructor.

Only one constructor is provided, which initializes the ChemicalSystem with
all the information read from the GEM input files.

@param Solut points to the Solution object for this system
@param GEMfilename is the name of the file holding GEM input data
@param GEMdbrname is the name of the GEM data bridge file
@param Interfacefilename is the name of the file containing information about how
            to relate GEM phases to microstructure phases
*/
ChemicalSystem (Solution *Solut,
                const string &GEMfilename,
                const string &GEMdbrname,
                const string &Interfacefilename);
    
/**
@brief Copy constructor.

@param obj is the ChemicalSysteme object to copy
*/
ChemicalSystem (const ChemicalSystem &obj);
    
/**
@brief Destructor (does nothing for now).
*/
~ChemicalSystem ();
    
/**
@brief Master function for parsing an input file in XML format.

@param docname is the name of the XML document to parse
*/
void parseDoc (const string &docname);

/**
@brief Parse input about a microstructure phase from an XML document.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param numentries is the number of entries in the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parsePhase (xmlDocPtr doc,
                 xmlNodePtr cur,
                 int numentries,
                 PhaseData &phasedata);

/**
@brief Parse input about a GEM CSD phase from an XML document.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parseGEMphasedata (xmlDocPtr doc,
                        xmlNodePtr cur,
                        PhaseData &phasedata);

/**
@brief Parse input about how to render a phase in an image.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parseDisplaydata (xmlDocPtr doc,
                       xmlNodePtr cur,
                       PhaseData &phasedata);

/**
@brief Parse input about dissolved impurities within a phase.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parseImpuritydata (xmlDocPtr doc,
                        xmlNodePtr cur,
                        PhaseData &phasedata);

/**
@brief Parse input about interfaces associated with a phase.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parseInterfacedata (xmlDocPtr doc,
                         xmlNodePtr cur,
                         PhaseData &phasedata);

/**
@brief Parse input about affinity for one phase to grow on another.

@param doc points to the XML file
@param cur points to the current location within the XML file
@param phasedata holds the structure of collected phase data from the document
*/
void parseAffinitydata (xmlDocPtr doc,
                        xmlNodePtr cur,
                        PhaseData &phasedata);
    

/**
@brief Set the total number of possible microstructure phases in the system.

@note NOT USED.

@param val is the total number of possible microstructure phases (non-negative)
*/
void setMicphasenum (const unsigned int val)
{
    micphasenum_ = val;
}

/**
@brief Get the total number of possible microstructure phases in the system.

@return the total number of possible microstructure phases (non-negative)
*/
unsigned int getMicphasenum () const
{
    return micphasenum_;
}

/**
@brief Set the number of possible dissolved impurities in the system.

@note NOT USED.

@param val is the total number of possible impurities (non-negative)
*/
void setMicimpuritynum (const unsigned int val)
{
    micimpuritynum_ = val;
}

/**
@brief get the number of possible dissolved impurities in the system.

@return the total number of possible impurities (non-negative)
*/
unsigned int getMicimpuritynum() const { return micimpuritynum_; }

/**
@brief Set the number of independent components (ICs).

@note NOT USED.

@param val is the total number of independent components (non-negative)
*/
void setICnum (const unsigned int val)
{
    ICnum_ = val;
}

/**
@brief Get the number of independent components (ICs).

@return the number of independent components (non-negative)
*/
unsigned int getICnum() const
{
    return ICnum_;
}

/**
@brief Set the number of dependent components (DCs).

@note NOT USED.

@param val is the total number of dependent components (non-negative)
*/
void setDCnum(const unsigned int val)
{
    DCnum_ = val;
}

/**
@brief Get the number of dependent components (DCs).

@return the number of dependent components (non-negative)
*/
unsigned int getDCnum () const
{
    return DCnum_;
}

/**
@brief Set the number of phases in the GEM CSD.

@note NOT USED.

@param val is the total number of phases recognized in the GEM CSD (non-negative)
*/
void setPhasenum (const unsigned int val)
{
    phasenum_ = val;
}

/**
@brief Get the number of phases in the GEM CSD.

@return the number of phases in the GEM CSD (non-negative)
*/
unsigned int getPhasenum () const
{
    return phasenum_;
}

/**
@brief Set the number of solution phases in the GEM CSD.

@note NOT USED.

@param val is the total number of solution phases recognized in the GEM CSD (non-negative)
*/
void setSolutionphasenum (const unsigned int val)
{
    solutionphasenum_ = val;
}

/**
@brief Get the number of solution phases in the GEM CSD.

@note Used only in this class's copy constructor.

@return the number of solution phases in the GEM CSD (non-negative)
*/
unsigned int getSolutionphasenum () const
{
    return phasenum_;
}

/**
@brief Set the name of a microstructure phase.

@note NOT USED.

@param idx is the phase index number (non-negative)
@param str is the name to assign to the phase
*/
void setMicphasename (const unsigned int idx,
                      const string &str)
{
    if (idx >= micphasename_.size()) {
        throw EOBException("ChemicalSystem","setMicphasename","micphasename_",
                            micphasename_.size(),idx);
    }
    micphasename_[idx] = str;
    return;
}

/**
@brief Get the name of a microstructure phase.

@param idx is the phase index number (non-negative)
@return the name of the phase
*/
string &getMicphasename (const unsigned int idx)
{
    if (idx >= micphasename_.size()) {
        throw EOBException("ChemicalSystem","getMicphasename","micphasename_",
                            micphasename_.size(),idx);
    }
    return (string &)micphasename_[idx];
}

/**
@brief Get the list of all microstructure phase names.

@return the vector of microstructure phase names
*/
vector<string> getMicphasename () const
{
    return micphasename_;
}

/**
@brief Set the name of an independent component (IC).

@note NOT USED.

@param idx is the IC index number (non-negative)
@param str is the name to assign to the IC
*/
void setICname (const unsigned int idx,
                const string &str)
{
    if (idx >= ICname_.size()) {
        throw EOBException("ChemicalSystem","setICname","ICname_",
                            ICname_.size(),idx);
    }
    ICname_[idx] = str;
    return;
}

/**
@brief Get the name of an independent component (IC).

@param idx is the IC index number (non-negative)
@return the name of the IC 
*/
string &getICname (const unsigned int idx)
{
    if (idx >= ICname_.size()) {
        throw EOBException("ChemicalSystem","getICname","ICname_",
                            ICname_.size(),idx);
    }
    return (string &)ICname_[idx];
}

/**
@brief Get the list of all independent component (IC) names.

@return the vector of IC names
*/
vector<string> getICname () const
{
    return ICname_;
}

/**
@brief Set the name of a dependent component (DC).

@note NOT USED.

@param idx is the DC index number (non-negative)
@param str is the name to assign to the DC
*/
void setDCname (const unsigned int idx,
                const string &str)
{
    if (idx >= DCname_.size()) {
        throw EOBException("ChemicalSystem","setDCname","DCname_",
                            DCname_.size(),idx);
    }
    DCname_[idx] = str;
    return;
}

/**
@brief Get the name of a dependent component (DC).

@param idx is the DC index number (non-negative)
@return the name of the DC 
*/
string &getDCname (const unsigned int idx)
{
    if (idx >= DCname_.size()) {
        throw EOBException("ChemicalSystem","getDCname","DCname_",
                            DCname_.size(),idx);
    }
    return (string &)DCname_[idx];
}

/**
@brief Get the list of all dependent component (DC) names.

@return the vector of DC names
*/
vector<string> getDCname () const
{
    return DCname_;
}

/**
@brief Set the name of a phase in the GEM CSD.

@note NOT USED.

@param idx is the GEM phase index number (non-negative)
@param str is the name to assign to the phase
*/
void setPhasename (const unsigned int idx,
                   const string &str)
{
    if (idx >= phasename_.size()) {
        throw EOBException("ChemicalSystem","setPhasename","phasename_",
                            phasename_.size(),idx);
    }
    phasename_[idx] = str;
    return;
}

/**
@brief Get the name of a phase in the GEM CSD.

@param idx is the GEM phase index number (non-negative)
@return the name of the GEM phase
*/
string &getPhasename (const unsigned int idx)
{
    if (idx >= phasename_.size()) {
        throw EOBException("ChemicalSystem","getPhasename","phasename_",
                            phasename_.size(),idx);
    }
    return (string &)phasename_[idx];
}

/**
@brief Get the list of all GEM CSD phase names.

@return the vector of GEM phase names
*/
vector<string> getPhasename () const
{
    return phasename_;
}

/**
@brief Set the integer id of a microstructure phase.

@note NOT USED.

@param idx is the vector element holding the phase id (non-negative)
@param val is the non-negative integer id to assign to the phase
*/
void setMicid (const unsigned int idx,
               const int val)
{
    if (idx >= micid_.size()) {
        throw EOBException("ChemicalSystem","setMicid","micid_",
                            micid_.size(),idx);
    }
    micid_[idx] = val;
    return;
}

/**
@brief Get the integer id of a microstructure phase by its position in the id vector.

@note NOT USED.

@param idx is the element of the vector holding the phase's id (non-negative)
@return the integer id stored at that element
*/
int getMicid (const unsigned int idx)
{
    if (idx >= micid_.size()) {
        throw EOBException("ChemicalSystem","getMicid","micid_",
                            micid_.size(),idx);
    }
    return micid_[idx];
}

/**
@brief Get the integer id of a microstructure phase by its name.

@param micname is the name of the microstructure phase
@return the integer id associated with that phase name
*/
int getMicid (const string &micname)
{
    string msg;
    map<string,int>::iterator p = micidlookup_.find(micname);
    if (p != micidlookup_.end()) {
        return p->second;
    } else {
        msg = "Could not find micid_ match to " + micname;
        throw EOBException("ChemicalSystem","getMicid",msg,micid_.size(),0);
    }
}

/**
@brief Get the list of all phase ids.

@note Used only in this class's copy constructor.

@return the vector holding phase id numbers
*/
vector<int> getMicid () const
{
    return micid_;
}

/**
@brief Get the integer id of an independent component (IC) by its name.

@param icname is the name of the IC
@return the integer id associated with that IC name (non-negative)
*/
unsigned int getICid (const string &icname)
{
    string msg;
    map<string,int>::iterator p = ICidlookup_.find(icname);
    if (p != ICidlookup_.end()) {
        return p->second;
    } else {
        msg = "Could not find ICidlookup_ match to " + icname;
        throw EOBException("ChemicalSystem","getICid",msg,ICidlookup_.size(),0);
    }
}

/**
@brief Get the integer id of a dependent component (DC) by its name.

@param dcname is the name of the DC
@return the integer id associated with that DC name
*/
unsigned int getDCid (const string &dcname)
{
    string msg;
    map<string,int>::iterator p = DCidlookup_.find(dcname);
    if (p != DCidlookup_.end()) {
        return p->second;
    } else {
        msg = "Could not find DCidlookup_ match to " + dcname;
        throw EOBException("ChemicalSystem","getDCid",msg,DCidlookup_.size(),0);
    }
}

/**
@brief Get the integer id of a GEM CSD phase by its name.

@param phasename is the name of the GEM phase
@return the integer id associated with that GEM phase name 
*/
unsigned int getPhaseid (const string &phasename)
{
    string msg;
    map<string,int>::iterator p = phaseidlookup_.find(phasename);
    if (p != phaseidlookup_.end()) {
        return p->second;
    } else {
        msg = "Could not find phaseidlookup_ match to " + phasename;
        throw EOBException("ChemicalSystem","getPhaseid",msg,phaseidlookup_.size(),0);
    }
}

/**
@brief Get the integer id of a GEM phase associated with a microstructure phase id.

@param i is integer id of the microstructure phase
@param idx is element position in the vector of associated GEM phases for that
        microstructure phase
@return the integer id  of the GEM phase stored at position idx in the list
*/
unsigned int getMic2phase (const int i,
                           const unsigned int idx)
{
    string msg;
    map<int,vector<int> >::iterator p = mic2phase_.find(i);
    if (p != mic2phase_.end()) {
        if (idx < (p->second).size()) {
            return (p->second)[idx];
        } else {
            msg = "mic2phase_";
            throw EOBException("ChemicalSystem","getMic2phase",msg,(p->second).size(),idx);
        }
    } else {
        msg = "Could not find mic2phase_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2phase",msg,mic2phase_.size(),0);
    }
}

/**
@brief Get the list of integer ids of GEM phases associated with a microstructure phase id.

@param i is integer id of the microstructure phase
@return the vector of integer ids of the GEM phases for that microstructure phase
*/
vector<int> getMic2phase (const int i)
{
    string msg;
    map<int,vector<int> >::iterator p = mic2phase_.find(i);
    if (p != mic2phase_.end()) {
        return p->second;
    } else {
        msg = "Could not find mic2phase_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2phase",msg,mic2phase_.size(),0);
    }
}

/**
@brief Get the integer id of a GEM DC associated with a microstructure phase id.

@param i is integer id of the microstructure phase
@param idx is element position in the vector of associated DCs for that
        microstructure phase
@return the integer id  of the DC stored at position idx in the list
*/
int getMic2DC (const int i,
               const unsigned int idx)
{
    string msg;
    map<int,vector<int> >::iterator p = mic2DC_.find(i);
    if (p != mic2DC_.end()) {
        if (idx < (p->second).size()) {
            return (p->second)[idx];
        } else {
            msg = "mic2DC_";
            throw EOBException("ChemicalSystem","getMic2DC",msg,(p->second).size(),idx);
        }
    } else {
        msg = "Could not find mic2DC_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2DC",msg,mic2DC_.size(),0);
    }
}

/**
@brief Get the list of integer ids of DCs associated with a microstructure phase id.

@note NOT USED.

@param i is integer id of the microstructure phase
@return the vector of integer ids of the DCs for that microstructure phase
*/
vector<int> getMic2DC (const int i)
{
    string msg;
    map<int,vector<int> >::iterator p = mic2DC_.find(i);
    if (p != mic2DC_.end()) {
        return p->second;
    } else {
        msg = "Could not find mic2DC_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2DC",msg,mic2DC_.size(),0);
    }
}

/**
@brief Get the integer id of a phase in the kinetic model that is associated.
with a given microstructure phase

@param i is integer id of the microstructure phase
@return the integer id  of the associated phase in the kinetic model
*/
int getMic2kinetic (const int i)
{
    string msg;
    map<int,int>::iterator p = mic2kinetic_.find(i);
    if (p != mic2kinetic_.end()) {
        return p->second;
    } else {
        msg = "Could not find mic2kinetic_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2kinetic",msg,mic2kinetic_.size(),0);
    }
}

/**
@brief Get the integer id of a microstructure phase that is associated.
with a given phase in the kinetic model

@param i is integer id of the phase in the kinetic model
@return the integer id  of the associated microstructure phase
*/
int getKinetic2mic (const int i)
{
    string msg;
    map<int,int>::iterator p = kinetic2mic_.find(i);
    if (p != kinetic2mic_.end()) {
        return p->second;
    } else {
        msg = "Could not find kinetic2mic_ match to index provided";
        throw EOBException("ChemicalSystem","getKinetic2mic",msg,kinetic2mic_.size(),0);
    }
}

/**
@brief Get the integer id of a phase under thermodynamic control that is associated.
with a given microstructure phase.

@note NOT USED.

@param i is integer id of the microstructure phase
@return the integer id of the thermodynamically controlled phase
*/
int getMic2thermo (const int i)
{
    string msg;
    map<int,int>::iterator p = mic2thermo_.find(i);
    if (p != mic2thermo_.end()) {
        return p->second;
    } else {
        msg = "Could not find mic2thermo_ match to index provided";
        throw EOBException("ChemicalSystem","getMic2thermo",msg,mic2thermo_.size(),0);
    }
}

/**
@brief Get the integer id of a microstructure phase that is associated.
with a given thermodynamically controlled phase.

@param i is integer id of the thermodynamically controlled phase
@return the integer id of the associated microstructure phase
*/
int getThermo2mic (const int i)
{
    string msg;
    map<int,int>::iterator p = thermo2mic_.find(i);
    if (p != thermo2mic_.end()) {
        return p->second;
    } else {
        msg = "Could not find thermo2mic_ match to index provided";
        throw EOBException("ChemicalSystem","getThermo2mic",msg,thermo2mic_.size(),0);
    }
}

/**
@brief Determine if a given microstructure phase is kinetically controlled.

@param idx is integer id of the microstructure phase
@return true if the phase is kinetically controlled, false otherwise
*/
bool isKineticphase (const int idx)
{
    map<int,int>::iterator p = mic2kinetic_.find(idx);
    if (p != mic2kinetic_.end()) return true;
    return false;
}

/**
@brief Set the integer id of a GEM DC associated with a microstructure phase id.

@note NOT USED.

@param micid is integer id of the microstructure phase
@param dcid is the integer id of the DC
*/
void setMic2DC (const int micid,
                const int dcid)
{
    map<int,vector<int> >::iterator p = mic2DC_.find(micid);
    if (p == mic2DC_.end()) {
        vector<int> dumvec;
        dumvec.clear();
        dumvec.push_back(dcid);
        mic2DC_.insert(make_pair(micid,dumvec));
    } else {
        bool found = false;
        for (register unsigned int i = 0; (i < (p->second).size()) && (!found); i++) {
            if ((p->second)[i] == dcid) found = true;
        }
        if (!found) (p->second).push_back(dcid);
    }
}

/**
@brief Set the integer id of a GEM phase associated with a microstructure phase id.

@note NOT USED.

@param micid is integer id of the microstructure phase
@param phaseid is the integer id of the GEM phase
*/
void setMic2phase (const int micid,
                   const int phaseid)
{
    map<int,vector<int> >::iterator p = mic2phase_.find(micid);
    if (p == mic2phase_.end()) {
        vector<int> dumvec;
        dumvec.clear();
        dumvec.push_back(phaseid);
        mic2phase_.insert(make_pair(micid,dumvec));
    } else {
        bool found = false;
        for (register unsigned int i = 0; (i < (p->second).size()) && (!found); i++) {
            if ((p->second)[i] == phaseid) found = true;
        }
        if (!found) (p->second).push_back(phaseid);
    }
}

/**
@brief Set the integer id of a kinetically controlled phase associated with a
microstructure phase id.

@param micid is integer id of the microstructure phase
@param kineticid is the integer id of the kinetically controlled phase
*/
void setMic2kinetic (const int micid,
                     const int kineticid)
{
    map<int,int>::iterator p = mic2kinetic_.find(micid);
    if (p == mic2kinetic_.end()) {
        mic2kinetic_.insert(make_pair(micid,kineticid));
        setKinetic2mic(kineticid,micid);
    } else {
        p->second = kineticid;
    }
}

/**
@brief Set the integer id of a microstructure phase associated with a
kinetically controlled phase id.

@param kineticid is the integer id of the kinetically controlled phase
@param micid is integer id of the microstructure phase
*/
void setKinetic2mic (const int kineticid,
                     const int micid)
{
    map<int,int>::iterator p = kinetic2mic_.find(kineticid);
    if (p == kinetic2mic_.end()) {
        kinetic2mic_.insert(make_pair(kineticid,micid));
        setMic2kinetic(micid,kineticid);
    } else {
        p->second = micid;
    }
}

/**
@brief Set the integer id of a thermodynamically controlled phase associated with a
microstructure phase id.

@param micid is integer id of the microstructure phase
@param thermoid is the integer id of the thermodynamically controlled phase
*/
void setMic2thermo (const int micid,
                    const int thermoid)
{
    map<int,int>::iterator p = mic2thermo_.find(micid);
    if (p == mic2thermo_.end()) {
        mic2thermo_.insert(make_pair(micid,thermoid));
        setThermo2mic(thermoid,micid);
    } else {
        p->second = thermoid;
    }
}

/**
@brief Set the integer id of a microstructure phase associated with a
thermodynamically controlled phase id.

@param thermoid is the integer id of the thermodynamically controlled phase
@param micid is integer id of the microstructure phase
*/
void setThermo2mic (const int thermoid,
                    const int micid)
{
    map<int,int>::iterator p = thermo2mic_.find(thermoid);
    if (p == thermo2mic_.end()) {
        thermo2mic_.insert(make_pair(thermoid,micid));
        setMic2thermo(micid,thermoid);
    } else {
        p->second = micid;
    }
}

/**
@brief Get the list of integer ids of GEM phases associated with a microstructure phase name.

@note NOT USED.

@param micname is the name of the microstructure phase
@return the vector of integer ids of the DCs for that microstructure phase
*/
vector<int> getMic2phase (const string &micname)
{
    int i = (int)(getMicid(micname));
    return (getMic2phase(i));
}

/**
@brief Get the map relating every microstructure phase to a list of its associated GEM phases.

@note Used only in this class's copy constructor.

@return the map relating every microstructure phase to a list of its associated GEM phases
*/
map<int,vector<int> > getMic2phase () const
{
    return mic2phase_;
}

/**
@brief Get the map relating every kinetically controlled phase to its microstructure phase.

@note Used only in this class's copy constructor.

@return the map relating every kinetically controlled phase to its microstructure phase
*/
map<int,int> getKinetic2mic () const
{
    return kinetic2mic_;
}

/**
@brief Get the map relating every microstructure phase to a kinetically controlled phase.

@note Used only in this class's copy constructor.

@return the map relating every microstructure phase to a kinetically controlled phase
*/
map<int,int> getMic2kinetic () const
{
    return mic2kinetic_;
}

/**
@brief Get the map relating every thermodynamically controlled phase to its microstructure phase.

@note Used only in this class's copy constructor.

@return the map relating every thermodynamically controlled phase to its microstructure phase
*/
map<int,int> getThermo2mic () const
{
    return thermo2mic_;
}

/**
@brief Get the map relating every microstructure phase to a thermodynamically controlled phase.

@note Used only in this class's copy constructor.

@return the map relating every microstructure phase to a thermodynamically controlled phase
*/
map<int,int> getMic2thermo () const
{
    return mic2thermo_;
}

/**
@brief Get the map relating every microstructure phase to a list of its DCs.

@note Used only in this class's copy constructor.

@return the map relating every microstructure phase to a list of its DCs
*/
map<int,vector<int> > getMic2DC () const
{
    return mic2DC_;
}

/**
@brief Add a phase id to the list of phases that are kinetically controlled.

@param idx is the phase id number to add to the kinetic control list
*/
void setKineticphase (const unsigned int idx)
{
    string msg;
    bool found = false;
    if (idx >= phasenum_) {
        msg = "kineticphase_";
        throw EOBException("ChemicalSystem","setKineticphase",msg,phasenum_,idx);
    }
    for (register unsigned int i = 0; (i < kineticphase_.size() && !found); i++) {
        if (kineticphase_[i] == idx) found = true;
    }
    if (!found) kineticphase_.push_back(idx);
    return;
}

/**
@brief Get the microstructure id of a kinetically controlled phase.

@note NOT USED.

@param idx is the element number in the vector of kinetic phases
@return the microstructdure phase id of that element
*/
int getKineticphase (const unsigned int idx)
{
    string msg;
    if (idx >= kineticphase_.size()) {
        msg = "kineticphase_";
        throw EOBException("ChemicalSystem","getKineticphase",msg,kineticphase_.size(),idx);
    }
    return kineticphase_[idx];
}

/**
@brief Add a phase id to the list of phases that are thermodynamically controlled.

@param idx is the phase id number to add to the thermodynamic control list
*/
void setThermophase (const unsigned int idx)
{
    string msg;
    bool found = false;
    if (idx >= phasenum_) {
        msg = "thermophase_";
        throw EOBException("ChemicalSystem","setThermophase",msg,phasenum_,idx);
    }
    for (register unsigned int i = 0; (i < thermophase_.size() && !found); i++) {
        if (thermophase_[i] == idx) found = true;
    }
    if (!found) thermophase_.push_back(idx);
    return;
}

/**
@brief Get the thermodynamic phase id of a microstructure phase.

@note NOT USED.

@param idx is the element number in the vector of thermodynamic phases
@return the microstructure phase id of that thermodynamically controlled phase
*/
int getThermophase (const unsigned int idx)
{
    string msg;
    if (idx >= thermophase_.size()) {
        msg = "thermophase_";
        throw EOBException("ChemicalSystem","getThermophase",msg,thermophase_.size(),idx);
    }
    return thermophase_[idx];
}

/**
@brief Get the whole list of kinetically controlled phases.

@note Used only in this class's copy constructor.

@return the vector of ids of all microstructure phases that are kinetically controlled
*/
vector<int> getKineticphase () const
{
    return kineticphase_;
}

/**
@brief Get the whole list of thermodynamically controlled phases.

@note Used only in this class's copy constructor.

@return the vector of ids of all microstructure phases that are thermodynamically controlled
*/
vector<int> getThermophase () const
{
    return thermophase_;
}

/**
@brief Set the integer phase id of the tricalcium silicate (alite) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note NOT USED.

@param val is the integer id of the tricalcium silicate phase
*/
void setC3sid (const int val)
{
    c3sid_ = val;
}

/**
@brief Get the integer phase id of the tricalcium silicate (alite) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note Used only in this class's copy constructor.

@return the integer id of the tricalcium silicate phase
*/
int getC3sid () const
{
    return c3sid_;
}

/**
@brief Set the integer phase id of the dicalcium silicate (belite) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note NOT USED.

@param val is the integer id of the dicalcium silicate phase
*/
void setC2sid (const int val)
{
    c2sid_ = val;
}

/**
@brief Get the integer phase id of the dicalcium silicate (belite) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note Used only in this class's copy constructor.

@return the integer id of the dicalcium silicate phase
*/
int getC2sid () const
{
    return c2sid_;
}

/**
@brief Set the integer phase id of the tricalcium aluminate (aluminate)  phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note NOT USED.

@param val is the integer id of the tricalcium aluminate phase
*/
void setC3aid (const int val)
{
    c3aid_ = val;
}

/**
@brief Get the integer phase id of the tricalcium aluminate (aluminate) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note Used only in this class's copy constructor.

@return the integer id of the tricalcium aluminate phase
*/
int getC3aid () const
{
    return c3aid_;
}

/**
@brief Set the integer phase id of the tetracalcium aluminoferrite (ferrite)  phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note NOT USED.

@param val is the integer id of the tetracalcium aluminoferrite phase
*/
void setC4afid (const int val)
{
    c4afid_ = val;
}

/**
@brief Get the integer phase id of the tetracalcium aluminoferrite (ferrite) phase.

Each clinker phase is especially important because it is kinetically controlled
and therefore is considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the four clinker phases.

@note Used only in this class's copy constructor.

@return the integer id of the tetracalcium aluminoferrite phase
*/
int getC4afid () const
{
    return c4afid_;
}

/**
@brief Set the integer phase id of the calcium sulfate dihydrate (gypsum) phase.

Soluble calcium sulfates are important because they are kinetically controlled
and therefore are considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the soluble calcium sulfates.

@note NOT USED.

@param val is the integer id of the calcium sulfate dihydrate phase
*/
void setGypsumid (const int val)
{
    gypsumid_ = val;
}

/**
@brief Get the integer phase id of the calcium sulfate dihydrate (gypsum) phase.

Soluble calcium sulfates are important because they are kinetically controlled
and therefore are considered to be outside the thermodynamic system.  For convenience,
we have special id tags for the soluble calcium sulfates.

@note Used only in this class's copy constructor.

@return the integer id of the calcium sulfate dihydrate phase
*/
int getGypsumid() const
{
    return gypsumid_;
}

/**
@brief Set the random growth tendency parameter for a microstructure phase.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note NOT USED.

@param idx is the microstructure phase id
@param gmval is the random growth parameter value to set for that phase
*/
void setRandomgrowth (const unsigned int idx,
                      double gmval)
{
    if (idx < randomgrowth_.size()) {
        if (gmval > 1.0) gmval = 1.0;
        if (gmval < 0.0) gmval = 0.0;
        randomgrowth_[idx] = gmval;
    }
}

/**
@brief Get the random growth tendency parameter for a microstructure phase.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@param idx is the microstructure phase id
@return the random growth parameter value for that phase
*/
double getRandomgrowth (const unsigned int idx)
{
    string msg;
    if (idx >= randomgrowth_.size()) {
        msg = "randomgrowth_";
        throw EOBException("ChemicalSystem","getRandomgrowth",msg,randomgrowth_.size(),idx);
    }
    return randomgrowth_[idx];
}

/**
@brief Get the whole random growth tendency parameters for all microstructure phases.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note Used only in this class's copy constructor.

@return the vector of random growth parameters for all microstructure phases
*/
vector<double> getRandomgrowth () const
{
    return randomgrowth_;
}

/**
@brief Set the list of all microstructure phases on which a given microstructure phase can grow.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.

@note NOT USED.

@param idx is the microstructure phase id
@param gtvec is the list of all microstructure phases that can template growth of that phase
*/
void setGrowthtemplate (const unsigned int idx,
                        vector<int> gtvec)
{
    string msg;
    if (idx >= growthtemplate_.size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","setGrowthtemplate",msg,growthtemplate_.size(),idx);
    }
    growthtemplate_[idx] = gtvec;
    return;
}

/**
@brief Set the list of all microstructure phases on which a given microstructure phase can grow.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.
This function seeks a requested template for growth of a phase, and sets the affinity
of the growing phase for that template phase.  A high affinity for a template mimics
a low energy barrier for heterogeneous nucleation.

@note NOT USED.

@param idx is the microstructure phase id
@param jdx is the sought template microstructure phase id
@param val is the affinity value to set for that template, an integer
*/
void setGrowthtemplate (const unsigned int idx,
                        const unsigned int jdx,
                        const int val)
{
    string msg;
    if (idx >= growthtemplate_.size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","setGrowthtemplate",msg,growthtemplate_.size(),idx);
    }
    if (jdx >= growthtemplate_[idx].size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","setGrowthtemplate",
                            msg,growthtemplate_[idx].size(),jdx);
    }
    growthtemplate_[idx][jdx] = val;
    return;
}

/**
@brief Get the list of all microstructure phases on which a given microstructure phase can grow.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.
This function seeks a requested template for growth of a phase, and sets the affinity
of the growing phase for that template phase.  A high affinity for a template mimics
a low energy barrier for heterogeneous nucleation.

@param idx is the microstructure phase id
@return the list of all templates for growth of that phase
*/
vector<int> getGrowthtemplate (const unsigned int idx)
{
    string msg;
    if (idx >= growthtemplate_.size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","getGrowthtemplate",msg,growthtemplate_.size(),idx);
    }
    return growthtemplate_[idx];
}

/**
@brief Get one of the growth templates for a microstructure phase.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.
This function seeks a requested template for growth of a phase, and sets the affinity
of the growing phase for that template phase.  A high affinity for a template mimics
a low energy barrier for heterogeneous nucleation.

@note NOT USED.

@param idx is the microstructure phase id
@param jdx is the element to query in the list of templates for phase idx
@return the integer id of the growth template phase at element jdx
*/
int getGrowthtemplate (const unsigned int idx,
                       const unsigned int jdx)
{
    string msg;
    if (idx >= growthtemplate_.size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","getGrowthtemplate",msg,growthtemplate_.size(),idx);
    }
    if (jdx >= growthtemplate_[idx].size()) {
        msg = "growthtemplate_";
        throw EOBException("ChemicalSystem","getGrowthtemplate",
                            msg,growthtemplate_[idx].size(),jdx);
    }
    return growthtemplate_[idx][jdx];
}

/**
@brief Get all the growth templates of all microstructure phases.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.
This function seeks a requested template for growth of a phase, and sets the affinity
of the growing phase for that template phase.  A high affinity for a template mimics
a low energy barrier for heterogeneous nucleation.

@note Used only in this class's copy constructor.

@return the whole growth template list for all phases
*/
vector<vector<int> > getGrowthtemplate () const
{
    return growthtemplate_;
}

/**
@brief Determine if a microstructure phase is a template for another phase's growth.

Any phase can grow on its own surface, and we also allow phases to grow on the
surface of other phases, mimicking hetergeneous nucleation and growth.
This function seeks a requested template for growth of a phase, and sets the affinity
of the growing phase for that template phase.  A high affinity for a template mimics
a low energy barrier for heterogeneous nucleation.

@param gphaseid is the integer id of the microstructure phase
@param gtmpid is the id of another phase that may or may not be a template
@return true if gtmpid is a template for growth of gphaseid, false otherwise
*/
bool isGrowthtemplate (const unsigned int gphaseid,
                       int gtmpid)
{
    bool answer = false;
    for (register unsigned int i = 0;
            (i < growthtemplate_[gphaseid].size()) && (!answer); i++) {
        if (growthtemplate_[gphaseid][i] == gtmpid) answer = true;
    }
    return answer;
}

/**
@brief Set the list of affinities for growth of a microstructure phase on its templates.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note NOT USED.

@param idx is the microstructure phase id
@param avec is the list of integer affinities for growth of the phase on its templates
*/
void setAffinity (const unsigned int idx,
                  vector<int> avec)
{
    string msg;
    if (idx >= affinity_.size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","setAffinity",msg,affinity_.size(),idx);
    }
    affinity_[idx] = avec;
    return;
}

/**
@brief Set the affinity for growth of a microstructure phase on one of its templates.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note NOT USED.

@param idx is the microstructure phase id
@param jdx is the element to access in the list of growth affinities
@param val is the integer affinity to assign to the phase at that element the list
*/
void setAffinity (const unsigned int idx,
                  const unsigned int jdx,
                  const int val) {
    string msg;
    if (idx >= affinity_.size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","setAffinity",msg,affinity_.size(),idx);
    }
    if (jdx >= affinity_[idx].size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","setAffinity",
                            msg,affinity_[idx].size(),jdx);
    }
    affinity_[idx][jdx] = val;
    return;
}

/**
@brief Get the list of affinities for growth of a microstructure phase on its templates.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note NOT USED.

@param idx is the microstructure phase id of which the affinities are sought
@return the list of integer affinities for all the templates for growth of phase idx
*/
vector<int> getAffinity (const unsigned int idx)
{
    string msg;
    if (idx >= affinity_.size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","getAffinity",msg,affinity_.size(),idx);
    }
    return affinity_[idx];
}

/**
@brief Get the affinitiy for growth of a microstructure phase on one of its templates.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@param idx is the microstructure phase id of which the affinity is sought
@param jdx is the element in the list of affinities being queried
@return the affinity for the template phase associated with list element jdx
*/
int getAffinity (const unsigned int idx,
                 const unsigned int jdx)
{
    string msg;
    if (idx >= affinity_.size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","getAffinity",msg,affinity_.size(),idx);
    }
    if (jdx >= affinity_[idx].size()) {
        msg = "affinity_";
        throw EOBException("ChemicalSystem","getAffinity",
                            msg,affinity_[idx].size(),jdx);
    }
    return affinity_[idx][jdx];
}

/**
@brief Get the list of affinities for growth of every microstructure phase.

A given phase, whether hydration product or product of chemical degradation,
will generally grow in a compact form by make the growth potential highest
at sites with low mean curvature and lowest at point with high mean curvature.
The growth sites are ordered from lowest to highest mean curvature when the
random growth parameter is set to zero.  Higer values of random growth parameter
cause more severe shuffling of the ordered list of growth sites.

@note Used only in this class's copy constructor.

@return the list of all integer affinities for all the templates for growth of all phases
*/
vector<vector<int> > getAffinity () const
{
    return affinity_;
}

/**
@brief Set the potassium impurity content for a given clinker phase, on an oxide basis.

@note NOT USED.

@param idx is the id of the microstructure phase
@param ival is the mass percentage of potassium to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
void setK2o (const unsigned int idx,
             double ival) {
    if (idx >= k2o_.size()) {
        throw EOBException("ChemicalSystem","setK2o","k2o_",k2o_.size(),idx);
    }
    k2o_[idx] = ival;
    return;
}

/**
@brief Set the sodium impurity content for a given phase, on an oxide basis.

@note NOT USED.

@param idx is the id of the microstructure phase
@param ival is the mass percentage of sodium to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
void setNa2o (const unsigned int idx,
              double ival)
{
    if (idx >= na2o_.size()) {
        throw EOBException("ChemicalSystem","setNa2o","na2o_",na2o_.size(),idx);
    }
    na2o_[idx] = ival;
    return;
}

/**
@brief Set the magnesium impurity content for a given phase, on an oxide basis.

@note NOT USED.

@param idx is the id of the microstructure phase
@param ival is the mass percentage of magnesium to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
void setMgo (const unsigned int idx,
             double ival)
{
    if (idx >= mgo_.size()) {
        throw EOBException("ChemicalSystem","setMgo","mgo_",mgo_.size(),idx);
    }
    mgo_[idx] = ival;
    return;
}

/**
@brief Set the sulfur impurity content for a given phase, on an oxide basis.

@note NOT USED.

@param idx is the id of the microstructure phase
@param ival is the mass percentage of sulfur to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
void setSo3 (const unsigned int idx,
             double ival)
{
    if (idx >= so3_.size()) {
        throw EOBException("ChemicalSystem","setSo3","so3_",so3_.size(),idx);
    }
    so3_[idx] = ival;
    return;
}

/**
@brief Get the potassium impurity content for a given clinker phase, on an oxide basis.

@param idx is the id of the microstructure phase
@return the mass percentage of potassium to assign to this phase on an oxide basis, with
        units of g per 100 g of the clinker phase
*/
double getK2o(const unsigned int idx) {
    if (idx >= k2o_.size()) {
        throw EOBException("ChemicalSystem","getK2o","k2o_",k2o_.size(),idx);
    }
    return k2o_[idx];
}

/**
@brief Get the sodium impurity content for a given clinker phase, on an oxide basis.

@param idx is the id of the microstructure phase
@return the mass percentage of sodium to assign to this phase on an oxide basis, with
       units of g per 100 g of the clinker phase
*/
double getNa2o (const unsigned int idx)
{
    if (idx >= na2o_.size()) {
        throw EOBException("ChemicalSystem","getNa2o","na2o_",na2o_.size(),idx);
    }
    return na2o_[idx];
}

/**
@brief Set the magnesium impurity content for a given clinker phase, on an oxide basis.

@param idx is the id of the microstructure phase
@return the mass percentage of magnesium to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
double getMgo (const unsigned int idx)
{
    if (idx >= mgo_.size()) {
        throw EOBException("ChemicalSystem","getMgo","mgo_",mgo_.size(),idx);
    }
    return mgo_[idx];
}

/**
@brief Set the sulfur impurity content for a given clinker phase, on an oxide basis.

@param idx is the id of the microstructure phase
@param ival is the mass percentage of sulfur to assign to this phase on an oxide basis,
       with units of g per 100 g of the clinker phase
*/
double getSo3 (const unsigned int idx)
{
    if (idx >= so3_.size()) {
        throw EOBException("ChemicalSystem","getSo3","so3_",so3_.size(),idx);
    }
    return so3_[idx];
}

/**
@brief Get the potassium impurity content for all clinker phases, on an oxide basis.

@note Used only in this class's copy constructor.

@return a list of the mass percentages of potassium in all clinker phases on an oxide basis,
       with units of g per 100 g of each clinker phase
*/
vector<double> getK2o () const
{
    return k2o_;
}

/**
@brief Get the sodium impurity content for all clinker phases, on an oxide basis.

@note Used only in this class's copy constructor.

@return a list of the mass percentages of sodium in all clinker phases on an oxide basis,
       with units of g per 100 g of each clinker phase
*/
vector<double> getNa2o () const
{
    return na2o_;
}

/**
@brief Get the magnesium impurity content for all clinker phases, on an oxide basis.

@note Used only in this class's copy constructor.

@return a list of the mass percentages of magnesium in all clinker phases on an oxide basis,
       with units of g per 100 g of each clinker phase
*/
vector<double> getMgo () const
{
    return mgo_;
}

/**
@brief Get the sulfur impurity content for all clinker phases, on an oxide basis.

@note Used only in this class's copy constructor.

@return a list of the mass percentages of sulfur in all clinker phases on an oxide basis,
       with units of g per 100 g of each clinker phase
*/
vector<double> getSo3 () const
{
    return so3_;
}

/**
@brief Set the internal porosity of a microstructure phases.

A few phases, mainly C-S-H gel, have finely dispersed porosity that is not resolved
at the microstructure scale, so these phases are given a property of their average
internal porosity on a scale of one micrometer.  This has important implications
when converting mass fractions of a phase to volume fractions within a microstructure,
because the GEM thermodynamic phase may or may not include such porosity in its
value of the molar volume and density of such phases.

@warning Each version of the GEMIPM library
needs to be checked to see what the convention ensures compatibility with that library.

@note NOT USED.

@param idx is the microstructure phase with internal porosity
@param pval is the value of the volume fraction of pores to assign at the scale of one
       micrometer
*/
void setPorosity (const unsigned int idx,
                  double pval)
{
    if (idx >= porosity_.size()) {
        throw EOBException("ChemicalSystem","setPorosity","porosity_",porosity_.size(),idx);
    }
    porosity_[idx] = pval;
    return;
}

/**
@brief Get the internal porosity of a microstructure phases.

A few phases, mainly C-S-H gel, have finely dispersed porosity that is not resolved
at the microstructure scale, so these phases are given a property of their average
internal porosity on a scale of one micrometer.  This has important implications
when converting mass fractions of a phase to volume fractions within a microstructure,
because the GEM thermodynamic phase may or may not include such porosity in its
value of the molar volume and density of such phases.

@warning Each version of the GEMIPM library
needs to be checked to see what the convention ensures compatibility with that library.

@param idx is the microstructure phase with internal porosity
@return the volume fraction of the phase occupied by pores at the scale of one micrometer
*/
double getPorosity (const unsigned int idx)
{
    if (idx >= porosity_.size()) {
        throw EOBException("ChemicalSystem","getPorosity","porosity_",porosity_.size(),idx);
    }
    return porosity_[idx];
}

/**
@brief Get the list of internal porosities of all microstructure phases.

A few phases, mainly C-S-H gel, have finely dispersed porosity that is not resolved
at the microstructure scale, so these phases are given a property of their average
internal porosity on a scale of one micrometer.  This has important implications
when converting mass fractions of a phase to volume fractions within a microstructure,
because the GEM thermodynamic phase may or may not include such porosity in its
value of the molar volume and density of such phases.

@warning Each version of the GEMIPM library
needs to be checked to see what the convention ensures compatibility with that library.

@note Used only in this class's copy constructor.

@return the list of porosities of all microstructure phases at the scale of one
micrometer
*/
vector<double> getPorosity() const
{
    return porosity_;
}

/**
@brief Set the list of all GEM CSD phases that are associated with a given microstructure phase.

@note NOT USED.

@param idx is the microstructure phase in question
@param mpvec is the vector of all GEM CSD phase ids associated with the microstructure phase
*/
void setMicphasemembers (const unsigned int idx,
                         vector<int> mpvec)
{
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        p->second = mpvec;
    } else {
        micphasemembers_.insert(make_pair(idx,mpvec));
    }
}

/**
@brief Set one of the GEM CSD phases that are associated with a given microstructure phase.

@note NOT USED.

@param idx is the microstructure phase in question
@param jdx is the element position in the list of all GEM CSD phases for microstructure phase idx
@param val is GEM CSD phase id to assign to this element in the list
*/
void setMicphasemembers(const unsigned int idx,
                        const unsigned int jdx,
                        const unsigned int val)
{
    string msg;
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        if (jdx < (p->second).size()) {
            (p->second)[jdx] = val;
        } else {
            throw EOBException("ChemicalSystem","setMicphasemembers","micphasemembers_",micphasemembers_.size(),jdx);
        }
    } else {
        msg = "Could not find micphasemembers_ match to index provided";
        throw EOBException("ChemicalSystem","setMicphasemembers",msg,micphasemembers_.size(),0);
    }
}

/**
@brief Get the list of all GEM CSD phases that are associated with a given microstructure phase id number.

@note NOT USED.

@param idx is the microstructure phase in question
@return the vector of all GEM CSD phase ids associated with the microstructure phase
*/
vector<int> getMicphasemembers (const unsigned int idx)
{
    string msg;
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find micphasemembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembers",msg,micphasemembers_.size(),0);
    }
}

/**
@brief Get the list of all GEM CSD phases that are associated with a given microstructure phase name.

@note NOT USED.

@param stre is the name of the microstructure phase in question
@return the vector of all GEM CSD phase ids associated with the microstructure phase
*/
vector<int> getMicphasemembers (const string &str)
{
    string msg;
    int idx = getMicid(str);
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find micphasemembers_ match to " + str;
        throw EOBException("ChemicalSystem","getMicphasemembers",msg,micphasemembers_.size(),0);
    }
}

/**
@brief Get one of the GEM CSD phases that are associated with a given microstructure phase id.

@note NOT USED.

@param idx is the microstructure phase in question
@param jdx is the element position in the list of all GEM CSD phases for microstructure phase idx
@return the GEM CSD phase id at element jdx in the list
*/
unsigned int getMicphasemembers(const unsigned int idx,
                                const unsigned int jdx)
{
    string msg;
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        if (jdx < (p->second).size()) {
            return (p->second)[jdx];
        } else {
            throw EOBException("ChemicalSystem","getMicphasemembers","micphasemembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micphasemembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembers",msg,micphasemembers_.size(),0);
    }
}

/**
@brief Get one of the GEM CSD phases that are associated with a given microstructure phase name.

@note NOT USED.

@param str is the name of the microstructure phase in question
@param jdx is the element position in the list of all GEM CSD phases for microstructure phase idx
@return the GEM CSD phase id at element jdx in the list
*/
unsigned int getMicphasemembers (const string &str,
                                 const unsigned int jdx)
{
    string msg;
    int idx = getMicid(str);
    map<int,vector<int> >::iterator p = micphasemembers_.find(idx);
    if (p != micphasemembers_.end()) {
        if (jdx < (p->second).size()) {
            return (p->second)[jdx];
        } else {
            throw EOBException("ChemicalSystem","getMicphasemembers","micphasemembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micphasemembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembers",msg,micphasemembers_.size(),0);
    }
}

/**
@brief Get the entire list of all GEM CSD phase associations for all microstructure phases.

@note Used only in this class's copy constructor.

@return the map of all GEM CSD associations for all microstructure phases
*/
map<int, vector<int> > getMicphasemembers () const
{
    return micphasemembers_;
}

/**
@brief Get the list of volume fractions of all the GEM CSD phases for a given microstructure phase.

@note NOT USED.

@param idx is the integer id of the microstructure phase
@return the vector of volume fractions of each GEM CSD phase for this microstructure phase
*/
vector<double> getMicphasemembersvolfrac (const unsigned int idx)
{
    string msg;
    map<int,vector<double> >::iterator p = micphasemembersvolfrac_.find(idx);
    if (p != micphasemembersvolfrac_.end()) {
        return p->second;
    } else {
        msg = "Could not find micphasemembersvolfrac_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembersvolfrac_",msg,micphasemembersvolfrac_.size(),0);
    }
}

/**
@brief Get the volume fraction of one of the GEM CSD phases associated with a given microstructure phase.

@note NOT USED.

@param idx is the microstructure phase in question
@param jdx is the element position in the list of all GEM CSD phases for microstructure phase idx
@return the volume fraction of that GEM CSD phase id at element jdx in the list
*/
double getMicphasemembersvolfrac (const unsigned int idx,
                                  const unsigned int jdx)
{
    string msg;
    map<int,vector<double> >::iterator p = micphasemembersvolfrac_.find(idx);
    if (p != micphasemembersvolfrac_.end()) {
        if (jdx < (p->second).size()) {
            return (p->second)[jdx];
        } else {
            throw EOBException("ChemicalSystem","getMicphasemembersvolfrac","micphasemembersvolfrac_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micphasemembersvolfrac_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembersvolfrac",msg,micphasemembersvolfrac_.size(),0);
    }
}


/**
@brief Set the list of all dependent component (DC) ids for a microstructure phase.

@note NOT USED.

@param idx is the microstructure phase in question
@param mpvec is the list of all DCs for that phase
*/
void setMicDCmembers (const unsigned int idx,
                      vector<int> mpvec)
{
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        p->second = mpvec;
    } else {
        micDCmembers_.insert(make_pair(idx,mpvec));
    }
}

/**
@brief Set one of the DC component ids for a given microstructure phase id.

@note NOT USED.

@param idx is the microstructure phase id in question
@param jdx is the element position in the list of all GEM CSD phases for microstructure phase idx
@param val is the DC component id to set at that position in the list
*/
void setMicDCmembers (const unsigned int idx,
                      const unsigned int jdx,
                      const unsigned int val)
{
    string msg;
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        if (jdx < (p->second).size()) {
            (p->second)[jdx] = val;
        } else {
            throw EOBException("ChemicalSystem","setMicDCmembers","micDCmembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","setMicDCmembers",msg,micDCmembers_.size(),0);
    }
}

/**
@brief Get the list of DC component ids associated with a given microstructure phase id.

@note NOT USED.

@param idx is the microstructure phase id in question
@return the list of all DC component ids for that phase
*/
vector<int> getMicDCmembers (const unsigned int idx)
{
    string msg;
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find micDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicDCmembers",msg,micDCmembers_.size(),0);
    }
}

/**
@brief Get the list of DC component ids associated with a given microstructure phase name.

@note NOT USED.

@param str is the name of the microstructure phase in question
@return the list of all DC component ids for that phase
*/
vector<int> getMicDCmembers (const string &str)
{
    string msg;
    int idx = getMicid(str);
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find micDCmembers_ match to " + str;
        throw EOBException("ChemicalSystem","getMicDCmembers",msg,micDCmembers_.size(),0);
    }
}

/**
@brief Get one of the DC component ids for a given microstructure phase id.

@note NOT USED.

@param idx is the microstructure phase id in question
@param jdx is the element position in the list of all DCs for microstructure phase idx
@return the DC component id to set at that position in the list
*/
unsigned int getMicDCmembers (const unsigned int idx,
                              const unsigned int jdx)
{
    string msg;
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        if (jdx < (p->second).size()) {
            return (p->second)[jdx];
        } else {
            throw EOBException("ChemicalSystem","getMicDCmembers","micDCmembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicDCmembers",msg,micDCmembers_.size(),0);
    }
}

/**
@brief Get one of the DC component ids for a given microstructure phase name.

@note NOT USED.

@param str is the microstructure phase name in question
@param jdx is the element position in the list of all DCs for microstructure phase idx
@return the DC component id to set at that position in the list
*/
unsigned int getMicDCmembers (const string &str,
                              const unsigned int jdx)
{
    string msg;
    int idx = getMicid(str);
    map<int,vector<int> >::iterator p = micDCmembers_.find(idx);
    if (p != micDCmembers_.end()) {
        if (jdx < (p->second).size()) {
            return (p->second)[jdx];
        } else {
            throw EOBException("ChemicalSystem","getMicDCmembers","micDCmembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicDCmembers",msg,micDCmembers_.size(),0);
    }
}


/**
@brief Get DC component ids for all microstructure phases at once.

@note Used only in this class's copy constructor.

@return the map of all DC component ids for every microstructure phase
*/
map<int, vector<int> > getMicDCmembers () const
{
    return micDCmembers_;
}

/**
@brief Get the list of DC component ids associated with a given GEM CSD phase id.

@param idx is the GEM phase id in question
@return the list of all DC component ids for that GEM phase
*/
vector<int> getPhaseDCmembers (const unsigned int idx)
{
    string msg;
    map<int,vector<int> >::iterator p = phaseDCmembers_.find(idx);
    if (p != phaseDCmembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find phaseDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","getPhaseDCmembers",msg,phaseDCmembers_.size(),0);
    }
} 

/**
@brief Get the list of DC component ids associated with a given GEM CSD phase name.

@note NOT USED.

@param str is the GEM phase name in question
@return the list of all DC component ids for that GEM phase
*/
vector<int> getPhaseDCmembers (const string &str)
{
    string msg;
    int pidx = getPhaseid(str);
    map<int,vector<int> >::iterator p = phaseDCmembers_.find(pidx);
    if (p != phaseDCmembers_.end()) {
        return p->second;
    } else {
        msg = "Could not find phaseDCmembers_ match to index provided";
        throw EOBException("ChemicalSystem","getPhaseDCmembers",msg,phaseDCmembers_.size(),0);
    }
}

/**
@brief Get one of the DC component ids for a given GEM CSD phase id.

@note NOT USED.

@param idx is the GEM phase id in question
@param jdx is the element position in the list of all DCs for GEM phase idx
@return the DC component id to set at that position in the list
*/
unsigned int getPhaseDCmembers (const unsigned int idx,
                                const unsigned int jdx)
{
    string msg;
    map<int,vector<int> >::iterator p = phaseDCmembers_.find(idx);
    if (p != phaseDCmembers_.end()) {
        if (jdx < (p->second).size()) {
           return (p->second)[jdx];
        } else {
           throw EOBException("ChemicalSystem","getPhaseDCmembers","phaseDCmembers_",(p->second).size(),jdx);
        }
    } else {
        msg = "Could not find micphasemembers_ match to index provided";
        throw EOBException("ChemicalSystem","getMicphasemembers",msg,phaseDCmembers_.size(),0);
    }
}

/**
@brief Get one of the DC component ids for a given GEM CSD phase name.

@note NOT USED.

@param str is the GEM phase name in question
@param jdx is the element position in the list of all DCs for GEM phase idx
@return the DC component id to set at that position in the list
*/
unsigned int getPhaseDCmembers (const string &str,
                                const unsigned int jdx);

/**
@brief Set the number of moles of a given independent component (IC) in the system.

@param idx is the GEM IC id to set
@param val is the number of moles to assign to that IC
*/
void setICmoles (const unsigned int idx,
                 const double val)
{
    if (idx >= ICnum_) {
        throw EOBException("ChemicalSystem","setICmoles","ICmoles_",ICnum_,idx);
    }
    ICmoles_[idx] = val;
    return;
}

/**
@brief Get the number of moles of every independent component (IC) in the system.

@note Used only in this class's copy constructor.

@return a pointer to the list of number of moles assigned to each IC
*/
double *getICmoles () const
{
    return ICmoles_;
}

/**
@brief Get the number of moles of an independent component (IC) in the system.

@param idx is the GEM IC id to query
@return the number of moles assigned to that IC
*/
double getICmoles (const unsigned int idx)
{
    if (idx >= ICnum_) {
        throw EOBException("ChemicalSystem","getICmoles","ICmoles_",ICnum_,idx);
    }
    return ICmoles_[idx];
}

/**
@brief Output the number of moles of every independent component (IC) in the system.

*/
void writeICmoles ()
{
    cout << endl;
    cout << "Vector of Independent Components:" << endl;
    for (register unsigned int i = 0; i < ICnum_; i++) {
        cout << "    " << ICname_[i] << ": " << ICmoles_[i] << " mol" << endl;
    }
    cout << endl;
    cout.flush();
    return;
}

/**
@brief Set the number of moles of a given dependent component (DC) in the system.

@param idx is the GEM DC id to set
@param val is the number of moles to assign to that DC
*/
void setDCmoles (const unsigned int idx, const double val)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","setDCmoles","DCmoles_",DCnum_,idx);
    }
    DCmoles_[idx] = val;
    return;
}

/**
@brief Get the number of moles of every dependent component (DC) in the system.

@note Used only in this class's copy constructor.

@return a pointer to the list of number of moles assigned to each DC
*/
double *getDCmoles () const
{
    return DCmoles_;
}

/**
@brief Get the number of moles of an dependent component (DC) in the system.

@param idx is the GEM DC id to query
@return the number of moles assigned to that DC
*/
double getDCmoles (const unsigned int idx)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","getDCmoles","DCmoles_",DCnum_,idx);
    }
    return DCmoles_[idx];
}

/**
@brief Output the number of moles of every independent component (IC) in the system.

*/
void writeDCmoles()
{
    cout << endl;
    cout << "Vector of Dependent Components:" << endl;
    for (register unsigned int i = 0; i < DCnum_; i++) {
        cout << "    " << DCname_[i] << ": " << DCmoles_[i] << " mol" << endl;
    }
    cout << endl;
    cout.flush();
}

/**
@brief Set the number of moles of a given GEM CSD phase in the system.

@note NOT USED.

@param idx is the GEM CSD phase id to set
@param val is the number of moles to assign to that GEM phase
*/
void setPhasemoles (const unsigned int idx,
                    const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setPhasemoles","phasemoles_",phasenum_,idx);
    }
    phasemoles_[idx] = val;
    return;
}

/**
@brief Get the number of moles of every GEM CSD phase in the system.

@note Used only in this class's copy constructor.

@return a pointer to the list of number of moles assigned to each GEM phase
*/
double *getPhasemoles () const
{
    return phasemoles_;
}

/**
@brief Get the number of moles of a GEM CSD phase in the system.

@param idx is the GEM phase id to query
@return the number of moles assigned to that GEM phase
*/
double getPhasemoles (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getPhasemoles","phasemoles_",phasenum_,idx);
    }
    return phasemoles_[idx];
}

/**
@brief Set the number of moles of all GEM CSD phase in the previous time step.

@note NOT USED.

This is just a simple copying of the `phasemoles_` vector to `ophasemoles_`.
*/
void setOphasemoles ()
{
    for (register int idx = 0; idx < phasenum_; idx++)
       ophasemoles_[idx] = phasemoles_[idx];
}

/**
@brief Set the number of moles of a given GEM CSD phase in the previous time step.

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the number of moles to assign to that GEM phase in the previous time step
*/
void setOphasemoles (const unsigned int idx,
                     const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setOphasemoles","ophasemoles_",phasenum_,idx);
    }
    ophasemoles_[idx] = val;
    return;
}

/**
@brief Get the number of moles of every GEM CSD phase in the previous time step.

@note Used only in this class's copy constructor.

@return a pointer to the list of number of moles assigned to each GEM phase in last time step
*/
double *getOphasemoles () const
{
    return ophasemoles_;
}

/**
@brief Get the number of moles of a GEM CSD phase in the previous time step.

@note NOT USED.

@param idx is the GEM phase id to query
@return the number of moles assigned to that GEM phase in the previous time step
*/
double getOphasemoles (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getOphasemoles","ophasemoles_",phasenum_,idx);
    }
    return ophasemoles_[idx];
}

/**
@brief Output the number of moles of every GEM phase in the system.

*/
void writePhasemoles ()
{
    cout << endl;
    cout << "Vector of Phases:" << endl;
    for (register unsigned int i = 0; i < phasenum_; i++) {
        cout << "    " << phasename_[i] << ": "
             << phasemoles_[i] << " mol" << endl;
    }
    cout << endl;
    cout.flush();
}

/**
@brief Set the mass of a given GEM CSD phase [g].

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the mass to assign to that GEM phase [g]
*/
void setPhasemass (const unsigned int idx,
                   const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setPhasemass","phasemass_",phasenum_,idx);
    }
    phasemass_[idx] = val;
    return;
}

/**
@brief Set the mass of all GEM CSD phases [g].

*/
void setPhasemass ()
{
    setOphasemass();
    for (register long int i = 0; i < phasenum_; i++) {
        phasemass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
    }
}

/**
@brief Get the mass of every GEM CSD phase in the system [g].

@note Used only in this class's copy constructor.

@return a pointer to the list of masses assigned to each GEM phase [g]
*/
double *getPhasemass () const
{
    return phasemass_;
}

/**
@brief Get the mass of a GEM CSD phase (by id) [g].

@param idx is the GEM phase id to query
@return the mass assigned to that GEM phase [g]
*/
double getPhasemass (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getPhasemass","phasemass_",phasenum_,idx);
   }
    return phasemass_[idx];
}

/**
@brief Get the mass of a GEM CSD phase (by name) [g].

@param str is the GEM phase name to query
@return the mass assigned to that GEM phase [g]
*/
double getPhasemass (const string &str)
{
    string msg;
    unsigned int idx = getPhaseid(str);
    if (idx >= phasenum_) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getPhasemass",msg,phasenum_,idx);
    }
    return phasemass_[idx];
}

/**
@brief Set the mass of a given GEM CSD phase in the previous time step [g].

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the mass to assign to that GEM phase in the previous time step [g]
*/
void setOphasemass (const unsigned int idx,
                    const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setOPhasemass","ophasemass_",phasenum_,idx);
    }
    ophasemass_[idx] = val;
    return;
}

/**
@brief Set the mass of all GEM CSD phases in the previous time step [g].

*/
void setOphasemass ()
{
    for (register long int i = 0; i < phasenum_; i++) {
        ophasemass_[i] = phasemass_[i];
    }
}

/**
@brief Get the mass of every GEM CSD phase in the system in the previous time step [g].

@note Used only in this class's copy constructor.

@return a pointer to the list of masses assigned to each GEM phase in the previous time step [g]
*/
double *getOphasemass () const
{
    return ophasemass_;
}

/**
@brief Get the mass of a GEM CSD phase (by id) in the previous time step [g].

@param idx is the GEM phase id to query
@return the mass assigned to that GEM phase in the previous time step [g]
*/
double getOphasemass (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setOPhasemass","ophasemass_",phasenum_,idx);
    }
    return ophasemass_[idx];
}

/**
@brief Get the mass of a GEM CSD phase (by name) in the previous time step.

@note NOT USED.

@param str is the GEM phase name to query
@return the mass assigned to that GEM phase in the previous time step
*/
double getOphasemass (const string &str)
{
    string msg;
    unsigned int idx = getPhaseid(str);
    if (idx >= phasenum_) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getOPhasemass",msg,phasenum_,0);
    }
    return ophasemass_[idx];
}

/**
@brief Set the volume of a given GEM CSD phase.

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the volume to assign to that GEM phase
*/
void setPhasevolume (const unsigned int idx,
                     const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setPhasevolume","phasevolume_",phasenum_,idx);
    }
    phasevolume_[idx] = val;
    return;
}

/**
@brief Set the volume of all GEM CSD phases.

*/
void setPhasevolume ()
{
    setOphasevolume();
    for (register long int i = 0; i < phasenum_; i++) {
        phasevolume_[i] = (double)(node_->Ph_Volume(i));
        //cout << "phasevolume_[" << i << "] = " <<phasevolume_[i] << endl;
    }
}

/**
@brief Get the volume of every GEM CSD phase in the system.

@note Used only in this class's copy constructor.

@return a pointer to the list of volumes assigned to each GEM phase
*/
double *getPhasevolume() const
{
    return phasevolume_;
}

/**
@brief Get the volume of a GEM CSD phase (by id).

@note NOT USED.

@param idx is the GEM phase id to query
@return the volume assigned to that GEM phase
*/
double getPhasevolume (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getPhasevolume","phasevolume_",phasenum_,idx);
    }
    return phasevolume_[idx];
}

/**
@brief Get the volume of a GEM CSD phase (by name).

@param str is the GEM phase name to query
@return the volume assigned to that GEM phase
*/
double getPhasevolume (const string &str)
{
    string msg;
    unsigned int idx = getPhaseid(str);
    if (idx >= phasenum_) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getPhasevolume",msg,phasenum_,idx);
    }
    return phasevolume_[idx];
}

/**
@brief Set the volume of a given GEM CSD phase in the previous time step.

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the volume to assign to that GEM phase in the previous time step
*/
void setOphasevolume (const unsigned int idx,
                      const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setOphasevolume","ophasevolume_",phasenum_,idx);
    }
    ophasevolume_[idx] = val;
    return;
}

/**
@brief Get the volume of every GEM CSD phase in the system in the previous time step.

@return a pointer to the list of volumes assigned to each GEM phase in the previous time step
*/
void setOphasevolume ()
{
    for (register long int i = 0; i < phasenum_; i++) {
        ophasevolume_[i] = phasevolume_[i];
    }
}

/**
@brief Get the volume of every GEM CSD phase in the system in the previous time step.

@note Used only in this class's copy constructor.

@return a pointer to the list of volumes assigned to each GEM phase in the previous time step
*/
double *getOphasevolume () const
{
    return ophasevolume_;
}

/**
@brief Get the volume of a GEM CSD phase (by id) in the previous time step.

@note NOT USED.

@param idx is the GEM phase id to query
@return the volume assigned to that GEM phase in the previous time step
*/
double getOphasevolume (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getOphasevolume","ophasevolume_",phasenum_,idx);
    }
    return ophasevolume_[idx];
}

/**
@brief Get the volume of a GEM CSD phase (by name) in the previous time step.

@note NOT USED.

@param str is the GEM phase name to query
@return the volume assigned to that GEM phase in the previous time step
*/
double getOphasevolume (const string &str)
{
    string msg;
    unsigned int idx = getPhaseid(str);
    if (idx >= phasenum_) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getOphasevolume",msg,phasenum_,idx);
    }
    return ophasevolume_[idx];
}

/**
@brief Set the moles of a GEM CSD phase identified as a solvent.

@note NOT USED.

@param idx is the GEM CSD solvent id
@param val is the moles to assign to that GEM phase
*/
void setCarrier (const unsigned int idx,
                 const double val)
{
    if (idx >= solutionphasenum_) {
        throw EOBException("ChemicalSystem","setCarrier","carrier_",solutionphasenum_,idx);
    }
    carrier_[idx] = val;
    return;
}

/**
@brief Get the moles of every GEM CSD phase identified as a solvent.

@note Used only in this class's copy constructor.

@return a pointer to the list of moles assigned to each GEM solvent phase
*/
double *getCarrier () const
{
    return carrier_;
}

/**
@brief Get the moles of a GEM CSD solvent phase (by id).

@note NOT USED.

@param idx is the GEM solvent phase name to query
@return the moles assigned to that GEM solvent phase
*/
double getCarrier (const unsigned int idx)
{
    if (idx >= solutionphasenum_) {
        throw EOBException("ChemicalSystem","getCarrier","carrier_",solutionphasenum_,idx);
    }
    return carrier_[idx];
}

/**
@brief Set the surface area of a GEM CSD phase (m<sup>2</sup>/kg).

@note NOT USED.

@param idx is the GEM CSD phase id
@param val is the surface area to assign to that GEM phase
*/
void setSurfacearea (const unsigned int idx,
                     const double val)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setSurfacearea",
                           "surfacearea_",phasenum_,idx);
    }
    surfacearea_[idx] = val;
}

/**
@brief Get the surface area of all GEM CSD phases (m<sup>2</sup>/kg).

@note Used only in this class's copy constructor.

@return a pointer to the list of all phase surface areas
*/
double *getSurfacearea () const
{
    return surfacearea_;
}

/**
@brief Get the surface area of a GEM CSD phase (m<sup>2</sup>/kg).

@note NOT USED.

@param idx is the GEM CSD phase id
@return the surface area assigned to that GEM phase
*/
double getSurfacearea (const unsigned int idx)
{
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","getSurfacearea",
                           "surfacearea_",phasenum_,idx);
    }
    return surfacearea_[idx];
}

/**
@brief Set the volume fraction of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@param val is the volume fraction to assign to that microstructure phase
*/
void setMicphasevolfrac (const unsigned int idx,
                         const double val)
{
    if (idx >= micphasevolfrac_.size()) {
        throw EOBException("ChemicalSystem","setMicphasevolfrac",
                           "micphasevolfrac_",micphasevolfrac_.size(),idx);
    }
    micphasevolfrac_[idx] = val;
    return;
}

/**
@brief Get the volume fractions of all microstructure phases.

@return a pointer to the list of volume fractions of every microstructure phase
*/
vector<double> getMicphasevolfrac () const
{
    return micphasevolfrac_;
}

/**
@brief Get the volume fraction of a microstructure phase (by id).

@param idx is the microstructure phase id
@return the volume fraction assigned to that microstructure phase
*/
double getMicphasevolfrac (const unsigned int idx)
{
    if (idx >= micphasevolfrac_.size()) {
        throw EOBException("ChemicalSystem","getMicphasevolfrac",
                           "micphasevolfrac_",micphasevolfrac_.size(),idx);
    }
    // cout << "        " << micphasename_[idx] << " volfrac = "
    //      << micphasevolfrac_[idx] << endl;
    return micphasevolfrac_[idx];
}

/**
@brief Set the volume of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@param val is the volume to assign to that microstructure phase
*/
void setMicphasevolume (const unsigned int idx,
                        const double val)
{
    if (idx >= micphasevolume_.size()) {
        throw EOBException("ChemicalSystem","setMicphasevolume",
                           "micphasevolume_",micphasevolume_.size(),idx);
    }
    micphasevolume_[idx] = val;
    return;
}

/**
@brief Get the volumes of all microstructure phases.

@note Used only in this class's copy constructor.

@return a pointer to the list of volumes of every microstructure phase
*/
vector<double> getMicphasevolume () const
{
    return micphasevolume_;
}

/**
@brief Get the volume of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@return the volume assigned to that microstructure phase
*/
double getMicphasevolume(const unsigned int idx)
{
    if (idx >= micphasevolume_.size()) {
        throw EOBException("ChemicalSystem","getMicphasevolume",
                           "micphasevolume_",micphasevolume_.size(),idx);
    }
    return micphasevolume_[idx];
}

/**
@brief Set the mass of a microstructure phase (by id).

@param idx is the microstructure phase id
@param val is the mass to assign to that microstructure phase
*/
void setMicphasemass (const unsigned int idx,
                      const double val)
{
    if (idx >= micphasemass_.size()) {
        throw EOBException("ChemicalSystem","setMicphasemass",
                           "micphasemass_",micphasemass_.size(),idx);
    }
    micphasemass_[idx] = val;
    double v0 = node_->DC_V0(getMic2DC(idx,0),P_,T_);
    double dcmm = getDCmolarmass(getMic2DC(idx,0));
    cout << "    " << micphasename_[idx] << ": v0 = "
         << v0 << ", dcmm = " << dcmm << endl;
    micphasevolume_[idx] = val * v0 / dcmm;
    return;
}

/**
@brief Get the masses of all microstructure phases.

@note Used only in this class's copy constructor.

@return a pointer to the list of masses of every microstructure phase
*/
vector<double> getMicphasemass () const
{
    return micphasemass_;
}

/**
@brief Get the mass of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@return the mass assigned to that microstructure phase
*/
double getMicphasemass (const unsigned int idx)
{
    if (idx >= micphasemass_.size()) {
        throw EOBException("ChemicalSystem","getMicphasemass",
                           "micphasemass_",micphasemass_.size(),idx);
    }
    return micphasemass_[idx];
}

/**
@brief Output all the microstructure phases.

@note NOT USED.
*/
void writeMicphases ()
{
    cout << endl;
    cout << "Microstructure phase quantities:" << endl;
    cout << "Name     Mass (g)     Volume (m3)     Volume Fraction" << endl;
    cout << "----     --------     -----------     ---------------" << endl;
    for (register unsigned int i = 1; i < micphasename_.size(); i++) {
        cout << micphasename_[i] << "     " << micphasemass_[i] << "     "
             << micphasevolume_[i] << "     " << micphasevolfrac_[i] << endl;
    }
    cout << "Void     0.0    " << micvoidvolume_ << " " << micvoidvolfrac_ << endl;
    cout << endl;
    cout.flush();
}

/**
@brief Set the mass dissolved of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@param val is the dissolved mass to assign to that microstructure phase
*/
void setMicphasemassdissolved (const unsigned int idx,
                               const double val)
{
    if (idx >= micphasemassdissolved_.size()) {
        throw EOBException("ChemicalSystem","setMicphasemassdissolved",
                           "micphasemassdissolved_",micphasemassdissolved_.size(),idx);
    }
    micphasemassdissolved_[idx] = val;
    return;
}

/**
@brief Get the dissolved masses of all microstructure phases.

@note Used only in this class's copy constructor.

@return a pointer to the list of dissolved masses of every microstructure phase
*/
vector<double> getMicphasemassdissolved () const
{
    return micphasemassdissolved_;
}

/**
@brief Get the dissolved mass of a microstructure phase (by id).

@note NOT USED.

@param idx is the microstructure phase id
@return the dissolved mass assigned to that microstructure phase
*/
double getMicphasemassdissolved (const unsigned int idx)
{
    if (idx >= micphasemassdissolved_.size()) {
        throw EOBException("ChemicalSystem","getMicphasemassdissolved",
                           "micphasemassdissolved_",micphasemassdissolved_.size(),idx);
    }
    return micphasemassdissolved_[idx];
}

/**
@brief Set the volume of void space in the microstructure.

@note NOT USED.

@param val is the volume of void space to assign to the microstructure
*/
void setMicvoidvolume (const double val)
{
    micvoidvolume_ = val;
}

/**
@brief Get the volume of void space in the microstructure.

@note Used only in this class's copy constructor.

@return the volume of void space in the microstructure
*/
double getMicvoidvolume () const
{
    return micvoidvolume_;
}

/**
@brief Set the volume fraction of void space in the microstructure.

@note NOT USED.

@param val is the volume fraction of void space to assign to the microstructure
*/
void setMicvoidvolfrac (const double val)
{
    micvoidvolfrac_ = val;
}

/**
@brief Get the volume fraction of void space in the microstructure.

@note Used only in this class's copy constructor.

@return the volume fraction of void space in the microstructure
*/
double getMicvoidvolfrac () const
{
    return micvoidvolfrac_;
}

/**
@brief Set the total microstructure volume.

@note NOT USED.

@param val is the total volume to assign to the microstructure
*/
void setMictotvolume (const double val)
{
    mictotvolume_ = val;
}

/**
@brief Get the total microstructure volume.

@note Used only in this class's copy constructor.

@return the total volume of the microstructure
*/
double getMictotvolume () const
{
    return mictotvolume_;
}

/**
@brief Set the initial total microstructure volume.

@note NOT USED.

@param val is the initial total volume to assign to the microstructure
*/
void setMictotinitvolume (const double val)
{
    mictotinitvolume_ = val;
}

/**
@brief Get the initial total microstructure volume.

@return the initial total volume of the microstructure
*/
double getMictotinitvolume () const
{
    return mictotinitvolume_;
}

/**
@brief Set the upper bound on the moles of a dependent component.

@param idx is the id of the DC component
@param val is the upper bound to set
*/
void setDCupperlimit (const unsigned int idx,
                      const double val)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","setDCupperlimit",
                           "DCupperlimit_",DCnum_,idx);
    }
    DCupperlimit_[idx] = val;
    return;
}

/**
@brief Get the upper bound on the moles of every dependent component (DC).

@note Used only in this class's copy constructor.

@return a pointer to the list of DC upper bounds
*/
double *getDCupperlimit () const
{
    return DCupperlimit_;
}

/**
@brief Get the upper bound on the moles of a dependent component (DC).

@note NOT USED.

@param idx is the id of the DC component
@return the upper bound on moles of that DC
*/
double getDCupperlimit (const unsigned int idx)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","getDCupperlimit",
                           "DCupperlimit_",DCnum_,idx);
    }
    return DCupperlimit_[idx];
}

/**
@brief Get the lower bound on the moles of a dependent component (DC).

@param idx is the id of the DC component
@param val is the lower bound to set for this DC
*/
void setDClowerlimit (const unsigned int idx,
                      const double val)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","setDClowerlimit",
                           "DClowerlimit_",DCnum_,idx);
    }
    DClowerlimit_[idx] = val;
    return;
}

/**
@brief Get the lower bound on the moles of a dependent component (DC).

@note Used only in this class's copy constructor.

@param idx is the id of the DC component
@return a pointer to the list of DC lower bounds
*/
double *getDClowerlimit () const
{
    return DClowerlimit_;
}

/**
@brief Get the lower bound on the moles of a dependent component (DC).

@note NOT USED.

@param idx is the id of the DC component
@return a pointer to the list of DC lower bounds
*/
double getDClowerlimit (const unsigned int idx)
{
    if (idx >= DCnum_) {
        throw EOBException("ChemicalSystem","getDClowerlimit",
                           "DClowerlimit_",DCnum_,idx);
    }
    return DClowerlimit_[idx];
}

/**
@brief Set the molar mass of a particular independent component (IC) [g/mol].

@note NOT USED.

@param idx is the id of the IC
@param val is the molar mass to assign to that IC [g/mol]
*/
void setICmolarmass (const unsigned int idx,
                     const double val)
{
    if (idx >= ICmolarmass_.size()) {
        throw EOBException("ChemicalSystem","setICmolarmass",
                           "ICmolarmass_",ICmolarmass_.size(),idx);
    }
    ICmolarmass_[idx] = val;
    return;
}

/**
@brief Get the molar mass of a particular independent component (IC), by id [g/mol].

@param idx is the id of the IC
@return the molar mass of that IC, [g/mol]
*/
double getICmolarmass (const unsigned int idx)
{
    if (idx >= ICmolarmass_.size()) {
        throw EOBException("ChemicalSystem","getICmolarmass",
                           "ICmolarmass_",ICmolarmass_.size(),idx);
    }
    return ICmolarmass_[idx];
}

/**
@brief Get the molar mass of a particular independent component (IC), by name [g/mol].

@param str is the name of the IC
@return the molar mass of that IC, [g/mol]
*/
double getICmolarmass (const string &str)
{
    string msg;
    int idx = getICid(str);
    if (idx >= ICmolarmass_.size()) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getICmolarmass",msg,ICmolarmass_.size(),0);
    }
    return ICmolarmass_[idx];
}

/**
@brief Get the molar masses of all independent components (IC) [g/mol].

@note Used only in this class's copy constructor.

@return the vector of IC molar masses [g/mol]
*/
vector<double> getICmolarmass () const
{
    return ICmolarmass_;
}

/**
@brief Set the molar mass of a particular dependent component (DC) [g/mol].

@note NOT USED.

@param idx is the id of the DC
@param val is the molar mass to assign to that DC [g/mol]
*/
void setDCmolarmass (const unsigned int idx,
                     const double val)
{
    if (idx >= DCmolarmass_.size()) {
        throw EOBException("ChemicalSystem","setDCmolarmass",
                           "DCmolarmass_",DCmolarmass_.size(),idx);
    }
    DCmolarmass_[idx] = val;
    return;
}

/**
@brief Get the molar mass of a particular dependent component (DC), by id [g/mol].

@param idx is the id of the DC
@return the molar mass of that DC, [g/mol]
*/
double getDCmolarmass (const unsigned int idx)
{
    if (idx >= DCmolarmass_.size()) {
        throw EOBException("ChemicalSystem","getDCmolarmass",
                           "DCmolarmass_",DCmolarmass_.size(),idx);
    }
    return DCmolarmass_[idx];
}

/**
@brief Get the molar mass of a particular dependent component (DC), by name [g/mol].

@param str is the name of the DC
@return the molar mass of that DC, [g/mol]
*/
double getDCmolarmass (const string &str)
{
    string msg;
    int idx = getDCid(str);
    if (idx >= DCmolarmass_.size()) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getDCmolarmass",msg,DCmolarmass_.size(),0);
    }
    return DCmolarmass_[idx];
}

/**
@brief Get the molar masses of all dependent components (DC) [g/mol].

@note Used only in this class's copy constructor.

@return the vector of DC molar masses [g/mol]
*/
vector<double> getDCmolarmass () const
{
    return DCmolarmass_;
}

/**
@brief Get the molar masses of all GEM CSD phases [g/mol].

@note NOT USED.

@return the vector of GEM phase molar masses [g/mol]
*/
vector<double> getPhasemolarmass () const
{
    return phasemolarmass_;
}

/**
@brief Set the molar masses of all GEM CSD phases [g/mol].

*/
void setPhasemolarmass ()
{
    double pmm;
    phasemolarmass_.clear();
    phasemolarmass_.resize(phasenum_,0.0);
    for (register unsigned int pidx = 0; pidx < vphasestoich_.size(); pidx++) {
        pmm = 0.0;
        for (register unsigned int icidx = 0; icidx < vphasestoich_[pidx].size(); icidx++) {
            pmm += ((vphasestoich_[pidx][icidx]) * getICmolarmass(icidx));
        }
        phasemolarmass_[pidx] = pmm;
    }
}

/**
@brief Get the molar mass of a particular GEM CSD phase, by id [g/mol].

@note NOT USED.

@param idx is the id of the GEM phase
@return the molar mass of that GEM phase [g/mol]
*/
double getPhasemolarmass (const unsigned int idx)
{
    if (idx >= phasemolarmass_.size()) {
        throw EOBException("ChemicalSystem","getPhasemolarmass",
                           "phasemolarmass_",phasemolarmass_.size(),idx);
    }
    return phasemolarmass_[idx];
}

/**
@brief Get the molar mass of a particular GEM CSD phase, by name [g/mol].

@note NOT USED.

@param str is the name of the GEM phase
@return the molar mass of that GEM phase [g/mol]
*/
double getPhasemolarmass (const string &str)
{
    string msg;
    unsigned int idx = getPhaseid(str);
    if (idx >= phasemolarmass_.size()) {
        msg = "Name " + str + " does not have a valid phase id";
        throw EOBException("ChemicalSystem","getPhasemolarmass",msg,phasemolarmass_.size(),idx);
    }
    return phasemolarmass_[idx];
}

/**
@brief Get the number of moles of a particular IC in a particular DC.

Dependent components (DC) comprise one or (usually) more independent components (IC).
For example, the DC component CO2(g) has one mole of C and 2 moles of O per mole
of the DC.

@param dcidx is the id of the DC being queried
@param icidx is the id of the IC that may be contained within this DC
@return the number of moles of the ID icidx contained within the DC dcidx
*/
double getDCstoich (const unsigned int dcidx,
                    const unsigned int icidx)
{
    if (dcidx >= DCstoich_.size()) {
        throw EOBException("ChemicalSystem","getDCstoich",
                           "DCstoich_",DCstoich_.size(),dcidx);
    }
    if(icidx >= DCstoich_.size()) {
        throw EOBException("ChemicalSystem","getDCstoich",
                           "DCstoich_",DCstoich_[dcidx].size(),icidx);
    }
    return DCstoich_[dcidx][icidx];
}

/**
@brief Get the number of moles all ICs in a particular DC.

Dependent components (DC) comprise one or (usually) more independent components (IC).
For example, the DC component CO2(g) has one mole of C and 2 moles of O per mole
of the DC.

@note NOT USED.

@param dcidx is the id of the DC being queried
@return the list of moles of each IC in the DC specified by dcidx
*/
vector<double> getDCstoich (const int dcidx)
{
    if (dcidx >= DCstoich_.size()) {
        throw EOBException("ChemicalSystem","getDCstoich",
                           "DCstoich_",DCstoich_.size(),dcidx);
    }
    return DCstoich_[dcidx];
}

/**
@brief Get the number of moles all ICs in all DCs.

Dependent components (DC) comprise one or (usually) more independent components (IC).
For example, the DC component CO2(g) has one mole of C and 2 moles of O per mole
of the DC.

@note Used only in this class's copy constructor.

@return the matrix of moles of each IC in every DCs.
*/
vector<vector<double> > getDCstoich () const
{
    return DCstoich_;
}

/**
@brief Get the number of moles all ICs in all GEM CSD phases.

GEM phases comprise one or (usually) more independent components (IC).
For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole of calcium,
one mole of sulfur, four moles of hydrogen and six moles of oxygen

This method makes calls to the GEM node to get the matrix of GEM phase
stoichiometries, and assigns the results into the 1D `phasestoich_` array within
the method itself.  For this reason, the function does not return anything.

@todo Consider renaming this method because it does not actually return anything
to the calling function
*/
void getGEMPhasestoich ();

/**
@brief Set the matrix of moles all ICs in all GEM CSD phases.

GEM phases comprise one or (usually) more independent components (IC).
For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole of calcium,
one mole of sulfur, four moles of hydrogen and six moles of oxygen

This particular function just calls the getGEMPhasestoich method, which sets all
the phase stoichiometries into the 1D `phasestoich_` array.

*/
void setPhasestoich ()
{
    getGEMPhasestoich();
    return;
}

/**
@brief Set the moles of one IC component of one GEM CSD phase.

GEM phases comprise one or (usually) more independent components (IC).
For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole of calcium,
one mole of sulfur, four moles of hydrogen and six moles of oxygen

@note NOT USED.

@param pidx is the GEM phase id to modify
@param icidx is the IC component index
@param val is the number of moles of IC icidx in GEM phase pidx
*/
void setPhasestoich (const unsigned int pidx,
                     const unsigned int icidx,
                     const double val)
{
    unsigned int idx = (pidx * ICnum_) + icidx;
    if (idx >= phasenum_) {
        throw EOBException("ChemicalSystem","setPhasestoich","phasestoich_",phasenum_,idx);
    }
    phasestoich_[idx] = val;
    return;
}

/**
@brief Renormalizes phase stoichiometries per mole of oxygen in the phase.

This is a poorly worded function because it does not return anything.
simply renormalizes the GEM phase stoichiometries by the number of moles
of oxygen in that phase.

@todo Find out why it is useful to normalize by moles of oxygen here but not for the 1D array
version
@todo Rename this function to something more descriptive
*/
void getGEMVphasestoich ();

/**
@brief Caller function for `getGEMVphasestoich`.

*/
void setVphasestoich ()
{
    getGEMVphasestoich();
    return;
}

/**
@brief Get the number of moles of an IC in a GEM phase.

This function uses the 1D array pointed to by the `phasestoich_` member
variable.  This variable lists the IC molar stoichiometry of phase one in the first
n elements (where n is the number of ICs), and of phase two in the second
n elements, etc., instead of storing as a 2D matrix

@note NOT USED.

@param pidx is the id of the GEM phase to query
@param icidx is the index of the IC being queried
@return the molar stoichiometry of the icidx-th IC in phase id pidx
*/
double getPhasestoich (const unsigned int pidx,
                      const unsigned int icidx)
{
    unsigned int index = (pidx * ICnum_) + icidx;
    if (index >= phasenum_) {
        throw EOBException("ChemicalSystem","getPhasestoich","phasestoich_",phasenum_,index);       
    }
    return phasestoich_[index];
}

/**
@brief Get the list of moles of all ICs in a GEM phase.

This function uses the 1D array pointed to by the `phasestoich_` member
variable.  This variable lists the IC molar stoichiometry of phase one in the first
n elements (where n is the number of ICs), and of phase two in the second
n elements, etc., instead of storing as a 2D matrix

@param pidx is the id of the GEM phase to query
@return a pointer to the list of IC moles in phase id pidx
*/
double *getPhasestoich(const unsigned int pidx)
{
    if (pidx >= phasenum_) {
        throw EOBException("ChemicalSystem","getPhasestoich","phasestoich_",phasenum_,pidx);       
    }
    double *p = phasestoich_;
    p += (pidx * ICnum_);
    return p;
}

/**
@brief Get the number of moles of an IC in a GEM phase by phase id.

This function uses the 2D matrix form of the molar stoichiometries, implemented
as a vector of vectors, `vphasestoich_`.  Each row is a phase and each column
is an IC.

@note NOT USED.

@param pidx is the id of the GEM phase to query
@param icidx is the index of the IC being queried
@return the molar stoichiometry of the icidx-th IC in phase id pidx
*/
double getVphasestoich(const unsigned int pidx,
                       const unsigned int icidx)
{
    if (pidx >= vphasestoich_.size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                           vphasestoich_.size(),pidx);       
    }
    if (icidx >= vphasestoich_[pidx].size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                           vphasestoich_[pidx].size(),icidx);
    }
    return ((double)(vphasestoich_[pidx][icidx]));
}

/**
@brief Get the number of moles of an IC in a GEM phase by phase name.

This function uses the 2D matrix form of the molar stoichiometries, implemented
as a vector of vectors, `vphasestoich_`.  Each row is a phase and each column
is an IC.

@param str is the name of the GEM phase to query
@param icidx is the index of the IC being queried
@return the molar stoichiometry of the icidx-th IC in phase id pidx
*/
double getVphasestoich (const string &str,
                        const unsigned int icidx)
{
    unsigned int pidx = getPhaseid(str);
    if (pidx >= vphasestoich_.size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                           vphasestoich_.size(),pidx);       
    }
    if (icidx >= vphasestoich_[pidx].size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                           vphasestoich_[pidx].size(),icidx);
    }
    return ((double)(vphasestoich_[pidx][icidx]));
}

/**
@brief Get the number of moles of each IC in a GEM phase by phase id.

This function uses the 2D matrix form of the molar stoichiometries, implemented
as a vector of vectors, `vphasestoich_`.  Each row is a phase and each column
is an IC.

@note NOT USED.

@param pidx is the id of the GEM phase to query
@return the list of moles of each IC in phase id pidx
*/
vector<double> getVphasestoich (const unsigned int pidx)
{
    if (pidx >= vphasestoich_.size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                            vphasestoich_.size(),pidx);       
    }
    return vphasestoich_[pidx];
}

/**
@brief Get the number of moles of each IC in a GEM phase by phase name.

This function uses the 2D matrix form of the molar stoichiometries, implemented
as a vector of vectors, `vphasestoich_`.  Each row is a phase and each column
is an IC.

@note NOT USED.

@param str is the name of the GEM phase to query
@return the list of moles of each IC in phase id pidx
*/
vector<double> getVphasestoich (const string &str)
{
    unsigned int pidx = getPhaseid(str);
    if (pidx >= vphasestoich_.size()) {
        throw EOBException("ChemicalSystem","getVphasestoich","vphasestoich_",
                           vphasestoich_.size(),pidx);       
    }
    return vphasestoich_[pidx];
}

/**
@brief Get the 1D array of IC stoichiometries of each GEM CSD phase.

This function uses the 1D array pointed to by the `phasestoich_` member
variable.  This variable lists the IC molar stoichiometry of phase one in the first
n elements (where n is the number of ICs), and of phase two in the second
n elements, etc., instead of storing as a 2D matrix

@return a pointer to the list of moles of each IC in every GEM phase
*/
double *getPhasestoich () const
{
    return phasestoich_;
}

/**
@brief Get the 2D matrix of IC stoichiometries of each GEM phase.

This function uses the 2D matrix form of the molar stoichiometries, implemented
as a vector of vectors, `vphasestoich_`.  Each row is a phase and each column
is an IC.

@note Used only in this class's copy constructor.

@return the 2D matrix of moles of each IC in each GEM phase
*/
vector<vector<double> > getVphasestoich () const
{
    return vphasestoich_;
}

/**
@brief Get the chemical potential of an independent component (IC) [mol/mol].

@note NOT USED.

@param icidx is the index of the IC being queried
@return the chemical potential of the IC
*/
double getICchempot (const unsigned int icidx)
{
    if (icidx >= ICnum_) {
        throw EOBException("ChemicalSystem","getICchempot","ICchempot_",ICnum_,icidx);       
    }
    return ICchempot_[icidx];
}

/**
@brief Get the list of chemical potentials of each independent component (IC) [mol/mol].

@note Used only in this class's copy constructor.

@return a pointer to the list of IC chemical potentials [mol/mol]
*/
double *getICchempot () const
{
    return ICchempot_;
}

/**
@brief Get the molal activity coefficient of a dependent component (DC).

@note NOT USED.

@param dcidx is the index of the DC being queried
@return the activity coefficient of the DC
*/
double getDCactivitycoeff (const unsigned int dcidx) {
    if (dcidx >= DCnum_) {
        throw EOBException("ChemicalSystem","getDCactivitycoeff","DCactivitycoeff_",
                            DCnum_,dcidx);       
    }
    return DCactivitycoeff_[dcidx];
}

/**
@brief Get the list of molal activity coefficients of each dependent component (DC).

@note Used only in this class's copy constructor.

@return a pointer to the list of DC molal activity coefficients
*/
double *getDCactivitycoeff () const
{
    return DCactivitycoeff_;
}

/**
@brief Get the error (residual) in the mass balance for an independent component (IC).

@note NOT USED.

@param icidx is the index of the IC being queried
@return the residual for the IC being queried
*/
double getICresiduals (const unsigned int icidx)
{
    if (icidx >= ICnum_) {
        throw EOBException("ChemicalSystem","getICresiduals","ICresiduals_",
                           ICnum_,icidx);       
    }
    return ICresiduals_[icidx];
}

/**
@brief Get the list of errors (residuals) in the mass balance for each independent component (IC).

@note Used only in this class's copy constructor.

@return a pointer to the list of IC residuals
*/
double *getICresiduals () const
{
    return ICresiduals_;
}

/**
@brief Get the chemical potential of a dependent component (DC) [mol/mol or J/mol].

@note NOT USED.

@param dcidx is the index of the DC being queried
@param norm is true if units will be in mol/mol, otherwise in J/mol
@return the chemical potential of the DC
*/
double getDCchempot(const long int dcidx, bool norm = false) {
    return (node_->Get_muDC(dcidx,norm));
}

/**
@brief Get the chemical activity of a dependent component (DC).

@note NOT USED.

@param dcidx is the index of the DC being queried
@return the chemical activity of the DC
*/
double getDCactivity (const long int dcidx)
{
    return (node_->Get_aDC(dcidx));
}

/**
@brief Get the concentration of a dependent component (DC).

The units of the returned value depend on the typd of DC being queried:

    - Aqueous species [molal]
    - Gas species [partial pressure]
    - Surface complexes [mol/m<sup>2</sup>]
    - Species in other phases [mole fraction]

@note NOT USED.

@param dcidx is the index of the DC being queried
@return the concentration of the DC in appropriate units
*/
double getDCconcentration (const long int dcidx)
{
    return (node_->Get_cDC(dcidx));
}

/**
@brief Set the red channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@param val is the value of the red channel to set
*/
void setRed (const unsigned int mpidx,
             const double rval)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","setRed","color_",color_.size(),mpidx);
    }
    color_[mpidx][0] = min(rval,COLORSATVAL);
    return;
}

/**
@brief Get the red channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@return the value of the red channel for this phase's color
*/
double getRed (const unsigned int mpidx)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","getRed","color_",color_.size(),mpidx);
    }
    return color_[mpidx][0];
}

/**
@brief Set the green channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@param val is the value of the green channel to set
*/
void setGreen (const unsigned int mpidx,
               const double rval)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","setGreen","color_",color_.size(),mpidx);
    }
    color_[mpidx][1] = min(rval,COLORSATVAL);
    return;
}
/**
@brief Get the green channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@return the value of the green channel for this phase's color
*/
double getGreen (const unsigned int mpidx)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","getGreen","color_",color_.size(),mpidx);
    }
    return color_[mpidx][1];
}

/**
@brief Set the blue channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@param val is the value of the blue channel to set
*/
void setBlue (const unsigned int mpidx,
              const double rval)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","setBlue","color_",color_.size(),mpidx);
    }
    color_[mpidx][2] = min(rval,COLORSATVAL);
    return;
}

/**
@brief Get the blue channel in the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@return the value of the blue channel for this phase's color
*/
double getBlue (const unsigned int mpidx)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","getBlue","color_",color_.size(),mpidx);
    }
    return color_[mpidx][2];
}

/**
@brief Set the rgb triplet for color of a microstructure phase.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@param cv is the vector of rgb values to set for that phase
*/
void setColor (const unsigned int mpidx,
               vector<double> cv)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","setColor","color_",color_.size(),mpidx);
    }
    color_[mpidx] = cv;
    return;
}

/**
@brief Get the rgb triplet for color of a microstructure phase.

@param mpidx is the index of the microstructure phase
@return the vector of rgb values defining this phase's color
*/
vector<double> getColor (const unsigned int mpidx)
{
    if (mpidx >= color_.size()) {
        throw EOBException("ChemicalSystem","getColor","color_",color_.size(),mpidx);
    }
    return color_[mpidx];
}

/**
@brief Set the rgb triplet for color of every microstructure phase.

@note NOT USED.

@param cv is the 2D matrix of rgb values to set each microstructure phase
*/
void setColor (vector<vector<double> > cv)
{
    color_ = cv;
}

/**
@brief Get the 2D matrix of rgb triplets for color of every microstructure phase.

@note Used only in this class's copy constructor.

@return the 2D matrix of rgb values of every microstructure phase
*/
vector<vector<double> > getColor () const
{
    return color_;
}

/**
@brief Set the grayscale value of a microstructure phase.

The grayscale value of each microstructure phase is in proportion to its relative
brightness in a backscattered electron image.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@param rval is the grayscale value to set for that microstructure phase
*/
void setGrayscale (const unsigned int mpidx,
                   const double rval)
{
    if (mpidx >= grayscale_.size()) {
        throw EOBException("ChemicalSystem","setGrayscale","grayscale_",grayscale_.size(),mpidx);
    }
    grayscale_[mpidx] = min(rval,COLORSATVAL);
    return;
}

/**
@brief Get the grayscale value of a microstructure phase.

The grayscale value of each microstructure phase is in proportion to its relative
brightness in a backscattered electron image.

@note NOT USED.

@param mpidx is the index of the microstructure phase
@return the grayscale value for the microstructure phase
*/
double getGrayscale (const unsigned int mpidx)
{
    if (mpidx >= grayscale_.size()) {
        throw EOBException("ChemicalSystem","getGrayscale","grayscale_",grayscale_.size(),mpidx);
    }
    return grayscale_[mpidx];
}

/**
@brief Set the grayscale value of every microstructure phase.

The grayscale value of each microstructure phase is in proportion to its relative
brightness in a backscattered electron image.

@note NOT USED.

@param gv is the list of grayscale values to set each microstructure phase
*/
void setGrayscale (vector<double> gv)
{
    grayscale_ = gv;
}

/**
@brief Get the list of grayscale values of every microstructure phase.

@note Used only in this class's copy constructor.

@return the list of grayscale values of every microstructure phase
*/
vector<double> getGrayscale () const
{
    return grayscale_;
}

/**
@brief Get the map of of the vector index of the microstructure phases by name.

The integer id of each microstructure phase is keyed to its name in this map,
so one can "look up" a phase id by the name of that phase.

@note Used only in this class's copy constructor.

@return the microstructure phase lookup map (look up by name)
*/
map<string,int> getMicidlookup () const
{
    return micidlookup_;
}

/**
@brief Get the map of of the vector index of the independent components (IC) by name.

The integer id of each independent component is keyed to its name in this map,
so one can "look up" a IC id by the name of that IC.

@note Used only in this class's copy constructor.

@return the IC lookup map (look up by name)
*/
map<string,int> getICidlookup () const
{
    return ICidlookup_;
}

/**
@brief Get the map of of the vector index of the dependent components (DC) by name.

The integer id of each dependent component is keyed to its name in this map,
so one can "look up" a DC id by the name of that DC.

@note Used only in this class's copy constructor.

@return the DC lookup map (look up by name)
*/
map<string,int> getDCidlookup () const
{
    return DCidlookup_;
}

/**
@brief Get the map of of the vector index of the GEM CSD phases by name.

The integer id of each phase defined in the GEM CSD is keyed to its name in this map,
so one can "look up" a GEM phase id by the name of that phase.

@note Used only in this class's copy constructor.

@return the DC lookup map (look up by name)
*/
map<string,int> getPhaseidlookup() const
{
return phaseidlookup_;
}

/**
@brief Get the class code of an independent component (IC).

@note NOT USED.

@param icidx is the index of the IC
@return the class code of the IC
*/
char getICclasscode (const unsigned int icidx)
{
    if (icidx >= ICclasscode_.size()) {
        throw EOBException("ChemicalSystem","getICclasscode","ICclasscode_",
                           ICclasscode_.size(),icidx);
    }
    return ICclasscode_[icidx];
}

/**
@brief Get the class codes of all independent components (IC).

@note Used only in this class's copy constructor.

@return the list of character class codes of the ICs
*/
vector<char> getICclasscode () const
{
    return ICclasscode_;
}

/**
@brief Get the class code of a dependent component (DC).

@param dcidx is the index of the DC
@return the class code of the DC
*/
char getDCclasscode (const unsigned int dcidx)
{
    if (dcidx >= DCclasscode_.size()) {
        throw EOBException("ChemicalSystem","getDCclasscode","DCclasscode_",
                           DCclasscode_.size(),dcidx);
    }
    return DCclasscode_[dcidx];
}

/**
@brief Get the class codes of all dependent components (DC).

@return the list of character class codes of the DCs
*/
vector<char> getDCclasscode () const
{
    return DCclasscode_;
}

/**
@brief Get the class code of a phase defined in the GEM CSD.

@param pidx is the index of a GEM phase
@return the class code of the GEM phase 
*/
char getPhaseclasscode (const unsigned int pidx)
{
    if (pidx >= phaseclasscode_.size()) {
        throw EOBException("ChemicalSystem","getPhaseclasscode","phaseclasscode_",
                           phaseclasscode_.size(),pidx);
    }
    return phaseclasscode_[pidx];
}

/**
@brief Get the class codes of all phases defined in the GEM CSD.

@return the list of character class codes of the GEM phases
*/
vector<char> getPhaseclasscode () const
{
    return phaseclasscode_;
}

/**
@brief Set the system temperature [K].

@note NOT USED.

@param val is the absolute temperature [K]
*/
void setT (const double val)
{
    T_ = val;
}

/**
@brief Get the system temperature [K].

@return the absolute temperature [K]
*/
double getT () const
{
    return T_;
}

/**
@brief Set the system pressure [Pa].

@note NOT USED.

@param val is the pressure [Pa]
*/
void setP (const double val)
{
    P_ = val;
}

/**
@brief Get the system pressure [Pa].

@return the pressure [Pa]
*/
double getP () const
{
    return P_;
}

/**
@brief Set the system volume [m<sup>3</sup>].

@note NOT USED.

@param val is the system volume [m<sup>3</sup>]
*/
void setVs (const double val)
{
    Vs_ = val;
}

/**
@brief Get the system volume [m<sup>3</sup>].

@note Used only in this class's copy constructor.

@return the system volume [m<sup>3</sup>]
*/
double getVs () const
{
    return Vs_;
}

/**
@brief Set the system mass [kg].

@note NOT USED.

@param val is the system mass [kg]
*/
void setMs (const double val)
{
    Ms_ = val;
}

/**
@brief Get the system mass [kg].

@note Used only in this class's copy constructor.

@return the system mass [kg]
*/
double getMs () const
{
    return Ms_;
}

/**
@brief Set the system Gibbs energy [J].

@note NOT USED.

@param val is the system Gibbs energy [J]
*/
void setGs (const double val)
{
    Gs_ = val;
}

/**
@brief Get the system Gibbs energy [J].

@note Used only in this class's copy constructor.

@return the system Gibbs energy [J]
*/
double getGs () const
{
    return Gs_;
}

/**
@brief Set the system enthalpy [J].

@note NOT USED.

@param val is the system enthalpy [J]
*/
void setHs (const double val)
{
    Hs_ = val;
}

/**
@brief Get the system enthalpy [J].

@note Used only in this class's copy constructor.

@return the system enthalpy [J]
*/
double getHs () const
{
    return Hs_;
}

/**
@brief Set the aqueous solution ionic strength [molal].

@note NOT USED.

@param isval is the aqueous solution ionic strength [molal]
*/
void setIonicstrength (const double isval)
{
    ionicstrength_ = isval;
}

/**
@brief Get the aqueous solution ionic strength [molal].

@note Used only in this class's copy constructor.

@return the aqueous solution ionic strength [molal]
*/
double getIonicstrength () const
{
    return ionicstrength_;
}

/**
@brief Set the aqueous solution pH.

@note NOT USED.

@param val is the aqueous solution pH
*/
void setPH (const double val)
{
    pH_ = val;
}

/**
@brief Get the aqueous solution pH.

@return the aqueous solution pH
*/
double getPH () const
{
    return pH_;
}

/**
@brief Set the aqueous solution pe.

@note NOT USED.

@param val is the aqueous solution pe
*/
void setPe (const double val)
{
    pe_ = val;
}

/**
@brief Get the aqueous solution pe.

@note Used only in this class's copy constructor.

@return the aqueous solution pe
*/
double getPe () const
{
    return pe_;
}

/**
@brief Set the aqueous solution Eh [V].

@note NOT USED.

@param val is the aqueous solution Eh [V]
*/
void setEh (const double val)
{
    Eh_ = val;
}

/**
@brief Get the aqueous solution Eh [V].

@note Used only in this class's copy constructor.

@return the aqueous solution Eh [V]
*/
double getEh () const
{
    return Eh_;
}

/**
@brief Set the handle identification for the GEM-IPM node doing the calculations.

@note NOT USED.

@param val is the node handle identification number
*/
void setNodehandle (const long int val)
{
    nodehandle_ = val;
}

/**
@brief Get the handle identification for the GEM-IPM node doing the calculations.

@note Used only in this class's copy constructor.

@return the node handle identification number
*/
long int getNodehandle () const
{
    return nodehandle_;
}

/**
@brief Set the status of the GEM-IPM node doing the calculations.

The node status can be:

    - 0 = No GEM re-calculation needed for this node
    - 1 = Need GEM calculation with simplex (automatic) initial approximation (IA)
    - 2 = OK after GEM calculation with simplex (automatic) IA
    - 3 = Bad (not fully trustful) result after GEM calculation with simplex (automatic) IA
    - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
    - 5 = Need GEM calculation with no-simplex (smart) IA
            the previous GEM solution (full DATABR lists only)
    - 6 = OK after GEM calculation with no-simplex (smart) IA
    - 7 = Bad (not fully trustful) result after GEM calculation with no-simplex (smart) IA
    - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
    - 9 = Terminal error has occurred in GEMIPM2K (e.g., memory corruption).
            Restart is required.

@note NOT USED.

@param val is the node status number to set
*/
void setNodestatus (const long int val)
{
    nodestatus_ = val;
}

/**
@brief Get the status of the GEM-IPM node doing the calculations.

The node status can be:

    - 0 = No GEM re-calculation needed for this node
    - 1 = Need GEM calculation with simplex (automatic) initial approximation (IA)
    - 2 = OK after GEM calculation with simplex (automatic) IA
    - 3 = Bad (not fully trustful) result after GEM calculation with simplex (automatic) IA
    - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
    - 5 = Need GEM calculation with no-simplex (smart) IA
            the previous GEM solution (full DATABR lists only)
    - 6 = OK after GEM calculation with no-simplex (smart) IA
    - 7 = Bad (not fully trustful) result after GEM calculation with no-simplex (smart) IA
    - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
    - 9 = Terminal error has occurred in GEMIPM2K (e.g., memory corruption).
            Restart is required.

@note Used only in this class's copy constructor.

@return the node status number
*/
long int getNodestatus () const
{
    return nodestatus_;
}

/**
@brief Set the number of iterations executed by the GEM-IPM calculation.

@note NOT USED.

@param val is the number of iterations
*/
void setIterdone (const long int val)
{
    iterdone_ = val;
}

/**
@brief Get the number of iterations executed by the GEM-IPM calculation.

@note Used only in this class's copy constructor.

@return the number of iterations executed by the GEM-IPM solver
*/
long int getIterdone () const
{
    return iterdone_;
}

/**
@brief Set the pointer to the GEM TNode object doing the GEM-IPM calculations.

@note NOT USED (possibly used in GEM3K library).

@param np is a pointer to the GEM TNode object doing the GEM-IPM calculations
*/
void setNode (TNode *np)
{
    node_ = np;
}

/**
@brief Get the pointer to the GEM TNode object doing the GEM-IPM calculations.

@return a pointer to the GEM TNode object doing the GEM-IPM calculations
*/
TNode *getNode ()
{
    return node_;
}

/**
@brief Set the simulated time at which to implement the sulfate attack calculations [days].

@param sattack_time is the time at which to start sulfate attack [days]
*/
void setSattack_time(const double sattack_time)
{
    sattack_time_ = sattack_time;
}

/**
@brief Get the simulated time at which to start sulfate attack [days].

@note NOT USED.

@return the simulated time at which to begin sulfate attack [days]
*/
double getSattack_time () const
{
    return sattack_time_;
}

/**
@brief Set the simulated time at which to implement the leaching calculations [days].

@param leach_time is the time at which to start leaching [days]
*/
void setLeach_time (const double leach_time)
{
    leach_time_ = leach_time;
}

/**
@brief Get the simulated time at which to start leaching [days].

@note NOT USED.

@return the simulated time at which to begin leaching [days]
*/
double getLeach_time () const
{
    return leach_time_;
}
   
/**
@brief Formatted writing of some microstructure phase data to a stream.

@note NOT USED.

@param stream is the output stream to which to direct output
*/
void writeDb (ostream &stream);

/**
@brief Formatted writing of the name, id, and internal porosity of a microstructure phase.

@param i is the index of the microstructure phase
@param stream is the output stream to which to direct output
*/
void writeMember (const unsigned int i,
                  ostream &stream);
    

/**
@brief Formatted writing of ChemicalSystem data to a stream.

The output stream is defined, opened, and closed within the function itself
*/
void writeChemSys ();

/**
@brief Formatted writing of ChemicalSystem data to a prescribed stream.

@param out is the output stream to which to direct output
*/
void writeChemSys (ostream &out);
    
/**
@brief Calculates new equilibrium state of the GEM system and relate to microstructure.

One of the main purposes of the [[ChemicalSystem]] class is to calculate changes
to the state from one time increment to another.  This is accomplished
by changing the chemical system definition to alter the quantitities
of one or more independent components (IC), and then call GEMIPM2K
to calculate the new assemblage of DCs and phases that minimizes the
Gibbs free energy subject to the IC quantities in the system.

As a result, we need methods that will allow one to change the quantities
of ICs and then to call GEMIPM2K functions that will calculate the new
state and return the values of the quantities of the DCs and phases.

The method returns the GEM-IPM TNode status, which can be:

    - 0 = No GEM re-calculation needed for this node
    - 1 = Need GEM calculation with simplex (automatic) initial approximation (IA)
    - 2 = OK after GEM calculation with simplex (automatic) IA
    - 3 = Bad (not fully trustful) result after GEM calculation with simplex (automatic) IA
    - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
    - 5 = Need GEM calculation with no-simplex (smart) IA
            the previous GEM solution (full DATABR lists only)
    - 6 = OK after GEM calculation with no-simplex (smart) IA
    - 7 = Bad (not fully trustful) result after GEM calculation with no-simplex (smart) IA
    - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
    - 9 = Terminal error has occurred in GEMIPM2K (e.g., memory corruption).
            Restart is required.

@param time is the simulated time associated with this state [days]
@param isfirst is true if this is the first state calculation, false otherwise
@return the node status handle
*/
int calculateState (double time,
                    bool isfirst);

/**
@brief Update the number of moles of each IC based on changes to a dependent component.

@note NOT USED.

@param dcid is the index of the DC
@param moles is the change in number of moles of the DC
*/
void DCimpact (int dcid,
               double moles)
{
    for (int i = 0; i < ICnum_; i++) {
        ICmoles_[i] += DCstoich_[dcid][i] * moles;
    }

    return;
}

/**
@brief Check whether growth of a given dependent component is allowed.

The change in moles of a particular dependent component (DC) could conceivably
cause the moles of one of the independent components to decrease to effectively
zero.  This is not allowed for numerical stability reasons, so this function
checks whether that would happen.

@note NOT USED.

@param dcid is the index of the DC that wants to change in number moles
@param moles is the proposed change in moles of that DC
@return true if the change in this DC would not cause an IC mole content to drop to near zero
*/
bool checkDCimpact (int dcid,
                    double moles)
{
    bool possible;
    possible = true;
    for (int i = 0; i < ICnum_; i++) {
        if ((ICmoles_[i] + DCstoich_[dcid][i] * moles) < 2.0e-17) {
            possible = false;
            cout << "The growth of this phase can cause one or more IC moles to be" 
                 << " lower than 2.0e-17, so this phase can not grow."
                 << endl;
            break;
        }    
    }
    return possible;
}

/**
@brief Set the vector of saturation indices of all GEM CSD phases.

*/
void setSI ()
{
    SI_.clear();
    double *Falp;
    Falp = (node_->ppmm())->Falp;
    
    for (int i = 0; i < phasenum_; i++) {
        cout << "logSI for " << phasename_[i] << " is: "
             << Falp[i] << endl;
        double si = pow(10,Falp[i]);
        SI_.push_back(si);
    }
    return;
}

/**
@brief Get the vector of saturation indices of all GEM CSD phases.

@return the vector of saturation indices of all GEM CSD phases
*/
vector<double> getSI ()
{
    return SI_;
}

/**
@brief Get the saturation index of a GEM CSD phase, by its id.

@param phaseid is the index of the GEM phase to query
@return the saturation index of the GEM phase
*/
double getSI (int phaseid)
{
    if (phaseid >= SI_.size()) {
        throw EOBException("ChemicalSystem","getSI","SI_",
                            SI_.size(),phaseid);
    }
    return SI_[phaseid];
}

/**
@brief Get the saturation index of a GEM CSD phase, by its name.

@param str is the name of the GEM phase to query
@return the saturation index of the GEM phase
*/
double getSI (const string &str)
{
    return SI_[getPhaseid(str)];
}

/**
@brief Get the list of IC moles in the aqueous solution.

@return the vector of moles of each IC in the aqueous solution
*/
vector<double> getSolution ();

};      // End of ChemicalSystem class

#endif
