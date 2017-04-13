/**
@file Solution.cc
@brief Define methods for the Solution class

*/
#include "Solution.h"
#include <iostream>

using namespace std;

Solution::Solution (const string &GEMfilename,
                    const string &GEMdbrname)
{
    pH_ = 7.0;
    pe_ = 1.0;
    Eh_ = 0.0;
    T_ = 0.0;
    P_ = 101325.0;
    Vs_ = Ms_ = 1.0;
    Gs_ = Hs_ = 0.0;
    ionicstrength_ = 0.0;
    nodestatus_ = nodehandle_ = iterdone_ = 0;

    ICnum_ = DCnum_ = phasenum_ = solutionphasenum_ = 0;
    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    SI_.clear();

    char *cGEMfilename = (char*)GEMfilename.c_str();
    char *cGEMdbrname = (char*)GEMdbrname.c_str();
    cout << "Trying to read chemical system definition file " << cGEMfilename << endl;

    ///
    /// A new GEM3K node is allocated, completely separate from the one
    /// already allocated and used by the simulation's Chemicalsystem object
    ///

    node_ = new TNode();

    ///
    /// Use GEM3K to open and process the chemical system definition (CSD) file
    ///

    long int gemflag = 0;
    gemflag = node_->GEM_init(cGEMfilename);
    if (gemflag == 1) {
        cout << "Bad return from GEM_init: " << GEMfilename << " missing or corrupt." << endl;
        exit(0);
    }
    if (gemflag == -1) {
        cout << "Bad return from GEM_init: internal memory allocation error." << endl;
        exit(0);
    }

    ///
    /// GEM3K has read the CSD file and now can be queried for class member values
    ///

    ICnum_ = (unsigned int)((node_->pCSD())->nIC);
    cout << "ICnum_ is: " << ICnum_ << endl;
    DCnum_ = (unsigned int)((node_->pCSD())->nDC);
    cout << "DCnum_ is: " << DCnum_ << endl;
    phasenum_ = (unsigned int)((node_->pCSD())->nPH);
    cout << "phasenum_ is: " << phasenum_ << endl;
    solutionphasenum_ = (unsigned int)((node_->pCSD())->nPS);
    cout << "solutionphasenum_ is: " << solutionphasenum_ << endl;

    ///
    /// Attempt to allocate memory for the various arrays
    ///

    string exmsg;
    try {
        exmsg = "ICmoles_";
        ICmoles_ = new double [ICnum_];
        exmsg = "ICresiduals_";
        ICresiduals_ = new double [ICnum_];
        exmsg = "ICchempot_";
        ICchempot_ = new double [ICnum_];
        exmsg = "DCmoles_";
        DCmoles_ = new double [DCnum_];
        exmsg = "DCactivitycoeff_";
        DCactivitycoeff_ = new double [DCnum_];
        exmsg = "DCupperlimit_";
        DCupperlimit_ = new double [DCnum_];
        exmsg = "DClowerlimit_";
        DClowerlimit_ = new double [DCnum_];
        exmsg = "phasemoles_";
        phasemoles_ = new double [phasenum_];
        exmsg = "phasemass_";
        phasemass_ = new double [phasenum_];
        exmsg = "phasevolume_";
        phasevolume_ = new double [phasenum_];
        exmsg = "surfacearea_";
        surfacearea_ = new double [phasenum_];
        exmsg = "carrier_";
        carrier_ = new double [solutionphasenum_];
        exmsg = "phasestoich_";
        phasestoich_ = new double [phasenum_ * ICnum_];
        exmsg = "solidstoich_";
        solidstoich_ = new double [ICnum_];
    }
    catch (bad_alloc& ba) {
        cout << endl << "Bad_alloc Exception Thrown:" << endl;
        cout << "    Details:" << endl;
        cout << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
        cout << "    Error in allocating memory for array " << exmsg << endl;
        cerr << endl << "Bad_alloc Exception Thrown:" << endl;
        cerr << "    Details:" << endl;
        cerr << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
        cerr << "    Error in allocating memory for array " << exmsg << endl;
        exit(0);
    }


    ///
    /// Initialize GEM
    ///

    (node_->pCNode())->NodeStatusCH = NEED_GEM_AIA;
    cout << "Solution::Constructor: Entering GEM_run with node status = "
         << nodestatus_ << endl;
    cout << "NodeStatusCH is: " << (node_->pCNode())->NodeStatusCH << endl;
    cout.flush();

    ///
    /// Attempt to run GEM
    ///

    try {
        nodestatus_ = node_->GEM_run(false);
        cout << "Solution::Constructor: Exited GEM_run with node status = "
             << nodestatus_ << endl;
        cout.flush();
        if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
            cout << "ERROR: Call to GEM_run failed..." << endl;
            exit(1);
        }
        cout << "Solution::Constructor: Entering GEM_restore_MT..." << endl;
        cout.flush();

        node_->GEM_restore_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,&ICmoles_[0],
          &DCupperlimit_[0],&DClowerlimit_[0],&surfacearea_[0]);

        cout << "Done!" << endl;
        cout << "T_ is: " << T_ << endl;
        cout << "Solution::Constructor: Entering GEM_to_MT..." << endl;
        cout.flush();
        node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
            Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
            &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&phasemoles_[0],
            &phasevolume_[0],&phasemass_[0],&phasestoich_[0],
            &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
        cout << "Done!" << endl;
    }
    catch (GEMException e) {
        e.printException();
        exit(0);
    }
    catch (long int flag) {
        cout << endl << "GEM Exception Thrown:" << endl;
        cout << "Offending Function ChemicalSystem::ChemicalSystem" << endl;
        cerr << endl << "GEM Exception Thrown:" << endl;
        cerr << "Offending Function ChemicalSystem::ChemicalSystem" << endl;
        if (flag < 0) {
            cout << "Details:  Unkown internal error in TNode::GEM_read_dbr" << endl;
            cerr << "Details:  Unkown internal error in TNode::GEM_read_dbr" << endl;
        } else {
            cout << "Details:  Check ipmlog.txt file for error details" << endl;
            cerr << "Details:  Check ipmlog.txt file for error details" << endl;
        }
        exit(0);
    }

    cout << "Solution::Constructor: Entering GEM_read_dbr..." << endl;
    cout.flush();

    nodestatus_ =  node_->GEM_read_dbr(cGEMdbrname);
    if (nodestatus_ != 0) {
        cout << "ERROR: Call to GEM_read_dbr failed..." << endl;
        exit(1);
    }
    cout << "Done!" << endl;
    cout.flush();

    ///
    /// Transfer phase names from GEM3K to member variables
    ///

    string string1;
    for (int i = 0; i < ICnum_; i++) {
        string1.assign((node_->pCSD())->ICNL[i]);
        ICname_.push_back(string1);
    }
 
    string1.clear();  // This command may be unnecessary
    for (int i = 0; i < DCnum_; i++) {
        string1.assign((node_->pCSD())->DCNL[i]);
        DCname_.push_back(string1);
    }
 
    string1.clear();  // This command may be unnecessary
    for (int i = 0; i < phasenum_; i++) {
        string1.assign((node_->pCSD())->PHNL[i]);
        phasename_.push_back(string1);
    }
}

Solution::~Solution ()
{
    ///
    /// Clear out the vectors
    ///

    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    SI_.clear();

    ///
    /// Delete previously allocated memory
    ///

    delete[]solidstoich_;
    delete[]phasestoich_;
    delete[]carrier_;
    delete[]surfacearea_;
    delete[]phasevolume_;
    delete[]phasemass_;
    delete[]phasemoles_;
    delete[]DClowerlimit_;
    delete[]DCupperlimit_;
    delete[]DCactivitycoeff_;
    delete[]DCmoles_;
    delete[]ICchempot_;
    delete[]ICresiduals_;
    delete[]ICmoles_;

    delete node_;
}

void Solution::calculateState (bool isfirst)
{
    int status = 0;

    ///
    /// Load GEM data to the GEM3K library
    ///

    nodestatus_ = NEED_GEM_SIA;
    cout << "    Going into Solution::calculateState::GEM_from_MT..." << endl;
    cout.flush();
 
    node_->GEM_from_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,
           ICmoles_,DCupperlimit_,DClowerlimit_,surfacearea_,
           DCmoles_,DCactivitycoeff_);
    cout << "Done!" << endl;
  
    ///
    /// Execute GEM calculation
    ///

    cout << "    Going into Solution::calculateState::GEM_set_MT..." << endl;
    cout.flush();
    if (isfirst) {
        nodestatus_ = node_->GEM_run(false);
    } else {
        nodestatus_ = node_->GEM_run(true);
    }

    cout << "Done! nodestatus is " << nodestatus_ << endl;
    cout.flush();

    ///
    /// Get the GEM data back from the GEM3K library, assuming that it ran
    /// without error.  Check for the errors first.
    ///

    if (nodestatus_ == ERR_GEM_AIA || nodestatus_ == ERR_GEM_SIA
       || nodestatus_ == T_ERROR_GEM) {

        switch (nodestatus_) {
            case ERR_GEM_AIA:
              cout << "Flag nodestatus_ = ERR_GEM_AIA";
              break;
            case ERR_GEM_SIA:
              cout << "Flag nodestatus_ = ERR_GEM_SIA";
              break;
            case T_ERROR_GEM:
              cout << "Flag nodestatus_ = T_ERROR_GEM";
              break;
        } 

    } else if (nodestatus_ == BAD_GEM_AIA || nodestatus_ == BAD_GEM_SIA) {
  
        switch (nodestatus_) {
            case BAD_GEM_AIA:
              cout << "Flag nodestatus_ = BAD_GEM_AIA";
              break;
            case BAD_GEM_SIA:
              cout << "Flag nodestatus_ = BAD_GEM_SIA";
              break;
        }

    } else {

        cout << "    Going into Solution::calculateState::GEM_to_MT...";
        cout.flush();
        node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
                Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
                &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&phasemoles_[0],
                &phasevolume_[0],&phasemass_[0],&phasestoich_[0],
                &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
        cout << "Done!" << endl;
        cout << "after GEM_to_MT...Ms_ = " << Ms_ << endl;
        cout.flush();    
    }

  ///
  /// The GEM calculation did not precipitate or dissolve solid, so it
  /// stores the saturation index of each solid phase.  This now needs to
  /// be assigned to the corresponding array in the Solution object.
  ///

  setSI();

  return;

}

double Solution::calculateCrystrain (double ettrSI,
                                     double porevolfrac,
                                     double Kp,
                                     double Ks)
{
    crystrain_ = 0.0;

    /*
    calculateState(false);
    */
 
    ///
    /// Crystallization pressure only exists if the solution is supersaturated,
    /// which means that \f$\beta > 1\f$.
    ///

    if (ettrSI > 1.0) {

        ///
        /// Estimate the hydrostatic pressure in the pore solution as 1 atmosphere
        /// Unit of pressure for this calculation is MPa
        ///

        double pl = 0.101;

        ///
        /// The assumed largest radius of capillary pores offering entrance
        /// into gel porosity of C-S-H or other nanoporous component.
        /// 500 nm is chosen because it is half the size of a typical lattice
        /// site in THAMES (although the lattice resolution can be varied in
        /// the future if desired, so we may want to revisit this assumption).
        /// Unit of length for this calculation is millimeters.
        ///

        double r = 5.0e-4;

        ///
        /// Thickness of liquid film separating the crystallizing solid and
        /// the pore walls, assumed to be 1 nm.
        /// Unit of length for this calculation is millimeters.
        ///

        double delta = 1.0e-6;

        ///
        /// Crystal-liquid surface energy, assumed to be 100 mJ/m<sup>2</sup>.
        /// Unit of surface energy for this calculation is N/mm.
        ///

        double gamma = 1.0e-4; // N/mm
  
        ///
        /// Stress-free molar volume of the growing crystal, assumed to be ettringite
        /// Units of molar volume for this calcualtion is mm<sup>3</sup>/mol.
        ///
        /// @note This could be loaded up directly from the GEM CSD, rather than
        ///       hardwiring it into the code here.
        ///

        double Vc = 7.070e5;

        ///
        /// The ideal gas constant, with units of (N mm)/(mol K)
        ///

        double Rg = 8.314e3; // gas constant; N.mm/mol.K
  
        cout << "SI for ettringite is: " << ettrSI << endl; 

        ///
        /// Calculate the crystal mean curvature in equilibrium with the
        /// solution with this saturation index (Thompson-Freundlich effect)
        ///

        double kcr = Rg * T_ * log(ettrSI) / (Vc * gamma);

        ///
        /// If the portion of the crystal near the wall is a hemispherical cap,
        /// then the mean curvature is 2/r, where r is the radius of curvature
        /// of the crystal.  Therefore, the smallest pore within which the crystal
        /// can fit is delta larger than this.
        ///

        double rcr = (2.0 / kcr) + delta;
  
        ///
        /// Crystallization pressure associated with this pore size
        ///

        double pa = 2.0 * gamma * (1.0 / (rcr - delta) - 1.0 / (r - delta));    

        ///
        /// Strain can be retrieved from the stress via the effeictive elastic
        /// of the porous medium (poromechanics assumption)
        ///

        crystrain_ = (1.0 / (3.0 * Kp) - 1.0 / (3.0 * Ks)) * (porevolfrac * pa + pl);

        cout << "crystrain is: " << crystrain_ << endl;
  
    }

    return crystrain_;
}
