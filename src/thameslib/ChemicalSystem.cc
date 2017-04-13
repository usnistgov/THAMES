/**
@file ChemicalSystem.cc
@brief Method definitions for the ChemicalSystem base class
*/

#include "ChemicalSystem.h"

ChemicalSystem::ChemicalSystem (Solution *Solut,
                                const string &GEMfilename,
                                const string &GEMdbrname,
                                const string &Interfacefilename) : solut_(Solut)
{
  unsigned int i,j;
  double *amat;
  string exmsg;
  long int gemflag = 0;

  double *icmolarmass,*dcmolarmass;
  char *cc;
  register unsigned int ii,jj;
  int k;
  bool found = false;
    
  ///  The constructor initializes all the members to default values,
  ///  then launches the initial thermodynamic calculation, and sets
  ///  up the correspondences between GEM CSD phases and microstructure
  ///  phases.
  ///
  ///  All members are initialized to default values first, and all
  ///  vectors and maps are cleared.  The thermodynamic variables
  ///  are set to be consistent with neutral water at STP
  ///
 
  micphasenum_ = phasenum_ = solutionphasenum_ = 0;
  micimpuritynum_ = 4;
  micphasename_.clear();
  micid_.clear();
  DCname_.clear();
  phasename_.clear();
  ICnum_ = DCnum_ = phasenum_ = 0;
  micphasemembers_.clear();
  micphasemembersvolfrac_.clear();
  micDCmembers_.clear();
  randomgrowth_.clear();
  growthtemplate_.clear();
  affinity_.clear();
  porosity_.clear();
  k2o_.clear();
  na2o_.clear();
  mgo_.clear();
  so3_.clear();
  color_.clear();
  vphasestoich_.clear();
  micidlookup_.clear();
  ICidlookup_.clear();
  DCidlookup_.clear();
  phaseidlookup_.clear();
  mic2phase_.clear();
  mic2DC_.clear();
  mic2kinetic_.clear();
  kinetic2mic_.clear();
  mic2thermo_.clear();
  thermo2mic_.clear();
  kineticphase_.clear();
  thermophase_.clear();
  DCstoich_.clear();
  ICclasscode_.clear();
  DCclasscode_.clear();
  phaseclasscode_.clear();
  ICmolarmass_.clear();
  DCmolarmass_.clear();
  phasemolarmass_.clear();
  pH_ = 7.0;
  pe_ = 1.0;
  Eh_ = 0.0;
  T_ = 298.0;
  P_ = 101325.0;
  Vs_ = Ms_ = 1.0;
  Gs_ = Ms_ = 0.0;
  nodestatus_ = nodehandle_ = iterdone_ = 0;
  sattack_time_ = 1.0e10;
  leach_time_ = 1.0e10; 
  ICname_.clear();
  DCname_.clear();
  phasename_.clear();
  ICidlookup_.clear();
  DCidlookup_.clear();
  phaseidlookup_.clear();
  micphasevolfrac_.clear();
  micphasevolume_.clear();
  micphasemass_.clear();
  micphasemassdissolved_.clear();
  SI_.clear();
    
  node_ = new TNode();
   
  ///
  /// Initialize the thermodynamic system for both hydrates and solution 
  /// in order to initialize phasevolume_ 
  ///

  char *cGEMfilename = (char*)GEMfilename.c_str();
  char *cGEMdbrname = (char*)GEMdbrname.c_str();
  cout << "Trying to read chemical system definition file " << cGEMfilename << endl;
  try {
    gemflag = node_->GEM_init(cGEMfilename);
    if (gemflag == 1) {
        exmsg = "Bad return from GEM_init: " + GEMfilename + " missing or corrupt";
        throw GEMException("ChemicalSystem","ChemicalSystem",exmsg);
    }
    if (gemflag == -1) {
        exmsg = "Bad return from GEM_init: internal memory allocation error";
        throw GEMException("ChemicalSystem","ChemicalSystem",exmsg);
    }
  }
  catch (GEMException e) {
    e.printException();
    exit(0);
  }
    
  /// 
  /// Determine the number of possible ICs, DCs, and phases from the
  /// GEM CSD input that was read by GEM-IPM during initialization
  ///
  
  ICnum_ = (unsigned int)((node_->pCSD())->nIC);
  DCnum_ = (unsigned int)((node_->pCSD())->nDC);
  phasenum_ = (unsigned int)((node_->pCSD())->nPH);
  solutionphasenum_ = (unsigned int)((node_->pCSD())->nPS);
    
  ///
  /// Knowing the dimensions, allocate the memory for all the arrays that
  /// must be created to store thermodynamic calculation results and communicate
  /// them to the microstructure
  ///

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
    exmsg = "phasemoles_";
    phasemoles_ = new double [phasenum_];
    exmsg = "solutphasemoles_";
    solutphasemoles_ = new double [phasenum_];
    exmsg = "ophasemoles_";
    ophasemoles_ = new double [phasenum_];
    exmsg = "phasevolume_";
    phasevolume_ = new double [phasenum_];
    exmsg = "solutphasevolume_";
    solutphasevolume_ = new double [phasenum_];
    exmsg = "ophasevolume_";
    ophasevolume_ = new double [phasenum_];
    exmsg = "phasemass_";
    phasemass_ = new double [phasenum_];
    exmsg = "solutphasemass_";
    solutphasemass_ = new double [phasenum_];
    exmsg = "ophasemass_";
    ophasemass_ = new double [phasenum_];
    exmsg = "surfacearea_";
    surfacearea_ = new double [phasenum_];
    exmsg = "carrier_";
    carrier_ = new double [solutionphasenum_];
    exmsg = "DCupperlimit_";
    DCupperlimit_ = new double [DCnum_];
    exmsg = "DClowerlimit_";
    DClowerlimit_ = new double [DCnum_];
    exmsg = "phasestoich_";
    phasestoich_ = new double [phasenum_ * ICnum_];
    exmsg = "solidstoich_";
    solidstoich_ = new double [ICnum_];
    exmsg = "solutphasestoich_";
    solutphasestoich_ = new double [phasenum_ * ICnum_];
    exmsg = "solutsolidstoich_";
    solutsolidstoich_ = new double [ICnum_];
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
  /// Ask to run with automatic initial approximation.
  ///
  /// This starts the thermodynamic calculation and returns the results, including
  /// the ionic strength, pH, IC chemical potentials, DC moles, phase moles, phase
  /// volumes, and other results of the calculation.  All of these parameters
  /// are loaded into the THAMES vectors that keep track of these things, since,
  /// they were passed to the GEM calculation by reference.
  /// 
  /// Refer to GEM-IPM documentation for details about different ways to run the
  /// thermodynamic calculations
  ///
 
  try {
    (node_->pCNode())->NodeStatusCH = NEED_GEM_AIA;
    cout << "ChemicalSystem::Constructor: Entering GEM_run with node status = "
         << nodestatus_ << endl;
    cout.flush();
    nodestatus_ = node_->GEM_run(false);
    cout << "ChemicalSystem::Constructor: Exited GEM_run with node status = "
         << nodestatus_ << endl;
    cout.flush();
    if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
        exmsg = "ERROR:    Call to GEM_run failed...";
        throw GEMException("ChemicalSystem","ChemicalSystem",exmsg);
    }
    cout << "ChemicalSystem::Constructor: Entering GEM_restore_MT... " << endl;
    cout.flush();
    node_->GEM_restore_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,&ICmoles_[0],
        &DCupperlimit_[0],&DClowerlimit_[0],&surfacearea_[0]);
    cout << "Done!" << endl;
    cout << "ChemicalSystem::Constructor: Entering GEM_to_MT... " << endl;
    cout.flush();
    node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
        Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
        &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&phasemoles_[0],
        &phasevolume_[0],&phasemass_[0],&phasestoich_[0],
        &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
    cout << "Done!" << endl;

    /*
    // reset IC moles according to prievous time step
    ifstream in10("initicmoles.dat");
    double icbuff;
    for (int i = 0; i < ICnum_; i++) {
      in10 >> icbuff;
      ICmoles_[i] = icbuff;
    }
    in10.close();

    nodestatus_ = NEED_GEM_AIA;
    node_->GEM_from_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,
          ICmoles_,DCupperlimit_,DClowerlimit_,surfacearea_,
          DCmoles_,DCactivitycoeff_);
    nodestatus_ = node_->GEM_run(false);
    node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
        Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
        &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&phasemoles_[0],
        &phasevolume_[0],&phasemass_[0],&phasestoich_[0],
        &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
    cout << "Done!" << endl;
    */  
	
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

	
  /// The results of the thermodynamic calculation are now known, and
  /// the constructor can cast them into appropriate units and set up
  /// the data structure to make correspondences between GEM and microstructure
  ///
  /// Convert all IC and DC molar masses from kg/mol to g/mol
  ///
 
  ICmolarmass_.resize(ICnum_,0.0);
  icmolarmass = (node_->pCSD())->ICmm;
  for (i = 0; i < ICnum_; i++) {
    ICmolarmass_[i] = (1000.0 * (double)(*icmolarmass));  // converts to g per mol
    icmolarmass++;
  }
  DCmolarmass_.resize(DCnum_,0.0);
  dcmolarmass = (node_->pCSD())->DCmm;
  for (i = 0; i < DCnum_; i++) {
    DCmolarmass_[i] = (1000.0 * (double)(*dcmolarmass));  // converts to g per mol
    dcmolarmass++;
  }
  string string1;
  for (i = 0; i < ICnum_; i++) {
    string1.assign((node_->pCSD())->ICNL[i]);
    ICname_.push_back(string1);
    ICidlookup_.insert(make_pair(string1,i));
  }
  for (i = 0; i < DCnum_; i++) {
    string1.assign((node_->pCSD())->DCNL[i]);
    DCname_.push_back(string1);
    DCidlookup_.insert(make_pair(string1,i));
  }
  for (i = 0; i < phasenum_; i++) {
    string1.assign((node_->pCSD())->PHNL[i]);
    phasename_.push_back(string1);
    phaseidlookup_.insert(make_pair(string1,i));
  }

  /*
  cout << "To initialize phasevolume_ and phasemass_, set DCupperlimit to be normal: " 
       << endl;
  for (int i = 0; i < DCnum_; i++) {
    cout << DCname_[i] << ": " << DCupperlimit_[i] << endl;
  }
  */

  ///
  /// Set up the stoichiometry matrix for dependent components (DCs) in terms
  /// of independent components (ICs).  This is the GEM CSD A matrix
  ///
 
  vector<double> scplaceholder;
  scplaceholder.clear();
  scplaceholder.resize(ICnum_,0);
  DCstoich_.resize(DCnum_,scplaceholder);
  amat = (node_->pCSD())->A;
  for (i = 0; i < DCnum_; i++) {
    for (j = 0; j < ICnum_; j++) {
      DCstoich_[i][j] = (double)(*amat);
      amat++;
    } 
  }

  ///
  /// Set up the stoichiometry and molar masses of the GEM CSD phases
  ///

  setPhasestoich();
  setVphasestoich();
  setPhasemass();
  setPhasevolume();
  setPhasemolarmass();

  ///
  /// Set up the class codes for ICs, DCs, and phases, based on the type of
  /// component they are.  Refer to the documentation for these individual members
  /// for more detailed information about allowable values of the class codes
  ///

  ICclasscode_.resize(ICnum_,' ');
  cc = (node_->pCSD())->ccIC;
  for (i = 0; i < ICnum_; i++) {
    ICclasscode_[i] = *cc;
    cc++;
  }

  DCclasscode_.resize(DCnum_,' ');
  cc = (node_->pCSD())->ccDC;
  for (i = 0; i < DCnum_; i++) {
    DCclasscode_[i] = *cc;
    cc++;
  }

  phaseclasscode_.resize(phasenum_,' ');
  cc = (node_->pCSD())->ccPH;
  for (i = 0; i < phasenum_; i++) {
    phaseclasscode_[i] = *cc;
    cc++;
  }

  ///
  /// Begin parsing the chemistry input XML file
  ///

  string msg;
  string xmlext = ".xml";
  size_t foundxml = Interfacefilename.find(xmlext);
  try {
    if (foundxml != string::npos) {
      parseDoc(Interfacefilename);
        
      micphasevolfrac_.resize(micphasenum_,0.0);
      micphasevolume_.resize(micphasenum_,0.0);
      micphasemass_.resize(micphasenum_,0.0);
      cout << " Setting micphasemass size to " << micphasenum_ << endl;
      cout.flush();
      micphasemassdissolved_.resize(micphasenum_,0.0);
    } else {
      msg = "Not an XML file";
      throw FileException("ChemicalSystem","ChemicalSystem",Interfacefilename,msg);
    }
  }
  catch (FileException e) {
    e.printException();
    exit(0);
  }
    
  ///
  /// Set up the main map that correlates microstructure phases with GEM CSD phases
  ///

  mictotinitvolume_ = 0.0;
  for (register unsigned int i = 0; i < micphasenum_; i++) {
    mic2phase_.insert(make_pair((int)i,micphasemembers_[i]));
    mic2DC_.insert(make_pair((int)i,micDCmembers_[i]));
  }

  ///
  /// Set up the vector of saturation indices for each GEM phase.
  /// This determines the driving force for growth and is also used in calculations
  /// of the crystallization pressure during external sulfate attack

  setSI();
  vector<double> SIforsystem = getSI();

  /*
  cout << "print out SI for each phase based on whole system: " << endl;
  for (int j = 0; j < phasenum_; j++) {
    cout << phasename_[j] << ": " << SIforsystem[j] << endl;
  }
  */

  ///
  /// Set up all the information for the composition of the aqueous solution
  ///

  vector<double> solutionICmoles = getSolution();

  /*
  cout << "print out moles in solution: " << endl;
  */
 
  solut_->setICmoles(solutionICmoles);
  solut_->calculateState(true);
  vector<double> solutionSI = solut_->getSI();

  /*
  cout << "print out SI for each phase based on solution: " << endl;
  for (int j = 0; j < phasenum_; j++) {
    cout << phasename_[j] << ": " << solutionSI[j] << endl;
  }
  cout << "print out concentration in solution: " << endl;
  for (int j = 0; j < DCnum_; j++) {
    char dd = getDCclasscode(j);
    if (dd == 'S' || dd == 'T' || dd == 'W') {
      cout << (solut_->getNode())->Get_cDC(j) << endl;
    }
  }
  */
}

vector<double> ChemicalSystem::getSolution ()
{
  ///
  /// Get IC moles for solution, which fully characterizes the
  /// composition of the aqueous solution
  ///
 
  double watermass = DCmoles_[getDCid("H2O@")] 
                 * DCmolarmass_[getDCid("H2O@")];

  cout << "water mass at the end of hydration is: " << watermass << endl;
  vector<double> tempicmoles;
  tempicmoles.clear();
  tempicmoles.resize(ICnum_,0.0);
  for (register unsigned int i = 0; i < DCnum_; i++) {
    char cc = getDCclasscode(i);
    if (cc == 'S' || cc == 'T') {
      double moles = node_->Get_cDC(i) * watermass * 1.0e-3;
      for(int j = 0; j < (ICnum_ - 1); j++) {
        tempicmoles[j] += moles * DCstoich_[i][j];
      }
    }
  }

  // Treat H2O separately

  double watermoles = DCmoles_[getDCid("H2O@")];
  for (int j = 0; j < ICnum_; j++) {
    if (ICname_[j] == "H") tempicmoles[j] += watermoles * 2;
    if (ICname_[j] == "O") tempicmoles[j] += watermoles;
  }

  return tempicmoles;  

}

void ChemicalSystem::parseDoc (const string &docname)
{
    string msg;
    PhaseData phasedata;
    xmlDocPtr doc;
    xmlChar *key;
    xmlNodePtr cur;
    cout.flush();
    doc = xmlParseFile(docname.c_str());

    /// Check if the xml file is valid

    string rxcsd = xsd_files_path;
    rxcsd+="/chemistry.xsd";
    cout << "Chemistry xsd file is at " << rxcsd << endl;
    if(!is_xml_valid(doc,rxcsd.c_str())) {
        cout << "Chemistry xml is NOT valid" <<endl;
        cout.flush();
    } else {
        cout << "Chemistry xml IS valid" << endl;
        cout.flush();
    }

    if (doc == NULL ) {
        msg = "XML file not parsed successfully";
        throw FileException("ChemicalSystem","parseDoc",docname,msg);
    }

    cur = xmlDocGetRootElement(doc);

    if (cur == NULL) {
        msg = "XML file is empty";
        xmlFreeDoc(doc);
        throw FileException("ChemicalSystem","parseDoc",docname,msg);
    }

    ///
    /// Go through the XML file one tag at a time
    ///
    /// The interface file contains information about each microstructure
    /// phase that is defined, including the list of GEM CSD phases that
    /// are to be associated with that phase, the phase's internal porosity,
    /// dissolved impurities, and visualization properties.
    ///

    cur = cur->xmlChildrenNode;
    int testnumentries;
    while (cur != NULL) {
        cout << "Key name = " << cur->name << endl;
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"numentries"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testnumentries,st);
            cout << "Num entries in interface xml file is " << testnumentries << endl;
            xmlFree(key);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {
            cout << "Preparing to parse a phase from interface xml file..." << endl;
            parsePhase(doc, cur, testnumentries, phasedata);
            cout << "Done with phase parse." << endl;
        }
        cur = cur->next;
    }

    xmlFreeDoc(doc);
    return;
}

void ChemicalSystem::parsePhase (xmlDocPtr doc,
                                 xmlNodePtr cur,
                                 int numentries,
                                 PhaseData &phasedata)
{
    xmlChar *key;
    int proposedgemphaseid,proposedgemdcid;

    cout << "In function ChemicalSystem::parsePhase" << endl;
    phasedata.gtmplt.clear();
    phasedata.atmpvec.clear();
    phasedata.atmpvec.resize(numentries,0);
    phasedata.gemphaseid.clear();
    phasedata.gemdcid.clear();
    phasedata.gemphasename.clear();
    phasedata.gemdcname.clear();
    phasedata.gemphaseid.clear();
    phasedata.gemdcid.clear();

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        cout << "    Key name = " << cur->name << endl;
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"id"))) {
            cout << "    Trying to get phase id" << endl;
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.id,st);
            cout << "    Phase id = " << phasedata.id << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
            cout << "    Trying to get phase name" << endl;
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            phasedata.thamesname = st;
            cout << "    Phase name = " << phasedata.thamesname << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gemphase_data"))) {
            cout << "    Parsing gemphase_data..." << endl;
            parseGEMphasedata(doc, cur, phasedata);
            cout << "    Done parsing gemphase_data." << endl;
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"porosity"))) {
            cout << "    Trying to get porosity" << endl;
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.porosity,st);
            cout << "    Phase porosity = " << phasedata.porosity << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"display_data"))) {
            cout << "    Parsing display data..." << endl;
            parseDisplaydata(doc, cur, phasedata);
            cout << "    Done parsing display data." << endl;
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"impurity_data"))) {
            cout << "    Parsing impurity data..." << endl;
            parseImpuritydata(doc, cur, phasedata);
            cout << "    Done parsing impurity data." << endl;
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"interface_data"))) {
            cout << "    Parsing interface data..." << endl;
            parseInterfacedata(doc, cur, phasedata);
            cout << "    Done parsing interface data." << endl;
        }
        cur = cur->next;
    }

    micphasename_.push_back(phasedata.thamesname);
    micid_.push_back(phasedata.id);
    micidlookup_.insert(make_pair(phasedata.thamesname,phasedata.id));
    randomgrowth_.push_back(phasedata.randomgrowth);
    growthtemplate_.push_back(phasedata.gtmplt);
    affinity_.push_back(phasedata.atmpvec);
    porosity_.push_back(phasedata.porosity);
    grayscale_.push_back(phasedata.gray);
    color_.push_back(phasedata.colors);
    k2o_.push_back(phasedata.k2o);
    na2o_.push_back(phasedata.na2o);
    mgo_.push_back(phasedata.mgo);
    so3_.push_back(phasedata.so3);
    // if (phasedata.gemphaseid.size() == 0) phasedata.gemphaseid.push_back(0);
    micphasemembers_.insert(make_pair(phasedata.id,phasedata.gemphaseid));
    // if (phasedata.gemdcid.size() == 0) phasedata.gemdcid.push_back(0);
    micDCmembers_.insert(make_pair(phasedata.id,phasedata.gemdcid));
    micphasenum_++;

    return;
}

void ChemicalSystem::parseGEMphasedata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    int gemphaseid = 0;
    phasedata.gemphasedcmembers.clear();

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemphasename"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phasedata.gemphasename.push_back((char *)key);
            cout << "    GEM phase name = " << (char *)key << endl;
            gemphaseid = getPhaseid((char *)key);
            cout << "    GEM phase id = " << gemphaseid << endl;
            phasedata.gemphaseid.push_back(gemphaseid);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemdcname"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phasedata.gemdcname.push_back((char *)key);
            cout << "    GEM DC name = " << (char *)key << endl;
            int dcid = getDCid((char *)key);
            phasedata.gemdcid.push_back(dcid);
            cout << "    GEM DC id = " << dcid << endl;
            phasedata.gemphasedcmembers.push_back(dcid);
            xmlFree(key);
        }
        cur = cur->next;
    } 
    phaseDCmembers_.insert(make_pair(gemphaseid,phasedata.gemphasedcmembers));   
}

void ChemicalSystem::parseDisplaydata (xmlDocPtr doc,
                                       xmlNodePtr cur,
                                       PhaseData &phasedata)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;
    double red,green,blue;

    red = green = blue = 0.0;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"red"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(red,st);
            cout << "        red = " << red << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"green"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(green,st);
            cout << "        green = " << green << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"blue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(blue,st);
            cout << "        blue = " << blue << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gray"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.gray,st);
            cout << "        gray = " << phasedata.gray << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }

    phasedata.colors.clear();
    phasedata.colors.push_back(red);
    phasedata.colors.push_back(green);
    phasedata.colors.push_back(blue);

    return;
}

void ChemicalSystem::parseImpuritydata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.k2o,st);
            cout << "        k2ocoeff = " << phasedata.k2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"na2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.na2o,st);
            cout << "        na2ocoeff = " << phasedata.na2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"mgocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.mgo,st);
            cout << "        mgocoeff = " << phasedata.mgo << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"so3coeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.so3,st);
            cout << "        so3coeff = " << phasedata.so3 << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseInterfacedata (xmlDocPtr doc,
                                         xmlNodePtr cur,
                                         PhaseData &phasedata)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int testtemplate;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"randomgrowth"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.randomgrowth,st);
            cout << "        randomgrowth = " << phasedata.randomgrowth << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"growthtemplate"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testtemplate,st);
            phasedata.gtmplt.push_back(testtemplate);
            cout << "        growth template = " << testtemplate << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinity"))) {
            cout << "        Parsing affinity data..." << endl;
            parseAffinitydata(doc,cur,phasedata);
            cout << "        Done parsing affinity data." << endl;
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseAffinitydata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int testaftyid,testaftyval;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityphaseid"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testaftyid,st);
            cout << "            affinity id = " << testaftyid << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityvalue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testaftyval,st);
            cout << "            affinity value = " << testaftyval << endl;
            xmlFree(key);
        }
        cur = cur->next;
    }
    phasedata.atmpvec[testaftyid] = ((int)testaftyval);
    return;
}

ChemicalSystem::ChemicalSystem (const ChemicalSystem &obj)
{

    ///
    /// This is a straightforward copy constructor
    /// Each member and vector is copied into the newly constructed object
    ///

    micphasenum_ = obj.getMicphasenum();
    micimpuritynum_ = obj.getMicimpuritynum();
    ICnum_ = obj.getICnum();
    DCnum_ = obj.getDCnum();
    phasenum_ = obj.getPhasenum();
    solutionphasenum_ = obj.getSolutionphasenum();
    micphasename_ = obj.getMicphasename();
    ICname_ = obj.getICname();
    DCname_ = obj.getDCname();
    phasename_ = obj.getPhasename();
    micid_ = obj.getMicid();
    c3sid_ = obj.getC3sid();
    c2sid_ = obj.getC2sid();
    c3aid_ = obj.getC3aid();
    c4afid_ = obj.getC4afid();
    gypsumid_ = obj.getGypsumid();
    randomgrowth_ = obj.getRandomgrowth();
    ICmoles_ = obj.getICmoles();
    DCmoles_ = obj.getDCmoles();
    ICmolarmass_ = obj.getICmolarmass();
    DCmolarmass_ = obj.getDCmolarmass();
    growthtemplate_ = obj.getGrowthtemplate();
    affinity_ = obj.getAffinity();
    micphasemembers_ = obj.getMicphasemembers();
    micDCmembers_ = obj.getMicDCmembers();
    porosity_ = obj.getPorosity();
    k2o_ = obj.getK2o();
    na2o_ = obj.getNa2o();
    mgo_ = obj.getMgo();
    so3_ = obj.getSo3();
    grayscale_ = obj.getGrayscale();
    color_ = obj.getColor();
    micidlookup_ = obj.getMicidlookup();
    ICidlookup_ = obj.getICidlookup();
    DCidlookup_ = obj.getDCidlookup();
    phaseidlookup_ = obj.getPhaseidlookup();
    mic2phase_ = obj.getMic2phase();
    mic2DC_ = obj.getMic2DC();
    mic2kinetic_ = obj.getMic2kinetic();
    kinetic2mic_ = obj.getKinetic2mic();
	mic2thermo_ = obj.getMic2thermo();
	thermo2mic_ = obj.getThermo2mic();
    kineticphase_ = obj.getKineticphase();
	thermophase_ = obj.getThermophase();
    ICclasscode_ = obj.getICclasscode();
    DCclasscode_ = obj.getDCclasscode();
    phaseclasscode_ = obj.getPhaseclasscode();
    DCstoich_ = obj.getDCstoich();
    phasestoich_ = obj.getPhasestoich();
    vphasestoich_ = obj.getVphasestoich();
    ICresiduals_ = obj.getICresiduals();
    ICchempot_ = obj.getICchempot();
    DCactivitycoeff_ = obj.getDCactivitycoeff();
    phasemoles_ = obj.getPhasemoles();
    ophasemoles_ = obj.getOphasemoles();
    phasemass_ = obj.getPhasemass();
    ophasemass_ = obj.getOphasemass();
    phasevolume_ = obj.getPhasevolume();
    ophasevolume_ = obj.getOphasevolume();
    carrier_ = obj.getCarrier();
    surfacearea_ = obj.getSurfacearea();
    DClowerlimit_ = obj.getDClowerlimit();
    DCupperlimit_ = obj.getDCupperlimit();
    /*
    node_ = obj.getNode();
    */
    T_ = obj.getT();
    P_ = obj.getP();
    Vs_ = obj.getVs();
    Ms_ = obj.getMs();
    pH_ = obj.getPH();
    pe_ = obj.getPe();
    Eh_ = obj.getEh();
    ionicstrength_ = obj.getIonicstrength();
    Gs_ = obj.getGs();
    Hs_ = obj.getHs();
    nodehandle_ = obj.getNodehandle();
    nodestatus_ = obj.getNodestatus();
    iterdone_ = obj.getIterdone();
    micphasevolfrac_ = obj.getMicphasevolfrac();
    micphasevolume_ = obj.getMicphasevolume();
    mictotvolume_ = obj.getMictotvolume();
    mictotinitvolume_ = obj.getMictotinitvolume();
    micphasemass_ = obj.getMicphasemass();
    micphasemassdissolved_ = obj.getMicphasemassdissolved();
    micvoidvolume_ = obj.getMicvoidvolume();
    micvoidvolfrac_ = obj.getMicvoidvolfrac();
}

ChemicalSystem::~ChemicalSystem ()
{
    ///
    /// Clear out the maps
    ///

    micidlookup_.clear();
    DCidlookup_.clear();
    ICidlookup_.clear();
    phaseidlookup_.clear();
    mic2phase_.clear();
    mic2DC_.clear();
    mic2kinetic_.clear();
    kinetic2mic_.clear();
	mic2thermo_.clear();
	thermo2mic_.clear();

    ///
    /// Clear out the vectors
    ///

    micphasename_.clear();
    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    micid_.clear();
    randomgrowth_.clear();
    DCstoich_.clear();
    growthtemplate_.clear();
    affinity_.clear();
    micphasemembers_.clear();
    micphasemembersvolfrac_.clear();
    micphasemass_.clear();
    micphasemassdissolved_.clear();
    micDCmembers_.clear();
    porosity_.clear();
    k2o_.clear();
    na2o_.clear();
    mgo_.clear();
    so3_.clear();
    grayscale_.clear();
    color_.clear();
    ICclasscode_.clear();
    DCclasscode_.clear();
    phaseclasscode_.clear();
    kineticphase_.clear();
	thermophase_.clear();
    vphasestoich_.clear();

    ///
    /// Free up the dynamically allocated memory
    ///

    delete[]DClowerlimit_;
    delete[]DCupperlimit_;
    delete[]surfacearea_;
    delete[]ophasemass_;
    delete[]ophasevolume_;
    delete[]ophasemoles_;
    delete[]phasemass_;
    delete[]phasevolume_;
    delete[]carrier_;
    delete[]phasemoles_;
    delete[]DCactivitycoeff_;
    delete[]DCmoles_;
    delete[]ICchempot_;
    delete[]ICresiduals_;
    delete[]ICmoles_;
    delete[]phasestoich_;

    delete node_;
}

void ChemicalSystem::getGEMPhasestoich ()
{
    double *arout = new double[ICnum_];
    for (long int i = 0; i < phasenum_; i++) {
        arout = node_->Ph_BC(i,arout);
        for (register unsigned int j = 0; j < ICnum_; j++) {
            phasestoich_[(i * ICnum_) + j] = arout[j];
        }
    }
    delete[] arout;
}

void ChemicalSystem::getGEMVphasestoich ()
{
    double minval = 0.0;
    vphasestoich_.clear();
    vector<double> vplace;
    vplace.clear();
    vplace.resize(ICnum_,0.0);
    vphasestoich_.resize(phasenum_,vplace);
    int indexval,oval;
    for (register unsigned int i = 0; i < phasenum_; i++) {
        if (phasename_[i] == "aq_gen") {

            ///
            /// Normalize to one mole of oxygen
            ///

            oval = (i * ICnum_) + getICid("O");
            if (phasestoich_[oval] > 0.0) {
                for (register unsigned int j = 0; j < ICnum_; j++) {
                    indexval = (i * ICnum_) + j;
                    vphasestoich_[i][j] = (phasestoich_[indexval]/phasestoich_[oval]);
                }
            }
        } else {
            for (register unsigned int j = 0; j < ICnum_; j++) {
                indexval = (i * ICnum_) + j;
                vphasestoich_[i][j] = (phasestoich_[indexval]/minval);
            }
        }
    }
}

void ChemicalSystem::writeDb (ostream &stream)
{
    register unsigned int i;

    ///
    /// Make the header
    ///

    stream << "--------------------------------------------------------" << endl;
    stream << "CONTENTS OF PHASE DATABASE:" << endl;
    stream <<  endl;

    ///
    /// Format one line at a time
    ///

    for (i = 0; i < micphasenum_; i++) {
        writeMember(i,stream);
    }

    stream << endl;
    stream << "--------------------------------------------------------" << endl;
    stream << endl;
}

void ChemicalSystem::writeMember (const unsigned int i,
                                  ostream &stream)
{

    unsigned int idnum;
    if (i >= micphasenum_) {
        throw EOBException("ChemicalSystem","writeMember","micphasename_",micphasenum_,i);
    }

    ///
    /// Format the output for one phase
    ///

    stream << "------------------------------------------------------" << endl;
    stream << "DATA FOR MATERIAL " << i << ":" << endl;
    stream << "       Name = " << micphasename_[i] << endl;
    stream << "         Id = " << micid_[i] << endl;
    stream << "   Porosity = " << porosity_[i] << endl;
    stream << "------------------------------------------------------" << endl;
}

void ChemicalSystem::writeChemSys ()
{
    register unsigned int j;

    ///
    /// First we will list details for the ICs
    ///

    string CSfilename("chemsys.report");
    ofstream out(CSfilename.c_str());
    out << "Report on the Material Database" << endl;
    out << "-------------------------------" << endl << endl;
    out << "List of Independent Components:" << endl << endl;
    for(register unsigned int i = 0; i < ICnum_; i++) {
        out << i << ")            Name: " << ICname_[i] << endl;
        out << "        classcode: " << ICclasscode_[i] << endl;
        out << "       molar mass: " << ICmolarmass_[i] << endl << endl;
    }

    out << "List of Dependent Components:" << endl << endl;
    for(register unsigned int i = 0; i < DCnum_; i++) {
        out << i << ")            Name: " << DCname_[i] << endl;
        out << "        classcode: " << DCclasscode_[i] << endl;
        out << "       molar mass: " << DCmolarmass_[i] << endl << endl;
    }

    out << "List of Phases:" << endl << endl;
    for(register unsigned int i = 0; i < phasenum_; i++) {
        out << i << ")            Name: " << phasename_[i] << endl;
        out << "        classcode: " << phaseclasscode_[i] << endl;
    }

    out << "List of Microstructure Phases:" << endl << endl;
    for(register unsigned int i = 0; i < micphasenum_; i++) {
        out << i << ")       Name: " << micphasename_[i] << endl;
        out << "               id: " << micid_[i] << endl;
        out << "    random growth: " << randomgrowth_[i] << endl;
        for (j = 0; j < affinity_[i].size(); j++) {
            out << "        affinity to " << j << ": " << affinity_[i][j] << endl;
        }
        cout << "the size of growthtemplate_[" << i << "] = "
             << growthtemplate_[i].size() << endl;
        for (j = 0; j < growthtemplate_[i].size(); j++) {
            out << "        growthtemplate: " << growthtemplate_[i][j] << endl;
        }
        out << "         porosity: " << porosity_[i] << endl;
        out << "              k2o: " << k2o_[i] << endl;
        out << "             na2o: " << na2o_[i] << endl;
        out << "              mgo: " << mgo_[i] << endl;
        out << "              so3: " << so3_[i] << endl;
    }

    out.close();
    return;
}

void ChemicalSystem::writeChemSys (ostream &out)
{
    register unsigned int j;

    ///
    /// First we will list details for the ICs
    ///

    out << "Report on the Material Database" << endl;
    out << "-------------------------------" << endl << endl;
    out << "List of Independent Components:" << endl << endl;
    for(register unsigned int i = 0; i < ICnum_; i++) {
        out << i << ")            Name: " << ICname_[i] << endl;
        out << "        classcode: " << ICclasscode_[i] << endl;
        out << "       molar mass: " << ICmolarmass_[i] << endl << endl;
    }

    out << "List of Dependent Components:" << endl << endl;
    for(register unsigned int i = 0; i < DCnum_; i++) {
        out << i << ")            Name: " << DCname_[i] << endl;
        out << "        classcode: " << DCclasscode_[i] << endl;
        out << "       molar mass: " << DCmolarmass_[i] << endl << endl;
    }

    out << "List of Phases:" << endl << endl;
    for(register unsigned int i = 0; i < phasenum_; i++) {
        out << i << ")            Name: " << phasename_[i] << endl;
        out << "        classcode: " << phaseclasscode_[i] << endl;
    }

    out << "List of Microstructure Phases:" << endl << endl;
    for(register unsigned int i = 0; i < micphasenum_; i++) {
        out << i << ")       Name: " << micphasename_[i] << endl;
        out << "               id: " << micid_[i] << endl;
        out << "    random growth: " << randomgrowth_[i] << endl;
        for (j = 0; j < affinity_[i].size(); j++) {
            out << "        affinity to " << j << ": " << affinity_[i][j] << endl;
        }
        for (j = 0; j < growthtemplate_[i].size(); j++) {
            out << "        growthtemplate: " << growthtemplate_[i][j] << endl;
        }
        out << "         porosity: " << porosity_[i] << endl;
        out << "              k2o: " << k2o_[i] << endl;
        out << "             na2o: " << na2o_[i] << endl;
        out << "              mgo: " << mgo_[i] << endl;
        out << "              so3: " << so3_[i] << endl;
    }

    return;
}

int ChemicalSystem::calculateState (double time,
                                    bool isfirst = false)
{
    int status = 0;
    string msg;
 
    /*
    isfirst = true; 
    */
      
    vector<double> oDCmoles;
    oDCmoles.clear();
    oDCmoles.resize(DCnum_,0.0);
  
    for (int i = 0; i < DCnum_; i++) {
      oDCmoles[i] = DCmoles_[i];
    }
 
    nodestatus_ = NEED_GEM_SIA;
    cout << "    Going into ChemicalSystem::calculateState::GEM_from_MT... " << endl;
    cout.flush();
    node_->GEM_from_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,
          ICmoles_,DCupperlimit_,DClowerlimit_,surfacearea_,
          DCmoles_,DCactivitycoeff_);
    cout << "Done!" << endl;
    
    cout << "    Going into ChemicalSystem::calculateState::GEM_set_MT... ";
    cout.flush();
    node_->GEM_set_MT(time,1.0);
    cout << "Done!" << endl
         << "    Going into ChemicalSystem::calculateState::GEM_run(true)... ";
    cout.flush();
    writeICmoles();

    /*
    writeDCmoles();
    */

    if (isfirst) {
        cout << "going into GEM_run(false)" << endl;
        nodestatus_ = node_->GEM_run(false);
    } else {
        cout << "going into GEM_run(true)" << endl;
        nodestatus_ = node_->GEM_run(true);
    }

    cout << "Done!  nodestatus is " << nodestatus_ << endl;
    cout.flush();
    
    if (nodestatus_ == ERR_GEM_AIA || nodestatus_ == ERR_GEM_SIA
        || nodestatus_ == T_ERROR_GEM) {
        switch (nodestatus_) {
            case ERR_GEM_AIA:
                msg = "Flag nodestatus_ = ERR_GEM_AIA";
                break;
            case ERR_GEM_SIA:
                msg = "Flag nodestatus_ = ERR_GEM_SIA";
                break;
            case T_ERROR_GEM:
                msg = "Flag nodestatus_ = T_ERROR_GEM";
                break;
        }
        throw GEMException("ChemicalSystem","calculateState",msg);
    } else {
        if (nodestatus_ == BAD_GEM_AIA || nodestatus_ == BAD_GEM_SIA) {
            switch (nodestatus_) {
                case BAD_GEM_AIA:
                    msg = "Flag nodestatus_ = BAD_GEM_AIA";
                    break;
                case BAD_GEM_SIA:
                    msg = "Flag nodestatus_ = BAD_GEM_SIA";
                    break;
            }
            throw GEMException("ChemicalSystem","calculateState",msg);
        } else {
            cout << "    Going into ChemicalSystem::calculateState::GEM_to_MT... ";
            cout.flush();
            node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
                    Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
                    &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&solutphasemoles_[0],
                    &solutphasevolume_[0],&solutphasemass_[0],&solutphasestoich_[0],
                    &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
            cout << "Done!" << endl;
            cout << "after GEM_to_MT...Ms_ = " << Ms_ << endl;
            cout.flush();
      }
    }

    /*
    writePhasemoles();
    */

    mictotvolume_ = 0.0;
    setPhasestoich();
    setVphasestoich();
    setPhasemass();
    setPhasevolume();
    setPhasemolarmass();
  
    for (register unsigned int i = 1; i < micphasenum_; i++) {
        cout << "Setting micphase amounts for " << i << " = " << micphasename_[i] << endl;
        cout.flush();
        if (!isKineticphase(i)) {
            micphasemass_[i] = micphasevolume_[i] = 0.0;
            for (register unsigned int j = 0; j < micphasemembers_[i].size(); j++) {
                cout << "    Is a THERMO phase composed of "
                     << phasename_[micphasemembers_[i][j]]
                     << " having mass = " << phasemass_[micphasemembers_[i][j]]
                     << " and volume = " << phasevolume_[micphasemembers_[i][j]] << endl;
                cout.flush();
                mictotvolume_ += phasevolume_[micphasemembers_[i][j]];
                micphasemass_[i] += phasemass_[micphasemembers_[i][j]];
                micphasevolume_[i] += phasevolume_[micphasemembers_[i][j]];
            }
        } else if (isKineticphase(i)) {
            cout << "    Is a KINETIC phase composed of "
                 << phasename_[micphasemembers_[i][0]]
                 << "  having mass = " << micphasemass_[i]
                 << "  and volume = " << micphasevolume_[i] << endl;
            cout.flush();

            ///
            /// micphasemass and micphasevolume need to be set elsewhere???
            ///

            mictotvolume_ += micphasevolume_[i];
        } else {
            micphasevolume_[i] = mictotinitvolume_ * micphasevolfrac_[i];
            cout << "micphasevolfrac_[DAMAGE] is: " << micphasevolfrac_[i] << endl;
            cout << "    Is a DAMAGE phase having volume = "
                 << micphasevolume_[i] << endl;
            cout.flush();
        }
    }

    if (isfirst) mictotinitvolume_ = mictotvolume_;
    cout << "isfirst = " << isfirst << endl;

    // Now we take care of the void volume separately
    cout << "now use water to keep total volume constantly." << endl;
  
    if (mictotinitvolume_ > mictotvolume_) {
        micphasevolume_[1] += mictotinitvolume_ - mictotvolume_;
        double water_molarv, water_molesincr;
        for (register int i = 0; i < micphasenum_; i++) {
            if (micphasename_[i] == "H2O") {
                water_molarv = node_->DC_V0(getMic2DC(i,0), P_, T_);
                water_molesincr = (mictotinitvolume_ - mictotvolume_) / water_molarv;
                cout << "water_molarv = " << water_molarv << endl;
                cout << "volume increase of water is: "
                     << (mictotinitvolume_ - mictotvolume_) << endl;
                cout << "water_molesincr = " << water_molesincr << endl;
            }
        }
        for (register int i = 0; i < ICnum_; i++) {
            if (ICname_[i] == "H") ICmoles_[i] += water_molesincr * 2.0;
            if (ICname_[i] == "O") ICmoles_[i] += water_molesincr;
        }
        mictotvolume_ += mictotinitvolume_ - mictotvolume_;
    }

    ///
    /// Manually adjust volume for CSH
    ///

    int CSHID = getMicid("CSH");
    int CSHgemphaseid = micphasemembers_[CSHID][0];
    vector<int> CSHgemDCid = getPhaseDCmembers(CSHgemphaseid);
    double CSH_corrV = 0.0;
    for (int i = 0; i < CSHgemDCid.size();i++) {
        double v = 0.0;
        if (CSHgemDCid[i] == 79) {
            cout << "correcting the volume for " << DCname_[79] << endl;
            v = DCmoles_[79] * 0.0001128;
        } else {
            v = DCmoles_[CSHgemDCid[i]] * node_->DC_V0(CSHgemDCid[i],P_,T_);
        }
        CSH_corrV += v;
    }
    cout << "before correction, the volume of CSH is: " << micphasevolume_[CSHID]
         << endl;
    cout << "after correction, the volume of CSH is: " << CSH_corrV << endl;
    micphasevolume_[1] = micphasevolume_[1] - (CSH_corrV - micphasevolume_[CSHID]);
    micphasevolume_[CSHID] = CSH_corrV;  

    ///
    /// End of manual adjustment  
    ///

    micphasevolume_[0] = 0.0;
    cout << "total volume of the microstructure is: " << mictotvolume_ << endl;
    if (mictotvolume_ > 0.0) {
        for (register unsigned int i = 0; i < micphasenum_; i++) {
            /*
            micphasevolfrac_[i] = micphasevolume_[i] / mictotvolume_;
            */
            micphasevolfrac_[i] = micphasevolume_[i] / mictotinitvolume_;
            cout << "    " << micphasename_[i] << "  V = " << micphasevolume_[i] << ", Vt = "
                 << mictotvolume_ << ", and volume fraction = " << micphasevolfrac_[i] << endl;
            cout.flush();
        }
    } else {
        throw DataException("ChemicalSystem","calculateState","mictotvolume_ is NOT positive");
    }
  
    setPhasestoich();
   
    ///
    /// Calculate driving force for ettringite growth
    ///

    double *soluticmoles;
    soluticmoles = solut_->getICmoles();
    for (int i = 0; i < DCnum_; i++) {
        char cc;
        cc = getDCclasscode(i);
        if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
            double dissolveDCmoles = oDCmoles[i] - DCmoles_[i];
            if (dissolveDCmoles > 0.0) {
                for (int j = 0; j < (ICnum_ - 1); j++) {
                    soluticmoles[j] += dissolveDCmoles * DCstoich_[i][j];
                }
            }
        }
    }
    for (int i = 0; i < ICnum_; i++) {
        solut_->setICmoles(i, soluticmoles[i]);
    }
    solut_->calculateState(true);

    ///
    /// Update solution
    ///

    vector<double> solutICmoles = getSolution();
    cout << "Now update solution IC moles...";
    for (int i = 0; i < ICnum_; i++) {
        solut_->setICmoles(i, solutICmoles[i]);
    } 
    cout << "Done." << endl;

    /*
    solut_->calculateState(false);
    */

    return status;
}
