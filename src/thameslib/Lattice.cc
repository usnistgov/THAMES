/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"
#include "Interface.h"

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut)
: siteneighbors_(18),chemsys_(cs),solut_(solut)
{
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP;   // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;       // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES;     // in micrometers (see global.h)
  site_.clear();
}

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut,
                  const string &fname)
: siteneighbors_(18),chemsys_(cs),solut_(solut)
{
  register unsigned int i,j,k;
  register unsigned long int ii;
  string buff;
  int xn,yn,zn;
  unsigned int idn;
  unsigned int pid;
  string msg;
    
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP;   // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;       // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES;     // in micrometers (see global.h)
  site_.clear();
    
  ///
  /// Open the microstructure input file and process it
  ///

  ifstream in(fname.c_str());
  try {
    if (!in) {
        throw FileException("Lattice","Lattice",fname,"Could not open");
    }
  }
  catch (FileException fex) {
    fex.printException();
    exit(1);
  }
    
  in >> buff;
  if (buff == VERSIONSTRING) {
    in >> version_;
    in >> buff;
    if (buff == XSIZESTRING) {
        in >> xdim_;
        in >> buff;
        in >> ydim_;
        in >> buff;
        in >> zdim_;
    } else if (buff == IMGSIZESTRING) {
        ydim_ = zdim_ = xdim_;
    }
    in >> buff;
    if (buff == IMGRESSTRING) {
        double testres;
        in >> testres;
        setResolution(testres);
    }
  } else {

    ///
    /// Image file generated prior to VCCTL Version 3.0.
    /// Allow backward compatibility by defaulting system
    /// size to 100 and resolution to 1.0 micrometers
    ///

    version_ = "2.0";
    double testres = 1.0;
    double firsttime = true;
    setResolution(testres);
    xdim_ = ydim_ = zdim_ = 100;
  }

  ///
  /// Print out the microstructure size and characteristics
  ///

  cout << "Read microstructure file header..." << endl;
  cout.flush();
  cout << "    Version = " << version_ << endl;
  cout.flush();
  cout << "    xdim_ = " << xdim_ << endl;
  cout.flush();
  cout << "    ydim_ = " << ydim_ << endl;
  cout.flush();
  cout << "    zdim_ = " << zdim_ << endl;
  cout.flush();
  numsites_ = (unsigned long int)(xdim_ * ydim_ * zdim_);
  cout << "    numsites_ = " << numsites_ << endl;
  cout.flush();
  cout << "    resolution_ = " << resolution_ << endl;
  cout.flush();  
  
  ///
  /// Allocate a random number generator object and seed it
  ///

  try { rg_ = new RanGen(); }
    catch (bad_alloc &ba) {
        cout << "Lattice constructor failed when allocating rg_";
        cout.flush();
        exit(1);
    }

    rg_->setSeed(-142234);
    cout << "Checking whether I can get a seed " << rg_->getSeed() << endl;
    cout.flush();
    cout << "Checking value of random number " << ", " << rg_->Ran3() << endl;
    cout.flush();

    count_.clear();
    count_.resize(chemsys_->getMicphasenum(),0);

    ettrSI_ = 0.0;
    expansion_.clear();
    expansion_coordin_.clear();

    waterchange_ = 0.0;

    ///
    /// Add the Site objects to the `site_` array, defining the structure
    /// of the lattice.  Not yet assigning a phase to the site because
    /// the file has not yet been read
    ///

    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
                addSite(i,j,k);
            }
        }
    }

    ///
    /// Set up the array of neighbors for each site
    ///

    for (ii = 0; ii < numsites_; ii++) {

      for (j = 0; j < (NUM_NEAREST_NEIGHBORS
                + NUM_SECONDNEAREST_NEIGHBORS); j++) {
        ///
        /// Store up to 2nd nearest neighbors
        ///

        switch (j) {

            case 0: // West neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ();
                break;
            case 1: // East neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ();
                break;
            case 2: // South neighbor
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                xn = (int)site_[ii].getX();
                zn = (int)site_[ii].getZ();
                break;
            case 3: // North neighbor
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                xn = (int)site_[ii].getX();
                zn = (int)site_[ii].getZ();
                break;
            case 4: // Down neighbor
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY();
                break;
            case 5: // Up neighbor
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY();
                break;
            case 6: // Southwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ();
                break;
            case 7: // Northwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ();
                break;
            case 8: // Northeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ();
                break;
            case 9: // Southeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ();
                break;
            case 10: // Downwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 11: // Downnorth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 12: // Downeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 13: // Downsouth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 14: // Upwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 15: // Upnorth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 16: // Upeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 17: // Upsouth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
        }

        idn = (unsigned int)(xn + (xdim_ * yn) +
                ((xdim_ * ydim_) * zn));
        site_[ii].setNb(j,&site_[idn]);
      }
    }

    ///
    /// This loop reads the phase for each site and assigns it,
    /// also updating the count of the phase.
    ///

    for (i = 0; i < numsites_; i++) {
      in >> pid;
      site_[i].setPhaseId(pid);
      count_[pid] += 1.0;
    }
    
    ///
    /// Done with the input microstructure file, so close it.
    ///

    in.close();
 
 
    ///
    /// With the phase counts known, calculate phase volume fractions
    ///

    volumefraction_.clear();
    volumefraction_.resize(chemsys_->getMicphasenum(),0.0);
    cout << "Calculating Volume Fractions now..." << endl;
    try {
      if (site_.size() > 0) {
        for (ii = 0; ii < chemsys_->getMicphasenum(); ii++) {
            volumefraction_[ii] = ((double)count_[ii])/((double)site_.size());
            if (volumefraction_[ii] > 0.0) cout << "***Volume fraction["
                << ii << "] = " << volumefraction_[ii] << endl;
        }
      } else {
        msg = "Divide by zero error:  site_.size() = 0";
        throw FloatException("Lattice","Lattice",msg);
      }
      cout << "...Done!" << endl;
    }
    catch (FloatException fex) {
      fex.printException();
      exit(1);
    }
    
    cout << "Calculating weighted mean curvatures now..." << endl;
    for (ii = 0; ii < site_.size(); ii++) {
      site_[ii].calcWmc();
    }
    cout << "...Done!" << endl;
}

Lattice::~Lattice ()
{
    interface_.clear();
    site_.clear();
    delete rg_;
}

void Lattice::addSite (const unsigned int x,
                       const unsigned int y,
                       const unsigned int z)
{
    string msg;
    
    try {
      if ( x >= xdim_ || x < 0 ) {
        throw EOBException("Lattice","addSite","site_",xdim_,x);
      } else if ( y >= ydim_ || y < 0 ) {
        throw EOBException("Lattice","addSite","site_",ydim_,y);
      } else if ( z >= zdim_ || z < 0 ) {
        throw EOBException("Lattice","addSite","site_",zdim_,z);
      }
    }
    catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    
    ///
    /// Why is numphase assigned here?  It is a local variable and is not
    /// used.
    ///

    unsigned int numphase = chemsys_->getMicphasenum();

    site_.push_back(Site(x,y,z,xdim_,ydim_,zdim_,siteneighbors_,chemsys_));
}

void Lattice::findInterfaces ()
{
    register unsigned int i,kk;
    register unsigned long int k;
    vector<Site *> gsite,dsite;
    
    ///
    /// An interface must have at least one adjacent site that is water or void
    ///

    interface_.clear();
    for (i = 0; i < chemsys_->getMicphasenum(); i++) { 
      cout << "  Database item " << i << ": ";
      cout.flush();
      if (i != WATERID && i != VOIDID) {
        cout << "Probing for interface... ";
        gsite.clear();
        dsite.clear();
        for (k = 0; k < site_.size(); k++) {
          if (site_[k].getWmc() > 0) {
            if ((site_[k].getPhaseId() == i)) {
              dsite.push_back(&site_[k]);
              site_[k].setDissolutionSite(i);
              for (kk = 0; kk < site_[k].nbSize(2); kk++) {
                if ((site_[k].nb(kk))->getPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
            } else if (chemsys_->isGrowthtemplate(i,site_[k].getPhaseId())) {
              for (kk = 0; kk < site_[k].nbSize(1); kk++) {
                if ((site_[k].nb(kk))->getPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
            }
          }
        }

        cout << " Done! " << dsite.size()
            << " dissolution sites and " << gsite.size()
            << " growth sites" << endl;

        if ((gsite.size() == 0) && (dsite.size() == 0)) {
          cout << "Testing phase " << i << " for nucleation ";

          ///
          /// We are dealing with a phase that may need to
          /// homogeneously nucleate.  Therefore, identify the
          /// eligible nucleation sites for that phase
          ///

          double thresh = (0.5 / pow(resolution_,3.0));
          cout << "Thresh = " << thresh << endl;
          cout.flush();
          double g = 0.0;
          for (k = 0; k < site_.size(); k++) {
            if ((site_[k].getPhaseId() == WATERID)) {
              g = rg_->Ran3();
              if (g < thresh) {
                gsite.push_back(&site_[k]);
                site_[k].setGrowthSite(i);
              }
            }
          }
        }
        
        cout << "Trying to add a water interface for phase " << i << "... ";
        cout.flush();
        interface_.push_back(Interface(chemsys_,rg_,gsite,dsite,i));   
        cout << "Done!" << endl;
        cout.flush();
      } else {
        cout << "Trying to add a regular interface for phase " << i << "... ";
        cout.flush();
        interface_.push_back(Interface(rg_));   
      }
    }

    return;
}

long int Lattice::growPhase (unsigned int phaseid,
                             long int numtoadd)
{
    unsigned int pid;
    register unsigned long int i,j;
    double dwmcval;
    Site *ste,*stenb;

    try {
        if (numtoadd == 0) return 0;
        if (phaseid >= interface_.size()) {
            throw EOBException("Lattice","growPhase","interface_",
                               interface_.size(),phaseid);
        }
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }

    vector<Isite> isite = interface_[phaseid].getGrowthSites();
    cout << "size of interface_[" << phaseid << "] is " << isite.size() << endl; 
    long int numleft,numchange = 0;

    ///
    /// We need to go through the interface list in
    /// normal order, which is the reverse order that
    /// we use for dissolution.
    ///

    numleft = numtoadd;
    cout << "-->Phase " << phaseid
        << " needs to grow at " << numtoadd
        << " sites" << endl;

    while ((numleft > 0) && (isite.size() >= 1)) {
      for (i = 0; (numleft > 0) && (i < isite.size()); i++) {
        ste = &site_[isite[i].getId()];
        pid = ste->getPhaseId();
     
        if (pid == WATERID) {
          removeGrowthSite(ste,phaseid);
          /*
          vector<unsigned int> plist = ste->getGrowthPhases();
          for (int ii = 0; ii < plist.size(); ii++) {
           removeGrowthSite(ste, plist[ii]);
          } 
          */

          ///
          /// Weighted mean curvature (wmc) is changed by the difference
          /// between the growing phase's porosity and the template's porosity.
          ///
          /// @todo Determine why the calculation works this way.
          ///

          dwmcval = chemsys_->getPorosity(phaseid)
                  - chemsys_->getPorosity(pid);
          setPhaseId(ste,phaseid);
          count_[phaseid] ++;
          ste->dWmc(dwmcval);

          ///
          /// Now that the site has been added, it is eligible for dissolution
          /// later on, so we add it to the list of dissolution sites.
          ///

          if (ste->getWmc() > 0.0) {
            addDissolutionSite(ste,phaseid);
          }

          ///
          /// Update the wmc of the neighboring sites.  This can be done
          /// because the wmc is originally calculated within a box around each
          /// site, so any time the id of a site within that box changes, it
          /// will change the wmc of the site at the box's center.
          ///

          for (j = 0; j < ste->nbSize(1); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
            if (stenb->getPhaseId() == WATERID) {
              addGrowthSite(stenb,phaseid);
            } else if (stenb->getWmc() <= 0.0) {
              removeDissolutionSite(stenb,stenb->getPhaseId());
            }
          }

          for (j = ste->nbSize(1); j < ste->nbSize(2); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }

          numleft--;
          numchange++;
        } else {
          removeGrowthSite(ste,phaseid);
        }
      }
  
      isite = interface_[phaseid].getGrowthSites();
    }
    
    return(numchange);
}


long int Lattice::dissolvePhase (unsigned int phaseid,
                                 long int numtotake)
{
    unsigned int pid;
    double dwmcval;
    Site *ste,*stenb;

    if (numtotake == 0) return 0;
    try {
        if (phaseid >= interface_.size()) {
            throw EOBException("Lattice","dissolvePhase",
                           "interface_",interface_.size(),phaseid);
        }
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
    
    register unsigned long int i;
    vector<Isite> isite = interface_[phaseid].getDissolutionSites();
    cout << "size of interface_[" << phaseid << "] is " << isite.size() << endl; 
    try {
      for (i = 0; i < isite.size(); i++) {
        if (site_.at(isite[i].getId()).getPhaseId() != phaseid) {
          cout << "Uh-oh... interface " << phaseid
               << " is corrupted with phase "
               << site_.at(isite[i].getId()).getPhaseId() << endl;
          cout << "Offending site is " << isite[i].getId() << endl;
        }
      }
    }
    catch (out_of_range &oor) {
      EOBException ex("Lattice","dissolvePhase","site_",site_.size(),i);
      ex.printException();
      exit(1);
    }
    cout << "-->Phase " << phaseid
         << " needs to dissolve at " << numtotake
         << " sites" << endl;
    
    long int numleft = numtotake;
    long int numchange = 0;
    while ((numleft > 0) && (isite.size() > 1)) {
  
      try {
        for (i = isite.size() - 1; (i > 0) && (numleft > 0); i--) {
          ste = &site_.at(isite[i].getId());
          pid = ste->getPhaseId();
          removeDissolutionSite(ste,pid);
      
          setPhaseId(ste,WATERID);
          count_[WATERID] += 1.0;
      
          ///
          /// Weighted mean curvature (wmc) is changed by the difference
          /// between the growing phase's porosity and the template's porosity.
          ///
          /// @todo Determine why the calculation works this way.
          ///

          dwmcval = chemsys_->getPorosity(WATERID)
                  - chemsys_->getPorosity(pid);
          ste->dWmc(dwmcval);
      
          register unsigned long int j;
          for (j = 0; j < ste->nbSize(1); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
            if ((stenb->getPhaseId() != WATERID)
                && (stenb->getPhaseId() != VOIDID)) {
              int nbpid;
              nbpid = stenb->getPhaseId();

              ///
              /// Now that the site has been dissolved, it is eligible for growth
              /// later on, so we add it to the list of growth sites.
              ///

              addDissolutionSite(stenb,nbpid);
              addGrowthSite(ste,nbpid);
    
              vector<int> nbgrowthtemp;
              nbgrowthtemp.clear();
              nbgrowthtemp = chemsys_->getGrowthtemplate(nbpid);
              for (int ii = 0; ii < nbgrowthtemp.size(); ii++) {
                if (nbgrowthtemp[ii] != WATERID) {
                  addGrowthSite(ste,nbgrowthtemp[ii]);
                }
              }
            }
          }

          ///
          /// Update the wmc of the neighboring sites.  This can be done
          /// because the wmc is originally calculated within a box around each
          /// site, so any time the id of a site within that box changes, it
          /// will change the wmc of the site at the box's center.
          ///

          for (j = ste->nbSize(1); j < ste->nbSize(2); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }

          numleft--;
          numchange++;
        }
      }
      catch (out_of_range &oor) {
        EOBException ex("Lattice","dissolvePhase","site_",site_.size(),i);
        ex.printException();
        exit(1);
      }
    }
    
    isite = interface_[phaseid].getDissolutionSites();

    return(numchange);
}


void Lattice::addDissolutionSite (Site *ste,
                                  unsigned int pid)
{
    try {
        interface_.at(pid).addDissolutionSite(ste);
        vector<unsigned int> plist = ste->getGrowthPhases();
        for (register unsigned int i = 0; i < plist.size(); i++) {
            interface_.at(plist[i]).removeGrowthSite(ste);
        }
        ste->setDissolutionSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","addDissolutionSite","interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
    return;
}

void Lattice::addGrowthSite (Site *ste,
                             unsigned int pid)
{
    try {
        interface_.at(pid).addGrowthSite(ste);
        ste->setGrowthSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","addGrowthSite","interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

void Lattice::removeDissolutionSite (Site *ste,
                                     unsigned int pid)
{
    try {
        interface_.at(pid).removeDissolutionSite(ste,false);
        ste->removeDissolutionSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","removeDissolutionSite",
                        "interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

void Lattice::removeGrowthSite (Site *ste,
                                unsigned int pid)
{
    try {
        interface_.at(pid).removeGrowthSite(ste);
        ste->removeGrowthSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","removeGrowthSite",
                        "interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

long int Lattice::emptyPorosity (long int numsites)
{
    register unsigned long int i,j;
    long int numemptied = 0;
    unsigned int cntpore,cntmax;
    bool placed;
    vector<Isite> emptysite;
    Isite isiteobj(0,0);

    if (numsites == 0) return (0);
    if (numsites < 0) {
        cout << "Going into Lattice::fillPorosity(" << -numsites << ")... " << endl;
        cout.flush();
        numemptied = fillPorosity(-numsites);
        return(-numemptied);
    }
    emptysite.clear();
    /*
    emptysite.resize(numsites,isiteobj);
    */
    vector<Isite>::iterator p = emptysite.begin();
    vector<Isite>::iterator start,end;
    cntmax = 0;

    
    ///
    /// Finding all potential VOID sites.
    ///
    /// @todo Consider removing some of the standard output, or setting a flag for it.
    ///

    cout << "Finding all potential void sites ... ";
    for (i = 0; i < site_.size(); i++) {
      if (site_[i].getPhaseId() == WATERID) {
        cntpore = countBox(7,i);
        if (cntpore > cntmax) cntmax = cntpore;
        isiteobj.setId(i);
        isiteobj.setAffinity(-cntpore);
        emptysite.push_back(isiteobj);
      }
    }
    cout << "OK, found " << emptysite.size()
         << " potential void sites." << endl;
    
    ///
    /// Sort the emptysite vector by the affinity as defined here.
    ///

    cout << "Now sorting the void sites... ";
    start = emptysite.begin();
    end = emptysite.end();
    sort(start,end,affinitySort);
    cout << "OK!" << endl;
    
    ///
    /// We want to empty the sites with the largest pore count
    ///

    if (emptysite.size() < numsites) numsites = emptysite.size();
    for (i = 0; i < numsites; i++) {
        setPhaseId(emptysite[i].getId(),VOIDID);
	count_[VOIDID] += 1.0;
        numemptied++;
    }

    emptysite.clear();
    
    return(numemptied);
}

long int Lattice::fillPorosity (long int numsites)
{
    register unsigned long int i,j;
    long int numfilled = 0;
    unsigned int cntpore,cntmin;
    bool placed;
    vector<Isite> fillsite;
    Isite isiteobj(0,0);

    cout << "In fillPorosity (" << numsites << ")" << endl;
    cout.flush();
    if (numsites == 0) return (0);
    if (numsites < 0) {
        numfilled = emptyPorosity(-numsites);
        return(numfilled);
    }
    fillsite.clear();
    /*
    fillsite.resize(numsites,isiteobj);
    */
    vector<Isite>::iterator p = fillsite.begin();
    vector<Isite>::iterator start,end;
    cntmin = 10000000;

    ///
    /// Finding all potential WATERID sites
    ///

    cout << "Finding all potential growth sites ... ";
    for (i = 0; i < site_.size(); i++) {
        if (site_[i].getPhaseId() == VOIDID) {
            cntpore = countBox(7,i);
            if (cntpore < cntmin) cntmin = cntpore;
            isiteobj.setId(i);
            isiteobj.setAffinity(cntpore);
            fillsite.push_back(isiteobj);
        }
    }

    cout << "OK, found " << fillsite.size()
         << " potential water sites." << endl;
    cout << "Now sorting the water sites... ";

    ///
    /// Sort the fillsite vector by the affinity as defined here
    ///

    start = fillsite.begin();
    end = fillsite.end();
    sort(start,end,affinitySort);
    cout << "OK!" << endl;

    ///
    /// We want to fill the sites with the smallest pore count
    ///

    if (fillsite.size() < numsites) numsites = fillsite.size();
    cout << "Setting the water sites in order from smallest to largest... ";
    cout.flush();
    for (i = 0; i < numsites; i++) {
        cout << "Setting site " << i << " with id "
             << fillsite[i].getId() << " and pore count "
             << -fillsite[i].getAffinity() << "... ";
        cout.flush();
        setPhaseId(fillsite[i].getId(),WATERID);
	count_[WATERID] += 1.0;
        numfilled++;
        cout << "Done.  Numfilled = " << numfilled << endl;
        cout.flush();
    }
    cout << "Done." << endl;
    cout.flush();

    fillsite.clear();
    return(numfilled);
}

int Lattice::countBox (int boxsize,
                       unsigned long int siteid)
{
    string msg;
    int boxhalf = boxsize / 2;
    int nfound = 0;
    register int ix,iy,iz,hx,hy,hz;

    try {
        int qxlo = site_[siteid].getX() - boxhalf;
        int qxhi = site_[siteid].getX() + boxhalf;
        int qylo = site_[siteid].getY() - boxhalf;
        int qyhi = site_[siteid].getY() + boxhalf;
        int qzlo = site_[siteid].getZ() - boxhalf;
        int qzhi = site_[siteid].getZ() + boxhalf;
    
        ///
        /// Make sure that the box bounds are all sensible given the
        /// boundary conditions
        ///

        qxlo += checkBC(qxlo,xdim_);
        qxhi += checkBC(qxhi,xdim_);
        qylo += checkBC(qylo,ydim_);
        qyhi += checkBC(qyhi,ydim_);
        qzlo += checkBC(qzlo,zdim_);
        qzhi += checkBC(qzhi,zdim_);

        ///
        /// Count the number of water or void sites in the box
        ///

        for (iz = qzlo; iz <= qzhi; iz++) {
            hz = iz + checkBC(iz,zdim_);
            for (iy = qylo; iy <= qyhi; iy++) {
                hy = iy + checkBC(iy,ydim_);
                for (ix = qxlo; ix <= qxhi; ix++) {
                    hx = ix + checkBC(ix,xdim_);
                    if (site_.at(getIndex(hx,hy,hz)).getPhaseId() == WATERID
                            || site_.at(getIndex(hx,hy,hz)).getPhaseId() == VOIDID) {
                        nfound++;
                    }
                }
            }
        }
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","countBox","site_",site_.size(),getIndex(hx,hy,hz));
        ex.printException();
        exit(1);
    }
    return(nfound);
}


void Lattice::setResolution (const double res)
{
    ///
    /// Make sure that resolution is a valid value
    ///

    try {
        string msg;
        if (res <= 0.001) {
            cout << endl;
            msg = "Lattice resolution <= 0.001";
            throw DataException("Lattice","setResolution",msg);
        }

        cout << "Changing lattice resolution from ";
        cout << resolution_ << " to " << res << endl;
        cout.flush();
        resolution_ = res;
    }
    catch (DataException dex) {
        dex.printException();
        exit(1);
    }
    return;
}

vector<unsigned long int> Lattice::getNeighborhood (const unsigned long int sitenum,
                                                    const int size)
{
    int xp,yp,zp;
    double dist;

    vector<unsigned long int> nh;
    nh.clear();

    unsigned int xc = site_[sitenum].getX();
    unsigned int yc = site_[sitenum].getY();
    unsigned int zc = site_[sitenum].getZ();

    if (size == 0) {
        nh.push_back(getIndex(xc,yc,zc));
        return nh;
    }
    
    for (register int k = -size; k <= size; k++) {
        zp = zc + k;
        for (register int j = -size; j <= size; j++) {
            yp = yc + j;
            for (register int i = -size; i <= size; i++) {
                xp = xc + i;
                dist = (double)((xc - xp) * (xc - xp));
                dist += (double)((yc - yp) * (yc - yp));
                dist += (double)((zc - zp) * (zc - zp));
                /*
                if ((sqrt(dist) - 0.5) <= (double)size) {
                    nh.push_back(getIndex(xp,yp,zp));
                }
                */
                nh.push_back(getIndex(xp,yp,zp));
            }
        } 
    }

    return nh;
}

unsigned long int Lattice::getIndex (int ix,
                                     int iy,
                                     int iz) const
{
   if (ix < 0) {
       if (BC != 1) {
           ix += (int)xdim_;
        } else {
           ix = 0;
        }
   } else if (ix >= (int)xdim_) {
       if (BC != 1) {
           ix -= (int)xdim_;
        } else {
           ix = (int)xdim_ - 1;
        }
   }

   if (iy < 0) {
       if (BC != 2) {
           iy += (int)ydim_;
        } else {
           iy = 0;
        }
   } else if (iy >= (int)ydim_) {
       if (BC != 2) {
           iy -= (int)ydim_;
        } else {
           iy = (int)ydim_ - 1;
        }
   }

   if (iz < 0) {
       if (BC != 3) {
           iz += (int)zdim_;
        } else {
           iz = 0;
        }
   } else if (iz >= (int)zdim_) {
       if (BC != 3) {
           iz -= (int)zdim_;
        } else {
           iz = (int)zdim_ - 1;
        }
   }

   return (unsigned long int)(ix + (xdim_ * iy)
           + ((xdim_ * ydim_) * iz)); 
}

void Lattice::changeMicrostructure (double time,
                                    bool isfirst)
{
    register unsigned int i,ii;
    long int numadded, numadded_actual;
    unsigned int tpid;
    unsigned long int cursites,newsites;
    long int tnetsites;
    double td,tvol,tmass;
    vector<double> vol,vfrac_next;
    vector<long int> netsites;
    vector<unsigned int> pid;
    vector<string> phasenames;

    ///
    /// @todo This function is very large; consider breaking it into small pieces
    /// for ease of maintenance and readability.
    ///
    /// Assign time of this state to the lattice time_ variable
    ///

    time_ = time;

    ///
    /// The next block, if uncommented, reads prior expansion
    /// strain data from a file and loads it into the appropriate
    /// class members.
    ///

    /*
    expansion_.clear();
    expansion_coordin_.clear();
    string fexpansion = jobroot_ + "_exp.dat";
    string fexpansioncoor = jobroot_ + "_exp_coor.dat";
    ifstream in(fexpansion.c_str());
    if (!in) {
      cout << "can't open file " << fexpansion << ", so exit program." << endl;
      exit(1);
    } else {
      while (!in.eof()) {
        int index;
        vector<double> expval;
        expval.clear();
        expval.resize(3,0.0);
        in >> index;
        //cout << "index = " << index << endl;
        in >> expval[0];
        in >> expval[1];
        in >> expval[2];
        //cout << "expval[0]: " << expval[0] << " expval[1]: " << expval[1] 
             //<< " expval[2]: " << expval[2] << endl;
        expansion_.insert(make_pair(index,expval));
        site_[index].setExpansionStrain(expval[0]);
      }
      in.close();
    }
    cout << "open fexpansioncoor file" << endl; 
    ifstream in1(fexpansioncoor.c_str());
    if (!in1) {
      cout << "can't open file " << fexpansioncoor << ", so exit program." << endl;
      exit(1);
    } else {
      while (!in1.eof()) {
        int index;
        vector<int> coordin;
        coordin.clear();
        coordin.resize(3,0);
        in1 >> index;
        in1 >> coordin[0];
        in1 >> coordin[1];
        in1 >> coordin[2];
        expansion_coordin_.insert(make_pair(index,coordin));
      }
      in1.close();
    }
    */

    cout << "expansion_.size() = " << expansion_.size() << " at the beginning of " 
         << time << endl;

    ///
    /// Zero out the amount of water to add to the microstructure.
    ///
    /// @todo Find out why; this class variable is not used again.
    ///

    waterchange_ = 0.0;

    ///
    /// Next, determine how many sites of each phase need to be assigned
    /// The ChemicalSystem object calculates the new volume fractions, so we
    /// just load them up here from the ChemicalSystem object.
    ///

    vfrac_next = chemsys_->getMicphasevolfrac();
    phasenames = chemsys_->getMicphasename();

    cout << "In Lattice::changeMicrostructure... new volume fractions are..." << endl;
    try {
        for (i = 0; i < vfrac_next.size(); i++) {
            cout << phasenames.at(i) << " = " << vfrac_next.at(i) << endl;
        }
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure","phasenames",
                        phasenames.size(),i);
        ex.printException();
        exit(1);
    }
    vol.clear();
    tvol = 0.0;

    ///
    /// Calculate number of sites of each phase in next state
    ///

    cout << "Calculating volume of each phase to be added..." << endl;
    try {
      for (i = 0; i < vfrac_next.size(); i++) {
        if (vfrac_next.at(i) > 0.0) {
          cout << "****Volume fraction[" << phasenames.at(i)
               << "] in next state should be = " << vfrac_next.at(i);
          cout << ", or "
               << (long int)((double)(numsites_ * vfrac_next.at(i)))
               << " sites" << endl;
        }
      }
    }
    catch (out_of_range &oor) {
      EOBException ex("Lattice","changeMicrostructure","phasenames",
                    phasenames.size(),i);
      ex.printException();
      exit(1);
    }

    ///
    /// The next block is executed only if there will eventually be some
    /// sulfate attack during this simulation.
    ///

    if (sattack_time_ != 1.0e10) {
        
        ///
        ///  Normalize to get volume fractions and compute number
        ///  of sites of each phase needed
        ///
        ///  @todo Find out why we need to do all of this just because
        ///  there will eventually be sulfate attack.  Why not wait until
        ///  sulfate attack simulation actually starts?
        ///

        netsites.clear();
        netsites.resize(chemsys_->getMicphasenum(),0);
        pid.clear();
        pid.resize(chemsys_->getMicphasenum(),0);

        try {
            int CSHID = 0, ETTRID = 0;
            int MONOID = 0, AFMCID = 0, HYDROTID = 0;
            CSHID = chemsys_->getMicid("CSH");
            ETTRID = chemsys_->getMicid("ETTR");
            MONOID = chemsys_->getMicid("MONOSULPH");
            AFMCID = chemsys_->getMicid("AFMC");
            HYDROTID = chemsys_->getMicid("HYDROTALC");
            for (i = 0; i < vfrac_next.size(); i++) {
                cursites = (unsigned long int)count_.at(i);
                newsites = (unsigned long int)((numsites_
                    * vfrac_next.at(i)) + 0.5);
                tnetsites = (i != WATERID) ? (newsites - cursites) : 0;
                netsites.at(i) = tnetsites;
                pid.at(i) = i;
                if (i == ETTRID && isfirst) {
                  netsites.at(i) = 0;
                  count_.at(i) = newsites;
                }        
                if (netsites.at(i) != 0) cout << "****netsites["
                                              << phasenames.at(i)
                                              << "] in this state = "
                                              << netsites.at(i)
                                              << endl;
              }
            
              ///
              /// Next block gets executed only if we are now simulating
              /// sulfate attack.
              ///

              if (time_ >= sattack_time_) {
                  cout << "start to transform aluminum phases into ETTR at time_ = " 
                       << time_ << endl;

                  ///
                  /// The relevant stress-free molar volume ratios for sulfate
                  /// attack phase transformations.
                  ///

                  double ETTRMONO = 2.288;
                  double ETTRAFMC = 2.699;
                  double ETTRHYDROT = 3.211;
                  vector<int> numchanged;
                  numchanged.clear();
                  numchanged.resize(2,0);
                  if ((netsites.at(MONOID) < 0) && (netsites.at(ETTRID) > 0)) {
                      numchanged = transform(MONOID,
                                             netsites.at(MONOID),
                                             ETTRID,
                                             netsites.at(ETTRID),
                                             ETTRMONO);

                      cout << "numchanged[0] is " << numchanged[0] << endl;
                      cout << "numchanged[1] is " << numchanged[1] << endl;
                      netsites.at(MONOID) += numchanged[0];
                      netsites.at(ETTRID) -= numchanged[1];
                      cout << "netsites.at(" << MONOID << ") is: "
                           << netsites.at(MONOID) << endl;
                      cout << "netsites.at(" << ETTRID << ") is: "
                           << netsites.at(ETTRID) << endl;
                  }

                  if ((netsites.at(AFMCID) < 0) && (netsites.at(ETTRID) > 0)) {
                      numchanged = transform(AFMCID,
                                             netsites.at(AFMCID),
                                             ETTRID,
                                             netsites.at(ETTRID),
                                             ETTRAFMC);
                      netsites.at(AFMCID) += numchanged[0];
                      netsites.at(ETTRID) -= numchanged[1];
                  }

                  if ((netsites.at(HYDROTID) < 0) && (netsites.at(ETTRID) > 0)) {
                      numchanged = transform(HYDROTID,
                                             netsites.at(HYDROTID),
                                             ETTRID,
                                             netsites.at(ETTRID),
                                             ETTRHYDROT);
                      netsites.at(HYDROTID) += numchanged[0];
                      netsites.at(ETTRID) -= numchanged[1];
                  }
              }
        }

        catch (out_of_range &oor) {
            EOBException ex("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                    phasenames.size(),i);
            ex.printException();
            exit(1);
        }

    } else {
        
        ///
        /// Sulfate attack will NEVER be done during this simulation.
        ///  Normalize to get volume fractions and compute number
        ///  of sites of each phase needed.
        ///

        netsites.clear();
        netsites.resize(chemsys_->getMicphasenum(),0);
        pid.clear();
        pid.resize(chemsys_->getMicphasenum(),0);

        try {
            for (i = 0; i < vfrac_next.size(); i++) {
                cursites = (long int)count_.at(i);
                newsites = (unsigned long int)((numsites_
                             * vfrac_next.at(i)) + 0.5);
                tnetsites = (i != WATERID) ? (newsites - cursites) : 0;
                netsites.at(i) = tnetsites;
                pid.at(i) = i;
                if (netsites.at(i) != 0) cout << "***netsites["
                                              << phasenames.at(i)
                                              << "] in this state = "
                                              << netsites.at(i)
                                              << endl;
            }
        }
        catch (out_of_range &oor) {
            EOBException ex("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                     phasenames.size(),i);
            ex.printException();
            exit(1);
        }
    }
    
    cout << "Sorting non-pore sites..." << endl;

    ///
    /// Sort netsites in ascending order, except we will handle
    /// void space differently.
    ///

    try {
        for (i = WATERID; i < netsites.size(); i++) {
            for (ii = i; ii < netsites.size(); ii++) {
                if (netsites[ii] < netsites[i]) {
                    tnetsites = netsites[i];
                    netsites[i] = netsites[ii];
                    netsites[ii] = tnetsites;
                    tpid = pid.at(i);
                    pid.at(i) = pid.at(ii);
                    pid.at(ii) = tpid;
                }
            }
        }
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure","pid",
                     pid.size(),ii);
        ex.printException();
        exit(1);
    }

    Interface ifc;
    vector<Isite> gs,ds;
    cout << "Getting change vectors for non-pore phases..." << endl;

    ///
    /// The next loop starts at 2 because we exclude void and water phases
    ///
    /// @todo Consider making the starting index more general
    ///

    for (i = 2; i < interface_.size(); i++) {
        ifc = interface_[i];
        gs = ifc.getGrowthSites();
        ds = ifc.getDissolutionSites();
    }

    cout << "Switching phase id values..." << endl;
    try {
        for (i = WATERID; i < netsites.size(); i++) {
            numadded = 0;
            numadded_actual = 0;
	  
	        if (netsites[i] < 0) {
                cout << "Going into dissolve_phase now... pid = "
                     << pid.at(i) << endl;
                cout.flush();
                numadded = dissolvePhase(pid.at(i),-netsites[i]);
                numadded_actual += numadded;
                cout << "...Done with dissolve_phase for phase "
                     << pid.at(i) << "!" << endl;
                cout.flush();
            } else if (netsites[i] > 0) {
                cout << "Going into grow_phase now... pid = "
                     << pid.at(i) << endl;
            
                cout.flush();
 
                numadded = growPhase(pid.at(i),netsites[i]);

                numadded_actual += numadded;
                int diff = netsites[i] - numadded;
                while (diff > 0) {
                    gs = interface_[pid.at(i)].getGrowthSites();
                    ds = interface_[pid.at(i)].getDissolutionSites();
                    cout << "gs.size() = " << gs.size() << endl;
                    cout << "ds.size() = " << ds.size() << endl;
                    if (gs.size() == 0) {
                        cout << "Phase " << pid.at(i) << " needs to be nucleated. " << endl;
                        int nuclei = diff;
	                    double thresh = (nuclei / (volumefraction_[WATERID] * site_.size()));
	                    cout << "Thresh = " << thresh << endl;
	                    cout.flush();
                        if (thresh < 1) {
    	                    double g = 0.0;
    	                    for (register int k = 0; k < site_.size(); k++) {
    	                        if (site_[k].getPhaseId() == WATERID) {
    	                            g = rg_->Ran3();
    	                            if (g < thresh) addGrowthSite(&site_[k],pid.at(i));
    	                        }
    	                    }
                        } else {
                            cout << "There is no room to grow, so exit the program." << endl;
                            exit(1);
                        }
	                }
                    numadded = growPhase(pid.at(i),diff);
    	            numadded_actual += numadded;
                    diff = diff - numadded;
	                cout << "diff = " << diff << endl;
	            }
    	        cout << "...Done with grow_phase for phase "
                     << pid.at(i) << "!" << endl;
                cout.flush();
            }
        
            if (numadded_actual*numadded_actual != netsites[i]*netsites[i]) {
                cout << "WARNING: Needed to switch on "
                     << netsites[i] << " of phase " << pid.at(i) << endl;
                cout << "         But actually did "
                     << numadded << " switches" << endl;
            }
        }
    }

    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure","pid",
                    pid.size(),ii);
        ex.printException();
        exit(1);
    }
    cout << "time_ in the Lattice.nw is: " << time_ << endl;

    cout << "Now emptying the necessary pore volume:  "
         << netsites[VOIDID] << " sites affected." << endl;
    long int numempty = emptyPorosity(netsites[VOIDID]);
    cout << "Number actually emptied was:  " << numempty << endl;

    ///
    /// Report on target and actual mass fractions
    ///

    cout << "*******************************" << endl;
    try {
        for (i = 0; i < vfrac_next.size(); i++) {
            volumefraction_.at(i) = ((double)(count_.at(i)))/((double)(site_.size()));
            cout << "Phase " << i << " Target volume fraction was "
                 << vfrac_next[i] << " and actual is "
                 << volumefraction_.at(i) << endl;
        }
    }

    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure",
                    "volumefraction_ or count_",volumefraction_.size(),i);
        ex.printException();
        exit(1);
    }

    cout << "*******************************" << endl;
    long int numfill = fillPorosity(1000000);

    /// 
    /// Report on target and actual mass fractions
    ///

    cout << "*******************************" << endl;
    try {
        for (i = 0; i < vfrac_next.size(); i++) {
            volumefraction_.at(i) = ((double)(count_.at(i)))/((double)(site_.size()));
            cout << "Phase " << i << " Target volume fraction was "
                 << vfrac_next[i] << " and actual is "
                 << volumefraction_.at(i) << endl;
        }
    }

    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure",
                    "volumefraction_ or count_",volumefraction_.size(),i);
        ex.printException();
        exit(1);
    }
    cout << "*******************************" << endl;
    
    for (i = 0; i < interface_.size(); i++) {
      if (i != VOIDID && i != WATERID) {
        interface_[i].sortGrowthSites(site_,i);
        interface_[i].sortDissolutionSites(site_,i);
      }
    }

    /*
    //update expansion and expansion_coordin files
    cout << "in Lattice, expansion_.size() = " << expansion_.size() << endl;
    ofstream out(fexpansion.c_str());
    for (map<int, vector<double> >::iterator it = expansion_.begin();
      it != expansion_.end(); it++) {
      int index = it -> first;
      vector<double> expval = it -> second;
      out << index << " " << expval[0] << " " << expval[1] << " " << expval[2] 
          << endl;
    }
    out.close();

    ofstream out1(fexpansioncoor.c_str());    
    for (map<int, vector<int> >::iterator it = expansion_coordin_.begin();
      it != expansion_coordin_.end(); it++) {
      int index = it -> first;
      vector<int> coor = it -> second;
      out1 << index << " " << coor[0] << " " << coor[1] << " " << coor[2]
           << endl;
    }
    out1.close();
    */

    ///
    ///  This is a local variable and the value is never used.
    ///
    ///  @todo Why not eliminate this line completely?
    ///

    double surfa = getSurfaceArea(chemsys_->getMicid("CSH"));
    
    return;
}

void Lattice::writeLattice (const string &root)
{
    register unsigned int i,j,k;
    string ofname(root);
    ostringstream ostr1,ostr2;
    ostr1 << (int)(time_ * 100.0);	// ten-thousandths of a second
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofname = ofname + "." + timestr + "." + tempstr + ".img";

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofname,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    // Write image header information first

    out << VERSIONSTRING << " " << version_ << endl;
    out << XSIZESTRING << " " << xdim_ << endl;
    out << YSIZESTRING << " " << ydim_ << endl;
    out << ZSIZESTRING << " " << zdim_ << endl;
    out << IMGRESSTRING << " " << resolution_ << endl;
 
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
               int index = getIndex(i,j,k);
               out << site_[index].getPhaseId() << endl;
            }
        }
    }

    out.close();
    
    ofname = root;
    ofname = ofname + "." + timestr + "." + tempstr + ".img.damage";

    ofstream out1(ofname.c_str());
    try {
        if (!out1.is_open()) {
            throw FileException("Lattice","writeLattice",ofname,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    // Write image header information first

    out1 << VERSIONSTRING << " " << version_ << endl;
    out1 << XSIZESTRING << " " << xdim_ << endl;
    out1 << YSIZESTRING << " " << ydim_ << endl;
    out1 << ZSIZESTRING << " " << zdim_ << endl;
    out1 << IMGRESSTRING << " " << resolution_ << endl;
 
    int DAMAGEID = chemsys_->getMicid("DAMAGE"); 
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
               int index = getIndex(i,j,k);
               if (site_[index].IsDamage()) {
                 out1 << DAMAGEID << endl;
               } else {
                 out1 << site_[index].getPhaseId() << endl;
               }
            }
        }
    }

    out1.close();
}

void Lattice::writeDamageLattice (const string &root)
{
    register unsigned int i,j,k;
    string ofname(root);
    ostringstream ostr1,ostr2;
    ostr1 << (int)(time_ * 100.0);	// ten-thousandths of a second
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofname = ofname + "." + timestr + "." + tempstr + ".img";

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofname,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    // Write image header information first

    out << VERSIONSTRING << " " << version_ << endl;
    out << XSIZESTRING << " " << xdim_ << endl;
    out << YSIZESTRING << " " << ydim_ << endl;
    out << ZSIZESTRING << " " << zdim_ << endl;
    out << IMGRESSTRING << " " << resolution_ << endl;
  
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
               int index = getIndex(i,j,k);
               if (site_[index].IsDamage()) {
                   out << "1" << endl;
               } else {
                   out << "0" << endl;
               }
            }
        }
    }

    out.close();
}

void Lattice::writeLatticePNG (const string &root)
{
    register unsigned int i,j,k;
    string ofname(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(zdim_,dumvec);
    image.resize(zdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    ostringstream ostr1,ostr2;
    ostr1 << (int)(time_ * 100.0);	// hundredths of an hour
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofname = ofname + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofname,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    ///
    /// Write PPM header for full color image
    ///

    out << "P3" << endl;
    out << ydim_ << " " << zdim_ << endl;
    out << COLORSATVAL << endl;

    unsigned int slice = xdim_/2;
    unsigned int nd,ixx,valout;
    long unsigned int sitenum;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           done = false;
           nd = 0;
           sitenum = getIndex(slice,j,k);
           ixx = slice;
           do {
               sitenum = getIndex(ixx,j,k);
               if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                   done = true;
               } else {
                   nd++;
                   ixx++;
                   if (ixx >= xdim_) ixx -= xdim_;
               }
           } while (!done);
           sitenum = getIndex(ixx,j,k);
           image[j][k] = site_[sitenum].getPhaseId();
           dshade[j][k] = 1.0;
           /*
           dshade[j][k] = 0.1 * (10.0 - ((double)nd));
           */
        }
    }          

    double red,green,blue;
    vector<double> colors;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           colors = chemsys_->getColor(image[j][k]);
           red = dshade[j][k]*colors[0] + 0.5;
           green = dshade[j][k]*colors[1] + 0.5;
           blue = dshade[j][k]*colors[2] + 0.5;
           out << (int)(red);
           out << " " << (int)(green);
           out << " " << (int)(blue) << endl;
        }
    }

    out.close();

    ///
    /// Execute system call to convert PPM to PNG.
    ///
    /// @warning This relies on installation of ImageMagick
    ///

    buff = "convert " + ofname + " " + ofpngname;
    system(buff.c_str());
    return;
}

void Lattice::writeDamageLatticePNG (const string &root)
{
    register unsigned int i,j,k;
    string ofname(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(zdim_,dumvec);
    image.resize(zdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    ostringstream ostr1,ostr2;
    ostr1 << (int)(time_ * 100.0);	// hundredths of an hour
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofname = ofname + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofname,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    ///
    /// Write PPM header for full color image
    ///

    out << "P3" << endl;
    out << ydim_ << " " << zdim_ << endl;
    out << COLORSATVAL << endl;

    unsigned int slice = xdim_/2;
    unsigned int nd,ixx,valout;
    long unsigned int sitenum;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           done = false;
           nd = 0;
           sitenum = getIndex(slice,j,k);
           ixx = slice;
           do {
               sitenum = getIndex(ixx,j,k);
               if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                   done = true;
               } else {
                   nd++;
                   ixx++;
                   if (ixx >= xdim_) ixx -= xdim_;
               }
           } while (!done);
           sitenum = getIndex(ixx,j,k);
           if (site_[sitenum].IsDamage()) {
               image[j][k] = 1;
           } else {
               image[j][k] = 5;
           }
           dshade[j][k] = 1.0;
           /*
           dshade[j][k] = 0.1 * (10.0 - ((double)nd));
           */
        }
    }          

    double red,green,blue;
    vector<double> colors;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           colors = chemsys_->getColor(image[j][k]);
           red = dshade[j][k]*colors[0] + 0.5;
           green = dshade[j][k]*colors[1] + 0.5;
           blue = dshade[j][k]*colors[2] + 0.5;
           out << (int)(red);
           out << " " << (int)(green);
           out << " " << (int)(blue) << endl;
        }
    }

    out.close();

    ///
    /// Execute system call to convert PPM to PNG.
    ///
    /// @warning This relies on installation of ImageMagick
    ///

    buff = "convert " + ofname + " " + ofpngname;
    system(buff.c_str());
    return;
}

void Lattice::makeMovie (const string &root)
{
    register unsigned int i,j,k;
    string ofname(root);
    string ofbasename(root);
    string ofgifname(root);
    string ofgifbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(zdim_,dumvec);
    image.resize(zdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    string buff;
    ostringstream ostr1,ostr2,ostr3;
    ostr1 << (int)(time_ * 100.0);	// hundredths of an hour
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());

    ///
    /// Loop over number of slices in the x direction, making one image per slice
    /// and appending it to the end of the master file.
    ///

    for (i = 10; i < xdim_; i++) {

        ///
        /// Open the output file.
        ///

        ostr3.clear();
        ostr3 << (int)(i);	// x slice number
        string istr(ostr3.str());
        ofname = ofbasename + "." + timestr + "."
            + tempstr + "." + istr + ".ppm";
        ofgifname = ofgifbasename + "." + timestr
            + "." + tempstr + "." + istr + ".gif";

        ofstream out(ofname.c_str());
        if (!out.is_open()) {
            throw FileException("Lattice","makeMovie",ofname,"Could not open");
        }

        ///
        /// Write PPM header for full color image.
        ///

        out << "P3" << endl;
        out << ydim_ << " " << zdim_ << endl;
        out << COLORSATVAL << endl;

        unsigned int slice = i;
        unsigned int nd,ixx,valout;
        long unsigned int sitenum;
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
               done = false;
               nd = 0;
               sitenum = getIndex(slice,j,k);
               ixx = slice;
               do {
                   sitenum = getIndex(ixx,j,k);
                   if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                       done = true;
                   } else {
                       nd++;
                       ixx++;
                       if (ixx >= xdim_) ixx -= xdim_;
                   }
               } while (!done);
               sitenum = getIndex(ixx,j,k);
               image[j][k] = site_[sitenum].getPhaseId();
               dshade[j][k] = 0.1 * (10.0 - ((double)nd));
            }
        }          

        double red,green,blue;
        vector<double> colors;
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
               colors = chemsys_->getColor(image[j][k]);
               red = dshade[j][k]*colors[0] + 0.5;
               green = dshade[j][k]*colors[1] + 0.5;
               blue = dshade[j][k]*colors[2] + 0.5;
               out << (int)(red);
               out << " " << (int)(green);
               out << " " << (int)(blue) << endl;
            }
        }
        out.close();


        ///
        /// Execute system call to convert PPM to GIF.
        ///
        /// @warning This relies on installation of ImageMagick
        ///

        buff = "convert " + ofname + " " + ofgifname;
        system(buff.c_str());
    }

    ///
    /// Execute system call to convert GIF frames to animated GIF.
    ///
    /// @warning This relies on installation of gifsicle
    ///

    buff = "gifsicle --delay=10 "
        + ofgifbasename + "*.gif > "
        + ofgifbasename + ".movie.gif";
    system(buff.c_str());
}

vector<int> Lattice::transform (int alphaseid,
                                int netsitesAlphaseid,
                                int ettrid,
                                int netsitesEttrid,
                                double volumeratio)
{
    ///
    /// @todo Consider breaking this method into smaller pieces
    ///

    int expindex;

    vector<double> expval;
    expval.clear();
    expval.resize(3,0.0);

    vector<int> coordin;
    coordin.clear();
    coordin.resize(3,0);

    vector<Isite> diss;
    diss = interface_[alphaseid].getDissolutionSites();
    cout << "The number of sites of phase " << alphaseid << " to dissolve is: "
         << diss.size() << endl;

    Site *ste;

    ///
    /// Construct the unordered list of sites to dissolve based only
    /// on the phase id and whether or not the total number of sites to dissolve
    /// has been reached.
    ///
    /// @remark Is this task biased by site position?  Has the list of sites been randomized?
    ///

    if (diss.size() < (-netsitesAlphaseid)) {
        for (int ii = 0; ii < numsites_; ii++) {
            ste = &site_[ii];
            if (ste->getPhaseId() == alphaseid) {
                addDissolutionSite(ste,alphaseid);
            }
        }
    }

    diss = interface_[alphaseid].getDissolutionSites();
    cout << "New size of diss of phase " << alphaseid << " is: " << diss.size() << endl;

    double ettrform = 0.0;
    int numtransform = 0;
    int max = (int) volumeratio;
    
    int CSHID = chemsys_->getMicid("CSH");
    int ETTRID = chemsys_->getMicid("ETTR");
    int HYDROTID = chemsys_->getMicid("HYDROTALC");
    cout << "CSHID = " << CSHID << "  ETTRID = " << ETTRID
         << " HYDROTID = " << HYDROTID << endl;    

    for (int ii = (diss.size() - 1); ii > 0 && ettrform < netsitesEttrid  
             && numtransform < (-netsitesAlphaseid); ii--) {
        ste = &site_[diss[ii].getId()];

        expval.clear();
        
        vector<Site *> cshneighbor, waterneighbor;
        cshneighbor.clear();
        waterneighbor.clear();
        for (int j = 0; j < ste->nbSize(1); j++) {
            Site *stenb;
            stenb = ste->nb(j);
            if (stenb->getPhaseId() == WATERID) {
              waterneighbor.push_back(stenb);
            } else if (stenb->getPhaseId() == CSHID) {
              cshneighbor.push_back(stenb);
            }
        }
        cout << "Having " << waterneighbor.size() << " water pixels and "
             << cshneighbor.size() << " CSH pixels in the neighborhood." << endl;             

        if ((waterneighbor.size() + 1) <= max) { // count the site itself

            /// Expansion should occur
            ///
            /// 1. Take the subvolume centered on aluminate site that will dissolve
            ///

            string fname(jobroot_ + "_alsubvol.dat");
            vector<unsigned long int> alnb = writeSubVolume(fname, ste, 1);
            ste->setDamage();  
            int numWater, numCSH;
            numWater = numCSH = 0;
            for (int nb = 0; nb < alnb.size(); nb++) {
              Site *alstenb = &site_[alnb[nb]];
              if (alstenb->getPhaseId() == WATERID) {
                numWater++;
              } else if (alstenb->getPhaseId() == CSHID) {
                numCSH++;
                alstenb->setDamage();
              } else if (alstenb->getPhaseId() == ETTRID || alstenb->getPhaseId() == HYDROTID) {
                alstenb->setDamage();
              }
            }

            ///
            /// 2. Calculate the effective bulk modulus of this subvolume
            ///

            double subbulk = FEsolver_->getBulkModulus(fname);
            cout << "subbulk = " << subbulk << " GPa." << endl;   
            subbulk = subbulk * 1.0e3; // convert GPa to MPa
            double subsolidbulk = subbulk;
            subsolidbulk *= ((1 + (double)numWater / 27.0) / (1 - (double)numWater / 27.0));
            
            ///
            /// @remark This seems like a double conversion.  Hasn't the conversion been done?
            ///

            subsolidbulk = subsolidbulk * 1.0e3; // convert GPa to MPa         

            ///
            /// 3. Calculate crystallization strain in this sub volume;
            ///    porevolfrac is the volume fraction of pore space occupied by crystal
            ///

            double porevolfrac = 0.0;
            if (numWater != 0 || numCSH != 0) {
              porevolfrac = (double)(volumeratio) / (numWater + numCSH * 0.25);
            } else {
              porevolfrac = 1.0;
            } 
            cout << "ettrSI_ = " << ettrSI_ << endl;
            double exp = solut_->calculateCrystrain(ettrSI_,
                                                    porevolfrac,
                                                    subbulk,
                                                    subsolidbulk);

            ///
            /// 4. Apply expansion strain on each voxel in this sub volume
            ///

            applyExp(alnb, exp);

            setPhaseId(ste, ettrid);

            ///
            /// The Al-bearing phase has dissolved, so remove it from the
            /// list of dissolution sites of this phase
            ///

            removeDissolutionSite(ste, alphaseid);
            numtransform++;
            ettrform++;
            count_.at(ettrid)++;
            
            for (int i = 0; i < waterneighbor.size(); i++) {
                setPhaseId(waterneighbor[i], ettrid);

                ///
                /// Ettringite has grown here, so remove this site from the
                /// list of growth sites of this phase
                ///

                removeGrowthSite(waterneighbor[i], ettrid);
                count_.at(ettrid)++;
                ettrform++;

                ///
                /// Weighted mean curvature (wmc) is changed by the difference
                /// between the growing phase's porosity and the template's porosity.
                ///
                /// @todo Determine why the calculation works this way.
                ///

                double dwmcval = chemsys_->getPorosity(ettrid)
                                - chemsys_->getPorosity(WATERID);
                for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
                    Site *nb = waterneighbor[i]->nb(j);
                    nb->dWmc(dwmcval);
                }
            }

            dWaterchange(volumeratio - (waterneighbor.size() + 1));
            ettrform += (volumeratio - (waterneighbor.size() + 1));
            count_.at(ettrid) += (volumeratio - (waterneighbor.size() + 1));

        } else {

            ///
            /// Expansion should not occur because there is sufficient
            /// free space for local ettringite growth.
            ///

            setPhaseId(ste, ettrid);

            ///
            /// The Al-bearing phase has dissolved, so remove it from the
            /// list of dissolution sites of this phase
            ///

            removeDissolutionSite(ste, alphaseid);
            numtransform++;
            ettrform++;
            count_.at(ettrid)++;
            double thresh = 0.0, g = 0.0;
            thresh = volumeratio - max;
            g = rg_->Ran3();

            int upperindex = (g < thresh) ? max : max - 1;
            for (int i = 0; i < upperindex; i++) {
                setPhaseId(waterneighbor[i], ettrid);
                removeGrowthSite(waterneighbor[i], ettrid);
                ettrform++;
                count_.at(ettrid)++;

                ///
                /// Weighted mean curvature (wmc) is changed by the difference
                /// between the growing phase's porosity and the template's porosity.
                ///
                /// @todo Determine why the calculation works this way.
                ///

                double dwmcval = chemsys_->getPorosity(ettrid)
                                 - chemsys_->getPorosity(WATERID);
                for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
                    Site *nb = waterneighbor[i]->nb(j);
                    nb->dWmc(dwmcval);
                }
            }

        }
    }       // End of loop over all ettringite sites to form

    /*
    netsites.at(alphaseid) += (long int) numtransform;
    netsites.at(ettrid) -= (long int) ettrform;	
    */

    cout << "The number of aluminum phase " << alphaseid
         << " transformed into ETTR is: " << numtransform << endl;

    vector<int> numchanged;
    numchanged.clear();
    numchanged.resize(2,0);
    numchanged[0] = numtransform;
    numchanged[1] = (long int)ettrform;
    cout << "The number of ettrform is: " << ettrform << endl;

    return numchanged;
}

vector<unsigned long int> Lattice::writeSubVolume (string fname,
                                                   Site *centerste,
                                                   int size)
{
    ofstream out(fname.c_str());
    
    out << "Version: 7.0" << endl;
    out << "X_Size: 3" << endl;
    out << "Y_Size: 3" << endl;
    out << "Z_Size: 3" << endl;
    out << "Image_Resolution: 1" << endl;

    vector<unsigned long int> alnb = getNeighborhood(centerste->getId(),size);

    for (int j = 0; j < alnb.size(); j++) {
      int phaseid = site_[alnb[j]].getPhaseId();
      out << phaseid << endl;
    }
    out.close();

    return alnb;
}

void Lattice::applyExp (vector<unsigned long int> alnb,
                        double exp)
{
    Site *ste;
    if (exp > 0.0) {
      for (int i = 0; i < alnb.size(); i++) {
        ste = &site_[alnb[i]];
        ste->setExpansionStrain(exp);

        map<int,vector<double> >::iterator p = expansion_.find(ste->getId());
        if (p != expansion_.end()) {
          for (int j = 0; j < (p->second).size(); j++) {
            (p->second)[j] = ste->getExpansionStrain();
          }
        } else {
          vector<double> expval;
          expval.clear();
          expval.resize(3,ste->getExpansionStrain());
          expansion_.insert(make_pair(ste->getId(),expval));
        }
      
        map<int,vector<int> >::iterator pp = expansion_coordin_.find(ste->getId());
        if (pp == expansion_coordin_.end()) {
          vector<int> coordin;
          coordin.clear();
          coordin.resize(3,0);
          coordin[0] = ste->getX();
          coordin[1] = ste->getY();
          coordin[2] = ste->getZ();
          expansion_coordin_.insert(make_pair(ste->getId(),coordin));
        }
      }
    }
  
    return;
}

double Lattice::getSurfaceArea (int phaseid)
{
    double surface1 = 0.0, surface2 = 0.0;
    Site *ste, *stenb;
    vector<Isite> isite = interface_[phaseid].getDissolutionSites();

    ///
    /// Method 1: Surface area is equal to the wmc (prop to volume of surrounding pores)
    ///

    for (int i = 0; i < isite.size(); i++) {
        ste = &site_[isite[i].getId()];
        surface1 += ste->getWmc();
    }

    ///
    /// Method 2: Surface area is related to interior porosity of neighbor sites
    ///

    for (int i = 0; i < site_.size(); i++) {
        ste = &site_[i];
        if (ste->getPhaseId() == phaseid) {
            for (int j = 0; j < ste->nbSize(1); j++) {
                stenb = ste->nb(j);
                surface2 += chemsys_->getPorosity(stenb->getPhaseId());
            }
        }
    }    

    cout << "surface area of phase " << phaseid << " calculated by method 1 is: "
         << surface1 << endl;
    cout << "surface area of phase " << phaseid << " calculated by method 2 is: "
         << surface2 << endl;

    ///
    /// Use Method 2
    ///

    surfacearea_ = surface2;

    return surfacearea_;
}
