/**
@file Interface.cc
@brief Definition of methods for the Interface class.

*/
#include "Interface.h"
#include "ChemicalSystem.h"
#include "RanGen.h"

bool cmp (const Site *s1,
          const Site *s2)
{
    return s1->getWmc() < s2->getWmc();
}

bool affinitySort (const Isite s1,
                   const Isite s2)
{
    return s1.getAffinity() > s2.getAffinity();
}

Interface::Interface ()
{
    phaseid_ = 0;
    growth_sites_.clear();
    dissolution_sites_.clear();
}

Interface::Interface (RanGen *rg)
{
    phaseid_ = 0;
    growth_sites_.clear();
    dissolution_sites_.clear();
    rg_ = rg;
}

Interface::Interface (ChemicalSystem *csys,
                      RanGen *rg,
                      vector<Site *> gv,
                      vector<Site *> dv,
                      unsigned int pid)
{
    register unsigned long int j;
    register unsigned int i;
    int afty;

    vector<Site *>::iterator beginLocation,endLocation;
    vector<Isite>::iterator start,end;

    rg_ = rg;
    phaseid_ = pid;
    chemsys_ = csys;

    dissolution_sites_.clear();
    growth_sites_.clear();

    ///
    /// Eliminate duplicate values from the growth site vector
    ///

    if (gv.size() > 0) {
        sort(gv.begin(),gv.end());
        i = 1;
        beginLocation = gv.begin();
        endLocation = unique(gv.begin(),gv.end());
        gv.erase(endLocation,gv.end());
    }


    ///
    /// Now do the same thing to the dissolution site vector
    ///

    if (dv.size() > 0) {
        sort(dv.begin(),dv.end());
        i = 1;
        beginLocation = dv.begin();
        endLocation = unique(dv.begin(),dv.end());
        dv.erase(endLocation,dv.end());
    }

    ///
    /// Now sort the growth sites according to the affinity
    ///

    for (j = 0; j < gv.size(); j++) {
        afty = 0;
        for (i = 0; i < gv[j]->nbSize(2); i++) {
            afty += chemsys_->getAffinity(pid,gv[j]->nb(i)->getPhaseId());
        }

        ///
        /// Add to the list of Isites.  An Isite is an object consisting
        /// of a pointer to a site and an affinity value
        ///

        growth_sites_.push_back(Isite(gv[j]->getId(),afty));
    }

    if (growth_sites_.size() > 0) {
        start = growth_sites_.begin();
        end = growth_sites_.end();

        ///
        /// The sort is built in to the STL for vectors or portions of vectors,
        /// as long as you give it a comparison function, which in this case is the
        /// `affinitySort` function already defined
        ///
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the sites are perfectly sorted in descending order of affinity.
    /// Now shuffle the sorted sites according to the random growth factor.
    /// Each phase has a certain amount of randomness to its growth, which can be
    /// accessed through the `getRandomgrowth` method of the ChemicalSystem object.
    ///

    unsigned long int site1,site2;
    unsigned long int numgsites = growth_sites_.size();

    for (j = 0; j < (chemsys_->getRandomgrowth(pid) * numgsites); j++) {

       // Choose two sites at random and switch their places
       site1 = (unsigned long int)(rg_->Ran3() * numgsites);
       site2 = (unsigned long int)(rg_->Ran3() * numgsites);
       swap(growth_sites_[site1],growth_sites_[site2]);
    }

    ///
    /// Now sort the dissolution sites according to the affinity
    ///

    try {
        for (j = 0; j < dv.size(); j++) {
            afty = 0;
            for (i = 0; i < dv[j]->nbSize(2); i++) {
                afty += chemsys_->getAffinity(pid,dv[j]->nb(i)->getPhaseId());
            }
            dissolution_sites_.push_back(Isite(dv[j]->getId(),afty));
        }
    }
    catch (EOBException e) {
        e.printException();
        exit(0);
    }

    if (dissolution_sites_.size() > 0) {
        start = dissolution_sites_.begin();
        end = dissolution_sites_.end();
        sort(start,end,affinitySort);
    }

    numgsites = dissolution_sites_.size();
    
    ///
    /// The dissolution sites are perfectly ordered by affinity, just like the
    /// growth sites were.  Now add the randomness factor for this phase.
    ///

    try {
        for (j = 0; j < (chemsys_->getRandomgrowth(pid) * numgsites); j++) {

           // Choose two sites at random and switch their places
           site1 = (unsigned long int)(rg_->Ran3() * numgsites);
           site2 = (unsigned long int)(rg_->Ran2() * numgsites);
           swap(dissolution_sites_[site1],dissolution_sites_[site2]);
        }
    }
    catch (EOBException e) {
        e.printException();
        exit(0);
    }
}       // End of constructors

Interface::~Interface ()
{
    growth_sites_.clear();
    dissolution_sites_.clear();
}

bool Interface::addGrowthSite (Site *loc)
{
    bool answer = false;
    bool found = false;
    register unsigned int i;
    vector<Isite>::iterator p,q,start,end;
    start = growth_sites_.begin();
    end = growth_sites_.end();

    ///
    /// See if site is already present
    ///

    for (i = 0; (i < growth_sites_.size()) && (!found); i++) {
        if (loc->getId() == growth_sites_[i].getId()) found = true;
    }

    ///
    /// Add the site only if it was not found already
    ///

    if (!found) {
        int afty = 0;
        for (i = 0; i < loc->nbSize(2); i++) {
            afty += chemsys_->getAffinity(phaseid_,loc->nb(i)->getPhaseId());
        }
        Isite tisite(loc->getId(),afty);
        q = lower_bound(start,end,tisite,affinitySort);
        growth_sites_.insert(q,tisite);
        answer = true;
    }

    return answer;
}

bool Interface::addDissolutionSite (Site *loc)
{
    bool answer = false;
    bool found = false;
    register unsigned int i;
    vector<Isite>::iterator p,q,start,end;
    start = dissolution_sites_.begin();
    end = dissolution_sites_.end();

    ///
    /// See if site is already present
    ///

    for (i = 0; (i < dissolution_sites_.size()) && (!found); i++) {
        if (loc->getId() == dissolution_sites_[i].getId()) found = true;
    }

    ///
    /// Add the site only if it was not found already
    ///

    if (!found) {
        int afty = 0;
        for (i = 0; i < loc->nbSize(2); i++) {
            afty += chemsys_->getAffinity(phaseid_,loc->nb(i)->getPhaseId());
        }
        Isite tisite(loc->getId(),afty);
        q = lower_bound(start,end,tisite,affinitySort);
        dissolution_sites_.insert(q,tisite);
        answer = true;
    }

    return answer;
}

bool Interface::sortGrowthSites (vector<Site> &ste,
                                 unsigned int pid)
{
    register unsigned int i,j;
    int afty;
    Site gs;

    ///
    /// The list of growth sites already exists or has been constructed, so we
    /// only need to update the affinities for each site
    ///

    for (j = 0; j < growth_sites_.size(); j++) {
        afty = 0;
        gs = ste[growth_sites_[j].getId()];
        for (i = 0; i < gs.nbSize(2); i++) {
            afty += chemsys_->getAffinity(pid,gs.nb(i)->getPhaseId());
        }
        growth_sites_[j].setAffinity(afty);
    }

    ///
    /// Now we sort the list of growth sites in descending order of affinity.
    /// Remember that `affinitySort` is the comparison function, already defined
    /// in this class, that must be passed to the STL sort function.
    ///

    if (growth_sites_.size() > 0) {
        vector<Isite>::iterator start,end;
        start = growth_sites_.begin();
        end = growth_sites_.end();
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the growth sites are perfectly sorted in descending
    /// order of affinity.  Next, we shuffle this sorting somewhat depending
    /// on how much randomization of growth sites is indicated for this particular
    /// phase.  The amount of randomness is obtained from the `getRandomgrowth` method
    /// of the ChemicalSystem object for this simulation
    ///

    unsigned long int site1,site2;
    unsigned long int numgsites = growth_sites_.size();
    for (j = 0; j < (chemsys_->getRandomgrowth(pid) * numgsites); j++) {

       ///
       /// Choose two sites at random and switch their places
       ///
 
       site1 = (unsigned long int)(rg_->Ran3() * numgsites);
       site2 = (unsigned long int)(rg_->Ran3() * numgsites);
       swap(growth_sites_[site1],growth_sites_[site2]);
    }

    return true;  // successful sorting
}

bool Interface::sortDissolutionSites (vector<Site> &ste,
                                      unsigned int pid)
{
    register unsigned int i,j;
    int afty;
    Site ds;

    ///
    /// The list of dissolution sites already exists or has been constructed, so we
    /// only need to update the affinities for each site
    ///

    for (j = 0; j < dissolution_sites_.size(); j++) {
       afty = 0;
       ds = ste[dissolution_sites_[j].getId()];
       for (i = 0; i < ds.nbSize(2); i++) {
           afty += chemsys_->getAffinity(pid,ds.nb(i)->getPhaseId());
       }
       dissolution_sites_[j].setAffinity(afty);
    }

    ///
    /// Now we sort the list of dissolution sites in descending order of affinity.
    /// Remember that `affinitySort` is the comparison function, already defined
    /// in this class, that must be passed to the STL sort function.
    ///

    if (dissolution_sites_.size() > 0) {
        vector<Isite>::iterator start,end;
        start = dissolution_sites_.begin();
        end = dissolution_sites_.end();
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the dissolution sites are perfectly sorted in descending
    /// order of affinity.  Next, we shuffle this sorting somewhat depending
    /// on how much randomization of growth sites is indicated for this particular
    /// phase.  The amount of randomness is obtained from the `getRandomgrowth` method
    /// of the ChemicalSystem object for this simulation
    ///

    unsigned long int site1,site2;
    unsigned long int numdsites = dissolution_sites_.size();
    for (j = 0; j < (chemsys_->getRandomgrowth(pid) * numdsites); j++) {

       // Choose two sites at random and switch their places
       site1 = (unsigned long int)(rg_->Ran3() * numdsites);
       site2 = (unsigned long int)(rg_->Ran3() * numdsites);
       swap(dissolution_sites_[site1],dissolution_sites_[site2]);
    }

    return true;
}

bool Interface::removeGrowthSite (Site *loc)
{
    bool found = false;
    vector<Isite>::iterator p = growth_sites_.begin();
    while ((p != growth_sites_.end()) && (!found)) {
        if (p->getId() == loc->getId()) {
            growth_sites_.erase(p);
            found = true;
        }
        p++;
    }

    return found;
}

bool Interface::removeDissolutionSite (Site *loc,
                                       bool verbose = false)
{
     if (verbose) {
        cout << "Removing dissolution site " << loc->getId()
             << ", size is " << dissolution_sites_.size()
             << endl;
        cout.flush();
        bool found = false;
        cout << "Trying to declare iterator to Isite vector... ";
        cout.flush();
        vector<Isite>::iterator p;
        cout << "Success!" << endl;
        cout.flush();
        cout << "Trying to set it to beginning of dissolution_sites_ ... ";
        cout.flush();
        p = dissolution_sites_.begin();
        cout << "Success!" << endl;
        cout.flush();
        for (register int i = dissolution_sites_.size() - 1; (i >= 0 && (!found)); i--) {
            if (dissolution_sites_[i].getId() == loc->getId()) {
                p += i;
                dissolution_sites_.erase(p);
                found = true;
            }
        }

        return found;
    } else {
        bool found = false;
        vector<Isite>::iterator p;
        p = dissolution_sites_.begin();
        for (register int i = dissolution_sites_.size() - 1; (i >= 0 && (!found)); i--) {
            if (dissolution_sites_[i].getId() == loc->getId()) {
                p += i;
                dissolution_sites_.erase(p);
                found = true;
            }
        }

        return found;
    }
}
