/**
@file Site.cc
@brief Method definitions for the Site class.

*/

#include "Site.h"
#include "ChemicalSystem.h"

Site::Site ()
{
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  strfreevolume_ = truevolume_ = 0.0;
  damage_ = false;
  expstrain_ = 0.0;
}

Site::Site (unsigned int xp,
            unsigned int yp,
            unsigned int zp,
            unsigned int xs,
            unsigned int ys,
            unsigned int zs,
            unsigned int neigh,
            ChemicalSystem *csys)
{    
    x_ = y_ = z_ = 0;
    id_ = 0;
    nb_.clear();
    strfreevolume_ = truevolume_ = 0.0;
    damage_ = false;
    x_ = xp;
    y_ = yp;
    z_ = zp;

    dissolution_.clear();
    growth_.clear();

    id_ = (unsigned long int)(x_ + (xs * y_) + ((xs * ys) * z_));

    if (neigh > 0) nb_.resize(neigh,0);

    chemsys_ = csys;

    strfreevolume_ = truevolume_ = 1.0;

    expstrain_ = 0.0;
}

void Site::calcWmc(void)
{
    wmc_ = chemsys_->getPorosity(getPhaseId());
    for (register unsigned int i = 0; i < nb_.size(); i++) {
        wmc_ += chemsys_->getPorosity(nb_[i]->getPhaseId());
    }
}
