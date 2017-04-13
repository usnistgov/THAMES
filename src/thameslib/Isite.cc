/**
@file Isite.cc
@brief Method definitions for the Isite class.

@remark These could all be placed in the Isite.h file as inline functions.
@todo Place in the Isite.h file as inline functions, and delete this file.
*/
#include "Isite.h"

Isite::Isite ()
{
    affinity_ = 0;
    id_ = 0;
}

Isite::Isite (unsigned long int idval,
              int aftyval)
{
    id_ = idval;
    affinity_ = aftyval;
}

Isite::Isite (const Isite &obj)
{
    id_ = obj.id_;
    affinity_ = obj.affinity_;
}
