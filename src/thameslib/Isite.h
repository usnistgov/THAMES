/**
@file Isite.h
@brief Declare a class that keeps track of the chemical potential of sites

This class instantiates to objects that have

    - an id that defines a corresponding `Site` object, 
    - an <i>affinity</i> value, which is a qualitative chemical potential variable,
        similar in some ways to weighted mean curvature of an interface
*/

#ifndef ISITEH
#define ISITEH

#include "Site.h"
#include "Site.h"

using namespace std;

/**
@class Declaration of the Isite class.

*/
class Isite {

private:
    
    unsigned long int id_;      /**< The id of the corresponding Site */
    int affinity_;              /**< The affinity for growth of a phase at the site */

public:
    
/**
@brief Default constructor initializes members to zero.

@note NOT USED.
*/
Isite (); 

/**
@brief Overloaded constructor sets the members to prescribed values at construction time.

@param idval is the id of the corresponding Site object
@param aftyval is the prescribed value of the affinity to set
*/
Isite (unsigned long int idval,
       int aftyval);

/**
@brief Copy constructor.

@param The Isite object to copy
*/
Isite (const Isite &obj);
    
/**
@brief Get the id number of the corresponding Site object.

@return the id number of the corresponding Site object
*/
unsigned long int getId () const
{
    return id_;
}

/**
@brief Set the id number of the corresponding Site object.

@todo Maybe the argument should be declared const

@param idval is the id number of the corresponding Site object
*/
void setId (unsigned long int idval)
{
    id_ = idval;
}
    
/**
@brief Get the growth affinity of the corresponding Site object.

@return the growth affinity of the corresponding Site object
*/
int getAffinity () const
{
    return affinity_;
}

/**
@brief Set the growth affinity of the corresponding Site object.

@todo Maybe the argument should be declared const

@param num is the growth affinity of the corresponding Site object
*/
void setAffinity (int num)
{
    affinity_ = num;
}

};  // End of the Isite class
#endif
