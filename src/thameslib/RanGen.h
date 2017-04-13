/**
@file RanGen.h
@brief Declares and defines the random number generator object.

*/
#ifndef RANGENH
#define RANGENH

#include "ran3.h"
#include "ran2.h"
#include <vector>
#include <stdexcept>
#include "Exceptions.h"

/**
@class RanGen
@brief Declares and defines the portable random number generator object.
*/
class RanGen {

private:

static int* seed_;      /**< The random number seed, a negative integer */

public:

/**
@brief Set the random number seed.

@param s is the random number seed to set
*/
static void setSeed (int s)
{
    *seed_ = ((s < 0) ? s : -s);
}

/**
@brief Get the random number seed.

@param s is the random number seed to set
*/
static int getSeed ()
{
    return *seed_;
}

/**
@brief Interface to the ran3 portable random number generator (see ran3.h)

@return A uniform random number on [0,1]
*/
static double Ran3 ()
{
    return ran3(seed_);
}

/**
@brief Interface to the ran2 portable random number generator (see ran2.h)

@return A uniform random number on [0,1]
*/
static double Ran2 ()
{
    return ran2(seed_);
}

/**
@brief Randomize the elements of a vector.

@param v is a reference ot the vector
@param ntimes is the number of shuffling iterations to execute
*/
template<class X> static void shuffle(vector<X> &v, int ntimes=1)
{
    X tmp;
    int j;
    register int ii;
    register unsigned int i;
         
    for (ii = 0; ii < ntimes; ii++) {
        for (i = 0; i < v.size(); i++) {
            j = (int)(Ran3() * v.size());
            tmp = v[i];
            v[i] = v[j];
            v[j] = tmp;
        }
    }
}

/**
@brief Default constructor.

The default constructor attempts to seed itself using the C++ rand function.

@note NOT USED.
*/
RanGen () {
   try {
       seed_ = new int;
       *seed_ = rand();
   }
   catch (bad_alloc &ba) { throw ba; }
}

/**
@brief Overloaded constructor.

This constructor is given a seed to use to initialize itself.  It allocates
new memory to hold the integer seed.

@param s is the initial random number seed (negative integer)
*/
RanGen (int s)
{
   try {
       seed_ = new int;
       *seed_ = ((s < 0) ? s : -s);
   }
   catch (bad_alloc &ba) { throw ba; }
}

/**
@brief Destructor deletes the allocated memory for the seed.

*/
~RanGen ()
{
   delete seed_;
}

};      // End of the RanGen class
#endif
