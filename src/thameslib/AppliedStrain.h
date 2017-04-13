/**
@file AppliedStrain.h
@brief Declaration of the AppliedStrain derived class.

Solves the linear elastic state of the finite element mesh.
*/

#ifndef APPLIEDSTRAIN_H
#define APPLIEDSTRAIN_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "ElasticModel.h"

using namespace std;
 
/**
@class AppliedStrain
@brief Solves the linear elastic state of the finite element mesh.

@section intro_sec Introduction
This program solves the linear elastic equations in a                    
random linear elastic material, subject to an applied macroscopic strain,
using the finite element method.  Each pixel in the 3D digital          
image is a cubic tri-linear finite element,  having its own              
elastic moduli tensor. Periodic boundary conditions are maintained.      
In the comments below, notes are given for sections of code that  
the user might have to change for his particular problem. Therefore the  
user is encouraged to search for this string.                            

@subsection p_and_v Problem and Variable Definitions

The problem being solved is the minimization of the energy               
@f[
E = \frac{1}{2} u \cdot A \cdot u + b \cdot u + C
@f]
where <i>A</i> is the Hessian matrix composed of the         
stiffness matrices (dk) for each pixel/element, <i>b</i> is a constant vector   
and <i>C</i> is a constant that are determined by the applied strain and        
the periodic boundary conditions, and <i>u</i> is a vector of                   
all the displacements. The solution                                      
method used is the conjugate gradient relaxation algorithm.              
Other variables are:
    - `gb` is the gradient, <i>Au</i> + <i>b</i>
    - `h` is an  auxiliary variable used in the conjugate gradient algorithm (in dembx)
    - `dk(n,i,j)` is the stiffness matrix of the n'th phase
    - `cmod(n,i,j) is the elastic moduli tensor of the n'th phase
    - `pix` is a vector that gives  the phase label of each pixel
    - `ib` is a matrix that gives the labels of the 27 (counting itself)
       neighbors of a given node
    - `prob` is the volume   fractions of the various phases
    - `strxx`, `stryy`, `strzz`, `strxz`, `stryz`, and `strxy` are the six Voigt volume
        averaged total stresses
    - `sxx`, `syy`, `szz`, `sxz`, `syz`, and `sxy` are the six Voigt volume averaged total strains.

@subsection dims Dimensions
The vectors `u`, `gb`, `b`, and `h` are dimensioned to be the system size,         
@f$n_s = n_x n_y n_z@f$, with three components, where the digital image of the       
microstructure considered is a rectangular paralleliped with
dimenions <i>n<sub>x</sub></i>,<i>n<sub>y</sub></i>, <i>n<sub>z</sub></i>.
The arrays `pix` and `ib` are are also dimensioned to the system   
size.  The array `ib` has 27 components, for the 27 neighbors of a node.   
Note that the program is set up at present to have at most 50            
different phases.  This can easily be changed, simply by changing        
the dimensions of `dk`, `prob`, and `cmod`. The parameter `nphase` gives the     
number of phases being considered in the problem.                        

@warning Read the manual before using this program!!         
Manual available at NISTIR 6269 from NTIS or at:                         
http://ciks.cbt.nist.gov/~garbocz/manual/man.html              
*/

class AppliedStrain : public ElasticModel {

protected:

double exx_; /**< xx diagonal component of applied strain */
double eyy_; /**< yy diagonal component of applied strain */
double ezz_; /**< zz diagonal component of applied strain */
double exz_; /**< xz off-diagonal component of applied strain */
double eyz_; /**< yz off-diagonal component of applied strain */
double exy_; /**< xy off-diagonal component of applied strain */

public:

/**
@brief Default constructor.

The default constructor is the only constructor provided for this class.
It simply initializes the class variables to zero.

@param nx is the x dimension of the mesh
@param ny is the y dimension of the mesh
@param nz is the z dimension of the mesh
@param dim is the total number of elements
@param nphase is the number of phases
@param npoints is the number of microstructures (usually 1)
*/
AppliedStrain (int nx,
               int ny,
               int nz,
               int dim,
               int nphase,
               int npoints);
	
/**
@brief Destructor.

*/
~AppliedStrain ()
{
    exx_ = 0.0;
}

/**
@brief Set up the stiffness matrix for the problem.

@param nx is the x dimension of the mesh
@param ny is the y dimension of the mesh
@param nz is the z dimension of the mesh
@param ns is the total number of elements
@param nphase is the number of phases
*/
void femat (int nx,
            int ny,
            int nz,
            int ns,
            int nphase);

/**
@brief Conjugate gradient minimization of elastic energy.

@param ns is the total number of elements
@param gg is the energy gradient
@param ldemb is the maximum number of conjugate gradient steps to use during this call
@param kkk is the number of times the function has been called
@return number of conjugate gradient steps used in this call
*/
long int dembx (int ns,
                double gg,
                int ldemb,
                int kkk);

/**
@brief Computes the total elastic energy of the finite element mesh.

@param nx is the x dimension of the mesh
@param ny is the y dimension of the mesh
@param nz is the z dimension of the mesh
@param ns is the total number of elements
@return the energy of the finite element mesh
*/
double energy (int nx,
             int ny,
             int nz,
             int ns);
	
/**
@brief Computes the stress and strain components.

@param nx is the x dimension of the mesh
@param ny is the y dimension of the mesh
@param nz is the z dimension of the mesh
@param ns is the total number of elements
*/
void stress (int nx,
             int ny,
             int nz,
             int ns);
	
/**
@brief Controls the relaxation process, including calls to dembx, energy, and stress
functions.

@param kmax is the maximum number of times dembx will be called
*/
void relax (int kmax);
	
/**
@brief Main block, including the initialization of the stiffness matrix
and the energy relaxation sequence.

@param fname is the file name for the microstructure image
@param exx is the xx component of the applied strain
@param eyy is the yy component of the applied strain
@param ezz is the zz component of the applied strain
@param exz is the xz component of the applied strain
@param eyz is the yz component of the applied strain
@param exy is the xy component of the applied strain
*/
void Calc (string fname,
           double exx,
           double eyy,
           double ezz,
           double exz,
           double eyz,
           double exy);

/**
@brief Calculates the effective bulk modulus of the relaxed mesh.

@param fname is the file name for the microstructure image
@return the bulk modulus (GPa)
*/
double getBulkModulus (string fname);

};      // End of AppliedStrain class

#endif
