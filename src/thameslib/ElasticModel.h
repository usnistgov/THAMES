/**
@file ElasticModel.h
@brief Declaration of the ElasticModel base class for finite element calculations.

@section Background  
This program solves the linear elastic equations in a
random linear elastic material, subject to an applied macroscopic strain,
using the finite element method.  Each pixel in the 3D digital
image is a cubic tri-linear finite element,  having its own
elastic moduli tensor. Periodic boundary conditions are maintained.
In the comments below, (USER) means that this is a section of code that
the user might have to change for his particular problem. Therefore the
user is encouraged to search for this string.

@section Problem and Variable Definitions
The problem being solved is the minimization of the energy
\f$1/2 u \cdot A \cdot u + b \cdot u + C\f$, where <i>A</i> is
the Hessian matrix composed of the stiffness matrices (`dk`) for each
pixel/element, <i>b</i> is a constant vector and <i>C</i> is a constant
that are determined by the applied strain and
the periodic boundary conditions, and <i>u</i> is a vector of
all the displacements. The solution
method used is the conjugate gradient relaxation algorithm.
Other variables are:  `gb` is the gradient = Au+b, `h` is an
auxiliary variable used in the conjugate gradient algorithm (in dembx),
`dk`(n,i,j) is the stiffness matrix of the n'th phase, `cmod`(n,i,j) is
the elastic moduli tensor of the n'th phase, `pix` is a vector that gives
the phase label of each pixel, `ib` is a matrix that gives the labels of
the 27 (counting itself) neighbors of a given node, `prob` is the volume
fractions of the various phases,
`strxx`, `stryy`, `strzz`, `strxz`, `stryz`, and `strxy` are the six Voigt
volume averaged total stresses, and
`sxx`, `syy`, `szz`, `sxz`, `syz`, and `sxy` are the six Voigt
volume averaged total strains.

@section Dimensions
The vectors `u`,`gb`,`b`, and `h` are dimensioned to be the system size,
`ns`=nx*ny*nz, with three components, where the digital image of the
microstructure considered is a rectangular paralleliped, nx x ny x nz
in size.  The arrays `pix` and `ib` are are also dimensioned to the system
size.  The array ib has 27 components, for the 27 neighbors of a node.
Note that the program is set up at present to have at most 50
different phases.  This can easily be changed, simply by changing
the dimensions of `dk`, `prob`, and `cmod`. The parameter `nphase` gives the
number of phases being considered in the problem.
All arrays are passed between subroutines using simple common statements.

Manual available at NISTIR 6269 from NTIS or at:
http://ciks.cbt.nist.gov/~garbocz/manual/man.html

@todo Devise better and more descriptive function names

@author Edward J. Garboczi (original Fortran algorithms)
@author Pan Feng, Jeffrey W. Bullard (C/C++ versions)
@author Jeffrey W. Bullard (Documentation)

@warning Read the manual before using this program!!!
*/

#ifndef ELASTICMODEL_H
#define ELASTICMODEL_H

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "StrainEnergy.h"

using namespace std;

/**
@class Declare the ElasticModel base class for finite element calculations
*/
class ElasticModel {

protected:

    int nx_;                            /**< Number of voxels in x dimension */
    int ny_;                            /**< Number of voxels in y dimension */
    int nz_;                            /**< Number of voxels in z dimension */
    int ns_;                            /**< Number of total voxels */
    int nphase_;                        /**< Maximimum allowed number of phases */
    int npoints_;                       /**< Number of microstructures to process */
    string phasemod_fname_;             /**< File name saving the phase elastic moduli */
    vector<vector<double> > u_;         /**< 2D displacement field */
    vector<vector<double> > gb_;        /**< Energy gradient */
    vector<vector<double> > b_;         /**< Coefficient of linear displacement in energy */
    vector<vector<double> > h_;         /**< Auxiliary conjugate gradient variable */
    vector<vector<double> > Ah_;        /**< Local stiffness matrix */

    /**
    @brief Array of elastic moduli tensors of each phase.

    The first dimension of this 3D vector is the index of the phase, while the
    second and third are the first and second index of the elastic moduli tensor,
    written with engineering notation:

        - 1 = 11
        - 2 = 22
        - 3 = 33
        - 4 = 23 (and 32)
        - 5 = 13 (and 31)
        - 6 = 12 (and 21)

    This means that, for example, `cmod(4,1,4)` is the tensor component 1123 (or 1132)
    of the 4th phase in the system.
    */
    vector<vector<vector<double> > > cmod_;

    /**
    @brief Array of stiffness matrices for each phase.

    The stiffness matrix relates the local displacement field of an element (defined
    at its eight corners) to the local stress tensor components, and is the centeral
    material property used to compute the contributions to the quadratic term of
    the energy function.  The indices have the following meaning

        - First index is the phase index number
        - Second index can be 0 to 8, defining which corner of the cubic element
        - Third index can be 0 to 3, defining the components of the displacement at that corner
        - Fourth index can be 0 to 8, defining which corner of the cubic element
        - Fifth index can be 0 to 3, defining the components of the displacement at that corner
    */
    vector<vector<vector<vector<vector<double> > > > > dk_;

    /**
    @brief The bulk and shear moduli of each phase.

    This is a 2D array.  The first index varies over the number of phases.  The
    second index is either 0 (for bulk modulus) or 1 (for shear modulus), both
    in units of GPa.
    */
    vector<vector<double> > phasemod_;

    vector<double> prob_;                   /**< Array of phase volume fractions */

    /**
    @brief The neighbor table defining the mesh topology.

    The elements in the finite element mesh are indexed as a 1D array, not the
    usual 3D array.  The conversion from 3D indices (<i>x</i>,<i>y</i>,<i>z</i>)
    to a unique 1D index, <i>k</i>, is
    \f{equation}
        k = n_x n_y z + n_x y + x
    \f}
    where <i>n</i><sub>x</sub>, <i>n</i><sub>y</sub>, and <i>n</i><sub>z</sub>
    are the <i>x</i>, <i>y</i>, and <i>z</i> dimensions of the mesh, respectively.

    The neighor table stores the list of nearest-neighbor indices for each element.
    */
    vector<vector<int> > ib_;


    vector<int> pix_;                       /**< Phase index label at each element */
    vector<vector<double> > elestress_;     /**< Stress components at each element */
    vector<vector<double> > elestrain_;     /**< Strain components at each element */
    vector<double> strainengy_;             /**< Local strain energy in each element */
    vector<double> avgStrainengy_;          /**< Average strain energy of each phase */
	
    double strxx_;                          /**< Prescribed stress xx component */
    double stryy_;                          /**< Prescribed stress yy component */
    double strzz_;                          /**< Prescribed stress zz component */
    double strxz_;                          /**< Prescribed stress xz component */
    double stryz_;                          /**< Prescribed stress yz component */
    double strxy_;                          /**< Prescribed stress xy component */
    double sxx_;                            /**< Prescribed strain xx component */
    double syy_;                            /**< Prescribed strain yy component */
    double szz_;                            /**< Prescribed strain zz component */
    double sxz_;                            /**< Prescribed strain xz component */
    double syz_;                            /**< Prescribed strain yz component */
    double sxy_;                            /**< Prescribed strain xy component */
    double C_;                              /**< Constant term in the energy equation */
	
    /**
    @brief The stopping criterion determining convergence to solution.

    This the maximum allowed square of the energy gradient for which the whole system can be
    considered to have converged to the energy minimum.  This is the value to which
    the quantity `gg_` (the computed square of the energy gradient, `gb_`)is compared.
    In most cases, `gtest_` will be set equal to a value <i>abc</i> times the total
    number of elements in the mesh, so that when `gg_` is less than `gtest_`, the
    RMS gradient per element is less than \f$\sqrt{abc}\f$.
    */
    double gtest_;          

    double gg_;                 /**< The square of the gradient, `gb_`<sup>2</sup> */

 
	
public:

/**
@brief The only constructor provided with the base class.

@param nx is the number of elements in the x direction
@param ny is the number of elements in the y direction
@param nz is the number of elements in the z direction
@param dim is the total number of elements in the system (plus two?)
@param npoints is the number of microstructures to process
*/
ElasticModel (int nx,
              int ny,
              int nz,
              int dim,
              int nphase,
              int npoints);

/**
@brief Set up the elastic modulus variables.

@todo Change the name of this method to setElasticModuli.

@param phasemod_fname is the file name for saving the phase moduli data
@param nphase is the maximum allowed number of phases in the microstructure
*/
void ElasModul (string phasemod_fname,
                int nphase);


/**
@brief Construct the neighbor table for each element.

The neighbor table has 27 elements for each finite element, which is
the combined number of face neighbrs (6), edge neighbors (12), and
corner neighbrs (8), plus the element itself
*/
void BuildNeighbor ();

/**
@brief Set one of the elastic moduli components of a phase.

It is up to the user to know whether the components should be the
Young's modulus and Poisson's ratio, or the the bulk and shear moduli

@todo Implement bounds checking to avoid segmentation violations.

@todo Change the name of this method to setPhaseModuli.

@param phaseid is the index of the phase being set
@param i is 0 (Young's or bulk modulus) or 1 (Poisson's ratio or shear modulus)
@param val is the value to set for the component [GPa or dimensionless]
*/
void setPhasemod (int phaseid,
                  int i,
                  double val)
{
    phasemod_[phaseid][i] = val;
    return;
}


/**
@brief Get one of the elastic moduli components of a phase.

It is up to the user to know whether the components should be the
Young's modulus and Poisson's ratio, or the the bulk and shear moduli

@todo Implement bounds checking to avoid segmentation violations.

@todo Change the name of this method to getPhaseModuli.

@param phaseid is the index of the phase being set
@param i is 0 (Young's or bulk modulus) or 1 (Poisson's ratio or shear modulus)
@return the value of the modulus sought [GPa or dimensionless]
*/
double getPhasemod (int phaseid,
                    int i)
{
    return phasemod_[phaseid][i];
}
	
/**
@brief Read a microstructure and set up the phase assignments for each element.

@todo Devise better and more descriptive function names.

@todo Change the name of this method to setMicrostructure.

@param fname is the input file containing the microstructure data
@param nphase is the maximum allowed number of phases in the microstructure
*/
void ppixel (string fname,
             int nphase);

/**
@brief Determine the volume fractions of the different phases.

@todo Change the name of this method to setVolumeFractions.

@param ns is the total number of elements in the mesh
@param nphase is the maximum allowed number of phases in the microstructure
*/
void assig (int ns,
            int nphase)
{
    for (int i = 0; i < nphase; i++) {
      prob_[i] = 0.0;
    }

    for (int m = 0; m < ns; m++) {
      prob_[pix_[m]] += 1;
    }

    for (int i = 0; i < nphase; i++ ) {
      prob_[i] = prob_[i] / (float) ns;
    }

    return;
}

/**
@brief Determine average strain energy stored by each GEM dependent component (DC).

@todo Change the name of this method to getAverageStrainEnergy.

@remarks There is a variable used in this function, called `strainenergy`, which
is defined in the StrainEnergy.h header file, and it is unclear what this
variable does or why it is set the way it is.  It does not appear to be used
anywhere else in the program.
*/
void getAvgStrainengy ();

/**
@brief Get a component of the prescribed stress tensor.

@todo Put in bounds checking by exception handling.  Will there be performance
degradation by doing this?

@param i is the component index
@return the value of the i-th component of the prescribed stress tensor
*/
double getStress (int i)
{
    double val = 0.0;
    switch (i) {
        case 0: 
            val = strxx_;
            break;
        case 1:
            val = stryy_;
            break;
        case 2:
            val = strzz_;
            break;
        case 3:
            val = strxz_;
            break;
        case 4:
            val = stryz_;
            break;
        case 5:
            val = strxy_;
            break;
        default:
            cout << "i (" << i << ") is not a recognized stress component" << endl;
            break;
    }

    return val;
}

/**
@brief Get a component of the prescribed strain tensor.

@todo Put in bounds checking by exception handling.  Will there be performance
degradation by doing this?

@param i is the component index
@return the value of the i-th component of the prescribed strain tensor
*/
double getStrain (int i)
{
    double val = 0.0;
    switch (i) {
        case 0: 
            val = sxx_;
            break;
        case 1:
            val = syy_;
            break;
        case 2:
            val = szz_;
            break;
        case 3:
            val = sxz_;
            break;
        case 4:
            val = syz_;
            break;
        case 5:
            val = sxy_;
            break;
        default:
            cout << "i (" << i << ") is not a recognized strain component" << endl;
            break;
    }

    return val;
}

/**
@brief Get the stress component at a particular finite element.

@todo Put in bounds checking by exception handling.  Will there be performance
degradation by doing this?

@param m is the element index
@param index is the stres component to retrieve
@return the ith stress component at the element m [GPa]
*/
double getElestress(int m,int index)
{
    return elestress_[m][index];

}

/**
@brief Set up the stiffness matrices and parameters for computing energy.

This is a virtual method which does nothing in the base class, but which
will be customized for each of the derived classes.

@param nx is the number of voxels in the x direction
@param ny is the number of voxels in the y direction
@param nz is the number of voxels in the z direction
@param ns is the total number of finite elements
@param nphase is the maximum allowed number of phases in the microstructure
*/
virtual void femat (int nx,
                    int ny,
                    int nz,
                    int ns,
                    int nphase)
{
    cout << "virtual function 'femat' in base class." << endl;	
    return;
}

/**
@brief Perform the congugate gradient energy minimization.

This is a virtual method which does nothing in the base class, but
which will be customized for each of the derived classes.

@param ns is the total number of finite elements
@param gg is the square of the energy gradient
@param ldemb is the maximum number of conjugate gradient steps to perform
@param kkk is the number of times this function has already been called
@return the number of iterations performed during this call
*/
virtual long int dembx (int ns,
                        double gg,
                        int ldemb,
                        int kkk)
{
    cout << "virtual function 'dembx' in base class." << endl;
    return 0;
}

/**
@brief Compute the total mesh elastic strain energy

This is a virtual method which does nothing in the base class, but which
will be customized for each of the derived classes.

@todo Determine the energy units

@param nx is the number of voxels in the x direction
@param ny is the number of voxels in the y direction
@param nz is the number of voxels in the z direction
@param ns is the total number of finite elements
@return the energy [units?]
*/
virtual double energy (int nx,
                       int ny,
                       int nz,
                       int ns)
{
    cout << "virtual function 'energy' in base class." << endl;
    return 0.0;
}

/**
@brief Calculate the stress and strain fields.

This is a virtual method which does nothing in the base class, but which
will be customized for each of the derived classes.

@param nx is the number of voxels in the x direction
@param ny is the number of voxels in the y direction
@param nz is the number of voxels in the z direction
@param ns is the total number of finite elements
*/
virtual void stress (int nx,
                     int ny,
                     int nz,
                     int ns)
{
    cout << "virtual function 'stress' in base class." << endl;
    return;
}

/**
@brief Controlling method for energy minimization.

This is a virtual method which does nothing in the base class, but which
will be customized for each of the derived classes.

@param time is the simulation time [days]
@param kmax is the maximum number of times to call the conjugate gradient solver
*/
virtual void relax (double time,
                    int kmax)
{
    cout << "virtual function 'relax' in base class." << endl;
    return;
}

/**
@brief The master controlling function for the finite element method

This is a virtual method which does nothing in the base class, but which
will be customized for each of the derived classes.

@param time is the simulation time [days]
@param fname is the file name to store the results
@param exx is the input prescribed xx component of the strain
@param eyy is the input prescribed yy component of the strain
@param ezz is the input prescribed zz component of the strain
@param exz is the input prescribed xz component of the strain
@param eyz is the input prescribed yz component of the strain
@param exy is the input prescribed xy component of the strain
*/
virtual void Calc (double time,
                   string fname,
                   double exx,
                   double eyy,
                   double ezz,
                   double exz,
                   double eyz,
                   double exy)
{
    cout << "virtual function 'thrMic' in base class." << endl;
    return;
}

/**
@brief Create visualization of the stress field in a 2D microstructure slice.

The method can output full color image, but currently just uses different
intensities of cyan (green + blue channels equal, red zero).
The file is created as a portable pixel map (PPM) file and then manually
converted to a PNG file using ImageMagick's convert command

@todo Use exception handling to check the index requested

@todo Use the libpng library to avoid having to rely on system call to ImageMagick

@param root is the root name of the file to create
@param time is the simulation time [days]
@param index is the stress component to visualize (values 0 to 5)
*/
void writeStress (string &root,
                  double time,
                  int index);

/**
@brief Create visualization of the strain field in a 2D microstructure slice.

The file is created as a portable pixel map (PPM) file and then manually
converted to a PNG file using ImageMagick's convert command

@todo Use exception handling to check the index requested

@todo Use the libpng library to avoid having to rely on system call to ImageMagick

@param root is the root name of the file to create
@param time is the simulation time [days]
@param index is the strain component to visualize (values 0 to 5)
*/
void writeStrain (string &root,
                  double time,
                  int index);


/**
@brief Create data file storing the 3D displacement field

@param root is the root name of the file to create
@param time is the simulation time [days]
*/
void writeDisp (string &root,
                double time);
    
/**
@brief Create visualization of the strain energy in a 2D microstructure slice.

The file is created as a portable pixel map (PPM) file and then manually
converted to a PNG file using ImageMagick's convert command

@todo Use the libpng library to avoid having to rely on system call to ImageMagick

@param root is the root name of the file to create
@param time is the simulation time [days]
*/
void writeStrainEngy (string &root,
                      double time);

/**
@brief Set the name of the file containing the phase elastic moduli

@param phasemod_fname is the file name to set
*/
void setPhasemodfname (string phasemod_fname)
{
    phasemod_fname_ = phasemod_fname;
    return;
}

/**
@brief Get the name of the file containing the phase elastic moduli

@return the phase modulus input file name
*/
string getPhasemodfname (void)
{
    return phasemod_fname_;
}

};      // End of ElasticModel class

#endif
