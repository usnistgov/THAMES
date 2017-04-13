# INSTALL

## INTRODUCTION

This file describes how to install THAMES.  The instructions in this file
are for the most common use cases, and cover the command line tools.

For further information, or in case of problems, please contact the author,
Jeff Bullard (jeffrey.bullard@nist.gov).

The API documentation of THAMES is available in PDF and HTML formats.  It
is not bundled with the software but can be created during the installation
of the software.

For more information about THAMES, see the API documentation overview, or
the accompanying user guide (coming soon).  More information about THAMES and
its applications can also be found in the following references:

    * Bullard, J.W., Lothenbach, B., Stutzman, P.E., Snyder, K.A., Coupling thermodynamics and digital image models to simulate hydration and microstructure development of portland cement pastes, _Journal of Materials Research_ 26, (2011) 609-622.

    * Feng, P., Miao, C., Bullard, J.W., A model of phase stability, microstructure and properties during leaching of portland cement binders, _Cement and Concrete Composites_ 49, (2014) 9-19.

    * Li, X., Grasley, Z.C., Garboczi, E.J., Bullard, J.W., Modeling the apparent intrinsic viscoelastic relaxation of hydrating cement paste, _Cement and Concrete Composites_ 55, (2014) 322-330.

    * Feng, P., Garboczi, E.J., Miao, C., Bullard, J.W., Microstructural origins of cement paste degradation by external sulfate attack, _Construction and Building Materials_ 96, (2015) 391-403.

    * Li, X., Grasley, Z.C., Bullard, J.W., Garboczi, E.J., Computing the time evolution of the apparent viscoelastic/viscoplastic Poisson's ratio of hydrating cement paste, _Cement and Concrete Composites_, 56 (2015) 121-133.


## PREREQUISITES

To install THAMES, you need 'cmake'.  Creation of the API documentation also
requires 'doxygen'.

    * CMake (>= 2.6), the build system used by THAMES
        Required for building THAMES

    * Doxygen (>= 1.8.13), the API documentation software
        Required for creating the API documentation

    * LaTeX 2e, the document preparation system
        Required only for creating the PDF version of the API documentation


## BUILDING

The recommended way to configure THAMES is to do an out-of-source build,
which means that the original files and directories are left untouched.
Doing this makes the re-compiling and cleaning of the installation files
much simpler.

From the directory that contains this INSTALL file, execute the following
two commands:

    cd build
    cmake ..

The cmake command will assess your operating system's and C++ compilers
and will automatically create the Makefiles needed to build THAMES.  When
this step is complete, then you can build the libraries and executables
from the same directory by executing the command

    make                    (or nmake in a Windows command prompt)

All of the libraries and the executable will be in the build directory and
its subdirectories.  You can install the executable and libraries in the
top-level directories of your working directory by executing the command

    make install            (or nmake install in a Windows command prompt)

This will install the "thames" executable in the bin/ directory, and the
static libraries in the lib/ directory.

Finally, to create the documentation for the THAMES API, execute the command

    make doc                (or nmake doc in a Windows command prompt)



## UNINSTALLING

To uninstall everything except the original source files, simply delete
recursively everything in the build/ directory.
