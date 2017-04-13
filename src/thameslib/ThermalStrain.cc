/**
@file ThermalStrain.cc
@brief Defines methods for the ThermalStrain class.

*/
#include "ThermalStrain.h"
#include <iostream>

using namespace std;

ThermalStrain::ThermalStrain (int nx,
                              int ny,
                              int nz,
                              int dim,
                              int nphase,
                              int npoints)
    : ElasticModel(nx,ny,nz,dim,nphase,npoints) 
{ 
    cout << "constructor 'THERMAL3D' in derived class 'THERMAL3D'." << endl;

    isfirst_ = true;
		
    Y_ = 0.0;
    
    localgg_ = 0.0;    

    boxsize_ = 27;
    boxnum_ = (int)pow((double)((int)(boxsize_ / 2) * 2 + 1),(double)(3));
    localgtest_ = 1.0e-30 * boxnum_;
    cout << "localgtest_ is: " << localgtest_ << endl;    

    eigen_.clear();
    eigen_.resize(ns_);
    for (int n = 0; n < ns_; n++) {
        eigen_[n].resize(6,0.0);
    }

    b0_.clear();
    b1_.clear();
    b2_.clear();
    b3_.clear();
    b4_.clear();
    b5_.clear();
    b0_.resize(dim);
    b1_.resize(dim);
    b2_.resize(dim);
    b3_.resize(dim);
    b4_.resize(dim);
    b5_.resize(dim);
    for (int n = 0; n < dim; n++) {
        b0_[n].resize(3,0.0);
        b1_[n].resize(3,0.0);
        b2_[n].resize(3,0.0);
        b3_[n].resize(3,0.0);
        b4_[n].resize(3,0.0);
        b5_[n].resize(3,0.0);
    }

    exp_.clear();
	
    T_.clear();
    T_.resize(dim);
    for (int n = 0; n < dim; n++) {
      T_[n].resize(3,0.0);
    }

    tstrength_.clear();
    tstrength_.resize(nphase_, 5.0);

    ///
    /// Indexes for CSH, ettringite, and hydrotalcite are hard-wired into
    /// the code.
    ///
    /// @todo Generalize this so that is not fragile to changes in the indexing
    /// of phases.
    ///

    tstrength_[12] = 1.0; // CSH
    tstrength_[15] = 1.0; // ETTR
    tstrength_[17] = 1.0; // HYDROTALC
		
    zcon_.clear();
    zcon_.resize(2);
    for (int i = 0; i < 2; i++) {
      zcon_[i].resize(3);
      for (int k = 0; k < 3; k++) {
        zcon_[i][k].resize(2);
        for (int j = 0; j < 2; j++) {
          zcon_[i][k][j].resize(3,0.0);
        }
      }
    }
		
    ss_.clear();
    ss_.resize(ns_);
    for (int i = 0; i < ns_; i++) {
      ss_[i].resize(8);
      for (int j = 0; j < 8; j++) {
        ss_[i][j].resize(3,0.0);
      }
    }
}

void ThermalStrain::femat (int nx,
                           int ny,
                           int nz,
                           int ns,
                           int nphase,
                           int iskip) 
{
    double dndx[8],dndy[8],dndz[8];
    double g[3][3][3];
    double es[6][8][3],delta[8][3];
    int is[8];
    double x,y,z;
    int nxy = nx * ny;
	
	
    ///
    /// Generate dk, zcon, T and Y on first pass. After that they are constant,    
    /// since they are independent of the macrostrains. Only b gets upgraded       
    /// as the macrostrains change.                                                
    ///

    if (iskip == 0) {
	  
      ///
      /// Initialize stiffness matrices.                                             
      ///

      for (int m = 0; m < nphase; m++) {
        for (int l = 0; l < 3; l++) {
          for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 8; j++) {
              for (int i = 0; i < 8; i++) {
                dk_[m][i][k][j][l] = 0.0;
              }
            }
          }
        }
      }
	  
      ///
      /// Initialize zcon matrix (gives C term for arbitrary macrostrains)           
      ///

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          for (int mi = 0; mi < 3; mi++) {
            for (int mj = 0; mj < 3; mj++) {
              zcon_[i][mi][j][mj] = 0.0;
            }
          }
        }
      }  
	  
      ///
      /// @note An anisotropic elastic moduli tensor could be input at this point,   
      /// bypassing this part, which assumes isotropic elasticity, so that there    
      /// are only two independent numbers making up the elastic moduli tensor,    
      /// the bulk modulus K and the shear modulus G.                               
      ///

      ElasModul(phasemod_fname_,nphase_);
	   
      ///
      /// Set up Simpson's integration rule weight vector                          
      ///

      for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            int nm = 0;
            if (i == 1) nm++;
            if (j == 1) nm++;
            if (k == 1) nm++;
            g[i][j][k] = pow((double)4,(double)nm);
          }
        }
      }
	 
      ///
      /// Loop over the nphase kinds of pixels and Simpson's rule quadrature points  
      /// in order to compute the stiffness matrices. Stiffness matrices of          
      /// trilinear finite elements are quadratic in x, y, and z, so that Simpson's  
      /// rule quadrature gives exact results.                                       
      ///

      for (int ijk = 0; ijk < nphase; ijk++) {
        for (int k = 0; k < 3; k++) {
          for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
              x = (float)(i) / 2.0;
              y = (float)(j) / 2.0;
              z = (float)(k) / 2.0;

              ///
              /// dndx means the negative derivative with respect to x, of the shape         
              /// matrix N (see manual, Sec. 2.2), dndy and dndz are similar.                
              ///

              dndx[0] = (-(1.0 - y)) * (1.0 - z);
              dndx[1] = (1.0 - y) * (1.0 - z);
              dndx[2] = y * (1.0 - z);
              dndx[3] = (-y) * (1.0 - z);
              dndx[4] = (-(1.0 - y)) * z;
              dndx[5] = (1.0 - y) * z;
              dndx[6] = y * z;
              dndx[7] = (-y) * z;
              dndy[0] = (-(1.0 - x)) * (1.0 - z);
              dndy[1] = (-x) * (1.0 - z);
              dndy[2] = x * (1.0 - z);
              dndy[3] = (1.0 - x) * (1.0 - z);
              dndy[4] = (-(1.0 - x)) * z;
              dndy[5] = (-x) * z;
              dndy[6] = x * z;
              dndy[7] = (1.0 - x) * z;
              dndz[0] = (-(1.0 - x)) * (1.0 - y);
              dndz[1] = (-x) * (1.0 - y);
              dndz[2] = (-x) * y;
              dndz[3] = (-(1.0 - x)) * y;
              dndz[4] = (1.0 - x) * (1.0 - y);
              dndz[5] = x * (1.0 - y);
              dndz[6] = x * y;
              dndz[7] = (1.0 - x) * y;

              ///
              /// Now build strain matrix                                                  
              ///

              for (int n1 = 0; n1 < 6; n1++) {
                for (int n2 = 0; n2 < 8; n2++) {
                  for (int n3 = 0; n3 < 3; n3++) {
                    es[n1][n2][n3] = 0.0;
                  }
                }
              }
              for(int n = 0; n < 8; n++) {
                es[0][n][0] = dndx[n];
                es[1][n][1] = dndy[n];
                es[2][n][2] = dndz[n];
                es[3][n][0] = dndz[n];
                es[3][n][2] = dndx[n];
                es[4][n][1] = dndz[n];
                es[4][n][2] = dndy[n];
                es[5][n][0] = dndy[n];
                es[5][n][1] = dndx[n];	  
              }
              
              ///
              /// Matrix multiply to determine value at (x,y,z),
              /// multiply by proper weight, 
              /// and sum into dk, the stiffness matrix.                                    
              ///

              for (int mm = 0; mm < 3; mm++) {
                for (int nn = 0; nn < 3; nn++) {
                  for (int ii = 0; ii < 8; ii++) {
                    for (int jj = 0; jj < 8; jj++) {
                      double sum = 0.0;
                      for (int kk = 0; kk < 6; kk++) {
                        for (int ll = 0; ll < 6; ll++) {
                          sum += es[kk][ii][mm] * cmod_[ijk][kk][ll] 
                               * es[ll][jj][nn];
                        }
                      }
                      dk_[ijk][ii][mm][jj][nn] += g[i][j][k] * sum / 216.0;
                    }
                  }
                }
              }
            }
          }
        }
      }

      ///
      /// Now compute the ss matrices, which give the thermal strain terms for      
      /// the i'th phase, single pixel.                                             
      ///

      dndx[0] = -0.25;
      dndx[1] = 0.25;
      dndx[2] = 0.25;
      dndx[3] = -0.25;
      dndx[4] = -0.25;
      dndx[5] = 0.25;
      dndx[6] = 0.25;
      dndx[7] = -0.25;
      dndy[0] = -0.25;
      dndy[1] = -0.25;
      dndy[2] = 0.25;
      dndy[3] = 0.25;
      dndy[4] = -0.25;
      dndy[5] = -0.25;
      dndy[6] = 0.25;
      dndy[7] = 0.25;
      dndz[0] = -0.25;
      dndz[1] = -0.25;
      dndz[2] = -0.25;
      dndz[3] = -0.25;
      dndz[4] = 0.25;
      dndz[5] = 0.25;
      dndz[6] = 0.25;
      dndz[7] = 0.25;

      ///
      /// Build average strain matrix.                                           
      ///

      for (int n1 = 0; n1 < 6; n1++) {
        for (int n2 = 0; n2 < 8; n2++) {
          for (int n3 = 0; n3 < 3; n3++) {
            es[n1][n2][n3] = 0.0;
          }
        }
      }
	  
      for (int n = 0; n < 8; n++) {
        es[0][n][0] = dndx[n];
        es[1][n][1] = dndy[n];
        es[2][n][2] = dndz[n];
        es[3][n][0] = dndz[n];
        es[3][n][2] = dndx[n];
        es[4][n][1] = dndz[n];
        es[4][n][2] = dndy[n];
        es[5][n][0] = dndy[n];
        es[5][n][1] = dndx[n];
      }
	  
      for (int mmm = 0; mmm < ns; mmm++) {
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int nm = 0; nm < 6; nm++) {
              for (int n = 0; n < 6; n++) {
                sum += cmod_[pix_[mmm]][n][nm] * es[n][mm][nn] * eigen_[mmm][nm];
              }
            }
            ss_[mmm][mm][nn] = sum;
          }
        }
      }
	  
      ///
      /// Call method to generate `zcon_`                           
      /// `zcon_` is a (2,3) X (2,3) matrix.                                           
      ///

      constfunc(ns_,nx_,ny_,nz_);
	  
      ///
      /// Set up linear term, T_, for thermal energy. It does not depend on the  
      /// macrostrains or displacement, so there is no need to update it as the      
      /// macrostrains change. T_ is built up out of the ss_ matrices.               
      ///
	  
      int nss = ns + 2;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          T_[m][m3] = 0.0;
        }
      }
	  
      ///
      /// For all cases, the correspondence between 0-7 finite element node (corner)
      /// labels  and 0-26 neighbor labels is (see Table 4 in manual):                       
      ///
      ///   - 0 = ib_[m][26]
      ///   - 1 = ib_[m][2]
      ///   - 2 = ib_[m][1]
      ///   - 3 = ib_[m][0]
      ///   - 4 = ib_[m][25]
      ///   - 5 = ib_[m][18]
      ///   - 6 = ib_[m][17]
      ///   - 7 = ib_[m][16]                                                  
      ///

      is[0] = 26;    
      is[1] = 2;
      is[2] = 1;   
      is[3] = 0;   
      is[4] = 25;   
      is[5] = 18;    
      is[6] = 17;   
      is[7] = 16;  
	  
      ///
      /// Do all points, but no macrostrain terms                                    
      ///
      /// @note Factor of 2 on linear thermal term is cancelled                      
      /// by factor of 1/2 out in front of total energy term.                        
      ///

      for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
          for (int i = 0; i < nx; i++) {
            int m = nxy * k + nx * j + i;
            for (int mm = 0; mm < 8; mm++) {
              for (int nn = 0; nn < 3; nn++) {
                T_[ib_[m][is[mm]]][nn] -= ss_[m][mm][nn];
              }
            }
          }
        }
      }
	  
      ///
      /// Pick up and sum in all terms multiplying macrostrains. 
      ///

      for (int ipp = 0; ipp < 2; ipp++) {
        for (int jpp = 0; jpp < 3; jpp++) {
          double exx = 0.0;
          double eyy = 0.0;
          double ezz = 0.0;
          double exz = 0.0;
          double eyz = 0.0;
          double exy = 0.0;
          if ((ipp == 0) && (jpp == 0)) exx = 1.0;
          if ((ipp == 0) && (jpp == 1)) eyy = 1.0;
          if ((ipp == 0) && (jpp == 2)) ezz = 1.0;
          if ((ipp == 1) && (jpp == 0)) exz = 1.0;
          if ((ipp == 1) && (jpp == 1)) eyz = 1.0;
          if ((ipp == 1) && (jpp == 2)) exy = 1.0;
		  
          // x = nx - 1, face 
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
                delta[i8][0] = exx * (double)nx;
                delta[i8][1] = exy * (double)nx;
                delta[i8][2] = exz * (double)nx;
              }
            }
          }
          for (int j = 0; j < (ny - 1); j++) {
            for (int k = 0; k < (nz - 1); k++) {
              int m = nxy * k + nx * j + (nx - 1);
              for (int nn = 0; nn < 3; nn++) {
                for (int mm = 0; mm < 8; mm++) {
                  T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
                }
              }
            }
          }
		  
          // y = (ny - 1), face
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
                delta[i8][0] = exy * (double)ny;
                delta[i8][1] = eyy * (double)ny;
                delta[i8][2] = eyz * (double)ny;
              }
            }
          }
          for (int i = 0; i < (nx - 1); i++) {
            for (int k = 0; k < (nz - 1); k++) {
              int m = nxy * k + nx * (ny - 1) + i;
              for (int nn = 0; nn < 3; nn++) {
                for (int mm = 0; mm < 8; mm++) {
                  T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
                }
              }
            }
          }
	  
          // z = (nz - 1), face 
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
                delta[i8][0] = exz * (double)nz;
                delta[i8][1] = eyz * (double)nz;
                delta[i8][2] = ezz * (double)nz;
              }
            }
          }
          for (int i = 0; i < (nx - 1); i++) {
            for (int j = 0; j < (ny - 1); j++) {
              int m = nxy * (nz - 1) + nx * j + i;
              for (int nn = 0; nn < 3; nn++) {
                for (int mm = 0; mm < 8; mm++) {
                  T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
                }
              }
            }
          }

          // x = (nx - 1), y = (ny - 1), edge  
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if ((i8 == 1) || (i8 == 5)) {
                delta[i8][0] = exx * (double)nx;
                delta[i8][1] = exy * (double)nx;
                delta[i8][2] = exz * (double)nx;
              }
              if ((i8 == 3) || (i8 == 7)) {
                delta[i8][0] = exy * (double)ny;
                delta[i8][1] = eyy * (double)ny;
                delta[i8][2] = eyz * (double)ny;
              }
              if ((i8 == 2) || (i8 == 6)) {
                delta[i8][0] = exy * (double)ny + exx * (double)nx;
                delta[i8][1] = eyy * (double)ny + exy * (double)nx;
                delta[i8][2] = eyz * (double)ny + exz * (double)nx;
              }
            }
          }	
          for (int k = 0; k < (nz - 1); k++) {
            int m = nxy * k + nx * (ny - 1) + (nx - 1);
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }

          // x = (nx - 1), z = (nz - 1), edge 
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i3][i8] = 0.0;
              if ((i8 == 1) || (i8 == 2)) {
                delta[i8][0] = exx * (double)nx;
                delta[i8][1] = exy * (double)nx;
                delta[i8][2] = exz * (double)nx;
              }
              if ((i8 == 4) || (i8 == 7)) {
                delta[i8][0] = exz * (double)nz;
                delta[i8][1] = eyz * (double)nz;
                delta[i8][2] = ezz * (double)nz;
              }
              if ((i8 == 5) || (i8 == 6)) {
                delta[i8][0] = exz * (double)nz + exx * (double)nx;
                delta[i8][1] = eyz * (double)nz + exy * (double)nx;
                delta[i8][2] = ezz * (double)nz + exz * (double)nx; 
              }
            }
          }
          for (int j = 0; j < (ny - 1); j++) {
            int m = nxy * (nz - 1) + nx * j + (nx - 1);
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }

          // y = (ny - 1), z = (nz - 1), edge 
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if ((i8 == 4) || (i8 == 5)) {
                delta[i8][0] = exz * (double)nz;
                delta[i8][1] = eyz * (double)nz;
                delta[i8][2] = ezz * (double)nz;
              }
              if ((i8 == 2) || (i8 == 3)) {
                delta[i8][0] = exy * (double)ny;
                delta[i8][1] = eyy * (double)ny;
                delta[i8][2] = eyz * (double)ny;
              }
              if ((i8 == 6) || (i8 == 7)) {
                delta[i8][0] = exy * (double)ny + exz * (double)nz;
                delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
                delta[i8][2] = eyz * (double)ny + ezz * (double)nz; 
              }
            }
          }
          for (int i = 0; i < (nx - 1); i++) {
            int m = nxy * (nz - 1) + nx * (ny - 1) + i;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }

          // x = (nx - 1), y = (ny - 1), z = (nz - 1), corner 
          for (int i3 = 0; i3 < 3; i3++) {
            for (int i8 = 0; i8 < 8; i8++) {
              delta[i8][i3] = 0.0;
              if (i8 == 1) {
                delta[i8][0] = exx * (double)nx;
                delta[i8][1] = exy * (double)nx;
                delta[i8][2] = exz * (double)nx;
              }
              if (i8 == 3) {
                delta[i8][0] = exy * (double)ny;
                delta[i8][1] = eyy * (double)ny;
                delta[i8][2] = eyz * (double)ny;
              }
              if (i8 == 4) {
                delta[i8][0] = exz * (double)nz;
                delta[i8][1] = eyz * (double)nz;
                delta[i8][2] = ezz * (double)nz;
              }
              if (i8 == 7) {
                delta[i8][0] = exy * (double)ny + exz * (double)nz;
                delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
                delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
              }
              if (i8 == 5) {
                delta[i8][0] = exx * (double)nx + exz * (double)nz;
                delta[i8][1] = exy * (double)nx + eyz * (double)nz;
                delta[i8][2] = exz * (double)nx + ezz * (double)nz;
              }
              if (i8 == 2) {
                delta[i8][0] = exx * (double)nx + exy * (double)ny;
                delta[i8][1] = exy * (double)nx + eyy * (double)ny;
                delta[i8][2] = exz * (double)nx + eyz * (double)ny;
              }
              if (i8 == 6) {
                delta[i8][0] = exx * (double)nx + exy * (double)ny + exz * (double)nz;
                delta[i8][1] = exy * (double)nx + eyy * (double)ny + eyz * (double)nz;
                delta[i8][2] = exz * (double)nx + eyz * (double)ny + ezz * (double)nz;
              }
            }
          }
          int m = nxy * (nz - 1) + nx * (ny - 1) + (nx - 1);
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              T_[ns+ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
            }
          }		  
  
        }
      }
	  
      ///
      /// Compute Y, the 0.5 (eigenstrain<sub>i</sub>) C<sub>ij</sub> (eigenstrain<sub>j</sub>)
      /// energy, doesn't ever change with macrostrain or displacements.
      ///

      Y_ = 0.0;
      for (int m = 0; m < ns; m++) {
        for (int n = 0; n < 6; n++) {
          for (int nn = 0; nn < 6; nn++) {
            Y_ += 0.5 * eigen_[m][n] * cmod_[pix_[m]][n][nn] * eigen_[m][nn];
          }
        }
      }
    }

    ///
    /// Use auxiliary variables (exx, etc.) instead of u_[] variable,              
    /// for convenience, and to make the following code easier to read.            
    ///

    double exx = u_[ns][0];
    double eyy = u_[ns][1];
    double ezz = u_[ns][2];
    double exz = u_[ns+1][0];
    double eyz = u_[ns+1][1];
    double exy = u_[ns+1][2];

    ///
    /// Set up vector for linear term that comes from periodic boundary        
    /// conditions. Notation and conventions same as for `T_` term.                  
    /// This is done using the stiffness matrices, and the periodic terms in       
    /// the macrostrains. It is easier to set up `b_` this way than to analytically  
    /// write out all the terms involved.                                          
    ///

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns; m++) {
        b_[m][m3] = 0.0;
      }
    }
	
    ///
    /// For all cases, the correspondence between 0-7 finite element node (corner)
    /// labels  and 0-26 neighbor labels is (see Table 4 in manual):                       
    ///
    ///   - 0 = ib_[m][26]
    ///   - 1 = ib_[m][2]
    ///   - 2 = ib_[m][1]
    ///   - 3 = ib_[m][0]
    ///   - 4 = ib_[m][25]
    ///   - 5 = ib_[m][18]
    ///   - 6 = ib_[m][17]
    ///   - 7 = ib_[m][16]                                                  
    ///

    is[0] = 26;    
    is[1] = 2;
    is[2] = 1;   
    is[3] = 0;   
    is[4] = 25;   
    is[5] = 18;    
    is[6] = 17;   
    is[7] = 16;  
	
    C_ = 0.0;

    // x = nx - 1, face                                                       
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;	
          delta[i8][2] = exz * (double)nx;		  
        }
      }
    }
    for (int j = 0; j < (ny - 1); j++) {
      for (int k = 0; k < (nz - 1); k++) {
        int m = nxy * k + nx * j + (nx -1);
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }

    // y = ny - 1, face                                                           
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      for (int k = 0; k < (nz - 1); k++) {
        int m = nxy * k + nx * (ny - 1) + i;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }	
	
    // z = nz - 1, face                                              
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      for (int j = 0; j < (ny - 1); j++) {
        int m = nxy * (nz - 1) + nx * j + i;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }
	
    // x = nx - 1, y = ny - 1, edge                                      
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 5)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;	
          delta[i8][2] = exz * (double)nx;	
        }
        if ((i8 == 3) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if ((i8 == 2) || (i8 == 6)) {
          delta[i8][0] = exy * (double)ny + exx * (double)nx;
          delta[i8][1] = eyy * (double)ny + exy * (double)nx;	
          delta[i8][2] = eyz * (double)ny + exz * (double)nx;
        }
      }
    }
    for (int k = 0; k < (nz - 1); k++) {
      int m = nxy * k + nx * (ny - 1) + (nx -1);
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    // x = nx - 1, z = nz - 1, edge                                        
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;
          delta[i8][2] = exz * (double)nx;		  
        }
        if ((i8 == 4) || (i8 == 7)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if ((i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exz * (double)nz + exx * (double)nx;
          delta[i8][1] = eyz * (double)nz + exy * (double)nx;
          delta[i8][2] = ezz * (double)nz + exz * (double)nx;
        }
      }
    }
    for (int j = 0; j < (ny - 1); j++) {
      int m = nxy * (nz - 1) + nx * j + (nx - 1);
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    // y = ny - 1, z = nz - 1, edge                                            
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if ((i8 == 2) || (i8 == 3)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if ((i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny + exz * (double)nz;
          delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      int m = nxy * (nz - 1) + nx * (ny - 1) + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }

    // x = nx - 1, y = ny - 1, z = nz - 1, corner                                 
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if (i8 == 1) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;
          delta[i8][2] = exz * (double)nx;
        }
        if (i8 == 3) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if (i8 == 4) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if (i8 == 7) {
          delta[i8][0] = exy * (double)ny + exz * (double)nz; 
          delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
        }
        if (i8 == 5) {
          delta[i8][0] = exx * (double)nx + exz * (double)nz;
          delta[i8][1] = exy * (double)nx + eyz * (double)nz;
          delta[i8][2] = exz * (double)nx + ezz * (double)nz;
        }
        if (i8 == 2) {
          delta[i8][0] = exx * (double)nx + exy * (double)ny;
          delta[i8][1] = exy * (double)nx + eyy * (double)ny;
          delta[i8][2] = exz * (double)nx + eyz * (double)ny;
        }
        if (i8 == 6) {
          delta[i8][0] = exx * (double)nx + exy * (double)ny + exz * (double)nz;
          delta[i8][1] = exy * (double)nx + eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = exz * (double)nx + eyz * (double)ny + ezz * (double)nz;
        }
      }
    }
    int m = nxy * (nz - 1) + nx * (ny - 1) + (nx - 1);
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        double sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        b_[ib_[m][is[mm]]][nn] += sum;
      }
    }
	
    return;

}

void ThermalStrain::bgrad (int nx,
                           int ny,
                           int nz,
                           int ns,
                           double exx,
                           double eyy,
                           double ezz,
                           double exz,
                           double eyz,
                           double exy) 
{
    int is[8];
    int nxy = nx * ny;
    double delta[8][3];
	
    ///
    /// exx, eyy, ezz, exz, eyz, exy are the artificial macrostrains used to       
    /// get the gradient terms (appropriate combinations of 1's and 0's).          

    ///
    /// Set up vector for linear term.                                             
    ///

    is[0] = 26;    
    is[1] = 2;
    is[2] = 1;   
    is[3] = 0;   
    is[4] = 25;   
    is[5] = 18;    
    is[6] = 17;   
    is[7] = 16; 
	
    // x = nx - 1, face                                                         
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;	
          delta[i8][2] = exz * (double)nx;		  
        }
      }
    }
    for (int j = 0; j < (ny - 1); j++) {
      for (int k = 0; k < (nz - 1); k++) {
        int m = nxy * k + nx * j + (nx -1);
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
            if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
            if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
            if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
            if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
            if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }

    // y = ny - 1, face                                                           
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      for (int k = 0; k < (nz - 1); k++) {
        int m = nxy * k + nx * (ny - 1) + i;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
            if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
            if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
            if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
            if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
            if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }	
	
    // z = nz - 1, face                                                           
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      for (int j = 0; j < (ny - 1); j++) {
        int m = nxy * (nz - 1) + nx * j + i;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            double sum = 0.0;
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
              }
            }
            if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
            if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
            if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
            if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
            if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
            if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }

    // x = nx - 1, y = ny - 1, edge                                               
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 5)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;	
          delta[i8][2] = exz * (double)nx;	
        }
        if ((i8 == 3) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if ((i8 == 2) || (i8 == 6)) {
          delta[i8][0] = exy * (double)ny + exx * (double)nx;
          delta[i8][1] = eyy * (double)ny + exy * (double)nx;	
          delta[i8][2] = eyz * (double)ny + exz * (double)nx;
        }
      }
    }
    for (int k = 0; k < (nz - 1); k++) {
      int m = nxy * k + nx * (ny - 1) + (nx -1);
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    // x = nx - 1, z = nz - 1, edge                                           
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2)) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;
          delta[i8][2] = exz * (double)nx;		  
        }
        if ((i8 == 4) || (i8 == 7)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if ((i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exz * (double)nz + exx * (double)nx;
          delta[i8][1] = eyz * (double)nz + exy * (double)nx;
          delta[i8][2] = ezz * (double)nz + exz * (double)nx;
        }
      }
    }
    for (int j = 0; j < (ny - 1); j++) {
      int m = nxy * (nz - 1) + nx * j + (nx - 1);
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    // y = ny - 1, z = nz - 1, edge                                               
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5)) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if ((i8 == 2) || (i8 == 3)) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if ((i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy * (double)ny + exz * (double)nz;
          delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
        }
      }
    }
    for (int i = 0; i < (nx - 1); i++) {
      int m = nxy * (nz - 1) + nx * (ny - 1) + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          double sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    // x = nx - 1, y = ny - 1, z = nz - 1, corner                             
    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if (i8 == 1) {
          delta[i8][0] = exx * (double)nx;
          delta[i8][1] = exy * (double)nx;
          delta[i8][2] = exz * (double)nx;
        }
        if (i8 == 3) {
          delta[i8][0] = exy * (double)ny;
          delta[i8][1] = eyy * (double)ny;
          delta[i8][2] = eyz * (double)ny;
        }
        if (i8 == 4) {
          delta[i8][0] = exz * (double)nz;
          delta[i8][1] = eyz * (double)nz;
          delta[i8][2] = ezz * (double)nz;
        }
        if (i8 == 7) {
          delta[i8][0] = exy * (double)ny + exz * (double)nz; 
          delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
        }
        if (i8 == 5) {
          delta[i8][0] = exx * (double)nx + exz * (double)nz;
          delta[i8][1] = exy * (double)nx + eyz * (double)nz;
          delta[i8][2] = exz * (double)nx + ezz * (double)nz;
        }
        if (i8 == 2) {
          delta[i8][0] = exx * (double)nx + exy * (double)ny;
          delta[i8][1] = exy * (double)nx + eyy * (double)ny;
          delta[i8][2] = exz * (double)nx + eyz * (double)ny;
        }
        if (i8 == 6) {
          delta[i8][0] = exx * (double)nx + exy * (double)ny + exz * (double)nz;
          delta[i8][1] = exy * (double)nx + eyy * (double)ny + eyz * (double)nz;
          delta[i8][2] = exz * (double)nx + eyz * (double)ny + ezz * (double)nz;
        }
      }
    }
    int m = nxy * (nz - 1) + nx * (ny - 1) + (nx - 1);
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        double sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        if (exx == 1.0) b0_[ib_[m][is[mm]]][nn] += sum;
        if (eyy == 1.0) b1_[ib_[m][is[mm]]][nn] += sum;
        if (ezz == 1.0) b2_[ib_[m][is[mm]]][nn] += sum;
        if (exz == 1.0) b3_[ib_[m][is[mm]]][nn] += sum;
        if (eyz == 1.0) b4_[ib_[m][is[mm]]][nn] += sum;
        if (exy == 1.0) b5_[ib_[m][is[mm]]][nn] += sum;
      }
    }

    return;

}
	 
void ThermalStrain::constfunc (int ns,
                               int nx,
                               int ny,
                               int nz) 
{
    double delta[8][3];
    double pp[6][6],s[6][6];
    double econ;
	
    ///
    /// Routine to set up 6 x 6 matrix for energy term involving macrostrains     
    /// only, pulled out of femat.                                                 
    ///
    /// Idea is to evaluate the quadratic term in the macrostrains repeatedly for   
    /// choices of strain like exx = 1, exy = 1, all others = 0, build up 21       
    /// choices, then recombine to get matrix elements by themselves.              
    ///

    int nxy = nx * ny;
 
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        s[i][j] = 0.0;
        pp[i][j] = 0.0;
      }
    }
	
    for (int ii = 0; ii < 6; ii++) {
      for (int jj = ii; jj < 6; jj++) {
        econ = 0.0;
        double exx,eyy,ezz;
        double exz,eyz,exy;
        exx = eyy = ezz = 0.0;
        exz = eyz = exy = 0.0;
        if ((ii == 0) && (jj == 0)) exx = 1.0;
        if ((ii == 1) && (jj == 1)) eyy = 1.0;
        if ((ii == 2) && (jj == 2)) ezz = 1.0;
        if ((ii == 3) && (jj == 3)) exz = 1.0;
        if ((ii == 4) && (jj == 4)) eyz = 1.0;
        if ((ii == 5) && (jj == 5)) exy = 1.0;
        if ((ii == 0) && (jj == 1)) {
          exx = 1.0;
          eyy = 1.0;
        }
        if ((ii == 0) && (jj == 2)) {
          exx = 1.0;
          ezz = 1.0;
        }
        if ((ii == 0) && (jj == 3)) {
          exx = 1.0;
          exz = 1.0;
        }
        if ((ii == 0) && (jj == 4)) {
          exx = 1.0;
          eyz = 1.0;
        }
        if ((ii == 0) && (jj == 5)) {
          exx = 1.0;
          exy = 1.0;
        }
        if ((ii == 1) && (jj == 2)) {
          eyy = 1.0;
          ezz = 1.0;
        }
        if ((ii == 1) && (jj == 3)) {
          eyy = 1.0;
          exz = 1.0;
        }
        if ((ii == 1) && (jj == 4)) {
          eyy = 1.0;
          eyz = 1.0;
        }
        if ((ii == 1) && (jj == 5)) {
          eyy = 1.0;
          exy = 1.0;
        }
        if ((ii == 2) && (jj == 3)) {
          ezz = 1.0;
          exz = 1.0;
        }
        if ((ii == 2) && (jj == 4)) {
          ezz = 1.0;
          eyz = 1.0;
        }
        if ((ii == 2) && (jj == 5)) {
          ezz = 1.0;
          exy = 1.0;
        }
        if ((ii == 3) && (jj == 4)) {
          exz = 1.0;
          eyz = 1.0;
        }
        if ((ii == 3) && (jj == 5)) {
          exz = 1.0;
          exy = 1.0;
        }
        if ((ii == 4) && (jj == 5)) {
          eyz = 1.0;
          exy = 1.0;
        }
        // x = nx - 1, face                                             
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
              delta[i8][0] = exx * (double)nx;
              delta[i8][1] = exy * (double)nx;
              delta[i8][2] = exz * (double)nx;
            }
          }
        }
        for (int j = 0; j < (ny - 1); j++) {
          for (int k = 0; k < (nz - 1); k++) {
            int m = nxy * k + nx * j + (nx - 1);
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                for (int m3 = 0; m3 < 3; m3++) {
                  for (int m8 = 0; m8 < 8; m8++) {
                    econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                  }
                }
              }
            }
          }
        }

        // y = ny - 1, face                                                    
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exy * (double)ny;
              delta[i8][1] = eyy * (double)ny;
              delta[i8][2] = eyz * (double)ny;
            }
          }
        }
        for (int i = 0; i < (nx - 1); i++) {
          for (int k = 0; k < (nz - 1); k++) {
            int m = nxy * k + nx * (ny - 1) + i;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                for (int m3 = 0; m3 < 3; m3++) {
                  for (int m8 = 0; m8 < 8; m8++) {
                    econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                  }
                }
              }
            }
          }
        }
	
        // z = nz - 1, face                                        
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exz * (double)nz;
              delta[i8][1] = eyz * (double)nz;
              delta[i8][2] = ezz * (double)nz;
            }
          }
        }
        for (int i = 0; i < (nx - 1); i++) {
          for (int j = 0; j < (ny - 1); j++) {
            int m = nxy * (nz - 1) + nx * j + i;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                for (int m3 = 0; m3 < 3; m3++) {
                  for (int m8 = 0; m8 < 8; m8++) {
                    econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                  }
                }
              }
            }
          }
        }

        // x = nx - 1, y = ny - 1, edge                                    
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 5)) {
              delta[i8][0] = exx * (double)nx;
              delta[i8][1] = exy * (double)nx;	
              delta[i8][2] = exz * (double)nx;	
            }
            if ((i8 == 3) || (i8 == 7)) {
              delta[i8][0] = exy * (double)ny;
              delta[i8][1] = eyy * (double)ny;
              delta[i8][2] = eyz * (double)ny;
            }
            if ((i8 == 2) || (i8 == 6)) {
              delta[i8][0] = exy * (double)ny + exx * (double)nx;
              delta[i8][1] = eyy * (double)ny + exy * (double)nx;	
              delta[i8][2] = eyz * (double)ny + exz * (double)nx;
            }
          }
        }
        for (int k = 0; k < (nz - 1); k++) {
          int m = nxy * k + nx * (ny - 1) + (nx -1);
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                }
              }
            }
          }
        }

        // x = nx - 1, z = nz - 1, edge 
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 2)) {
              delta[i8][0] = exx * (double)nx;
              delta[i8][1] = exy * (double)nx;
              delta[i8][2] = exz * (double)nx;		  
            }
            if ((i8 == 4) || (i8 == 7)) {
              delta[i8][0] = exz * (double)nz;
              delta[i8][1] = eyz * (double)nz;
              delta[i8][2] = ezz * (double)nz;
            }
            if ((i8 == 5) || (i8 == 6)) {
              delta[i8][0] = exz * (double)nz + exx * (double)nx;
              delta[i8][1] = eyz * (double)nz + exy * (double)nx;
              delta[i8][2] = ezz * (double)nz + exz * (double)nx;
            }
          }
        }
        for (int j = 0; j < (ny - 1); j++) {
          int m = nxy * (nz - 1) + nx * j + (nx - 1);
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                }
              }
            }
          }
        }

        // y = ny - 1, z = nz - 1, edge                          
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 4) || (i8 == 5)) {
              delta[i8][0] = exz * (double)nz;
              delta[i8][1] = eyz * (double)nz;
              delta[i8][2] = ezz * (double)nz;
            }
            if ((i8 == 2) || (i8 == 3)) {
              delta[i8][0] = exy * (double)ny;
              delta[i8][1] = eyy * (double)ny;
              delta[i8][2] = eyz * (double)ny;
            }
            if ((i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exy * (double)ny + exz * (double)nz;
              delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
              delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
            }
          }
        }
        for (int i = 0; i < (nx - 1); i++) {
          int m = nxy * (nz - 1) + nx * (ny - 1) + i;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
                }
              }
            }
          }
        }

        // x = nx - 1, y = ny - 1, z = nz - 1, corner                            
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if (i8 == 1) {
              delta[i8][0] = exx * (double)nx;
              delta[i8][1] = exy * (double)nx;
              delta[i8][2] = exz * (double)nx;
            }
            if (i8 == 3) {
              delta[i8][0] = exy * (double)ny;
              delta[i8][1] = eyy * (double)ny;
              delta[i8][2] = eyz * (double)ny;
            }
            if (i8 == 4) {
              delta[i8][0] = exz * (double)nz;
              delta[i8][1] = eyz * (double)nz;
              delta[i8][2] = ezz * (double)nz;
            }
            if (i8 == 7) {
              delta[i8][0] = exy * (double)ny + exz * (double)nz; 
              delta[i8][1] = eyy * (double)ny + eyz * (double)nz;
              delta[i8][2] = eyz * (double)ny + ezz * (double)nz;
            }
            if (i8 == 5) {
              delta[i8][0] = exx * (double)nx + exz * (double)nz;
              delta[i8][1] = exy * (double)nx + eyz * (double)nz;
              delta[i8][2] = exz * (double)nx + ezz * (double)nz;
            }
            if (i8 == 2) {
              delta[i8][0] = exx * (double)nx + exy * (double)ny;
              delta[i8][1] = exy * (double)nx + eyy * (double)ny;
              delta[i8][2] = exz * (double)nx + eyz * (double)ny;
            }
            if (i8 == 6) {
              delta[i8][0] = exx * (double)nx + exy * (double)ny + exz * (double)nz;
              delta[i8][1] = exy * (double)nx + eyy * (double)ny + eyz * (double)nz;
              delta[i8][2] = exz * (double)nx + eyz * (double)ny + ezz * (double)nz;
            }
          }
        }
        int m = nxy * (nz - 1) + nx * (ny - 1) + (nx - 1);
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] * delta[mm][nn];
              }
            }
          }
        }
        pp[ii][jj] = econ * 2;
      }
    }

    for (int i = 0; i < 6; i++) {
      for (int j = i; j < 6; j++) {
        if (i == j) {
          s[i][j] = pp[i][j];
        } else {
          s[i][j] = pp[i][j] - pp[i][i] - pp[j][j];
        }
      }
    }
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        pp[i][j] = 0.5 * (s[i][j] + s[j][i]);
      }
    }

    ///
    /// Map pp[i][j] into zcon_[2][3][2][3]                                   
    ///

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int mi = 0; mi < 3; mi++) {
          for (int mj = 0; mj < 3; mj++) {
            int ii, jj;
            if (i == 0) ii = i + mi;
            if (i == 1) ii = i + mi + 2;
            if (j == 0) jj = j + mj;
            if (j == 1) jj = j + mj + 2;
            zcon_[i][mi][j][mj] = pp[ii][jj];
            cout << "i, mi, j, mj, zcon_[" << i << "][" 
                 << mi << "][" << j << "][" << mj << "]"
                 << i << " " << mi << " " << j << " "
                 << mj << " " << zcon_[i][mi][j][mj]
                 << endl;
          }
        }
      }
    }

    return;
}


double ThermalStrain::energy (int nx,
                            int ny,
                            int nz,
                            int ns) 
{
    int nss = ns + 2;
    double utot;
	
    ///
    /// Do global matrix multiply via small stiffness matrices, Ah_ = A * h_       
    /// The long statement below correctly brings in all the terms from the        
    /// global matrix A using only the small stiffness matrices.                   
    ///

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        gb_[m][m3] = 0.0;
      }
    }
	
    for (int j = 0; j < 3; j++) {
      for (int n = 0; n < 3; n++) {
        for (int m = 0; m < ns; m++) {
          gb_[m][j] += u_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n]
          + dk_[pix_[ib_[m][6]]][1][j][2][n] + dk_[pix_[ib_[m][24]]][4][j][7][n]
          + dk_[pix_[ib_[m][14]]][5][j][6][n]) + u_[ib_[m][1]][n] * (dk_[pix_[ib_[m][26]]][0][j][2][n]
          + dk_[pix_[ib_[m][24]]][4][j][6][n]) + u_[ib_[m][2]][n] * (dk_[pix_[ib_[m][26]]][0][j][1][n]
          + dk_[pix_[ib_[m][4]]][3][j][2][n] + dk_[pix_[ib_[m][12]]][7][j][6][n] 
          + dk_[pix_[ib_[m][24]]][4][j][5][n]) + u_[ib_[m][3]][n] * (dk_[pix_[ib_[m][4]]][3][j][1][n]
          + dk_[pix_[ib_[m][12]]][7][j][5][n]) + u_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] 
          + dk_[pix_[ib_[m][4]]][3][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][5][n] 
          + dk_[pix_[ib_[m][12]]][7][j][4][n]) + u_[ib_[m][5]][n] * (dk_[pix_[ib_[m][5]]][2][j][0][n] 
          + dk_[pix_[ib_[m][13]]][6][j][4][n]) + u_[ib_[m][6]][n] * (dk_[pix_[ib_[m][5]]][2][j][3][n] 
          + dk_[pix_[ib_[m][6]]][1][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][7][n] 
          + dk_[pix_[ib_[m][14]]][5][j][4][n]) + u_[ib_[m][7]][n] * (dk_[pix_[ib_[m][6]]][1][j][3][n] 
          + dk_[pix_[ib_[m][14]]][5][j][7][n]) + u_[ib_[m][8]][n] * (dk_[pix_[ib_[m][24]]][4][j][3][n]
          + dk_[pix_[ib_[m][14]]][5][j][2][n]) + u_[ib_[m][9]][n] * (dk_[pix_[ib_[m][24]]][4][j][2][n]) 
          + u_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] + dk_[pix_[ib_[m][24]]][4][j][1][n]) 
          + u_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n])
          + u_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][1][n]) 
          + u_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) 
          + u_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] + dk_[pix_[ib_[m][14]]][5][j][0][n])
          + u_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) 
          + u_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] + dk_[pix_[ib_[m][6]]][1][j][6][n]) 
          + u_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n])
          + u_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] + dk_[pix_[ib_[m][4]]][3][j][6][n]) 
          + u_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) 
          + u_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] + dk_[pix_[ib_[m][5]]][2][j][5][n])
          + u_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) 
          + u_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] + dk_[pix_[ib_[m][6]]][1][j][4][n]) 
          + u_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n])
          + u_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] + dk_[pix_[ib_[m][12]]][7][j][3][n] 
          + dk_[pix_[ib_[m][14]]][5][j][1][n] + dk_[pix_[ib_[m][24]]][4][j][0][n]) 
          + u_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] + dk_[pix_[ib_[m][4]]][3][j][7][n] 
          + dk_[pix_[ib_[m][26]]][0][j][4][n] + dk_[pix_[ib_[m][6]]][1][j][5][n]) 
          + u_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] + dk_[pix_[ib_[m][6]]][1][j][1][n] 
          + dk_[pix_[ib_[m][5]]][2][j][2][n] + dk_[pix_[ib_[m][4]]][3][j][3][n]
          + dk_[pix_[ib_[m][24]]][4][j][4][n] + dk_[pix_[ib_[m][14]]][5][j][5][n] 
          + dk_[pix_[ib_[m][13]]][6][j][6][n] + dk_[pix_[ib_[m][12]]][7][j][7][n]);
        }
      }
    }
	
    utot = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns; m++) {
        utot += 0.5 * u_[m][m3] * gb_[m][m3] + b_[m][m3] * u_[m][m3];

        ///
        /// This is the gradient of energy with respect to normal displacements.       
        ///

        gb_[m][m3] += b_[m][m3];
      }
    }
	
    ///
    /// Compute "constant" macrostrain energy term.                            
    ///

    C_ = 0.0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int mi = 0; mi < 3; mi++) {
          for (int mj = 0; mj < 3; mj++) {
            C_ += 0.5 * u_[ns+i][mi] * zcon_[i][mi][j][mj] * u_[ns+j][mj];
          }
        }
      }
    }

    utot += C_;

    ///
    /// Add in constant term from thermal energy, Y.
    ///

    utot += Y_;

    ///
    /// Add in linear term in thermal energy.
    ///

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        utot += T_[m][m3] * u_[m][m3];
      }
    }
	
    ///
    /// Compute the gradient with respect to macrostrains;
    /// put in piece from first derivative of zcon quadratic term.
    ///

    for (int i = 0; i < 2; i++) {
      for (int mi = 0; mi < 3; mi++) {
        double sum = 0.0;
        for (int j = 0; j < 2; j++) {
          for (int mj = 0; mj < 3; mj++) {
            sum += zcon_[i][mi][j][mj] * u_[ns+j][mj];
          }
        }
        gb_[ns+i][mi] = sum;
      }
    }

    ///
    /// Add in piece of gradient, for displacements as well as macrostrains,
    /// that come from linear term in thermal energy.
    ///

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        gb_[m][m3] += T_[m][m3];
      }
    }
	
    ///
    /// Generate part that comes from b*u term; do this by calling
    /// `b_` generation with appropriate macrostrain set to 1 to
    /// get that partial derivative. Just use bgrad (taken from femat),
    /// skip `dk_` and `zcon_` part.                                                   
    ///

    for (int ii = 0; ii < 6; ii++) {
      double sum = 0.0;
      if (ii == 0) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b0_[m][m3];
          }
        }
        gb_[ns][0] += sum;
      }
      if (ii == 1) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b1_[m][m3];
          }
        }
        gb_[ns][1] += sum;
      }
      if (ii == 2) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b2_[m][m3];
          }
        }
        gb_[ns][2] += sum;
      }
      if (ii == 3) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b3_[m][m3];
          }
        }
        gb_[ns+1][0] += sum;
      }
      if (ii == 4) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b4_[m][m3];
          }
        }
        gb_[ns+1][1] += sum;
      }
      if (ii == 5) {
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m = 0; m < ns; m++) {
            sum += u_[m][m3] * b5_[m][m3];
          }
        }
        gb_[ns+1][2] += sum;
      }
    }

    return utot;
	
}

long int ThermalStrain::dembx (int ns,
                               double gg,
                               int ldemb,
                               int kkk) 
{
    double lambda,gamma;
    double hAh,gglast;
    long int Lstep;
    int nss = ns + 2;
    
    ///
    /// Initialize the conjugate direction vector on first call to dembx only.     
    /// For calls to dembx after the first, we want to continue using the value    
    /// of h_ determined in the previous call. Of course, if npoints is greater    
    /// than 0, then this initialization step will be run each time a new          
    /// microstructure is used, as kkk will be reset to 0 every time the counter   
    /// micro is increased.                                                        
    ///

    if (kkk == 0) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          h_[m][m3] = gb_[m][m3];
        }
      }
    }
	
    ///
    /// Lstep counts the number of conjugate gradient steps taken                  
    /// in each call to dembx.                                                     
    ///

    Lstep = 0;
	
    for (int ijk = 0; ijk < ldemb; ijk++) {
      Lstep++;
  
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          Ah_[m][m3] = 0.0;
        }
      }

      ///
      /// Do global matrix multiply via small stiffness matrices, Ah_ = A * h_       
      /// The long statement below correctly brings in all the terms from the global 
      /// matrix A using only the small stiffness matrices.                          
      ///

      for (int j = 0; j < 3; j++) {
        for (int n = 0; n < 3; n++) {
          for (int m = 0; m < ns; m++) {
            Ah_[m][j] += h_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n]
               + dk_[pix_[ib_[m][6]]][1][j][2][n] + dk_[pix_[ib_[m][24]]][4][j][7][n]
               + dk_[pix_[ib_[m][14]]][5][j][6][n]) + h_[ib_[m][1]][n] * (dk_[pix_[ib_[m][26]]][0][j][2][n]
               + dk_[pix_[ib_[m][24]]][4][j][6][n]) + h_[ib_[m][2]][n] * (dk_[pix_[ib_[m][26]]][0][j][1][n]
               + dk_[pix_[ib_[m][4]]][3][j][2][n] + dk_[pix_[ib_[m][12]]][7][j][6][n] 
               + dk_[pix_[ib_[m][24]]][4][j][5][n]) + h_[ib_[m][3]][n] * (dk_[pix_[ib_[m][4]]][3][j][1][n]
               + dk_[pix_[ib_[m][12]]][7][j][5][n]) + h_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] 
               + dk_[pix_[ib_[m][4]]][3][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][5][n] 
               + dk_[pix_[ib_[m][12]]][7][j][4][n]) + h_[ib_[m][5]][n] * (dk_[pix_[ib_[m][5]]][2][j][0][n] 
               + dk_[pix_[ib_[m][13]]][6][j][4][n]) + h_[ib_[m][6]][n] * (dk_[pix_[ib_[m][5]]][2][j][3][n] 
               + dk_[pix_[ib_[m][6]]][1][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][7][n] 
               + dk_[pix_[ib_[m][14]]][5][j][4][n]) + h_[ib_[m][7]][n] * (dk_[pix_[ib_[m][6]]][1][j][3][n] 
               + dk_[pix_[ib_[m][14]]][5][j][7][n]) + h_[ib_[m][8]][n] * (dk_[pix_[ib_[m][24]]][4][j][3][n]
               + dk_[pix_[ib_[m][14]]][5][j][2][n]) + h_[ib_[m][9]][n] * (dk_[pix_[ib_[m][24]]][4][j][2][n]) 
               + h_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] + dk_[pix_[ib_[m][24]]][4][j][1][n]) 
               + h_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n]) + 
               + h_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][1][n]) 
               + h_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) 
               + h_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] + dk_[pix_[ib_[m][14]]][5][j][0][n])
               + h_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) 
               + h_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] + dk_[pix_[ib_[m][6]]][1][j][6][n]) 
               + h_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n])
               + h_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] + dk_[pix_[ib_[m][4]]][3][j][6][n]) 
               + h_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) 
               + h_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] + dk_[pix_[ib_[m][5]]][2][j][5][n])
               + h_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) 
               + h_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] + dk_[pix_[ib_[m][6]]][1][j][4][n]) 
               + h_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n])
               + h_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] + dk_[pix_[ib_[m][12]]][7][j][3][n] 
               + dk_[pix_[ib_[m][14]]][5][j][1][n] + dk_[pix_[ib_[m][24]]][4][j][0][n]) 
               + h_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] + dk_[pix_[ib_[m][4]]][3][j][7][n] 
               + dk_[pix_[ib_[m][26]]][0][j][4][n] + dk_[pix_[ib_[m][6]]][1][j][5][n]) 
               + h_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] + dk_[pix_[ib_[m][6]]][1][j][1][n] 
               + dk_[pix_[ib_[m][5]]][2][j][2][n] + dk_[pix_[ib_[m][4]]][3][j][3][n]
               + dk_[pix_[ib_[m][24]]][4][j][4][n] + dk_[pix_[ib_[m][14]]][5][j][5][n] 
               + dk_[pix_[ib_[m][13]]][6][j][6][n] + dk_[pix_[ib_[m][12]]][7][j][7][n]);
          }
        }
      }
	
      ///
      /// The above accurately gives the second derivative matrix with respect to    
      /// nodal displacements, but fails to give the 2nd derivative terms that       
      /// include the macrostrains [du d(strain) and d(strain)d(strain)].            
      /// Use repeated calls to bgrad to generate mixed 2nd derivatives terms,       
      /// plus use zcon_ in order to correct the matrix multiply and correctly bring 
      /// in macrostrain terms (see manual, Sec. 2.4).                               
      ///

      for (int ii = 0; ii < 6; ii++) {

        ///
        /// Fill in terms from matrix multiply right hand sides, 1 to ns_.  
        ///

        for (int m = 0; m < ns; m++) {
          for (int m1 = 0; m1 < 3; m1++) {
            if (ii == 0) Ah_[m][m1] += b0_[m][m1] * h_[ns][0];
            if (ii == 1) Ah_[m][m1] += b1_[m][m1] * h_[ns][1];
            if (ii == 2) Ah_[m][m1] += b2_[m][m1] * h_[ns][2];
            if (ii == 3) Ah_[m][m1] += b3_[m][m1] * h_[ns+1][0];
            if (ii == 4) Ah_[m][m1] += b4_[m][m1] * h_[ns+1][1];
            if (ii == 5) Ah_[m][m1] += b5_[m][m1] * h_[ns+1][2];
          }
        }

        // now do across bottom, 1 to ns_.     
        for (int m = 0; m < ns; m++) {
          if (ii == 0)
            Ah_[ns][0] += b0_[m][0] * h_[m][0] + b0_[m][1] * h_[m][1] + b0_[m][2] * h_[m][2];
          if (ii == 1)
            Ah_[ns][1] += b1_[m][0] * h_[m][0] + b1_[m][1] * h_[m][1] + b1_[m][2] * h_[m][2];
          if (ii == 2)
            Ah_[ns][2] += b2_[m][0] * h_[m][0] + b2_[m][1] * h_[m][1] + b2_[m][2] * h_[m][2];
          if (ii == 3)
            Ah_[ns+1][0] += b3_[m][0] * h_[m][0] + b3_[m][1] * h_[m][1] + b3_[m][2] * h_[m][2];
          if (ii == 4)
            Ah_[ns+1][1] += b4_[m][0] * h_[m][0] + b4_[m][1] * h_[m][1] + b4_[m][2] * h_[m][2];
          if (ii == 5)
            Ah_[ns+1][2] += b5_[m][0] * h_[m][0] + b5_[m][1] * h_[m][1] + b5_[m][2] * h_[m][2];
        }
        // now do right hand corner terms, ns to (ns+1).          
        for (int m = 0; m < 2; m++) {
          for (int m1 = 0; m1 < 3; m1++) {
            if (ii == 0)
              Ah_[ns][0] += zcon_[0][0][m][m1] * h_[ns+m][m1];
            if (ii == 1)
              Ah_[ns][1] += zcon_[0][1][m][m1] * h_[ns+m][m1];
            if (ii == 2)
              Ah_[ns][2] += zcon_[0][2][m][m1] * h_[ns+m][m1];
            if (ii == 3)
              Ah_[ns+1][0] += zcon_[1][0][m][m1] * h_[ns+m][m1];
            if (ii == 4)
              Ah_[ns+1][1] += zcon_[1][1][m][m1] * h_[ns+m][m1];
            if (ii == 5)
              Ah_[ns+1][2] += zcon_[1][2][m][m1] * h_[ns+m][m1];
          }
        }
      }
	  
      hAh = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          hAh += h_[m][m3] * Ah_[m][m3];
        }
      }

      lambda = gg_ / hAh;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          u_[m][m3] -= lambda * h_[m][m3];
          gb_[m][m3] -= lambda * Ah_[m][m3];
        }
      }

      /*
      string outfilename = "displacement.dat";
      ofstream out1(outfilename.c_str(),ios::app);
      out1 << "after one step of relaxation..." << endl;
      for (int k = 0; k < nz_; k++) {
        int m = nx_ * ny_ * k + nx_ * 50 + 50;
        out1 << u_[m][0] << endl;
      }
      out1 << "macrostrain:" << endl;
      out1 << u_[ns_][0] << endl;
      out1 << u_[ns_][1] << endl;
      out1 << u_[ns_][2] << endl;
      out1 << u_[ns_+1][0] << endl;
      out1 << u_[ns_+1][1] << endl;
      out1 << u_[ns_+1][2] << endl;
      */
	  
      gglast = gg_;
      gg_ = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          gg_ += gb_[m][m3] * gb_[m][m3];
        }
      }
      if (gg_ < gtest_) {
        return Lstep;
      } else {
        gamma = gg_ / gglast;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m= 0; m < nss; m++) {
            h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
          }
        }
      }
	  
    }
	
    return Lstep;
}

void ThermalStrain::relax (double time,
                           int kmax)
{
    int ldemb = 150,ltot = 0;
    int iskip;
    double utot;
    long int Lstep;

    iskip = 1;
    femat(nx_,ny_,nz_,ns_,nphase_,iskip);
    utot = energy(nx_,ny_,nz_,ns_);
    
    gg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < (ns_+2); m++) {
        gg_ += gb_[m][m3] * gb_[m][m3];
      }
    }
    cout << "In the global relax, gg_ is: " << gg_ << endl;
    
    /*
    string outfilename = "displacement.dat";
    ofstream out(outfilename.c_str(),ios::app);
    out << "displacement got from previous step:" << endl;
    for (int k = 0; k < nz_; k++) {
      int m = nx_ * ny_ * k + nx_ * 50 + 50;
      out << u_[m][0] << endl;
    }
    out << "macrostrain:" << endl;
    out << u_[ns_][0] << endl;
    out << u_[ns_][1] << endl;
    out << u_[ns_][2] << endl;
    out << u_[ns_+1][0] << endl;
    out << u_[ns_+1][1] << endl;
    out << u_[ns_+1][2] << endl;
    */
	
    for (int kkk = 0; kkk < kmax; kkk++) {
      ///
      /// Call dembx to implement conjugate gradient routine.             
      ///

      cout << "Going into dembx, call no. " << kkk << endl;
      Lstep = dembx(ns_,gg_,ldemb,kkk); 
      ltot += Lstep;

      ///
      /// Call energy to compute energy after dembx call. If gg_ < gtest_,          
      /// this will be the final energy. If gg_ is still larger than gtest_, then    
      /// this will give an intermediate energy with which to check how the          
      /// relaxation process is coming along. The call to energy does not change the 
      /// gradient or the value of gg_.                                              
      /// Need to first call femat to update the vector b_, as the value of the      
      /// components of b_ depend on the macrostrains.                               
      ///

      iskip = 1;
      femat(nx_,ny_,nz_,ns_,nphase_,iskip);
      utot = energy(nx_,ny_,nz_,ns_);
      cout << " energy = " << utot << " gg_ = " << gg_ 
           << " ltot = " << ltot << endl;
      cout.flush();
	  
      ///
      /// If relaxation process is finished, jump out of loop.                       
      ///

      if (gg_ >= gtest_) {

        ///
        /// Output stresses, strains, and macrostrains as an additional aid in         
        /// judging how well the relaxation process is proceeding.                     
        ///

        stress(nx_,ny_,nz_,ns_);
        cout << " stresses: xx, yy, zz, xz, yz, xy " << strxx_ << " "
             << stryy_ << " " << strzz_ << " " << strxz_ << " "
             << stryz_ << " " << strxy_ << endl;
        cout << " strains: xx, yy, zz, xz, yz, xy " << sxx_ << " "
             << syy_ << " " << szz_ << " " << sxz_ << " "
             << syz_ << " " << sxy_ << endl;
        
        cout << " macrostrains in same order "
             << u_[ns_][0] << " " << u_[ns_][1] << " " << u_[ns_][2] << " "
             << u_[ns_+1][0] << " " << u_[ns_+1][1] << " " << u_[ns_+1][2] << " "
             << endl;
        cout << "avg = " << (u_[ns_][0] + u_[ns_][1] + u_[ns_][2]) / 3.0 << endl;
      } else {
        break; 
      }
    }
    stress(nx_,ny_,nz_,ns_);
    cout << " stresses: xx, yy, zz, xz, yz, xy " << strxx_ << " "
         << stryy_ << " " << strzz_ << " " << strxz_ << " "
         << stryz_ << " " << strxy_ << endl;
    cout << " strains: xx, yy, zz, xz, yz, xy " << sxx_ << " "
         << syy_ << " " << szz_ << " " << sxz_ << " "
         << syz_ << " " << sxy_ << endl;
        
    cout << " macrostrains in same order "
         << u_[ns_][0] << " " << u_[ns_][1] << " " << u_[ns_][2] << " "
         << u_[ns_+1][0] << " " << u_[ns_+1][1] << " " << u_[ns_+1][2] << " "
         << endl;
    cout << "avg = " << (u_[ns_][0] + u_[ns_][1] + u_[ns_][2]) / 3.0 << endl;
    ostringstream ostr2;
    ostr2 << ltot;
    string step(ostr2.str());
    /*
    writeStress(time,step,0);
    writeStrainEngy(time,step); 
    */

    ofstream outstr("macrstrains.dat");
    outstr << u_[ns_][0];
    outstr.close();
   	
    return;
}

void ThermalStrain::stress (int nx,
                            int ny,
                            int nz,
                            int ns) 
{
    double dndx[8],dndy[8],dndz[8];
    double es[6][8][3];
    double uu[8][3];
    double str11,str22,str33,str13,str23,str12;
    double s11,s22,s33,s13,s23,s12;
    double exx,eyy,ezz,exz,eyz,exy;
    int nxy = nx * ny;
    exx = u_[ns][0];
    eyy = u_[ns][1];
    ezz = u_[ns][2];
    exz = u_[ns+1][0];
    eyz = u_[ns+1][1];
    exy = u_[ns+1][2];
	
    ///
    /// Set up single element strain matrix.
    ///

    dndx[0] = (-0.25);
    dndx[1] = 0.25;
    dndx[2] = 0.25;
    dndx[3] = (-0.25);
    dndx[4] = (-0.25);
    dndx[5] = 0.25;
    dndx[6] = 0.25;
    dndx[7] = (-0.25);
    dndy[0] = (-0.25);
    dndy[1] = (-0.25);
    dndy[2] = 0.25;
    dndy[3] = 0.25;
    dndy[4] = (-0.25);
    dndy[5] = (-0.25);
    dndy[6] = 0.25;
    dndy[7] = 0.25;
    dndz[0] = (-0.25);
    dndz[1] = (-0.25);
    dndz[2] = (-0.25);
    dndz[3] = (-0.25);
    dndz[4] = 0.25;
    dndz[5] = 0.25;
    dndz[6] = 0.25;
    dndz[7] = 0.25;
	
    ///
    /// Build average strain matrix, follows code in femat, 
    /// but for average strain over an element, not the strain at a point.
    ///

    for (int n1 = 0; n1 < 6; n1++) {
      for (int n2 = 0; n2 < 8; n2++) {
        for (int n3 = 0; n3 < 3; n3++) {
          es[n1][n2][n3] = 0.0;
        }
      }
    }
    for (int n = 0; n < 8; n++) {
      es[0][n][0] = dndx[n];
      es[1][n][1] = dndy[n];
      es[2][n][2] = dndz[n];
      es[3][n][0] = dndz[n];
      es[3][n][2] = dndx[n];
      es[4][n][1] = dndz[n];
      es[4][n][2] = dndy[n];
      es[5][n][0] = dndy[n];
      es[5][n][1] = dndx[n];
    }

    ///
    /// Compute average stresses and strains in each pixel.
    ///

    sxx_ = 0.0;
    syy_ = 0.0;
    szz_ = 0.0;
    sxz_ = 0.0;
    syz_ = 0.0;
    sxy_ = 0.0;
    strxx_ = 0.0;
    stryy_ = 0.0;
    strzz_ = 0.0;
    strxz_ = 0.0;
    stryz_ = 0.0;
    strxy_ = 0.0;

    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          int m = nxy * k + nx * j + i;

          ///
          /// Load in elements of 8-vector using periodic boundary conditions. 
          ///

          for (int mm = 0; mm < 3; mm++) {
            uu[0][mm] = u_[m][mm];
            uu[1][mm] = u_[ib_[m][2]][mm];
            uu[2][mm] = u_[ib_[m][1]][mm];
            uu[3][mm] = u_[ib_[m][0]][mm];
            uu[4][mm] = u_[ib_[m][25]][mm];
            uu[5][mm] = u_[ib_[m][18]][mm];
            uu[6][mm] = u_[ib_[m][17]][mm];
            uu[7][mm] = u_[ib_[m][16]][mm];
          }
 
          ///
          /// Correct for periodic boundary conditions, some displacements are wrong 
          /// for an element on a periodic boundary. Since they come from an opposite 
          /// face, need to put in applied strain to correct them. 
          ///

          if (i == (nx - 1)) {
            uu[1][0] += exx * nx;
            uu[1][1] += exy * nx;
            uu[1][2] += exz * nx;
            uu[2][0] += exx * nx;
            uu[2][1] += exy * nx;
            uu[2][2] += exz * nx;
            uu[5][0] += exx * nx;
            uu[5][1] += exy * nx;
            uu[5][2] += exz * nx;
            uu[6][0] += exx * nx;
            uu[6][1] += exy * nx;
            uu[6][2] += exz * nx;
          }
          if (j == (ny - 1)) {
            uu[2][0] += exy * ny;
            uu[2][1] += eyy * ny;
            uu[2][2] += eyz * ny;
            uu[3][0] += exy * ny;
            uu[3][1] += eyy * ny;
            uu[3][2] += eyz * ny;
            uu[6][0] += exy * ny;
            uu[6][1] += eyy * ny;
            uu[6][2] += eyz * ny;
            uu[7][0] += exy * ny;
            uu[7][1] += eyy * ny;
            uu[7][2] += eyz * ny;
          }
          if (k == (nz - 1)) {
            uu[4][0] += exz * nz;
            uu[4][1] += eyz * nz;
            uu[4][2] += ezz * nz;
            uu[5][0] += exz * nz;
            uu[5][1] += eyz * nz;
            uu[5][2] += ezz * nz;
            uu[6][0] += exz * nz;
            uu[6][1] += eyz * nz;
            uu[6][2] += ezz * nz;
            uu[7][0] += exz * nz;
            uu[7][1] += eyz * nz;
            uu[7][2] += ezz * nz;
          }

          ///
          /// Stresses and strains in an element. 
          ///

          str11 = 0.0;
          str22 = 0.0;
          str33 = 0.0;
          str13 = 0.0;
          str23 = 0.0;
          str12 = 0.0;
          s11 = 0.0;
          s22 = 0.0;
          s33 = 0.0;
          s13 = 0.0;
          s23 = 0.0;
          s12 = 0.0;

          ///
          /// Compute average stress and strain tensor in each element
          /// First put thermal strain-induced stresses into stress tensor       
          ///

          for (int n = 0; n < 6; n++) {
            str11 -= cmod_[pix_[m]][0][n] * eigen_[m][n];
            str22 -= cmod_[pix_[m]][1][n] * eigen_[m][n];
            str33 -= cmod_[pix_[m]][2][n] * eigen_[m][n];
            str13 -= cmod_[pix_[m]][3][n] * eigen_[m][n];
            str23 -= cmod_[pix_[m]][4][n] * eigen_[m][n];
            str12 -= cmod_[pix_[m]][5][n] * eigen_[m][n];
          }
          for (int n3 = 0; n3 < 3; n3++) {
            for (int n8 = 0; n8 < 8; n8++) {
              s11 += es[0][n8][n3] * uu[n8][n3];
              s22 += es[1][n8][n3] * uu[n8][n3];
              s33 += es[2][n8][n3] * uu[n8][n3];
              s13 += es[3][n8][n3] * uu[n8][n3];
              s23 += es[4][n8][n3] * uu[n8][n3];
              s12 += es[5][n8][n3] * uu[n8][n3];
              for (int n = 0; n < 6; n++) {

                ///
                /// Compute stresses in each element that include both non-thermal
                /// and thermal strains.
                ///

                str11 += cmod_[pix_[m]][0][n] * es[n][n8][n3] * uu[n8][n3];
                str22 += cmod_[pix_[m]][1][n] * es[n][n8][n3] * uu[n8][n3];
                str33 += cmod_[pix_[m]][2][n] * es[n][n8][n3] * uu[n8][n3];
                str13 += cmod_[pix_[m]][3][n] * es[n][n8][n3] * uu[n8][n3];
                str23 += cmod_[pix_[m]][4][n] * es[n][n8][n3] * uu[n8][n3];
                str12 += cmod_[pix_[m]][5][n] * es[n][n8][n3] * uu[n8][n3];
              }
            }
          }
          elestress_[m][0] = str11;
          elestress_[m][1] = str22;
          elestress_[m][2] = str33;
          elestress_[m][3] = str13;
          elestress_[m][4] = str23;
          elestress_[m][5] = str12;
          elestrain_[m][0] = s11;
          elestrain_[m][1] = s22;
          elestrain_[m][2] = s33;
          elestrain_[m][3] = s13;
          elestrain_[m][4] = s23;
          elestrain_[m][5] = s12;
          
          ///
          /// Compute the strain energy for each element.
          ///

          for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
              strainengy_[m] = 0.5 * elestrain_[m][i] * cmod_[pix_[m]][i][j]
                             * elestrain_[m][j];
            }
          }
          
          ///
          /// Calculate the strain energy for each GEM dependent component,
          /// which will be called in GEMS.
          ///

          getAvgStrainengy();

          ///
          /// Sum local stresses and strains into global stresses and strains.         
          ///

          strxx_ += str11;
          stryy_ += str22;
          strzz_ += str33;
          strxz_ += str13;
          stryz_ += str23;
          strxy_ += str12;
          sxx_ += s11;
          syy_ += s22;
          szz_ += s33;
          sxz_ += s13;
          syz_ += s23;
          sxy_ += s12;
        }
      }
    }
	
    ///
    /// Volume average global stresses and strain.                
    ///

    strxx_ = strxx_ / float(ns);
    stryy_ = stryy_ / float(ns);
    strzz_ = strzz_ / float(ns);
    strxz_ = strxz_ / float(ns);
    stryz_ = stryz_ / float(ns);
    strxy_ = strxy_ / float(ns);

    sxx_ = sxx_ / float(ns);
    syy_ = syy_ / float(ns);
    szz_ = szz_ / float(ns);
    sxz_ = sxz_ / float(ns);
    syz_ = syz_ / float(ns);
    sxy_ = sxy_ / float(ns);

    return;
	
}

void ThermalStrain::Calc (double time,
                          string fname,
                          double exx,
                          double eyy,
                          double ezz,
                          double exz,
                          double eyz,
                          double exy) 
{
    int kmax = 60;
    int iskip;
    double utot;
    
    ///
    /// Read in a microstructure in subroutine ppixel, and set up pix_[m] with     
    /// the appropriate phase assignments.                                         
    ///

    ppixel(fname,nphase_);

    ///
    /// Count and output the volume fractions of the different phases.             
    ///

    assig(ns_,nphase_);

    /*
    for (int i = 0; i < nphase_; i++) {
        cout << " Phase " << i << " bulk = " << phasemod_[i][0] << " shear = " 
             << phasemod_[i][1] << endl;
    }
    */
  
    for (int i = 0; i < nphase_; i++) {
      cout << " Volume fraction of phase " << i << "  is " << prob_[i] << endl;
    }
	 
    ///
    /// Output thermal strains for each phase.    
    ///

    /*
    cout << "Thermal Strains" << endl;
    for (int i = 0; i < ns_; i++) {
      bool flag = false;
      for (int j = 0; j < 6; j++) {
        if (eigen_[i][j] != 0.0)
          flag = true;
      }
      if (flag)
        cout << "eigen of site[" << i << "]: " << eigen_[i][0] << "   "
             << eigen_[i][1] << "   " << eigen_[i][2] << "   "
             << eigen_[i][3] << "   " << eigen_[i][4] << "   "
             << eigen_[i][5] << "   " << "phase is: "
             << pix_[i] << endl;
    }
    */
  
    if (isfirst_) {
      cout << " This is the first time Calc to be called, so u_ should be initialized."
           << endl;

      ///
      /// @note Set applied strains. Actual shear strain applied in is `exy_`, `exz_`,
      /// and `eyz_` as given in the statements below. The engineering shear strain,
      /// by which the shear modulus is usually defined, is twice these values.
      ///

      u_[ns_][0] = exx;
      u_[ns_][1] = eyy;
      u_[ns_][2] = ezz;
      u_[ns_+1][0] = exz;
      u_[ns_+1][1] = eyz;
      u_[ns_+1][2] = exy;
      cout << "Applied engineering strains" << endl;
      cout << " exx,   eyy,   ezz,   exz,   eyz,   exy" << endl;
      cout << exx << "   " << eyy << "   " << ezz << "   "
           << exz << "   " << eyz << "   " << exy << endl;
	  
      ///
      /// Apply homogeneous macroscopic strain as the initial condition             
      /// to displacement variables.                                                
      ///

      for (int k = 0; k < nz_; k++) {
        for (int j = 0; j < ny_; j++) {
          for (int i = 0; i < nx_; i++) {
            int m = nx_ * ny_ * k + nx_ * j + i;
            double x = (double)i;
            double y = (double)j;
            double z = (double)k;
            u_[m][0] = x * u_[ns_][0] + y * u_[ns_+1][2] + z * u_[ns_+1][0];
            u_[m][1] = x * u_[ns_+1][2] + y * u_[ns_][1] + z * u_[ns_+1][1];
            u_[m][2] = x * u_[ns_+1][0] + y * u_[ns_+1][1] + z * u_[ns_][2];
          }   
        }
      }
      if (time > 0) {
        cout << "read displacement.dat file to initialized u_" << endl;
        ifstream in("displacement.dat");
        if (!in) {
          cout << "hasn't find displacement.dat file." << endl;
        } else {
          double buff;
          for (int k = 0; k < nz_; k++) {
            for (int j = 0; j < ny_; j++) {
              for (int i = 0; i < nx_; i++) {
                int m = nx_ * ny_ * k + nx_ * j + i;
                in >> buff;
                u_[m][0] = buff;
                in >> buff;
                u_[m][1] = buff;
                in >> buff;
                u_[m][2] = buff;
              }
            }
          }
          in >> buff;
          u_[ns_][0] = buff;
          in >> buff;
          u_[ns_][1] = buff;
          in >> buff;
          u_[ns_][2] = buff;
          in >> buff;
          u_[ns_+1][0] = buff;
          in >> buff;
          u_[ns_+1][1] = buff;
          in >> buff;
          u_[ns_+1][2] = buff;
        }
      }
  
      isfirst_ = false;
    }
	  
    ///
    /// Set up the finite element stiffness matrices, the constant, `C_`,             
    /// the vector, `b_`, required for the energy. `b_` and `C_` depend on the           
    /// macrostrains.                                                              
    /// When they are updated, the values of `b_` and `C_` are updated too via calling 
    /// subroutine femat.                                                          
    /// Only compute the thermal strain terms the first time femat is calles,      
    /// (iskip = 0) as they are unaffected by later changes (iskip = 1) in         
    /// displacements and macrostrains.                                            
    /// Compute initial value of gradient `gb_` and `gg_ = gb_ * gb_`.                 
    ///

    iskip = 0;
    femat(nx_,ny_,nz_,ns_,nphase_,iskip); 
    utot = energy(nx_,ny_,nz_,ns_);
    gg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < (ns_+2); m++) {
        gg_ += gb_[m][m3] * gb_[m][m3];
      }
    }
    cout << "energy = " << utot << " gg_ = " << gg_ << endl;
    cout.flush();

    for (int m1 = 0; m1 < 3; m1++) {
      for (int m = 0; m < (ns_+2); m++) {
        b0_[m][m1] = 0.0;
        b1_[m][m1] = 0.0;
        b2_[m][m1] = 0.0;
        b3_[m][m1] = 0.0;
        b4_[m][m1] = 0.0;
        b5_[m][m1] = 0.0;
      }
    }

    ///
    /// Since `b0`, `b1`, `b2`, `b3`, `b4`, `b5` always keep constant,
    /// `bgrad` is only needed to be called for 6 times to get the values of them.
    ///

    for (int ii = 0; ii < 6; ii++) {
      double e11,e22,e33;
      double e13,e23,e12;
      e11 = e22 = e33 = 0.0;
      e13 = e23 = e12 = 0.0;
      if (ii == 0) {
        e11 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
      if (ii == 1) {
        e22 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
      if (ii == 2) {
        e33 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
      if (ii == 3) {
        e13 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
      if (ii == 4) {
        e23 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
      if (ii == 5) {
        e12 = 1.0;
        bgrad(nx_,ny_,nz_,ns_,e11,e22,e33,e13,e23,e12);
      }
    }

    ///
    /// This is the main relaxation loop.
    /// @note `kmax` is the maximum number of times dembx will be called, with       
    /// `ldemb` conjugate gradient steps performed during each call. The total       
    /// number of conjugate gradient steps allowed for a given elastic computation 
    /// is `kmax * ldemb`.
    ///

    for (int count = 0; count < 0; count++) {
      cout << "boxsize_ is: " << boxsize_ << endl;
      for (map<int, vector<int> >:: iterator it = exp_.begin();
           it != exp_.end(); it++) {
        localgg_ = 0.0;
        int halfbox = boxsize_ / 2;
        int expindex = it->first;
        vector<int> expcor = it->second;
        int xlo = expcor[0] - halfbox;
        int xhi = expcor[0] + halfbox;
        int ylo = expcor[1] - halfbox;
        int yhi = expcor[1] + halfbox;
        int zlo = expcor[2] - halfbox;
        int zhi = expcor[2] + halfbox;

        femat(nx_,ny_,nz_,ns_,nphase_,1);
        double tempenergy = 0.0;

        ///
        /// Call energy method to get the global gradient, `gb_`.
        ///

        tempenergy = energy(nx_,ny_,nz_,ns_);

        for (int m3 = 0; m3 < 3; m3++) {
          for (int k = zlo; k <= zhi; k++) {
            for (int j = ylo; j <= yhi; j++) {
              for (int i = xlo; i <= xhi; i++) {
                int x,y,z;
                x = i;
                y = j;
                z = k;
                if (x < 0) x += nx_;
                if (x >= nx_) x -= nx_;
                if (y < 0) y += ny_;
                if (y >= ny_) y -= ny_; 
                if (z < 0) z += nz_;
                if (z >= nz_) z -= nz_;
                int m = nx_ * ny_ * z + nx_ * y + x;
                localgg_ += gb_[m][m3] * gb_[m][m3];
              }
            }
          }
        }
        cout << "localgg_ is: " << localgg_ << endl;
        localRelax(boxsize_,expcor[0],expcor[1],expcor[2],expindex);
      }
    }                                                           
    relax(time,kmax);

    string outfilename = "displacement.dat";
    ofstream out(outfilename.c_str());
    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          int m = nx_ * ny_ * k + nx_ * j + i;
          out << u_[m][0] << " " << u_[m][1] << " " << u_[m][2] << endl;
        }
      }
    }
    out << u_[ns_][0] << " " << u_[ns_][1] << " " << u_[ns_][2] << endl;
    out << u_[ns_+1][0] << " " << u_[ns_+1][1] << " " << u_[ns_+1][2] << endl;
    out.close();

    return;  
}

void ThermalStrain::localRelax (int boxsize,
                                int x,
                                int y,
                                int z,
                                int index)
{
    cout << "in the localRelax..." << endl;
    cout << "x, y, z and index are: " << x << " " << y << " " << z << " " << index << endl;
    int localldemb = 60, ltot = 0;
    int iskip;
    double utot;
    long int Lstep;

    int halfbox = (int) (boxsize / 2);
    
    Lstep = localDembx(boxsize,x,y,z,localldemb,0);
    ltot += Lstep;
    cout << " localgg_ = " << localgg_ << " local relaxation steps ltot = "
         << ltot << endl;
    cout.flush();

    return;     
}



int ThermalStrain::localDembx (int boxsize,
                               int x,
                               int y,
                               int z,
                               int localldemb,
                               int kkk) 
{
    int nss = ns_ + 2;
    double lambda,gamma;
    double hAh,gglast;
    long int Lstep;
    int halfbox = (int) (boxsize / 2);
    bool flag = false; // determines whether bgrad function needs to be called or not
  
    int xlo,xhi,ylo,yhi,zlo,zhi;
    xlo = x - halfbox;
    xhi = x + halfbox;
    ylo = y - halfbox;
    yhi = y + halfbox;
    zlo = z - halfbox;
    zhi = z + halfbox;
	
    if (kkk == 0) {
      for (int k = zlo; k <= zhi; k++) {
        for (int j = ylo; j <= yhi; j++) {
          for (int i = xlo; i <= xhi; i++) {
            for (int m3 = 0; m3 < 3; m3++) {
              int x,y,z;
              x = i;
              y = j;
              z = k;
              if (x < 0) {x += nx_; flag = true;}
              if (x >= nx_) {x -= nx_; flag = true;}
              if (y < 0) {y += ny_; flag = true;}
              if (y >= ny_) {y -= ny_; flag = true;}
              if (z < 0) {z += nz_; flag = true;}
              if (z >= nz_) {z -= nz_; flag = true;}
              int m = nx_ * ny_ * z + nx_ * y + x;
              h_[m][m3] = gb_[m][m3];
            }
          }
        }
      }
    }
	
    ///
    /// Lstep counts the number of conjugate gradient steps taken
    /// in each call to dembx.
    ///

    Lstep = 0;
	
    for (int ijk = 0; ijk < localldemb; ijk++) {
      Lstep++;
  
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          Ah_[m][m3] = 0.0;
        }
      }

      ///
      /// Do global matrix multiply via small stiffness matrices, `Ah_ = A * h_`
      /// The long statement below correctly brings in all the terms from the global
      /// matrix A using only the small stiffness matrices.
      ///

      for (int j = 0; j < 3; j++) {
        for (int n = 0; n < 3; n++) {
          for (int kk = zlo; kk <= zhi; kk++) {
            for (int jj = ylo; jj <= yhi; jj++) {
              for (int ii = xlo; ii <= xhi; ii++) {
                int x,y,z;
                x = ii;
                y = jj;
                z = kk;
                if (x < 0) x += nx_;
                if (x >= nx_) x -= nx_;
                if (y < 0) y += ny_;
                if (y >= ny_) y -= ny_;
                if (z < 0) z += nz_;
                if (z >= nz_) z -= nz_;
                int m = nx_ * ny_ * z + nx_ * y + x;
                Ah_[m][j] += h_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n]
                   + dk_[pix_[ib_[m][6]]][1][j][2][n] + dk_[pix_[ib_[m][24]]][4][j][7][n]
                   + dk_[pix_[ib_[m][14]]][5][j][6][n]) + h_[ib_[m][1]][n] *
                    (dk_[pix_[ib_[m][26]]][0][j][2][n]
                   + dk_[pix_[ib_[m][24]]][4][j][6][n]) + h_[ib_[m][2]][n] *
                    (dk_[pix_[ib_[m][26]]][0][j][1][n]
                   + dk_[pix_[ib_[m][4]]][3][j][2][n] + dk_[pix_[ib_[m][12]]][7][j][6][n] 
                   + dk_[pix_[ib_[m][24]]][4][j][5][n]) + h_[ib_[m][3]][n] *
                    (dk_[pix_[ib_[m][4]]][3][j][1][n] + dk_[pix_[ib_[m][12]]][7][j][5][n])
                    + h_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] 
                   + dk_[pix_[ib_[m][4]]][3][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][5][n] 
                   + dk_[pix_[ib_[m][12]]][7][j][4][n]) + h_[ib_[m][5]][n] *
                    (dk_[pix_[ib_[m][5]]][2][j][0][n] 
                   + dk_[pix_[ib_[m][13]]][6][j][4][n]) + h_[ib_[m][6]][n] *
                    (dk_[pix_[ib_[m][5]]][2][j][3][n] 
                   + dk_[pix_[ib_[m][6]]][1][j][0][n] + dk_[pix_[ib_[m][13]]][6][j][7][n] 
                   + dk_[pix_[ib_[m][14]]][5][j][4][n]) + h_[ib_[m][7]][n] *
                    (dk_[pix_[ib_[m][6]]][1][j][3][n] 
                   + dk_[pix_[ib_[m][14]]][5][j][7][n]) + h_[ib_[m][8]][n] *
                    (dk_[pix_[ib_[m][24]]][4][j][3][n]
                   + dk_[pix_[ib_[m][14]]][5][j][2][n]) + h_[ib_[m][9]][n] *
                    (dk_[pix_[ib_[m][24]]][4][j][2][n]) 
                   + h_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] +
                           dk_[pix_[ib_[m][24]]][4][j][1][n]) 
                   + h_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n]) +
                   + h_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] +
                           dk_[pix_[ib_[m][13]]][6][j][1][n]) 
                   + h_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) 
                   + h_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] +
                           dk_[pix_[ib_[m][14]]][5][j][0][n])
                   + h_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) 
                   + h_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] +
                           dk_[pix_[ib_[m][6]]][1][j][6][n]) 
                   + h_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n])
                   + h_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] +
                           dk_[pix_[ib_[m][4]]][3][j][6][n]) 
                   + h_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) 
                   + h_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] +
                           dk_[pix_[ib_[m][5]]][2][j][5][n])
                   + h_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) 
                   + h_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] +
                           dk_[pix_[ib_[m][6]]][1][j][4][n]) 
                   + h_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n])
                   + h_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] +
                           dk_[pix_[ib_[m][12]]][7][j][3][n] 
                   + dk_[pix_[ib_[m][14]]][5][j][1][n] + dk_[pix_[ib_[m][24]]][4][j][0][n]) 
                   + h_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] +
                           dk_[pix_[ib_[m][4]]][3][j][7][n] 
                   + dk_[pix_[ib_[m][26]]][0][j][4][n] + dk_[pix_[ib_[m][6]]][1][j][5][n]) 
                   + h_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] +
                           dk_[pix_[ib_[m][6]]][1][j][1][n] 
                   + dk_[pix_[ib_[m][5]]][2][j][2][n] + dk_[pix_[ib_[m][4]]][3][j][3][n]
                   + dk_[pix_[ib_[m][24]]][4][j][4][n] + dk_[pix_[ib_[m][14]]][5][j][5][n] 
                   + dk_[pix_[ib_[m][13]]][6][j][6][n] + dk_[pix_[ib_[m][12]]][7][j][7][n]);
              }
            }
          }
        }
      }
	
      ///
      /// The above accurately gives the second derivative matrix with respect to    
      /// nodal displacements, but fails to give the 2nd derivative terms that       
      /// include the macrostrains [du d(strain) and d(strain)d(strain)].            
      /// Use repeated calls to bgrad to generate mixed 2nd derivatives terms,       
      /// plus use zcon_ in order to correct the matrix multiply and correctly bring 
      /// in macrostrain terms (see manual, Sec. 2.4).                               
      ///

      for (int ii = 0; ii < 6; ii++) {
        if (flag) {

          ///
          /// Fill in terms from matrix multiply right hand sides, 1 to ns_.  
          ///

          for (int m1 = 0; m1 < 3; m1++) {
            for (int k = zlo; k <= zhi; k++) {
              for (int j = ylo; j <= yhi; j++) {
                for (int i = xlo; i <= xhi; i++) {
                  int x,y,z;
                  x = i;
                  y = j;
                  z = k;
                  if (x < 0) x += nx_;
                  if (x >= nx_) x -= nx_;
                  if (y < 0) y += ny_;
                  if (y >= ny_) y -= ny_;
                  if (z < 0) z += nz_;
                  if (z >= nz_) z -= nz_;
                  int m = nx_ * ny_ * z + nx_ * y + x;
                  if (ii == 0) Ah_[m][m1] += b0_[m][m1] * h_[ns_][0];
                  if (ii == 1) Ah_[m][m1] += b1_[m][m1] * h_[ns_][1];
                  if (ii == 2) Ah_[m][m1] += b2_[m][m1] * h_[ns_][2];
                  if (ii == 3) Ah_[m][m1] += b3_[m][m1] * h_[ns_+1][0];
                  if (ii == 4) Ah_[m][m1] += b4_[m][m1] * h_[ns_+1][1];
                  if (ii == 5) Ah_[m][m1] += b5_[m][m1] * h_[ns_+1][2];
                }
              }
            }
          }
          // now do across bottom, 1 to ns_.     
          for (int k = zlo; k <= zhi; k++) {
            for (int j = ylo; j <= yhi; j++) {
              for (int i = xlo; i <= xhi; i++) {
                int x,y,z;
                x = i;
                y = j;
                z = k;
                if (x < 0) x += nx_;
                if (x >= nx_) x -= nx_;
                if (y < 0) y += ny_;
                if (y >= ny_) y -= ny_;
                if (z < 0) z += nz_;
                if (z >= nz_) z -= nz_;
                int m = nx_ * ny_ * z + nx_ * y + x;
                if (ii == 0)
                  Ah_[ns_][0] += b0_[m][0] * h_[m][0]
                                 + b0_[m][1] * h_[m][1]
                                 + b0_[m][2] * h_[m][2];
                if (ii == 1)
                  Ah_[ns_][1] += b1_[m][0] * h_[m][0]
                                 + b1_[m][1] * h_[m][1]
                                 + b1_[m][2] * h_[m][2];
                if (ii == 2)
                  Ah_[ns_][2] += b2_[m][0] * h_[m][0]
                                 + b2_[m][1] * h_[m][1]
                                 + b2_[m][2] * h_[m][2];
                if (ii == 3)
                  Ah_[ns_+1][0] += b3_[m][0] * h_[m][0]
                                   + b3_[m][1] * h_[m][1]
                                   + b3_[m][2] * h_[m][2];
                if (ii == 4)
                  Ah_[ns_+1][1] += b4_[m][0] * h_[m][0]
                                   + b4_[m][1] * h_[m][1]
                                   + b4_[m][2] * h_[m][2];
                if (ii == 5)
                  Ah_[ns_+1][2] += b5_[m][0] * h_[m][0]
                                   + b5_[m][1] * h_[m][1]
                                   + b5_[m][2] * h_[m][2];
              }
            }
          }
        }
        // now do right hand corner terms, ns to (ns+1).          
        for (int m = 0; m < 2; m++) {
          for (int m1 = 0; m1 < 3; m1++) {
            if (ii == 0)
              Ah_[ns_][0] += zcon_[0][0][m][m1] * h_[ns_+m][m1];
            if (ii == 1)
              Ah_[ns_][1] += zcon_[0][1][m][m1] * h_[ns_+m][m1];
            if (ii == 2)
              Ah_[ns_][2] += zcon_[0][2][m][m1] * h_[ns_+m][m1];
            if (ii == 3)
              Ah_[ns_+1][0] += zcon_[1][0][m][m1] * h_[ns_+m][m1];
            if (ii == 4)
              Ah_[ns_+1][1] += zcon_[1][1][m][m1] * h_[ns_+m][m1];
            if (ii == 5)
              Ah_[ns_+1][2] += zcon_[1][2][m][m1] * h_[ns_+m][m1];
          }
        }
      }
	  
      hAh = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int k = zlo; k <= zhi; k++) { 
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              int x,y,z;
              x = i;
              y = j;
              z = k;
              if (x < 0) x += nx_;
              if (x >= nx_) x -= nx_;
              if (y < 0) y += ny_;
              if (y >= ny_) y -= ny_;
              if (z < 0) z += nz_;
              if (z >= nz_) z -= nz_;
              int m = nx_ * ny_ * z + nx_ * y + x;
              hAh += h_[m][m3] * Ah_[m][m3];
            }
          }
        }
      }
	  
      lambda = localgg_ / hAh;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int k = zlo; k <= zhi; k++) {
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              int x,y,z;
              x = i;
              y = j;
              z = k;
              if (x < 0) x += nx_;
              if (x >= nx_) x -= nx_;
              if (y < 0) y += ny_;
              if (y >= ny_) y -= ny_;
              if (z < 0) z += nz_;
              if (z >= nz_) z -= nz_;
              int m = nx_ * ny_ * z + nx_ * y + x;
              u_[m][m3] -= lambda * h_[m][m3];
              gb_[m][m3] -= lambda * Ah_[m][m3];
            }
          }
        }
      }
	  
      gglast = localgg_;
      localgg_ = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int k = zlo; k <= zhi; k++) { 
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              int x,y,z;
              x = i;
              y = j;
              z = k;
              if (x < 0) x += nx_;
              if (x >= nx_) x -= nx_;
              if (y < 0) y += ny_;
              if (y >= ny_) y -= ny_;
              if (z < 0) z += nz_;
              if (z >= nz_) z -= nz_;
              int m = nx_ * ny_ * z + nx_ * y + x;
              localgg_ += gb_[m][m3] * gb_[m][m3];
            }
          }
        }
      }
      cout << "localgg_ is: " << localgg_ << endl;
      if (localgg_ < localgtest_) {
        return Lstep;
      } else {
        gamma = localgg_ / gglast;
        cout << "gamma = " << gamma << endl;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int k = zlo; k <= zhi; k++) { 
            for (int j = ylo; j <= yhi; j++) {
              for (int i = xlo; i <= xhi; i++) {
                int x,y,z;
                x = i;
                y = j;
                z = k;
                if (x < 0) x += nx_;
                if (x >= nx_) x -= nx_;
                if (y < 0) y += ny_;
                if (y >= ny_) y -= ny_;
                if (z < 0) z += nz_;
                if (z >= nz_) z -= nz_;
                int m = nx_ * ny_ * z + nx_ * y + x;
                h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
              }
            }
          }
        }
      }
	  
    }
	
    return Lstep;
}
