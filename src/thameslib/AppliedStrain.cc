/**
@file AppliedStrain.cc
@brief Method definitions for the AppliedStrain derived class.
*/
#include "AppliedStrain.h"
#include <iostream>

using namespace std;

AppliedStrain::AppliedStrain (int nx,
                              int ny,
                              int nz,
                              int dim,
                              int nphase,
                              int npoints) 
    : ElasticModel (nx,ny,nz,dim,nphase,npoints) 
{
    cout << "constructor AppliedStrain in derived class AppliedStrain." << endl;
    exx_ = eyy_ = ezz_ = 0.0;
    exz_ = eyz_ = exy_ = 0.0;
}

void AppliedStrain::femat (int nx,
                           int ny,
                           int nz,
                           int ns,
                           int nphase) 
{ 
    double dndx[8],dndy[8],dndz[8];
    double g[3][3][3];
    double es[6][8][3],delta[8][3];
    int is[8];
    double x,y,z;
    int nxy = nx * ny;

    ///
    /// (User) NOTE: complete elastic modulus matrix is used, so an anisotropic    
    /// matrix could be directly input at any point, since program is written to   
    /// use a general elastic moduli tensor, but is only explicitly implemented    
    /// for isotropic materials.                                                   
    ///

    ///
    /// Initialize stiffness matrices                                              
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
    /// Set up elastic moduli matrices for each kind of element                    
    ///

    ElasModul(phasemod_fname_,nphase_);

    ///
    /// Set up Simpson's integration rule weight vector         
    ///

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          int nm = 0;
          if (i == 1) nm += 1;
          if (j == 1) nm += 1;
          if (k == 1) nm += 1;
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
            x = ((float)(i)) / 2.0;
            y = ((float)(j)) / 2.0;
            z = ((float)(k)) / 2.0;

            ///
            /// dndx means the negative derivative, with respect to x, 
            /// of the shape matrix 
            /// N (see Manual, Sec. 2.2), dndy, and dndz are similar.                      
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
            /// Matrix multiply to determine value at (x,y,z), multiply by proper weight,  
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
    /// Set up vector for linear term, b, and constant term, C, in the elastic     
    /// energy. This is done using the stiffness matrices, and the periodic terms  
    /// in the applied strain that come in at the boundary pixels via the periodic  
    /// boundary conditions and the condition that an applied macroscopic strain   
    /// exists (see Sec. 2.2 in the manual). It is easier to set b up this way     
    /// than to analytically write out all the terms involved.                     
    ///

    ///
    /// Initialize b and C                                                         
    ///

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns; m++) {
        b_[m][m3] = 0.0;
      }
    }
    C_ = 0.0;
	
    ///
    /// For all cases, the correspondence between 0-7 finite element node labels   
    /// and 0-26 neighbor labels is (see Table 4 in manual):                       
    /// 0:ib_[m][26], 1:ib_[m][2],                                                 
    /// 2:ib_[m][1], 3:ib_[m][0],                                                  
    /// 4:ib_[m][25], 5:ib_[m][18],                                                
    /// 6:ib_[m][17], 7:ib_[m][16].                                                  
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
    /// x = nx - 1, face                                                           
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exx_ * (double)nx;
          delta[i8][1] = exy_ * (double)nx;	
          delta[i8][2] = exz_ * (double)nx;		  
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
                C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                    * delta[mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }	

    ///
    /// y = ny - 1, face                                                           
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy_ * (double)ny;
          delta[i8][1] = eyy_ * (double)ny;
          delta[i8][2] = eyz_ * (double)ny;
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
                C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                    * delta[mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }

    ///
    /// z = nz - 1, face                                                 
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exz_ * (double)nz;
          delta[i8][1] = eyz_ * (double)nz;
          delta[i8][2] = ezz_ * (double)nz;
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
                C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                    * delta[mm][nn];
              }
            }
            b_[ib_[m][is[mm]]][nn] += sum;
          }
        }
      }
    }

    ///
    /// x = nx - 1, y = ny - 1, edge                
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 5)) {
          delta[i8][0] = exx_ * (double)nx;
          delta[i8][1] = exy_ * (double)nx;	
          delta[i8][2] = exz_ * (double)nx;	
        }
        if ((i8 == 3) || (i8 == 7)) {
          delta[i8][0] = exy_ * (double)ny;
          delta[i8][1] = eyy_ * (double)ny;
          delta[i8][2] = eyz_ * (double)ny;
        }
        if ((i8 == 2) || (i8 == 6)) {
          delta[i8][0] = exy_ * (double)ny + exx_ * (double)nx;
          delta[i8][1] = eyy_ * (double)ny + exy_ * (double)nx;	
          delta[i8][2] = eyz_ * (double)ny + exz_ * (double)nx;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                 * delta[mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	
    ///
    /// x = nx - 1, z = nz - 1, edge                          
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 1) || (i8 == 2)) {
          delta[i8][0] = exx_ * (double)nx;
          delta[i8][1] = exy_ * (double)nx;
          delta[i8][2] = exz_ * (double)nx;		  
        }
        if ((i8 == 4) || (i8 == 7)) {
          delta[i8][0] = exz_ * (double)nz;
          delta[i8][1] = eyz_ * (double)nz;
          delta[i8][2] = ezz_ * (double)nz;
        }
        if ((i8 == 5) || (i8 == 6)) {
          delta[i8][0] = exz_ * (double)nz + exx_ * (double)nx;
          delta[i8][1] = eyz_ * (double)nz + exy_ * (double)nx;
          delta[i8][2] = ezz_ * (double)nz + exz_ * (double)nx;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                 * delta[mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }

    ///
    /// y = ny - 1, z = nz - 1, edge                                              
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if ((i8 == 4) || (i8 == 5)) {
          delta[i8][0] = exz_ * (double)nz;
          delta[i8][1] = eyz_ * (double)nz;
          delta[i8][2] = ezz_ * (double)nz;
        }
        if ((i8 == 2) || (i8 == 3)) {
          delta[i8][0] = exy_ * (double)ny;
          delta[i8][1] = eyy_ * (double)ny;
          delta[i8][2] = eyz_ * (double)ny;
        }
        if ((i8 == 6) || (i8 == 7)) {
          delta[i8][0] = exy_ * (double)ny + exz_ * (double)nz;
          delta[i8][1] = eyy_ * (double)ny + eyz_ * (double)nz;
          delta[i8][2] = eyz_ * (double)ny + ezz_ * (double)nz;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                 * delta[mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
	  
    ///
    /// x = nx - 1, y = ny - 1, z = nz - 1, corner                                 
    ///

    for (int i3 = 0; i3 < 3; i3++) {
      for (int i8 = 0; i8 < 8; i8++) {
        delta[i8][i3] = 0.0;
        if (i8 == 1) {
          delta[i8][0] = exx_ * (double)nx;
          delta[i8][1] = exy_ * (double)nx;
          delta[i8][2] = exz_ * (double)nx;
        }
        if (i8 == 3) {
          delta[i8][0] = exy_ * (double)ny;
          delta[i8][1] = eyy_ * (double)ny;
          delta[i8][2] = eyz_ * (double)ny;
        }
        if (i8 == 4) {
          delta[i8][0] = exz_ * (double)nz;
          delta[i8][1] = eyz_ * (double)nz;
          delta[i8][2] = ezz_ * (double)nz;
        }
        if (i8 == 7) {
          delta[i8][0] = exy_ * (double)ny + exz_ * (double)nz; 
          delta[i8][1] = eyy_ * (double)ny + eyz_ * (double)nz;
          delta[i8][2] = eyz_ * (double)ny + ezz_ * (double)nz;
        }
        if (i8 == 5) {
          delta[i8][0] = exx_ * (double)nx + exz_ * (double)nz;
          delta[i8][1] = exy_ * (double)nx + eyz_ * (double)nz;
          delta[i8][2] = exz_ * (double)nx + ezz_ * (double)nz;
        }
        if (i8 == 2) {
          delta[i8][0] = exx_ * (double)nx + exy_ * (double)ny;
          delta[i8][1] = exy_ * (double)nx + eyy_ * (double)ny;
          delta[i8][2] = exz_ * (double)nx + eyz_ * (double)ny;
        }
        if (i8 == 6) {
          delta[i8][0] = exx_ * (double)nx + exy_ * (double)ny + exz_ * (double)nz;
          delta[i8][1] = exy_ * (double)nx + eyy_ * (double)ny + eyz_ * (double)nz;
          delta[i8][2] = exz_ * (double)nx + eyz_ * (double)ny + ezz_ * (double)nz;
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
            C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn]
                * delta[mm][nn];
          }
        }
        b_[ib_[m][is[mm]]][nn] += sum;
      }
    }

    return;
}



double AppliedStrain::energy (int nx,
                            int ny,
                            int nz,
                            int ns) 
{
    double utot = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns; m++) {
        gb_[m][m3] = 0.0;
      }
    }

    ///
    /// Do global matrix multiply via small stiffness matrices, gb_ = A * u_      
    /// The long statement below correctly brings in all the terms from the       
    /// global matrix A using only the small stiffness matrices.                  
    ///

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
	
    utot = C_;
    for (int m3 = 0; m3 < 3; m3++){
      for (int m = 0; m < ns; m++) {
        utot += 0.5 * u_[m][m3] * gb_[m][m3] + b_[m][m3] * u_[m][m3];
        gb_[m][m3] += b_[m][m3];
      }
    }

    return utot;
}

long int AppliedStrain::dembx (int ns,
                               double gg,
                               int ldemb,
                               int kkk) 
{
    double lambda,gamma;
    double hAh,gglast;
    long int Lstep;
	
    if (kkk == 0) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          h_[m][m3] = gb_[m][m3];
        }
      }
    }

    ///
    /// Lstep counts the number of conjugate gradient steps taken in each      
    /// call to dembx.                                                             

    Lstep = 0;
	
    for (int ijk = 0; ijk < ldemb; ijk++) {
      Lstep += 1;
	  
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          Ah_[m][m3] = 0.0;
        }
      }

      ///
      /// Do global matrix multiply via small stiffness matrices, Ah_ = A_ * h_.     
      /// The long statement below correctly brings in all the terms from the global 
      /// matrix A_ using only the small stiffness matrices dk_.                     

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
	  
      hAh = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          hAh += h_[m][m3] * Ah_[m][m3];
        }
      }
      lambda = gg_ / hAh;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          u_[m][m3] = u_[m][m3] - lambda * h_[m][m3];
          gb_[m][m3] = gb_[m][m3] - lambda * Ah_[m][m3];
        }
      }
	  
      gglast = gg_;
      gg_ = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          gg_ += gb_[m][m3] * gb_[m][m3];
        }
      }
      if (gg_ < gtest_) return Lstep;
	  
      gamma = gg_ / gglast;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns; m++) {
          h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
        }
      }
	  
    }
    return Lstep;
}

void AppliedStrain::stress (int nx,
                            int ny,
                            int nz,
                            int ns) 
{
    int nxy = nx * ny;
    double dndx[8],dndy[8],dndz[8];
    double es[6][8][3];
    double uu[8][3];
    double str11,str22,str33,str13,str23,str12;
    double s11,s22,s33,s13,s23,s12;
    int m = 0;

    ///
    /// Set up single element strain matrix 
    /// dndx, dndy, and dndz are the components of the average strain matrix 
    /// in a pixel.                                                           
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
    /// Build averaged strain matrix, follows code in femat, but for average    
    /// strain over the pixel, not the strain at a point.                      
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
    /// Compute components of the average stress and strain tensors in each pixel.
    ///

    strxx_ = stryy_ = strzz_ = strxz_ = stryz_ = strxy_ = 0.0;
    sxx_ = syy_ = szz_ = sxz_ = syz_ = sxy_ = 0.0;

    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          m = k * nxy + j * nx + i;

          ///
          /// Load in elements of 8-vector using periodic boundary condition.     
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
          /// for a pixel on periodic boundary. Since they come from an opposite face,    
          /// need to put in applied strain to correct them.                             
          ///

          if (i == (nx - 1)) {
            uu[1][0] = uu[1][0] + exx_ * (double)nx;
            uu[1][1] = uu[1][1] + exy_ * (double)nx;
            uu[1][2] = uu[1][2] + exz_ * (double)nx;
            uu[2][0] = uu[2][0] + exx_ * (double)nx;
            uu[2][1] = uu[2][1] + exy_ * (double)nx;
            uu[2][2] = uu[2][2] + exz_ * (double)nx;
            uu[5][0] = uu[5][0] + exx_ * (double)nx;
            uu[5][1] = uu[5][1] + exy_ * (double)nx;
            uu[5][2] = uu[5][2] + exz_ * (double)nx;
            uu[6][0] = uu[6][0] + exx_ * (double)nx;
            uu[6][1] = uu[6][1] + exy_ * (double)nx;
            uu[6][2] = uu[6][2] + exz_ * (double)nx;
          }
          if (j == (ny - 1)) {
            uu[2][0] = uu[2][0] + exy_ * (double)ny;
            uu[2][1] = uu[2][1] + eyy_ * (double)ny;
            uu[2][2] = uu[2][2] + eyz_ * (double)ny;
            uu[3][0] = uu[3][0] + exy_ * (double)ny;
            uu[3][1] = uu[3][1] + eyy_ * (double)ny;
            uu[3][2] = uu[3][2] + eyz_ * (double)ny;
            uu[6][0] = uu[6][0] + exy_ * (double)ny;
            uu[6][1] = uu[6][1] + eyy_ * (double)ny;
            uu[6][2] = uu[6][2] + eyz_ * (double)ny;
            uu[7][0] = uu[7][0] + exy_ * (double)ny;
            uu[7][1] = uu[7][1] + eyy_ * (double)ny;
            uu[7][2] = uu[7][2] + eyz_ * (double)ny;
          }
          if (k == (nz - 1)) {
            uu[4][0] = uu[4][0] + exz_ * (double)nz;
            uu[4][1] = uu[4][1] + eyz_ * (double)nz;
            uu[4][2] = uu[4][2] + ezz_ * (double)nz;
            uu[5][0] = uu[5][0] + exz_ * (double)nz;
            uu[5][1] = uu[5][1] + eyz_ * (double)nz;
            uu[5][2] = uu[5][2] + ezz_ * (double)nz;
            uu[6][0] = uu[6][0] + exz_ * (double)nz;
            uu[6][1] = uu[6][1] + eyz_ * (double)nz;
            uu[6][2] = uu[6][2] + ezz_ * (double)nz;
            uu[7][0] = uu[7][0] + exz_ * (double)nz;
            uu[7][1] = uu[7][1] + eyz_ * (double)nz;
            uu[7][2] = uu[7][2] + ezz_ * (double)nz;
          }

          ///
          /// Local stresses and strains in a pixel
          ///

          str11 = str22 = str33 = str13 = str23 = str12 = 0.0;
          s11 = s22 = s33 = s13 = s23 = s12 = 0.0;

          for (int n3 = 0; n3 < 3; n3++){
            for (int n8 = 0; n8 < 8; n8++) {
              s11 += es[0][n8][n3] * uu[n8][n3];
              s22 += es[1][n8][n3] * uu[n8][n3];
              s33 += es[2][n8][n3] * uu[n8][n3];
              s13 += es[3][n8][n3] * uu[n8][n3];
              s23 += es[4][n8][n3] * uu[n8][n3];
              s12 += es[5][n8][n3] * uu[n8][n3];
              for ( int n = 0; n < 6; n++) {
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
          /// Compute the strain energy for each element
          ///

          for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
              strainengy_[m] = 0.5 * elestrain_[m][i] * cmod_[pix_[m]][i][j]
                             * elestrain_[m][j];
            }
          }
          
          ///
          /// Calculate the strain energy for each DC, which would be called in GEMS
          ///

          getAvgStrainengy();

          ///
          /// Sum local strains and stresses into global values     
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
    /// Volume average of global stresses and strains    
    ///

    strxx_ = strxx_ / (double)ns;
    stryy_ = stryy_ / (double)ns;
    strzz_ = strzz_ / (double)ns;
    strxz_ = strxz_ / (double)ns;
    stryz_ = stryz_ / (double)ns;
    strxy_ = strxy_ / (double)ns;
    sxx_ = sxx_ / (double)ns;
    syy_ = syy_ / (double)ns;
    szz_ = szz_ / (double)ns;
    sxz_ = sxz_ / (double)ns;
    syz_ = syz_ / (double)ns;
    sxy_ = sxy_ / (double)ns;

    return;
}


void AppliedStrain::relax (int kmax)
{

    ///
    /// RELAXATION LOOP          
    /// (USER) kmax is the maximum number of times dembx will be called, with  
    /// ldemb conjugate gradient steps performed during each call. The total  
    /// number of conjugate gradient steps allowed for a given elastic computation
    /// is kmax * ldemb.                                                          
    ///


    int ldemb = 50, ltot = 0;
    double utot;
    long int Lstep;

    ///
    /// Call energy to get initial energy and initial gradient.       
    ///

    utot = energy(nx_,ny_,nz_,ns_);

    ///
    /// gg_ is the norm squared of the gradient (gg_ = gb_ * gb_)     
    ///

    gg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        gg_ += (gb_[m][m3] * gb_[m][m3]);
      }
    }
    cout.flush();
	
    for (int kkk = 0; kkk < kmax; kkk++) {

      ///
      /// Call dembx to go into the conjugate gradient solver      
      ///

      Lstep = dembx(ns_, gg_, ldemb, kkk);
      ltot += Lstep;

      ///
      /// Call energy to compute energy after dembx call. If gg_ < gtest_, this will 
      /// be the final energy. If gg_ is still larger than gtest_, then this will    
      /// give an intermediate energy with which to check how the relaxation process 
      /// is coming along.                                                           
      ///

      utot = energy(nx_,ny_,nz_,ns_);

      /*
      cout << "Energy = " << utot << ", gg_ = " << gg_ << endl;
      cout << "Number of conjugate steps = " << ltot << endl;
      cout.flush();
      */

      ///
      /// If relaxation process is finished, jump out of loop.
      ///

      if (gg_ >= gtest_) { 

        ///
        /// If relaxation process will continue, compute an d output stresses and     
        /// strains as an additional aid to judge how the relaxation procedure is     
        /// progressing.                                                              
        ///

        stress(nx_,ny_,nz_,ns_);

        /*
        cout << " stresses: xx,    yy,    zz,    xz,    yz,    xy" << endl;
        cout << "         " << strxx_ << " " << stryy_ << " "
             << strzz_ << " " << strxz_ << " " << stryz_ << " "
             << strxy_ << endl;	
        cout << " strains: xx,    yy,    zz,    xz,    yz,    xy" << endl;
        cout << "         " << sxx_ << " " << syy_ << " "
             << szz_ << " " << sxz_ << " " << syz_ << " "
             << sxy_ << endl;
        */

      } else {

        break; 

      }
    }

    stress(nx_,ny_,nz_,ns_);

    /*
    cout << " stresses: xx,    yy,    zz,    xz,    yz,    xy" << endl;
    cout << "         " << strxx_ << " " << stryy_ << " "
         << strzz_ << " " << strxz_ << " " << stryz_ << " "
         << strxy_ << endl;	
    cout << " strains: xx,    yy,    zz,    xz,    yz,    xy" << endl;
    cout << "         " << sxx_ << " " << syy_ << " "
         << szz_ << " " << sxz_ << " " << syz_ << " "
         << sxy_ << endl;
    */
	
    return;
}


void AppliedStrain::Calc (string fname,
                          double exx,
                          double eyy,
                          double ezz,
                          double exz,
                          double eyz,
                          double exy)
{
    int kmax = 40;
   
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
	  
    for (int i = 0; i < nphase_; i++) {
      cout << " Volume fraction of phase " << i << "  is " << prob_[i] << endl;
    }
    */

    ///
    /// (USER) Set applied strains.                                                
    /// Actual shear strain applied in is exy_, exz_, and eyz_ as given in the     
    /// statements below. The engineering shear strain, by which the shear modulus 
    /// is usually defined, is twice these values.                                 
    ///

    exx_ = exx;
    eyy_ = eyy;
    ezz_ = ezz;
    exz_ = exz;
    eyz_ = eyz;
    exy_ = exy;

    /*
    cout << "Applied engineering strains" << endl;
    cout << " exx_,   eyy_,   ezz_,   exz_,   eyz_,   exy_" << endl;
    cout << exx_ << "   " << eyy_ << "   " << ezz_ << "   "
         << exz_ << "   " << eyz_ << "   " << exy_ << endl;
    */
		  
    ///
    /// Set up elastic modulus variables, finite element stiffness matrices,       
    /// the constant, C, and vector, b, required for computing the energy.         
    /// (USER) If anisotropic elastic moduli tensors are used, these need to be    
    /// input in subroutine femat.                                                 
    ///

    femat(nx_,ny_,nz_,ns_,nphase_);
	  
    ///
    /// Apply chosen strains as a homogeneous macroscopic strain                   
    /// as the initial condition.                                                  
    ///

    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          int m = nx_ * ny_ * k + nx_ * j + i;
	      double x = (double)i;
	      double y = (double)j;
	      double z = (double)k;
	      u_[m][0] = x * exx + y * exy + z * exz;
	      u_[m][1] = x * exy + y * eyy + z * eyz;
	      u_[m][2] = x * exz + y * eyz + z * ezz;
	    }  
      }
    }

    ///
    /// RELAXATION LOOP                                                            
    /// (USER) kmax is the maximum number of times dembx will be called, with       
    /// ldemb conjugate gradient steps performed during each call. The total       
    /// number of conjugate gradient steps allowed for a given elastic computation 
    /// is kmax * ldemb.                                                           
    ///

    relax(kmax);
      
    return;
}

double AppliedStrain::getBulkModulus (string fname)
{
    double bulk;
    double Stress, Strain;
    Stress = Strain = 0.0;

    Calc(fname,0.1,0.1,0.1,0.05,0.05,0.05);

    for (int i = 0; i < 3; i++) {
      Stress += getStress(i);
      Strain += getStrain(i);
    }
    bulk = Stress / 3.0 / Strain;

    return bulk;
}
