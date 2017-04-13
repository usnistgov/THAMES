/**
@file ElasticModel.cc
@brief Method definitions for the ElasticModel base class

*/
#include "ElasticModel.h"

ElasticModel::ElasticModel (int nx,
                            int ny,
                            int nz,
                            int dim,
                            int nphase,
                            int npoints) 
{
    ///
    /// Assign the dimensions of the finite element (FE) mesh
    ///

    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    ns_ = nx_ * ny_ * nz_;

    cout << "nx_ = " << nx_ << " ny_ = " << ny_
         << " nz_ = " << nz_ << " ns_ = " << ns_
         << endl;

    ///
    /// Initialize the prescribed stresses and strains
    ///

    strxx_ = stryy_ = strzz_ = 0.0;
    strxz_ = stryz_ = strxy_ = 0.0;
    sxx_ = syy_ = szz_ = 0.0;
    sxz_ = syz_ = sxy_ = 0.0;

    ///
    /// Initialize the stresses and strains (arrays) for each element
    ///

    elestress_.clear();
    elestrain_.clear();
    elestress_.resize(ns_);
    elestrain_.resize(ns_);

    for (int i = 0; i < ns_; i++) {
      elestress_[i].resize(6,0.0);
      elestrain_[i].resize(6,0.0);
    }
    
    ///
    /// Initialize the array of strain energies for each element
    ///

    strainengy_.clear();
    strainengy_.resize(ns_,0.0);

    ///
    /// `nphase_` is the number of phases being considered in the problem.           
    /// The values of `pix(m)` will run from 1 to `nphase_`.                           		 
    ///

    nphase_ = nphase;

    ///
    /// Establish the stopping criterion for convergence (proportional to
    /// the number of elements, `ns_`
    ///

    gtest_ = (1.0e-8) * ns_;
    gg_ = 0.0;

    ///
    /// The parameter `phasemod_[i][j]` is the `bulk[i][0]` and `shear[i][1]`
    /// moduli of the i'th phase. These can be input in terms of Young's moduli    
    /// E[i][0] and Poisson's ratio nu[i][1]. The program, then changes them to    
    /// bulk and shear moduli. For anisotropic elastic material, one can directly  
    /// input the elastic moduli tensor cmod in subroutine femat, and skip this    
    /// part.
    ///
    /// @warning If you wish to input in terms of bulk (0) and shear (1), then make   
    /// sure to comment out the following loop.                                    

    phasemod_.clear();
    phasemod_.resize(nphase_);
    for (int ijk = 0; ijk < nphase_; ijk++) {
      phasemod_[ijk].resize(2,0.0);
    }
	
    ///
    /// Initialize the neighbor table for each element
    ///

    ib_.clear();
    ib_.resize(dim);
    for (int m = 0; m < dim; m++) {
      ib_[m].resize(27,0);
    }

    BuildNeighbor();
	
    ///
    /// `npoints_` is the number of microstructres to use.                    
    ///

    npoints_ = npoints;
	
    ///
    /// Initialize `pix_[m]` and volume fraction of each phase `prob_[]`.       
    ///

    pix_.clear();
    pix_.resize(ns_,0);
	
    prob_.clear();
    prob_.resize(nphase_,0.0);
	
    ///
    /// Initialize displacement at each node `u_`.
    ///

    u_.clear();
    u_.resize(dim);
    for (int m = 0; m < dim; m++) {
      u_[m].resize(3,0.0);
    }
	
    ///
    /// Initialize elastic modulus variables, `cmod_`,                              
    /// finite element stiffness matrices, `dk_`, the energy constant, `C_`,
    /// and the linear term coefficient vector, `b_`, required for
    /// computing the energy.                                         

    cmod_.clear();
    cmod_.resize(nphase_);
    for (int ijk = 0; ijk < nphase_; ijk++) {
      cmod_[ijk].resize(6);
      for (int j = 0; j < 6; j++) {
        cmod_[ijk][j].resize(6,0.0);
      }
    }
	
    dk_.clear();
    dk_.resize(nphase_);
    for (int ijk = 0; ijk < nphase_; ijk++) {
      dk_[ijk].resize(8);
      for (int i = 0; i < 8; i++) {
        dk_[ijk][i].resize(3);
        for (int k = 0; k < 3; k++) {
          dk_[ijk][i][k].resize(8);
          for (int j = 0; j < 8; j++) {
            dk_[ijk][i][k][j].resize(3,0.0);
          }
        }
      }
    }
	
    b_.clear();
    b_.resize(dim);
    for (int m = 0; m < dim; m++) {
      b_[m].resize(3,0.0);
    }
	
    C_ = 0.0;	

    ///
    /// Initialize the gradient at each element, `gb_[ns_][3]`,
    /// the auxiliary conjugate gradient vector variable at each element, `h_[ns_][3]`,
    /// and the other conjugate gradient vector variable at each element,
    /// `Ah_[ns_][3]`.                  
    ///

    gb_.clear();
    h_.clear();
    Ah_.clear();

    gb_.resize(dim);
    h_.resize(dim);
    Ah_.resize(dim);
    for (int m = 0; m < dim; m++) {
      gb_[m].resize(3,0.0);
      h_[m].resize(3,0.0);
      Ah_[m].resize(3,0.0);
    }

    return;
	
}

void ElasticModel::BuildNeighbor () 
{
    cout << "now in the BuildNeighbor function." << endl;

    ///
    /// First construct the 27 neighbor table in terms of delta i, delta j, and    
    /// delta k information. (see Table 3 in manual)                               
    ///

    int in[27], jn[27], kn[27];
    in[0] = 0;
    in[1] = 1;
    in[2] = 1;
    in[3] = 1;
    in[4] = 0;
    in[5] = (-1);
    in[6] = (-1);
    in[7] = (-1);

    jn[0] = 1;
    jn[1] = 1;
    jn[2] = 0;
    jn[3] = (-1);
    jn[4] = (-1);
    jn[5] = (-1);
    jn[6] = 0;
    jn[7] = 1;
	
    for (int n = 0; n < 8; n++) {
      kn[n] = 0;
      kn[n+8] = (-1);
      kn[n+16] = 1;
      in[n+8] = in[n];
      in[n+16] = in[n];
      jn[n+8] = jn[n];
      jn[n+16] = jn[n];
    }

	
    in[24] = 0;
    in[25] = 0;
    in[26] = 0;
    jn[24] = 0;
    jn[25] = 0;
    jn[26] = 0;
    kn[24] = (-1);
    kn[25] = 1;
    kn[26] = 0;
	
    ///
    /// Now construct neighbor table according to 1D labels                       
    /// Matrix ib_[m][n] gives the 1-d label of the n'th neighbor (n=0,26) of the  
    /// node labelled m.
    ///
    /// This set of for loops also manually checks for wrapping due to
    /// periodic boundary conditions in all three directions
    ///
   
    int nxy = nx_ * ny_;
    int m, m1;
    int i1, j1, k1;
    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          m = nxy * k + nx_ * j + i;
          for (int n = 0; n < 27; n++) {
            i1 = i + in[n];
            j1 = j + jn[n];
            k1 = k + kn[n];
            if (i1 < 0) i1 += nx_;
            else if (i1 >= nx_) i1 -= nx_;
            if (j1 < 0) j1 += ny_;
            else if (j1 >= ny_) j1 -= ny_;
            if (k1 < 0) k1 += nz_;
            else if (k1 >= nz_) k1 -= nz_;
            m1 = nxy * k1 + nx_ * j1 + i1;
            ib_[m][n] = m1;
          }
        }
      }
    }

    cout << "after building the neighbors." << endl;
    return;
}

void ElasticModel::ElasModul (string phasemod_fname,
                              int nphase) 
{
    ///
    /// Open the file input stream with the specified name.
    /// This file holds the values of the Young's modulus [GPa] and
    /// Poisson's ratio of each phase
    ///

    ifstream in(phasemod_fname.c_str());
    if(!in) {
      cout << "can't open the file: " << phasemod_fname << endl;
      exit(1);
    } else {
      string buff;
      int phaseid;
      double elsmodul = 0.0;
      for (int i = 0; i < nphase; i++) {
        in >> phaseid;
        in >> buff;
        in >> elsmodul;
        phasemod_[phaseid][0] = elsmodul;
        in >> elsmodul;
        phasemod_[phaseid][1] = elsmodul;
      }
    }

    ///
    /// The program uses bulk modulus (0) and shear modulus (1), so transform   
    /// Young's modulus (0) and Poisoon's ratio (1) to them.
    ///

    for (int i = 0; i < nphase_; i++) {  
      double save = phasemod_[i][0];
      phasemod_[i][0] = (phasemod_[i][0]/3.0)/(1 - 2 * phasemod_[i][1]);
      phasemod_[i][1] = (save/2.0)/(1.0 + phasemod_[i][1]);
    }

    ///
    /// Set up the elastic modulus tensor for each phase, which uses
    /// engineering notation for the components as described in the
    /// documentation for the `cmod_` member.
    ///
    /// The ck and cmu matrices are used to multiply by the (scalar)
    /// bulk and shear moduli values, `phasemod_[0]` and `phasemod_[1]`,
    /// respectively.
    ///

    double ck[6][6],cmu[6][6];
    ck[0][0] = 1.0;
    ck[0][1] = 1.0;
    ck[0][2] = 1.0;
    ck[0][3] = 0.0;
    ck[0][4] = 0.0;
    ck[0][5] = 0.0;
    ck[1][0] = 1.0;
    ck[1][1] = 1.0;
    ck[1][2] = 1.0;
    ck[1][3] = 0.0;
    ck[1][4] = 0.0;
    ck[1][5] = 0.0;
    ck[2][0] = 1.0;
    ck[2][1] = 1.0;
    ck[2][2] = 1.0;
    ck[2][3] = 0.0;
    ck[2][4] = 0.0;
    ck[2][5] = 0.0;
    ck[3][0] = 0.0;
    ck[3][1] = 0.0;
    ck[3][2] = 0.0;
    ck[3][3] = 0.0;
    ck[3][4] = 0.0;
    ck[3][5] = 0.0;
    ck[4][0] = 0.0;
    ck[4][1] = 0.0;
    ck[4][2] = 0.0;
    ck[4][3] = 0.0;
    ck[4][4] = 0.0;
    ck[4][5] = 0.0;
    ck[5][0] = 0.0;
    ck[5][1] = 0.0;
    ck[5][2] = 0.0;
    ck[5][3] = 0.0;
    ck[5][4] = 0.0;
    ck[5][5] = 0.0;
	
    cmu[0][0] = 4.0/3.0;
    cmu[0][1] = (-2.0)/3.0;
    cmu[0][2] = (-2.0)/3.0;
    cmu[0][3] = 0.0;
    cmu[0][4] = 0.0;
    cmu[0][5] = 0.0;
    cmu[1][0] = (-2.0)/3.0;
    cmu[1][1] = 4.0/3.0;
    cmu[1][2] = (-2.0)/3.0;
    cmu[1][3] = 0.0;
    cmu[1][4] = 0.0;
    cmu[1][5] = 0.0;
    cmu[2][0] = (-2.0)/3.0;
    cmu[2][1] = (-2.0)/3.0;
    cmu[2][2] = 4.0/3.0;
    cmu[2][3] = 0.0;
    cmu[2][4] = 0.0;
    cmu[2][5] = 0.0;
    cmu[3][0] = 0.0;
    cmu[3][1] = 0.0;
    cmu[3][2] = 0.0;
    cmu[3][3] = 1.0;
    cmu[3][4] = 0.0;
    cmu[3][5] = 0.0;
    cmu[4][0] = 0.0;
    cmu[4][1] = 0.0;
    cmu[4][2] = 0.0;
    cmu[4][3] = 0.0;
    cmu[4][4] = 1.0;
    cmu[4][5] = 0.0;
    cmu[5][0] = 0.0;
    cmu[5][1] = 0.0;
    cmu[5][2] = 0.0;
    cmu[5][3] = 0.0;
    cmu[5][4] = 0.0;
    cmu[5][5] = 1.0;

    ///
    /// Construct `cmod_` tensor for each phase by matrix multiplication
    ///
    /// \f{equation}
    ///     c_{ij} = K c^k_{ij} + \mu c^{\mu}_{ij}
    /// \f}
    ///

    for (int ijk = 0; ijk < nphase; ijk++) {
      for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
          cmod_[ijk][i][j] = phasemod_[ijk][0] * ck[i][j] + 
                             phasemod_[ijk][1] * cmu[i][j];
        } 
      }
    }
    return;
	
}

void ElasticModel::ppixel (string fname,
                           int nphase)
{
    ///
    /// If you want to set up a test image inside the program, instead of
    /// reading it in from a file, this should be done inside this method
    ///

    string buff,version;
    double resolution;
    int m;

    ///
    /// Open and read the input file stream with the microstructure data
    ///
  
    ifstream in(fname.c_str());

    if (!in) {

      cout << "can't open the file: " << fname << endl;
      exit(1);

    } else {

      const string VERSIONSTRING("Version:");
      const string IMGSIZESTRING("Image_Size:");
      const string IMGRESSTRING("Image_Resolution:");
      const string XSIZESTRING("X_Size:");
      const string YSIZESTRING("Y_Size:");
      const string ZSIZESTRING("Z_Size:");
      in >> buff;
      if(buff == VERSIONSTRING) {
        in >> version;
        in >> buff;
        if(buff == XSIZESTRING) {
          in >> nx_;
          in >> buff;
          in >> ny_;
          in >> buff;
          in >> nz_;
        } else if (buff == IMGSIZESTRING) {
          in >> nx_;
          ny_ = nz_ = nx_;
        }
        in >> buff;
        if(buff == IMGRESSTRING) {
          in >> resolution;
        }
      } else {

        ///
        /// Image file does not have the expected header.
        /// Assume default size to 100 and resolution to 1.0
        /// micrometers.
        ///

        version = "2.0";
        resolution = 1.0;
        nx_ = ny_ = nz_ = 100;
      }

      ///
      /// Each line of the microstructure file contains the phase id
      /// to assign, and the microstructure file must be written in
      /// the same order as the finite elements are populated
      ///

      ns_ = nx_ * ny_ * nz_;
      int nxy = nx_ * ny_;
      for(int k = 0; k < nz_; k++) {
        for(int j = 0; j < ny_; j++) {
          for(int i = 0; i < nx_; i++) {
            m = nxy * k + nx_ * j + i;
            in >> pix_[m];
          }
        }
      }
		
      in.close();

      ///
      /// Check for wrong phase labels--less than 1 or greater than nphase_.    
      /// Note that nothing is done about the error; it is just reported
      /// to standard out.
      ///
      /// @todo See if it makes sense to do value checking as exception handling
      ///

      for (int m = 0; m < ns_; m++) {
        if (pix_[m] < 0) {
          cout << "Phase label in pix < 0 --- error at "
               << m << endl;
        } else if (pix_[m] >= nphase_) {
          cout << "Phase label in pix >= nphase_ --- error at "
               << m << endl;
        }
      }

    }

    return;
}

void ElasticModel::getAvgStrainengy ()
{
    ///
    /// Clear and initialize the array that will hold the
    /// volume averaged strain energy of each phase in the microstructure
    ///

    avgStrainengy_.clear();
    avgStrainengy_.resize(nphase_,0.0);

    ///
    /// Sum the strain energy in each phase by one sweep through the mesh
    ///

    for (int m = 0; m < ns_; m++) {
      avgStrainengy_[pix_[m]] += strainengy_[m];
    }
    
    ///
    /// Normalize by the volume of the phase
    ///

    for (int i = 0; i < nphase_; i++) {
      avgStrainengy_[i] = avgStrainengy_[i] / (float)(prob_[i] * ns_);
    }
   
    ///
    /// The variable `strainenergy` (not `strainenergy_`) is defined
    /// in the StrainEnergy.h header file.  It is used in the
    /// GEM-IPM node class to facilitate communication of phase
    /// strain energies back and forth between the two, because
    /// strain energy adds to the Gibbs energy of a phase and therefore
    /// should influence equilibrium calculations.
    ///
    /// The block below relates the COMPUTED average strain energy density,
    /// `avgStrainengy_`, with units of GJ/m3, to the strain
    /// energy PER MOLE in the corresponding GEM dependent components that
    /// make up that phase, `strainenergy`, in units of J/mol.
    /// To make the conversion, we must multiply avgStrainengy_ by the molar
    /// volume (units of m3/mol) and then multiply by e9 to convert GJ to J.
    ///
    /// To take the example of C3S, it is a microstructure phase (id = 2),
    /// with one associated GEM DC (also called C3S, with GEM id 115 and
    /// a molar volume of 7.318e-5 m3/mol.  Multiplying that molar volume by
    /// e9 to convert from GJ to J gives an overall converstion factor of
    /// 7.318e4 J m3 / GJ mol.
    ///
    /// @remarks It would be better if these index values
    /// and molar volumes were not hard-coded like this.
    ///
    /// @todo Generalize the indices and molar volumes
    ///

    strainenergy.clear();
    strainenergy.resize(156,0.0);
    
    ///
    /// DC component making up the H2O microstructure phase
    ///

    strainenergy[70] = avgStrainengy_[1] * (1.8068E4);  /* Water */

    ///
    /// DC component making up the C3S microstructure phase
    ///

    strainenergy[115] = avgStrainengy_[2] * (7.318E4);  /* C3S */

    ///
    /// DC component making up the C2S microstructure phase
    ///

    strainenergy[113] = avgStrainengy_[3] * (5.179E4);  /* Beta-C2S */

    ///
    /// DC component making up the C3A microstructure phase
    ///

    strainenergy[114] = avgStrainengy_[4] * (8.9217E4); /* C3A */

    ///
    /// DC component making up the C4AF microstructure phase
    ///

    strainenergy[116] = avgStrainengy_[5] * (1.30202E5); /* C4AF */

    ///
    /// DC components making up the GYPSUM microstructure phase
    ///

    strainenergy[134] = avgStrainengy_[8] * (7.469E4);   /* Gypsum */

    ///
    /// DC components making up the HEMIANH microstructure phase
    ///

    strainenergy[135] = avgStrainengy_[9] * (6.173E4);   /* Hemihydrate */
    strainenergy[133] = avgStrainengy_[9] * (4.594E4);   /* Anhydrite */

    ///
    /// DC components making up the CACO3 microstructure phase
    ///

    strainenergy[128] = avgStrainengy_[10] * (3.6934E4); /* Calcite */
    strainenergy[127] = avgStrainengy_[10] * (3.415E4);  /* Aragonite */

    strainenergy[132] = avgStrainengy_[11] * (3.306E4);  /* Portlandite */
	
    ///
    /// DC components making up the CSH microstructure phase
    ///

    strainenergy[80] = avgStrainengy_[12] * (5.87E4);    /* Tob-II */
    strainenergy[79] = avgStrainengy_[12] * (7.84E4);    /* Jennite */
	
    ///
    /// DC components making up the AFMC microstructure phase
    ///

    strainenergy[95] = avgStrainengy_[13] * (2.84515E5); /* Hemicarbonate */
    strainenergy[96] = avgStrainengy_[13] * (2.96472E5); /* Fe-hemicarbonate */
    strainenergy[97] = avgStrainengy_[13] * (2.61958E5); /* Monocarbonate */
    strainenergy[98] = avgStrainengy_[13] * (2.90190E5); /* Fe-monocarbonate */
	
    ///
    /// DC components making up the MONOSULPH microstructure phase
    ///

    strainenergy[91] = avgStrainengy_[14] * (3.0903E5);  /* Monosulfate */
    strainenergy[92] = avgStrainengy_[14] * (3.2114E5);  /* Fe-monosulfate */
    strainenergy[99] = avgStrainengy_[14] * (2.7398E5);  /* 1-C4AH13 */
    strainenergy[100] = avgStrainengy_[14] * (3.0903E5); /* 1-Monosulfate */
    strainenergy[101] = avgStrainengy_[14] * (2.7398E5); /* 2-C4AH13 */
    strainenergy[102] = avgStrainengy_[14] * (3.0903E5); /* 2-Monosulfate */
	
    ///
    /// DC components making up the ETTR microstructure phase
    ///

    strainenergy[87] = avgStrainengy_[15] * (7.0703E5);  /* Ettringite */
    strainenergy[88] = avgStrainengy_[15] * (7.1756E5);  /* Fe-ettringite */
    strainenergy[89] = avgStrainengy_[15] * (7.0703E5);  /* 1-Ettringite */
    strainenergy[90] = avgStrainengy_[15] * (7.1756E5);  /* 1-Fe-ettringite */
    strainenergy[103] = avgStrainengy_[15] * (6.504E5);  /* Tricarboaluminate */
    strainenergy[104] = avgStrainengy_[15] * (7.0703E5); /* 2-Ettringite */
    strainenergy[105] = avgStrainengy_[15] * (6.504E5);  /* 1-Tricarboaluminate */
    strainenergy[106] = avgStrainengy_[15] * (7.0703E5); /* 3-Ettringite */
    strainenergy[126] = avgStrainengy_[15] * (7.0703E5); /* 4-Ettringite */
	
    ///
    /// DC component making up the BRUCITE microstructure phase
    ///

    strainenergy[152] = avgStrainengy_[16] * (2.463E4);  /* Brucite */

    ///
    /// DC components making up the HYDROTALC microstructure phase
    ///

    strainenergy[107] = avgStrainengy_[17] * (2.202E5);  /* Hydrotalcite */
    strainenergy[108] = avgStrainengy_[17] * (2.324E5);  /* Fe-hydrotalcite */
    strainenergy[150] = avgStrainengy_[17] * (2.204E5);  /* CO3-hydrotalcite */
	
    ///
    /// DC components making up the AFM microstructure phase
    ///

    strainenergy[81] = avgStrainengy_[18] * (1.8386E5);  /* C2AH8 */
    strainenergy[82] = avgStrainengy_[18] * (1.9359E5);  /* C2FH8 */
    strainenergy[83] = avgStrainengy_[18] * (1.49702E5); /* C3AH6 */
    strainenergy[84] = avgStrainengy_[18] * (1.55287E5); /* C3FH6 */
    strainenergy[85] = avgStrainengy_[18] * (2.7398E5);  /* C4AH13 */
    strainenergy[86] = avgStrainengy_[18] * (2.8594E5);  /* C4FH13 */
    strainenergy[120] = avgStrainengy_[18] * (1.93985E5); /* CAH10 */
	
    ///
    /// DC component making up the LIME microstructure phase
    ///

    strainenergy[131] = avgStrainengy_[19] * (1.6764E4);  /* Free lime */
   
    return;
}

void ElasticModel::writeStress (string &root,
                                double time,
                                int index)
{
    if (index >= 0 && index < 6) {
      cout << "writing ppm file..." << endl;
      double min,max;
      min = max = 0.0;

      ///
      /// Create and initialize the local rgb vector
      ///
 
      vector<int> color;
      color.clear();
      color.resize(3,0);

      ///
      /// Specify the file name and open the output stream
      ///

      ostringstream ostr;
      ostr << (int) (time * 100.0);
      string timestr(ostr.str());
      string ofname(root);
      string ofpngname(root);
      if (index == 0) {
        ofname = ofname + "." + "stress-xx." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-xx." + timestr + ".png";
      } else if (index == 1) {
        ofname = ofname + "." + "stress-yy." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-yy." + timestr + ".png";
      } else if (index == 2) {
        ofname = ofname + "." + "stress-zz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-zz." + timestr + ".png";
      } else if (index == 3) {
        ofname = ofname + "." + "stress-xz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-xz." + timestr + ".png";
      } else if (index == 4) {
        ofname = ofname + "." + "stress-yz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-yz." + timestr + ".png";
      } else if (index == 5) {
        ofname = ofname + "." + "stress-xy." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "stress-xy." + timestr + ".png";
      }
      ofstream out(ofname.c_str());
      if (!out.is_open()) {
        cout << ofname << "could not open." << endl;
        exit(1);
      }

      ///
      /// Write the PPM header for a full color image
      /// Currently only uses cyan of different intensities
      ///

      out << "P3" << endl;
      out << nx_ << " " << nz_ << endl;
      out << "255" << endl;

      unsigned int slice = nx_/2;

      for (int j = 0; j < ny_; j++) {
        for (int k = 0; k < nz_; k++) {
          int m = nx_ * ny_ * k + nx_ * j + slice;
          if (min > elestress_[m][index]) min = elestress_[m][index];
          if (max < elestress_[m][index]) max = elestress_[m][index];
        }
      }
      cout << "minimum stress-" << index << " is: " << min << endl;
      cout << "maximum stress-" << index << " is: " << max << endl;
      for (int k = 0; k < nz_; k++) {
        for (int j = 0; j < nz_; j++) {
          int m = nx_ * ny_ * k + nx_ * j + slice;
          color[1] = (int) (((elestress_[m][index] - min) / (max - min)) * 255);
          color[2] = (int) (((elestress_[m][index] - min) / (max - min)) * 255);
          out << color[0] << " " << color[1] << " " << color[2] << endl;
        }
      }
    
      out.close();

      ///
      /// PPM file is finished.  Now convert to PNG using ImageMagick convert
      /// command via a system call (not recommended).
      ///

      string buff = "convert " + ofname + " " + ofpngname;
      system(buff.c_str());
      return;

    } else {

      cout << "index out of range. should be between 0 and 6." << endl;
      exit(1); 

    }
}

void ElasticModel::writeStrain (string &root,
                                double time,
                                int index)
{
    if (index >= 0 && index < 6) {
      cout << "writing ppm file..." << endl;
      double min,max;
      min = max = 0.0;

      ///
      /// Create and initialize the local rgb vector
      ///
 
      vector<int> color;
      color.clear();
      color.resize(3,0);

      ///
      /// Specify the file name and open the output stream
      ///

      ostringstream ostr;
      ostr << (int) (time * 100.0);
      string timestr(ostr.str());
      string ofname(root);
      string ofpngname(root);
      if (index == 0) {
        ofname = ofname + "." + "strain-xx." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-xx." + timestr + ".png";
      } else if (index == 1) {
        ofname = ofname + "." + "strain-yy." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-yy." + timestr + ".png";
      } else if (index == 2) {
        ofname = ofname + "." + "strain-zz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-zz." + timestr + ".png";
      } else if (index == 3) {
        ofname = ofname + "." + "strain-xz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-xz." + timestr + ".png";
      } else if (index == 4) {
        ofname = ofname + "." + "strain-yz." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-yz." + timestr + ".png";
      } else if (index == 5) {
        ofname = ofname + "." + "strain-xy." + timestr + ".ppm";
        ofpngname = ofpngname + "." + "strain-xy." + timestr + ".png";
      }
      ofstream out(ofname.c_str());
      if (!out.is_open()) {
        cout << ofname << "could not open." << endl;
        exit(1);
      }

      ///
      /// Write the PPM header for a full color image
      /// Currently only uses cyan of different intensities
      ///

      out << "P3" << endl;
      out << nx_ << " " << nz_ << endl;
      out << "255" << endl;

      unsigned int slice = nx_/2;

      for (int j = 0; j < ny_; j++) {
        for (int k = 0; k < nz_; k++) {
          int m = nx_ * ny_ * k + nx_ * j + slice;
          if (min > elestrain_[m][index]) min = elestrain_[m][index];
          if (max < elestrain_[m][index]) max = elestrain_[m][index];
        }
      }
      cout << "minimum strain-" << index << " is: " << min << endl;
      cout << "maximum strain-" << index << " is: " << max << endl;
      for (int k = 0; k < nz_; k++) {
        for (int j = 0; j < nz_; j++) {
          int m = nx_ * ny_ * k + nx_ * j + slice;
          color[1] = (int) (((elestrain_[m][index] - min) / (max - min)) * 255);
          color[2] = (int) (((elestrain_[m][index] - min) / (max - min)) * 255);
          out << color[0] << " " << color[1] << " " << color[2] << endl;
        }
      }
    
      out.close();

      ///
      /// PPM file is finished.  Now convert to PNG using ImageMagick convert
      /// command via a system call (not recommended).
      ///

      string buff = "convert " + ofname + " " + ofpngname;
      system(buff.c_str());
      return;

    } else {

      cout << "index out of range. should be between 0 to 5." << endl;
      exit(1); 

    }
}

void ElasticModel::writeDisp (string &root,
                              double time)
{
  cout << "writing displacement file..." << endl;

  ///
  /// Specify the file name and open the output stream
  ///

  ostringstream ostr;
  ostr << (int) (time * 100.0);
  string timestr(ostr.str());
  string ofname(root);
  ofname = ofname + "." + "disp." + timestr + ".dat";
  ofstream out(ofname.c_str());

  if (!out.is_open()) {
    cout << ofname << "could not open." << endl;
    exit(1);
  }

  for (int k = 0; k < nz_; k++) {
    for (int j = 0; j < ny_; j++) {
      for (int i = 0; i < nx_; i++) {
        int m = nx_ * ny_ * k + nx_ * j + i;
        out << u_[m][0] << "    " << u_[m][1] << "    " << u_[m][2] << endl;
      }
    }
  }

  out.close();

  return;
}

void ElasticModel::writeStrainEngy (string &root,
                                    double time)
{
    cout << "writing ppm file..." << endl;
    double min,max;
    min = max = 0.0;

    ///
    /// Create and initialize the local rgb vector
    ///
 
    vector<int> color;
    color.clear();
    color.resize(3,0);

    ///
    /// Specify the file name and open the output stream
    ///
   
    ostringstream ostr;
    ostr << (int) (time * 100.0);
    string timestr(ostr.str());
    string ofname(root);
    string ofpngname(root);
    ofname = ofname + "." + "strainengy." + timestr + ".ppm";
    ofpngname = ofpngname + "." + "strainengy." + timestr + ".png";
    ofstream out(ofname.c_str());
    if (!out.is_open()) {
      cout << ofname << "could not open." << endl;
      exit(1);
    }

    ///
    /// Write the PPM header for a full color image
    /// Currently only uses cyan of different intensities
    ///

    out << "P3" << endl;
    out << nx_ << " " << nz_ << endl;
    out << "255" << endl;

    unsigned int slice = nx_/2;
    for (int j = 0; j < ny_; j++) {
      for (int k = 0; k < nz_; k++) {
        int m = nx_ * ny_ * k + nx_ * j + slice;
        if (min > strainengy_[m]) min = strainengy_[m];
        if (max < strainengy_[m]) max = strainengy_[m];
      }
    }
    cout << "minimum strainengy is: " << min << endl;
    cout << "maximum strainengy is: " << max << endl;
    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < nz_; j++) {
        int m = nx_ * ny_ * k + nx_ * j + slice;
        color[1] = (int) (((strainengy_[m] - min) / (max - min)) * 255);
        color[2] = (int) (((strainengy_[m] - min) / (max - min)) * 255);
        out << color[0] << " " << color[1] << " " << color[2] << endl;
      }
    }
    
    out.close();

    ///
    /// PPM file is finished.  Now convert to PNG using ImageMagick convert
    /// command via a system call (not recommended).
    ///

    string buff = "convert " + ofname + " " + ofpngname;
    system(buff.c_str());
    return;
}
