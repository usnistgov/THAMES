/**
@file thames.cc

*/

#include "thames.h"
#include "version.h"

int* RanGen::seed_;

/**
@brief The main block for running THAMES.

@return 0 on successful completion, non-zero otherwise
*/
int main (void)
{
    //
    // Set up the strainenergy vector.  We allow no more than 156 phases,
    // but this can be changed below.
    //

    strainenergy.clear();
    strainenergy.resize(156,0.0); 
   
    int choice;
    string buff = "";
    ChemicalSystem *ChemSys;
    Solution *Solut;
    Lattice *Mic;
    ThermalStrain *ThermalStrainSolver;
    AppliedStrain *AppliedStrainSolver;

    KineticModel *KMod;
    
    //
    // Main menu where user decides what kind of simulation this will be.
    //

    cout << "Enter choice: " << endl;
    cout << "  " << QUIT_PROGRAM << ") Exit program " << endl;
    cout << "  " << HYDRATION << ") Hydration " << endl;
    cout << "  " << LEACHING << ") Leaching " << endl;
    cout << "  " << SULFATE_ATTACK << ") Sulfate attack " << endl;
    cin >> choice;
    cout << choice << endl;

    if (choice <= QUIT_PROGRAM || choice > SULFATE_ATTACK) {

        cout << "Exiting program now." << endl << endl;
        exit(1);

    } else {
        
        time_t lt = time('\0');
        struct tm *inittime;
        inittime = localtime(&lt);
        cout << asctime(inittime);
        clock_t starttime = clock();
        
        //
        // User must provide the name of the GEM chemical system definition (CSD) file
        // for the aqueous solution
        //

        cout << "What is the name of the GEM input file for solution? " << endl;
        cin >> buff;
        const string geminput_filename_solution(buff);
        cout << geminput_filename_solution << endl;

        //
        // User must provide the name of the GEM data bridge (DBR) file
        // for the aqueous solution
        //

        cout << "What is the name of the GEM DBR file for solution? " << endl;
        cin >> buff;
        const string geminput_dbrname_solution(buff);
        cout << geminput_dbrname_solution << endl;

        //
        // Create the aqueous solution object
        //

        try {
            Solut = new Solution(geminput_filename_solution,geminput_dbrname_solution);
        }
        catch (bad_alloc &ba) {
            cout << "Bad memory allocation in Solution constructor: " << ba.what() << endl;
        }

        //
        // User must provide the name of the GEM CSD for the whole system
        //

        cout << "What is the name of the GEM input file? " << endl;
        cin >> buff;
        const string geminput_filename(buff);
        cout << geminput_filename << endl;

        //
        // User must provide the name of the GEM data bridge (DBR) file
        // for the aqueous solution
        //

        cout << "What is the name of the GEM DBR file? " << endl;
        cin >> buff;
        const string geminput_dbrname(buff);
        cout << geminput_dbrname << endl;

        //
        // User must provide the name of the file specifying the microstructre
        // phase data
        //

        cout << "What is the name of the phase interface file? " << endl;
        cin >> buff;
        const string pi_filename(buff);
        const string cement_filename(buff);
        cout << pi_filename << endl;

        //
        // Create the ChemicalSystem object
        //

        try {
            ChemSys = new ChemicalSystem(Solut,geminput_filename,geminput_dbrname,pi_filename);
        }
        catch (bad_alloc &ba) {
            cout << "Bad memory allocation in ChemicalSystem constructor: "
                 << ba.what() << endl;
            delete Solut;
        }
        
        //
        // Create the random number generator for phase placement and for shuffling
        // ordered lists
        //

        RanGen rg(-2814357);
        
        //
        // User must specifiy the file containing the 3D microstructure itself
        //

        cout << "What is the name of the MICROSTRUCTURE file? " << endl;
        buff = "";
        cin >> buff;
        const string mic_filename(buff);
        cout << mic_filename << endl;

        //
        // Create the Lattice object to hold the microstructure
        //

        try {
            Mic = new Lattice(ChemSys, Solut, mic_filename);
            cout << "Lattice creation done... " << endl;
            cout << "X size of lattice is " << Mic->getXDim() << endl;
            cout << "Y size of lattice is " << Mic->getYDim() << endl;
            cout << "Z size of lattice is " << Mic->getZDim() << endl;
            cout << " " << endl;
            cout << "Total number of sites is " << Mic->getNumsites() << endl;
        }
        catch (bad_alloc &ba) {
            cout << "Bad memory allocation in Lattice constructor: "
                 << ba.what() << endl;
            delete ChemSys;
            delete Solut;
            exit(1);
        }

        if (choice == SULFATE_ATTACK) {
            
          //
          // This block is executed only if simulating external sulfate attack,
          // in which case we need information about the elastic moduli of the
          // constituent phases, and will need to include a finite element solver
          //

          cout << "What is the name of the elastic modulus file?" << endl;
          buff = "";
          cin >> buff;
          const string phasemod_fname(buff);
          cout << phasemod_fname << endl;

          //
          // Create the ThermalStrain FE solver, which handles phase transformation misfit
          //

          try {
              ThermalStrainSolver = new ThermalStrain(Mic->getXDim(),
                                                      Mic->getYDim(),
                                                      Mic->getZDim(),
                                                      (Mic->getNumsites() + 2),
                                                      ChemSys->getMicphasenum(),1);
              cout << "ThermalStrain object creation done... " << endl;
              ThermalStrainSolver->setPhasemodfname(phasemod_fname);
          }
          catch (bad_alloc &ba) {
              cout << "Bad memory allocation in ThermalStrain constructor: "
                   << ba.what() << endl;
              delete Mic;
              delete ChemSys;
              delete Solut;
              exit(1);
          }

          int nx, ny, nz;
          nx = ny = nz = 3;
          int ns = nx * ny * nz;

          //
          // Create the AppliedStrain FE solver, which handles applied external strains
          //

          try {
              AppliedStrainSolver = new AppliedStrain(nx,ny,nz,ns,ChemSys->getMicphasenum(),1);
              AppliedStrainSolver->setPhasemodfname(phasemod_fname);
          }
          catch (bad_alloc &ba) {
              cout << "Bad memory allocation in AppliedStrain constructor: "
                   << ba.what() << endl;
              delete ThermalStrainSolver;
              delete Mic;
              delete ChemSys;
              delete Solut;
              exit(1);
          }
        
          Mic->setFEsolver(AppliedStrainSolver);
        }

        string jobroot,par_filename,statfilename;
        Controller *Ctrl;
        cout << "About to enter KineticModel constructor" << endl;
        cout.flush();

        //
        // Create the KineticModel object
        //

        try {
            KMod = new KineticModel(ChemSys, Solut, Mic, cement_filename);
        }
        catch (bad_alloc &ba) {
            cout << "Bad memory allocation in KineticModel constructor: "
                 << ba.what() << endl;
            if (choice == SULFATE_ATTACK) {
                delete AppliedStrainSolver;
                delete ThermalStrainSolver;
            }
            delete Mic;
            delete ChemSys;
            delete ThermalStrainSolver;
            delete AppliedStrainSolver;
          exit(1);
        }

        cout << "Finished constructing KineticModel KMod" << endl;
        cout.flush();
        cout << "What is the name of the simulation parameter file? " << endl;
        buff = "";
        cin >> par_filename;
        cout << par_filename << endl;
 
        cout << "What is the root name of this job? " << endl;
        buff = "";
        cin >> jobroot;
        cout << jobroot << endl;
 
        statfilename = jobroot + ".stats";
        cout << "About to go into Controller constructor" << endl;
        cout.flush();

        //
        // Create the Controller object to direct flow of the program
        //

        try {
          Ctrl = new Controller(Mic,
                                KMod,
                                ChemSys,
                                Solut,
                                ThermalStrainSolver,
                                par_filename,
                                jobroot);
        }
        catch (bad_alloc &ba) {
          cout << "Bad memory allocation in Controller constructor: "
               << ba.what() << endl;
          delete KMod;
          if (choice == SULFATE_ATTACK) {
              delete AppliedStrainSolver;
              delete ThermalStrainSolver;
          }
          delete Mic;
          delete ChemSys;
          delete ThermalStrainSolver;
          delete AppliedStrainSolver;
          exit(1);
        }
        
        //
        // Write a formatted output of the simulation parameters for later reference
        //

        writeReport(jobroot,inittime,mic_filename,par_filename,
                    geminput_filename,ChemSys,Ctrl);

        //
        // Launch the main controller to run the simulation
        //

        cout << "Going into Controller::doCycle now" << endl;
        cout.flush();
        Ctrl->doCycle(statfilename,choice);
        
        // 
        // Simulation is finished.  Record and output the timing data.
        //

        time_t lt1 = time('\0');
        struct tm *inittime1;
        inittime1 = localtime(&lt1);
        cout << asctime(inittime1);
        clock_t endtime = clock();
        double elapsedtime = (double)(endtime - starttime)/CLOCKS_PER_SEC;
        double ltD = difftime(lt1,lt);
        cout << endl << "Total time = " << ltD << " seconds" << endl;
        cout << endl << "Total time with clock = " << elapsedtime
             << " seconds" << endl;
        
        //
        // Delete the dynamically allocated memory
        //

        cout << "About to delete Ctrl pointer... ";
        cout.flush();
        delete Ctrl;
        cout << "Done!" << endl;
        cout.flush();

        cout << "About to delete KMod pointer... ";
        cout.flush();
        delete KMod;
        cout << "Done!" << endl;

        cout << "About to delete Mic pointer... ";
        cout.flush();
        delete Mic;
        cout << "Done!" << endl;

        cout << "About to delete ChemSys pointer... ";
        cout.flush();
        delete ChemSys;
        cout << "Done!" << endl;
        cout.flush();
    }

    return 0;
}

void writeReport (const string &jobroot,
                  struct tm *itime,
                  const string &mfname,
                  const string &parfilename,
                  const string &csname,
                  ChemicalSystem *csys,
                  Controller *ctr)
{
    string statname = jobroot + ".stats";
    string jfilename = jobroot + ".report";
    ofstream out(jfilename.c_str(),ios::app);
    if (!out.is_open()) {
        cout << "WARNING:  Could not open report file" << endl;
        return;
    }

    //
    // Write the time the job was executed
    //

    out << "THAMES job " << jobroot;
    out << " initialized on " << asctime(itime) << endl;
    out << endl;
    out << "INPUT FILES USED:" << endl;
    out << "   Microstructure file name: " << mfname << endl;
    out << "        GEM input file name: " << csname << endl;
    out << endl;
    out << "OUTPUT FILES GENERATED:" << endl;
    out << "     Global phase fractions: " << statname << endl;
    out << endl;
    csys->writeChemSys();
    out << endl;

    out.close();
    return;
}
