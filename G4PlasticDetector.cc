//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// ********************************************************************
// G4PlasticDetector.cc
//
// Description: Main Geant4 simulation program for the plastic scintillator detector
// used by McMaster University to perform dosimetry measurements for the lens of the eye.
// ********************************************************************

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

// Simulation Files
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	// Choose the Random engine
  	G4Random::setTheEngine(new CLHEP::RanecuEngine);
     
  	// Construct the default run manager
	#ifdef G4MULTITHREADED
		G4cout << "Simulation is running on " << G4Threading::G4GetNumberOfCores() << " threads";
        G4MTRunManager* runManager = new G4MTRunManager;
        G4int nCores = G4Threading::G4GetNumberOfCores();
        runManager->SetNumberOfThreads(nCores);
	#else
  		G4RunManager* runManager = new G4RunManager;
	#endif

	//G4RunManager* runManager = new G4RunManager;

	// Set mandatory initialization classes
  	DetectorConstruction* detector = new DetectorConstruction();
  	runManager->SetUserInitialization(detector);
	  
	// Select a physics list
	PhysicsList* physicsList = new PhysicsList();
  	runManager->SetUserInitialization(physicsList);

	// Set user action classes
	runManager->SetUserInitialization(new ActionInitialization(detector)); 
  
	#ifdef G4VIS_USE
		// Initialize visualization
		G4VisManager* visManager = new G4VisExecutive;
		// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
		// G4VisManager* visManager = new G4VisExecutive("Quiet");
		visManager->Initialize();
	#endif

  	// Get the pointer to the User Interface manager
  	G4UImanager* UImanager = G4UImanager::GetUIpointer();

  	if (argc!=1)   // batch mode
    {
    	G4String command = "/control/execute ";
      	G4String fileName = argv[1];
      	UImanager->ApplyCommand(command+fileName);
    	}
  	else
    	{  // interactive mode : define UI session
		#ifdef G4UI_USE
      		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      		ui->SessionStart();
      		delete ui;
		#endif
    }

  	// Job termination
  	// Free the store: user actions, physics_list and detector_description are
  	// owned and deleted by the run manager, so they should not be deleted 
  	// in the main() program !

	#ifdef G4VIS_USE
  		delete visManager;
	#endif
  	delete runManager;

  	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
