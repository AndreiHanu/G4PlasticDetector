#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Run;
class DetectorConstruction;
class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
	public: 
		// Constructor
  		RunAction(DetectorConstruction* det, PrimaryGeneratorAction* primary=0);
  		// Destructor
  		virtual ~RunAction();

		// Methods
  		virtual G4Run* GenerateRun(); 

  		virtual void BeginOfRunAction(const G4Run*);
 		virtual void EndOfRunAction(const G4Run*);

	private:
		DetectorConstruction* detector;
		PrimaryGeneratorAction* particleGun;
		
		// Output File
		G4String outputFile_INFO;
		FILE* pFile_INFO;
};

#endif

