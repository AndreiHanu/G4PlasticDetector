#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class DetectorConstruction;
class PrimaryGeneratorAction;

class Run : public G4Run
{
	public:
		// Constructor
  		Run(DetectorConstruction* det, PrimaryGeneratorAction* primary=0);
  		// Destructor
  		virtual ~Run();
		
		// Methods
		virtual void RecordEvent(const G4Event*);
// 		virtual void Merge(const G4Run*);

	private:
		DetectorConstruction* detector;
		PrimaryGeneratorAction* particleGun;
		G4int ID_eDep;
		G4int ID_kinEGamma;
		G4int ID_kinEElectron;
};

#endif