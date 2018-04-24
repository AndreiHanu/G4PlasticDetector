#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Select output format for Analysis Manager
#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* primary):G4Run(),
detector(det), particleGun(primary)
{
	G4SDManager* SDMan = G4SDManager::GetSDMpointer(); 
    	
    ID_eDep = SDMan->GetCollectionID("PlasticDetector/eDep");
	ID_kinEGamma = SDMan->GetCollectionID("Source/kinEGamma");
	ID_kinEElectron = SDMan->GetCollectionID("Source/kinEElectron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{ 	
  	// Get hits collections
  	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  	if(!HCE) { 
    	G4ExceptionDescription msg; 
    	msg << "No hits collection of this event found.\n"; 
    	G4Exception("Run::RecordEvent()","Code001", JustWarning, msg); 
    	return; 
  	} 
  	
	// Zero out the variables
	G4double eDep = 0.;
	G4double kinEGamma = 0.;
	G4double kinEElectron = 0.;
	
	// Get the HitMaps for this event
	G4THitsMap<G4double>* event_eDep = (G4THitsMap<G4double>*)(HCE->GetHC(ID_eDep));
	G4THitsMap<G4double>* event_kinEGamma = (G4THitsMap<G4double>*)(HCE->GetHC(ID_kinEGamma));
	G4THitsMap<G4double>* event_kinEElectron = (G4THitsMap<G4double>*)(HCE->GetHC(ID_kinEElectron));
	
	std::map<G4int,G4double*>::iterator itr;

	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// Calculate the fluence for this event taking into account any angular biasing
	G4double fluence = 1/(3.14159*std::pow(detector->GetSourceRadius()/cm, 2)*(std::pow(std::sin(particleGun->GetGPS()->GetCurrentSource()->GetAngDist()->GetMaxTheta()), 2)-std::pow(std::sin(particleGun->GetGPS()->GetCurrentSource()->GetAngDist()->GetMinTheta()), 2)));
	//G4double fluence = 1.;
	
	// Get the total energy deposited for this event
	for (itr = event_eDep->GetMap()->begin(); itr != event_eDep->GetMap()->end(); itr++) {
		eDep += *(itr->second);
	}

	// Get the incident kinetic energy for gammas in this event
	for (itr = event_kinEGamma->GetMap()->begin(); itr != event_kinEGamma->GetMap()->end(); itr++) {
		kinEGamma = *(itr->second);
		if (kinEGamma > 0){ 
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Gamma)"), kinEGamma/keV, fluence);
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Gamma) Linear"), kinEGamma/keV, fluence);
			if (eDep > 0){
				analysisManager->FillH1(analysisManager->GetH1Id("Detector True Spectrum (Gamma)"), kinEGamma/keV);
				analysisManager->FillH1(analysisManager->GetH1Id("Detector True Spectrum (Gamma) Linear"), kinEGamma/keV);
				analysisManager->FillH2(analysisManager->GetH2Id("Energy Migration Matrix (Gamma)"), kinEGamma/keV, eDep/keV);
				analysisManager->FillH2(analysisManager->GetH2Id("Energy Migration Matrix (Gamma) Linear"), kinEGamma/keV, eDep/keV);
			}
		}
	}

	// Get the incident kinetic energy for electrons in this event
	for (itr = event_kinEElectron->GetMap()->begin(); itr != event_kinEElectron->GetMap()->end(); itr++) {
		kinEElectron = *(itr->second);
		if (kinEElectron > 0){
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Electron)"), kinEElectron/keV, fluence);
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Electron) Linear"), kinEElectron/keV, fluence);
			if (eDep > 0){
				analysisManager->FillH1(analysisManager->GetH1Id("Detector True Spectrum (Electron)"), kinEElectron/keV);
				analysisManager->FillH1(analysisManager->GetH1Id("Detector True Spectrum (Electron) Linear"), kinEElectron/keV);
				analysisManager->FillH2(analysisManager->GetH2Id("Energy Migration Matrix (Electron)"), kinEElectron/keV, eDep/keV);
				analysisManager->FillH2(analysisManager->GetH2Id("Energy Migration Matrix (Electron) Linear"), kinEElectron/keV, eDep/keV);
			}
		}
	}

	// Record events with non-zero deposited energy
	if (eDep > 0) {
		analysisManager->FillH1(analysisManager->GetH1Id("Detector Measured Spectrum"), eDep/keV);
		analysisManager->FillH1(analysisManager->GetH1Id("Detector Measured Spectrum Linear"), eDep/keV);
	}
	
	// Invoke base class method
  	G4Run::RecordEvent(event); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void Run::Merge(const G4Run* aRun)
{
   	const Run* localRun = static_cast<const Run*>(aRun);
   	
   	//  Invoke base class method
   	G4Run::Merge(aRun); 
} 
*/
