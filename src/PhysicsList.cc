#include "PhysicsList.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmConfigurator.hh"
#include "G4EmParameters.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList(),
	fEmPhysicsList(0),
  	fDecayPhysicsList(0),
	fRadioactiveDecayPhysicsList(0)
{	
	// Default cut value
  	SetDefaultCutValue(0.01*mm);
  	
  	// Verbosity level
  	SetVerboseLevel(1);

  	// Decay Physics is always defined
  	fDecayPhysicsList = new G4DecayPhysics();
	
	// Radioactive Decay Physics
	fRadioactiveDecayPhysicsList = new G4RadioactiveDecayPhysics();

  	// EM physics
  	fEmPhysicsList = new G4EmStandardPhysics();
  	// fEmPhysicsList = new G4EmStandardPhysics_option3();
  	// fEmPhysicsList = new G4EmPenelopePhysics();
  	
  	// EM physics in gas cavity (Photoabsorption Ionization Model)
  	// fEmName = G4String("pai_photon");
  	fEmName = G4String("pai");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
	delete fDecayPhysicsList;
	delete fRadioactiveDecayPhysicsList;
	delete fEmPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	fDecayPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
	// Transportation Physics
	AddTransportation();
	
	// Electromagnetic Physics
	fEmPhysicsList->ConstructProcess();
	AddPAIModel(fEmName);
	
	// Decay Physics
	fDecayPhysicsList->ConstructProcess();

	// Radioactive Decay Physics
	fRadioactiveDecayPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
	if (verboseLevel>1) {
    	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  	}

  	if (name == fEmName) {
    	return;

	} else if (name == "emstandard_opt1") {

    	fEmName = name;
    	delete fEmPhysicsList;
    	fEmPhysicsList = new G4EmStandardPhysics_option1();

  	} else if (name == "emstandard_opt2") {

    	fEmName = name;
    	delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option2();

	} else if (name == "emstandard_opt3") {

    	fEmName = name;
    	delete fEmPhysicsList;
    	fEmPhysicsList = new G4EmStandardPhysics_option3();

	} else if (name == "emstandard_opt4") {

		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option4();

	} else if (name == "emlivermore") {

		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmLivermorePhysics();

	} else if (name == "empenelope") {

		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmPenelopePhysics();
		
	} else if (name == "pai") {

    	fEmName = name;
    	AddPAIModel(name);

  	} else if (name == "pai_photon") { 

    	fEmName = name;
    	AddPAIModel(name);

	} else {

		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
        << " is not defined"
        << G4endl;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{

	// Special threshold for low energy physics
  	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(990*eV, 10*GeV);

	// Default production thresholds for the world volume
 	SetCutsWithDefault();
 	 
 	if (verboseLevel > 0) { DumpCutValuesTable(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPAIModel(const G4String& modname)
{
	auto theParticleIterator=GetParticleIterator();
	theParticleIterator->reset();

  	while ((*theParticleIterator)())
  	{
    	G4ParticleDefinition* particle = theParticleIterator->value();
    	G4String partname = particle->GetParticleName();
    	if(partname == "e-" || partname == "e+") {
      		NewPAIModel(particle, modname, "eIoni");
    	} else if(partname == "mu-" || partname == "mu+") {
      		NewPAIModel(particle, modname, "muIoni");
    	}
  	}
  	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::NewPAIModel(const G4ParticleDefinition* part, 
                              const G4String& modname,
                              const G4String& procname)
{
	G4EmConfigurator* fConfig = G4LossTableManager::Instance()->EmConfigurator();
	
	G4String partname = part->GetParticleName();
  	if(modname == "pai") {
    	G4PAIModel* pai = new G4PAIModel(part,"PAIModel");
    	fConfig->SetExtraEmModel(partname,procname,pai,"", 0.0,100.*TeV,pai);
  	} else if(modname == "pai_photon") {
    	G4PAIPhotModel* pai = new G4PAIPhotModel(part,"PAIPhotModel");
    	fConfig->SetExtraEmModel(partname,procname,pai,"", 0.0,100.*TeV,pai);
  	}
}
