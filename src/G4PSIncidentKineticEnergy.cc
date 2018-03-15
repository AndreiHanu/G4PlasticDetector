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
// G4PSIncidentKineticEnergy.cc
//
// Description: This is a custom primitive scorer class for scoring the
// kinetic energy of a particle entering a volume.
// ********************************************************************
#include "G4PSIncidentKineticEnergy.hh"

#include "G4StepStatus.hh"
#include "G4VProcess.hh"
#include "G4UnitsTable.hh"

G4PSIncidentKineticEnergy::G4PSIncidentKineticEnergy(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit("MeV");
}

G4PSIncidentKineticEnergy::G4PSIncidentKineticEnergy(G4String name, const G4String& unit, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit(unit);
}

G4PSIncidentKineticEnergy::~G4PSIncidentKineticEnergy()
{;}

G4bool G4PSIncidentKineticEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    // The kinetic energy of a source particle will always be the largest
    // kinetic energy ever recorded by this scorer. So instead of checking
    // if the particle is a primary followed by checking for a volume boundary,
    // it's better to compare the kinetic energies during each step.

    // Check if this is the primary particle
    if ( aStep->GetTrack()->GetParentID() != 0 ) {
        // It's a secondary particles, check if it came from radioactive decay
       if( aStep->GetTrack()->GetCreatorProcess()->GetProcessName() != "RadioactiveDecay") return FALSE;
    }

    // Kinetic energy of this particle at the starting point.
    G4double kineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();

    // Get the current stored value in the event map
    G4int index = GetIndex(aStep);
    G4double* mapValue = ((*EvtMap)[index]);

    // If mapValue exits (e.g not NULL ), compare it with the current kinetic energy.
    if ( mapValue && ( *mapValue > kineticEnergy ) ) return FALSE ;
    else EvtMap->set(index, kineticEnergy);

    return TRUE;
}

void G4PSIncidentKineticEnergy::Initialize(G4HCofThisEvent* HCE)
{
    EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if ( HCID < 0 ) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSIncidentKineticEnergy::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSIncidentKineticEnergy::clear(){
    EvtMap->clear();
}

void G4PSIncidentKineticEnergy::DrawAll()
{;}

void G4PSIncidentKineticEnergy::PrintAll()
{
    G4cout << " MultiFunctionalDet " << detector->GetName() << G4endl;
    G4cout << " PrimitiveScorer " << GetName() << G4endl;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl;
    std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
    for(; itr != EvtMap->GetMap()->end(); itr++) {
        G4cout << "  copy no.: " << itr->first
	      << " Kinetic Energy: " << *(itr->second)/GetUnitValue()
	      << " ["<<GetUnit()<<"]"
	      << G4endl;
    }
}


void G4PSIncidentKineticEnergy::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Energy");
}

