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
// kinetic energy of a particle entering or leaving a G4Sphere shape. 
//
// Surface is defined  at the inside of sphere.
// Direction                  -Rmin   +Rmax
//   0  IN || OUT            ->|<-     |
//   1  IN                   ->|       |
//   2  OUT                    |<-     |
//
// ********************************************************************
#include "G4PSIncidentKineticEnergy.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"

G4PSIncidentKineticEnergy::G4PSIncidentKineticEnergy(G4String name, G4int direction, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),EvtMap(0)
{
    SetUnit("MeV");
}

G4PSIncidentKineticEnergy::G4PSIncidentKineticEnergy(G4String name, G4int direction, const G4String& unit, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),EvtMap(0)
{
    SetUnit(unit);
}

G4PSIncidentKineticEnergy::~G4PSIncidentKineticEnergy()
{;}

G4bool G4PSIncidentKineticEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4StepPoint* preStep = aStep->GetPreStepPoint();
    G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
    G4VSolid* solid = physVol->GetLogicalVolume()->GetSolid();

    G4Sphere* sphereSolid = (G4Sphere*)(solid);

    G4int dirFlag = IsSelectedSurface(aStep,sphereSolid);
    if ( dirFlag > 0 ) {
        if ( fDirection == dirFlag ){
            // Kinetic energy of this particle at the starting point.
            G4double kineticEnergy = preStep->GetKineticEnergy();
            G4int index = GetIndex(aStep);
            EvtMap->add(index,kineticEnergy);
        }
    }

    return TRUE;
}

G4int G4PSIncidentKineticEnergy::IsSelectedSurface(G4Step* aStep, G4Sphere* sphereSolid)
{
    
    G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
    G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    
    if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
        // Entering Geometry
        G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
        G4ThreeVector localpos1 = theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
        G4double localR2 = localpos1.x()*localpos1.x()
                        +localpos1.y()*localpos1.y()
                        +localpos1.z()*localpos1.z();
        G4double InsideRadius = sphereSolid->GetInnerRadius();

        if ( localR2 > (InsideRadius-kCarTolerance)*(InsideRadius-kCarTolerance) && localR2 < (InsideRadius+kCarTolerance)*(InsideRadius+kCarTolerance)){
            return fCurrent_In;
        }
    }

    if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
        // Exiting Geometry
        G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
        G4ThreeVector localpos2 = theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
        G4double localR2 = localpos2.x()*localpos2.x()
                        +localpos2.y()*localpos2.y()
                        +localpos2.z()*localpos2.z();
        G4double InsideRadius = sphereSolid->GetInnerRadius();

        if ( localR2 > (InsideRadius-kCarTolerance)*(InsideRadius-kCarTolerance) && localR2 < (InsideRadius+kCarTolerance)*(InsideRadius+kCarTolerance)){
            return fCurrent_Out;
        }
    }

    return -1;
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

