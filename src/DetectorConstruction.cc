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
// G4SiDetector.cc
//
// Description: Definition of the Canberra PD450-15-500AM Passivated
// Implanted Planar Silicon (PIPS) detector used by McMaster University
// to perform dosimetry measurements for the lens of the eye.
//
// NOTE1: McMaster is actually using the ORTEC CR-020-450-500 detector
// but the two models are essentially identical in terms of response.
//
// ********************************************************************

#include "DetectorConstruction.hh"
#include <cmath>

// Units and constants
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Manager classes
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

// Store classes
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// Geometry classes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// Primitive geometry types
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"

// Boolean operations on volumes
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

// Regions
#include "G4Region.hh"

// Messenger classes
#include "G4GenericMessenger.hh"

// Scoring Components
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4PSIncidentKineticEnergy.hh"
#include "G4SDParticleFilter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), fCheckOverlaps(true),
WorldPhysical(0)
{	
	// Geometry Parameters (Default)

	// LaBr3(Ce) Crystal
	fLaBr3Diameter = 50.8*mm;
	fLaBr3Length = 50.8*mm;

	// Housing
	fHousingThickness = 0.5*mm;

	// Reflector
	// The LaBr3(Ce) crystal is wrapped in a material acting as a Lambertian reflector surface
	// which is designed to increase light collection efficiency 
	fReflectorThickness = 2.0*mm;

	// Light Guide
	fLightGuideDiameter = fLaBr3Diameter;
	fLightGuideThickness = 5.0*mm;

	// Rotation Angle
	rotX = 0.0*deg;		

	// Source Radius
	sourceRadius = 20.*cm;
			 
	// Define Materials
	DefineMaterials();
	
	// Define commands to control the geometry
   	DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    
    // NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	nistManager->SetVerbose(0);

	// Define elements
	G4Element* elLa = nistManager->FindOrBuildElement(57);
	G4Element* elBr = nistManager->FindOrBuildElement(35);

  	// Define materials from elements
  	G4Material* LaBr3 = new G4Material("LaBr3", 5.06*g/cm3, 2);
  	LaBr3->AddElement(elLa, 0.366875);
  	LaBr3->AddElement(elBr ,0.633124);
	LaBr3->GetIonisation()->SetMeanExcitationEnergy(454.5*eV);
	  
  	// Set the materials for the Geometry
	fMatWorld = nistManager->FindOrBuildMaterial("G4_Galactic");
    fMatHousing = nistManager->FindOrBuildMaterial("G4_Al");
	fMatReflector = nistManager->FindOrBuildMaterial("G4_TEFLON");
    fMatLaBr3 = LaBr3;
    fMatLightGuide = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
  	
  	// Print materials
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 	
	// Cleanup old geometry
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	G4PhysicalVolumeStore::GetInstance()->Clean();
  	G4LogicalVolumeStore::GetInstance()->Clean();
  	G4SolidStore::GetInstance()->Clean();	
  	
	////////////////////////////////////////////////////////////////////////
	// Construct The World Volume

	G4double world_X = 2*(sourceRadius + 1.*cm);
	G4double world_Y = world_X;
	G4double world_Z = world_X;
	
	G4Box* WorldSolid = new G4Box("World", world_X/2, world_Y/2, world_Z/2);
  
	WorldLogical = 
		new G4LogicalVolume(WorldSolid,						// The Solid
							fMatWorld,						// Material
							"World");						// Name
  
	WorldPhysical = 
		new G4PVPlacement(	0,								// Rotation
							G4ThreeVector(),				// Translation vector
							WorldLogical,					// Logical volume
							"World",						// Name
							0,								// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
							

	////////////////////////////////////////////////////////////////////////
	// Construct the Source sphere
	// Note: The actual radius of the Source solid will be slightly smaller (0.1 mm) than
	// specified in the macro files in order to allow tracking the incident kinetic energy
	// of particles.
	G4VSolid* SourceSolid = new G4Sphere("SourceSolid", 0., sourceRadius/2, 0., 360.0*degree, 0., 180.0*degree);

	SourceLogical = 
		new G4LogicalVolume(SourceSolid,						// The Solid
							fMatWorld,		    				// Material
							"SourceLogical");	     			// Name

	SourcePhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,0),
							SourceLogical,					// Logical volume
							"SourcePhysical",				// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
							
	////////////////////////////////////////////////////////////////////////
	// Construct the detector housing	
	
	G4VSolid* HousingSolid = new G4Tubs("Housing_Solid",
										0., 
										fLaBr3Diameter/2 + fReflectorThickness + fHousingThickness, 
										(fLaBr3Length + fReflectorThickness + fLightGuideThickness + fHousingThickness)/2,
										0.,
										360.*deg);
	
	HousingLogical = 
		new G4LogicalVolume(HousingSolid,					// The Solid
							fMatHousing,    				// Material
							"Housing_Logical");     		// Name

	// Rotation, Translation, and Transformation of the detector			
	G4RotationMatrix Housing_Rot = G4RotationMatrix();
	Housing_Rot.rotateX(rotX);
	
	G4ThreeVector Housing_Trans = G4ThreeVector(0,0,0);
			
	HousingPhysical = 
		new G4PVPlacement(	G4Transform3D(Housing_Rot,Housing_Trans),	// Translation
							HousingLogical,					// Logical volume
							"Housing_Physical",		        // Name
							SourceLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the reflector
	G4VSolid* ReflectorSolid = new G4Tubs("Reflector_Solid",
										0., 
										fLaBr3Diameter/2 + fReflectorThickness, 
										(fLaBr3Length + fReflectorThickness + fLightGuideThickness)/2,
										0.,
										360.*deg);

	ReflectorLogical = 
		new G4LogicalVolume(ReflectorSolid,					// The Solid
							fMatReflector,		    		// Material
							"Reflector_Logical");	     	// Name
	
	ReflectorPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,-fHousingThickness/2),
							ReflectorLogical,				// Logical volume
							"Reflector_Physical",			// Name
							HousingLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the LaBr3 crystal
	G4VSolid* LaBr3Solid = new G4Tubs("LaBr3_Solid",
										0., 
										fLaBr3Diameter/2, 
										fLaBr3Length/2,
										0.,
										360.*deg);
	LaBr3Logical = 
		new G4LogicalVolume(LaBr3Solid,						// The Solid
							fMatLaBr3,    					// Material
							"LaBr3_Logical");     			// Name

	LaBr3Physical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(fLightGuideThickness-fReflectorThickness)/2),
							LaBr3Logical,					// Logical volume
							"LaBr3Physical",				// Name
							ReflectorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	// Create a region for the LaBr3 crystal so we can apply the PAI model to it
	G4Region* regLaBr3 = new G4Region("Region_LaBr3");
  	LaBr3Logical->SetRegion(regLaBr3);
	regLaBr3->AddRootLogicalVolume(LaBr3Logical);
	
	////////////////////////////////////////////////////////////////////////
	// Construct the light guide
	G4VSolid* LightGuideSolid = new G4Tubs("LightGuide_Solid",
										0., 
										fLightGuideDiameter/2, 
										fLightGuideThickness/2,
										0.,
										360.*deg);

	LightGuideLogical = 
		new G4LogicalVolume(LightGuideSolid,				// The Solid
							fMatLightGuide,    				// Material
							"LightGuideLogical");     		// Name

	LightGuidePhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,-(fLaBr3Length + fReflectorThickness)/2),
							LightGuideLogical,				// Logical volume
							"LightGuidePhysical",			// Name
							ReflectorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
  	// Visualisation attributes
  	
  	// World Volume (White)
  	G4VisAttributes* Vis_World = new G4VisAttributes(G4Colour(1.,1.,1.,0.1));
  	Vis_World->SetForceWireframe(true);
  	WorldLogical->SetVisAttributes(Vis_World);

	// Source Volume (Light Yellow)
    G4VisAttributes* Vis_Source = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
    Vis_Source->SetForceWireframe(true);
    SourceLogical->SetVisAttributes(Vis_Source);

    // Housing Volume (Gray)
    G4VisAttributes* Vis_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,.2));
    Vis_Housing->SetForceWireframe(false);
    HousingLogical->SetVisAttributes(Vis_Housing);

	// Reflector (Cyan)
    G4VisAttributes* Vis_Reflector = new G4VisAttributes(G4Colour(0.,1.0,1.0,0.2));
    Vis_Reflector->SetForceWireframe(false);
    ReflectorLogical->SetVisAttributes(Vis_Reflector);

	// LaBr3 crystal (Magenta)
    G4VisAttributes* Vis_LaBr3 = new G4VisAttributes(G4Colour(1.,0.,1.,1.));
    Vis_LaBr3->SetForceWireframe(false);
    LaBr3Logical->SetVisAttributes(Vis_LaBr3);

	// Light Guide (Blue)
    G4VisAttributes* Vis_LightGuide = new G4VisAttributes(G4Colour(0.,0.,1.,0.3));
    Vis_LightGuide->SetForceWireframe(false);
    LightGuideLogical->SetVisAttributes(Vis_LightGuide);

	////////////////////////////////////////////////////////////////////////
	// Return world volume
	return WorldPhysical; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	G4String filterName, particleName;
  
  	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  	G4SDParticleFilter* electronFilter = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
	
	////////////////////////////////////////////////////////////////////////
	// Construct the Multi Functional Detector for the Si chip
	
	G4MultiFunctionalDetector* LaBr3Scorer = new G4MultiFunctionalDetector("LaBr3");
	G4SDManager::GetSDMpointer()->AddNewDetector(LaBr3Scorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	LaBr3Logical->SetSensitiveDetector(LaBr3Scorer);
 	
 	G4PSEnergyDeposit* eDep = new G4PSEnergyDeposit("eDep");
    LaBr3Scorer->RegisterPrimitive(eDep);
	
	G4MultiFunctionalDetector* SourceScorer = new G4MultiFunctionalDetector("Source");
	G4SDManager::GetSDMpointer()->AddNewDetector(SourceScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	WorldLogical->SetSensitiveDetector(SourceScorer);

	G4VPrimitiveScorer* kinEGamma = new G4PSIncidentKineticEnergy("kinEGamma");
	kinEGamma->SetFilter(gammaFilter);
    SourceScorer->RegisterPrimitive(kinEGamma);

	G4VPrimitiveScorer* kinEElectron = new G4PSIncidentKineticEnergy("kinEElectron");
	kinEElectron->SetFilter(electronFilter);
    SourceScorer->RegisterPrimitive(kinEElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorAngle(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetDetectorAngle()", "G4LaBr3Detector", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		rotX = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceRadius(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetSourceRadius()", "G4LaBr3Detector", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		sourceRadius = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorAngle()
{
	// Return the detector angle
	return rotX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetSourceRadius()
{
	// Return the detector angle
	return sourceRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/G4LaBr3Detector/", "Geometry control");

    // Detector Angle Command
    G4GenericMessenger::Command& DetectorAngleCmd
      = fMessenger->DeclareMethodWithUnit("DetectorAngle","deg",
                                  &DetectorConstruction::SetDetectorAngle, 
                                  "Set the angle of the detector within the world volume.");
    DetectorAngleCmd.SetParameterName("angle", true);
    DetectorAngleCmd.SetDefaultValue("0.0");

	// Source Radius Command
	G4GenericMessenger::Command& SourceRadiusCmd
      = fMessenger->DeclareMethodWithUnit("SourceRadius","cm",
                                  &DetectorConstruction::SetSourceRadius, 
                                  "Set the radius of the source volume within the world volume.");
    SourceRadiusCmd.SetParameterName("radius", true);
    SourceRadiusCmd.SetDefaultValue("20.0");
}
