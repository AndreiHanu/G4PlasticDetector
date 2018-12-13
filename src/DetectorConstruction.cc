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
// Description: Definition of the Eljen Technologies M550-20x8-1 (P/N: 4467-2233)
// plastic scintillator detector (EJ-204) that is used by McMaster University to
// perform dosimetry measurements for the lens of eye.
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
#include "G4Polycone.hh"

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
#include "G4PSSphereSurfaceCurrent3D.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSSphereSurfaceFlux3D.hh"
#include "G4PSIncidentKineticEnergy.hh"
#include "G4SDParticleFilter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), fCheckOverlaps(true),
WorldPhysical(0)
{	
	// Rotation Angle
	rotX = 0.0*deg;		

	// Source Radius
	sourceRadius = 25.*cm;
			 
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
	  
  	// Set the materials for the Geometry
	fMatWorld = nistManager->FindOrBuildMaterial("G4_Galactic");
	//fMatWorld = nistManager->FindOrBuildMaterial("G4_AIR");
	fMatDetHousing = nistManager->FindOrBuildMaterial("G4_Al");
	fMatEntranceWindow = nistManager->FindOrBuildMaterial("G4_MYLAR");
	fMatEJ204 = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	fMatLightGuide = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	fMatDetInterior = nistManager->FindOrBuildMaterial("G4_AIR");
	fMatPMT = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
	fMatPMTInterior = nistManager->FindOrBuildMaterial("G4_Galactic");
  	
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
	// Note: The actual radius of the Source solid will be slightly smaller (1 mm) than
	// specified in the macro files in order to allow tracking the incident kinetic energy
	// of particles.
	// NOTE: Make sure the outer radius is smaller than the source radius otherwise you'll
	// get "track stuck" warnings during navigation.
	G4Sphere* SourceSolid = new G4Sphere("SourceSolid", sourceRadius - 1.5*mm, sourceRadius - 0.5*mm, 0., 360.0*degree, 0., 180.0*degree);
	SetSourceInnerRadius(SourceSolid->GetInnerRadius());

	SourceLogical = 
		new G4LogicalVolume(SourceSolid,					// The Solid
							fMatWorld,		    			// Material
							"SourceLogical");	     		// Name

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
	// numZPlanes
	// zPlane
	// rInner
	// rOuter
	G4double startPhi = 0.0*deg;
   	G4double endPhi = 360.0*deg;
   	G4int nrRZ = 8;
   	G4double zPlane[]={0, 11.05*mm, 11.05*mm, 37.77*mm, 37.77*mm, 58.34*mm, 58.34*mm, 240.28*mm};
   	G4double rInner[]={0, 0, 0, 0, 0, 0, 0, 0};
	G4double rOuter[]={76.2/2*mm, 76.2/2*mm, 54.36/2*mm, 54.36/2*mm, 76.2/2*mm, 76.2/2*mm, 60.40/2*mm, 60.40/2*mm};

   	G4VSolid* DetHousingSolid = new G4Polycone("Detector_Housing_Solid",
	   											startPhi,
												endPhi,
												nrRZ,
												zPlane,
												rInner,
												rOuter);
	
	DetHousingLogical = 
		new G4LogicalVolume(DetHousingSolid,				// The Solid
							fMatDetHousing,    				// Material
							"Detector_Housing_Logical");	// Name

	// Rotation, Translation, and Transformation of the detector			
	G4RotationMatrix Housing_Rot = G4RotationMatrix();
	Housing_Rot.rotateX(rotX);
	
	G4ThreeVector Housing_Trans = G4ThreeVector(0,0,0);
			
	DetHousingPhysical = 
		new G4PVPlacement(	G4Transform3D(Housing_Rot,Housing_Trans),	// Translation
							DetHousingLogical,					// Logical volume
							"Detector_Housing_Physical",	// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the detector interior (air)	
	// numZPlanes
	// zPlane
	// rInner
	// rOuter
	G4double startPhi2 = 0.0*deg;
   	G4double endPhi2 = 360.0*deg; 
	G4int nrRZ2 = 6;
   	G4double zPlane2[]={0, 4.7*mm, 4.7*mm, 44.12*mm, 44.12*mm, 236.28*mm};
   	G4double rInner2[]={0, 0, 0, 0, 0, 0};
	G4double rOuter2[]={55.245/2*mm, 55.245/2*mm, 50.8/2*mm, 50.8/2*mm, 58.4/2*mm, 58.4/2*mm};
   	

   	G4VSolid* DetInteriorSolid = new G4Polycone("Detector_Interior_Solid",
	   											startPhi2,
												endPhi2,
												nrRZ2,
												zPlane2,
												rInner2,
												rOuter2);
	
	DetInteriorLogical = 
		new G4LogicalVolume(DetInteriorSolid,				// The Solid
							fMatDetInterior,    			// Material
							"Detector_Interior_Logical");	// Name
			
	DetInteriorPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,0),
							DetInteriorLogical,				// Logical volume
							"Detector_Interior_Physical",	// Name
							DetHousingLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the Entrance Window
   	G4VSolid* EntranceWindowSolid = new G4Tubs("Entrance_Window_Solid",
	   										0.,
											50./2*mm,
											8.69/2*um,
											0.,
											360.*deg);

	EntranceWindowLogical = 
		new G4LogicalVolume(EntranceWindowSolid,			// The Solid
							fMatEntranceWindow,    			// Material
							"Entrance_Window_Logical");		// Name

	EntranceWindowPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,4.7*mm - 8.69/2*um),
							EntranceWindowLogical,			// Logical volume
							"Entrance_Window_Physical",		// Name
							DetInteriorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps); 				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the EJ-204 Scintillator
   	G4VSolid* ScintillatorSolid = new G4Tubs("Scintillator_Interior_Solid",
	   										0.,
											50./2*mm,
											20./2*mm,
											0.,
											360.*deg);

	ScintillatorLogical = 
		new G4LogicalVolume(ScintillatorSolid,					// The Solid
							fMatEJ204,    						// Material
							"Scintillator_Interior_Logical");	// Name

	ScintillatorPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,4.7*mm + 20./2*mm),
							ScintillatorLogical,			// Logical volume
							"Scintillator_Interior_Physical",		// Name
							DetInteriorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps); 				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the light guide
   	G4VSolid* LightGuideSolid = new G4Tubs("Light Guide_Solid",
	   										0.,
											50./2*mm,
											20./2*mm,
											0.,
											360.*deg);

	LightGuideLogical = 
		new G4LogicalVolume(LightGuideSolid,				// The Solid
							fMatLightGuide,    				// Material
							"Light_Guide_Logical");			// Name

	LightGuidePhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,4.7*mm + 20.*mm + 20./2*mm),
							LightGuideLogical,				// Logical volume
							"Light_Guide_Physical",			// Name
							DetInteriorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps); 				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the PMT (glass)
   	G4VSolid* PMTSolid = new G4Tubs("PMT_Solid",
	   								0.,
									52./2*mm,
									112./2*mm,
									0.,
									360.*deg);

	PMTLogical = 
		new G4LogicalVolume(PMTSolid,						// The Solid
							fMatPMT,    					// Material
							"PMT_Logical");					// Name

	PMTPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,4.7*mm + 20.*mm + 20.*mm + 112./2*mm),
							PMTLogical,						// Logical volume
							"PMT_Physical",					// Name
							DetInteriorLogical,				// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps); 				// Overlap Check
						
	////////////////////////////////////////////////////////////////////////
	// Construct the PMT Interior (air)
   	G4VSolid* PMTInteriorSolid = new G4Tubs("PMT_Interior_Solid",
	   								0.,
									46./2*mm,
									106./2*mm,
									0.,
									360.*deg);

	PMTInteriorLogical = 
		new G4LogicalVolume(PMTInteriorSolid,				// The Solid
							fMatPMTInterior,    			// Material
							"PMT_Interior_Logical");		// Name

	PMTInteriorPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,0),
							PMTInteriorLogical,				// Logical volume
							"PMT_Interior_Physical",		// Name
							PMTLogical,						// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps); 				// Overlap Check

	////////////////////////////////////////////////////////////////////////
  	// Visualisation attributes
  	
  	// World Volume (White)
  	G4VisAttributes* Vis_World = new G4VisAttributes(G4Colour(0.,0.,0.,0.));
  	Vis_World->SetForceWireframe(false);
  	//WorldLogical->SetVisAttributes(Vis_World);
	WorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
	  
	// Source Volume (Light Yellow)
    G4VisAttributes* Vis_Source = new G4VisAttributes(G4Colour(1.,1.,1.,0.1));
    Vis_Source->SetForceWireframe(false);
    SourceLogical->SetVisAttributes(Vis_Source);

    // Housing Volume (Gray)
    G4VisAttributes* Vis_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.4));
    Vis_Housing->SetForceWireframe(false);
    DetHousingLogical->SetVisAttributes(Vis_Housing);

	// Housing Interior (Cyan)
    G4VisAttributes* Vis_Housing_Interior = new G4VisAttributes(G4Colour(0.,1.0,0.,0.1));
    Vis_Housing_Interior->SetForceWireframe(false);
    DetInteriorLogical->SetVisAttributes(Vis_Housing_Interior);

	// Entrance Window
    G4VisAttributes* Vis_EntranceWindow_Interior = new G4VisAttributes(G4Colour(1.,0.,1.0,0.4));
    Vis_EntranceWindow_Interior->SetForceWireframe(false);
    EntranceWindowLogical->SetVisAttributes(Vis_EntranceWindow_Interior);

	// Scintillator Interior (Cyan)
    G4VisAttributes* Vis_Scintillator_Interior = new G4VisAttributes(G4Colour(0.,1.0,1.0,0.4));
    Vis_Scintillator_Interior->SetForceWireframe(false);
    ScintillatorLogical->SetVisAttributes(Vis_Scintillator_Interior);

	// Light Guide (Red)
    G4VisAttributes* Vis_LightGuide = new G4VisAttributes(G4Colour(1.,0.,0.,0.3));
    Vis_LightGuide->SetForceWireframe(false);
    LightGuideLogical->SetVisAttributes(Vis_LightGuide);
	
	// PMT Glass (Magenta)
    G4VisAttributes* Vis_PMT = new G4VisAttributes(G4Colour(1.,0.,1.,0.5));
    Vis_PMT->SetForceWireframe(false);
    PMTLogical->SetVisAttributes(Vis_PMT);

	// PMT Interior (White)
    G4VisAttributes* Vis_PMT_Interior = new G4VisAttributes(G4Colour(1.,1.,1.,0.5));
    Vis_PMT_Interior->SetForceWireframe(false);
    PMTInteriorLogical->SetVisAttributes(Vis_PMT_Interior);

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
	// Construct the Multi Functional Detector
	G4MultiFunctionalDetector* PlasticDetectorScorer = new G4MultiFunctionalDetector("PlasticDetector");
	G4SDManager::GetSDMpointer()->AddNewDetector(PlasticDetectorScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	ScintillatorLogical->SetSensitiveDetector(PlasticDetectorScorer);
 	
 	G4PSEnergyDeposit* eDep = new G4PSEnergyDeposit("eDep");
    PlasticDetectorScorer->RegisterPrimitive(eDep);
	
	G4MultiFunctionalDetector* SourceScorer = new G4MultiFunctionalDetector("Source");
	G4SDManager::GetSDMpointer()->AddNewDetector(SourceScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	SourceLogical->SetSensitiveDetector(SourceScorer);

	G4VPrimitiveScorer* kinEGamma = new G4PSIncidentKineticEnergy("kinEGamma", fCurrent_Out);
	kinEGamma->SetFilter(gammaFilter);
    SourceScorer->RegisterPrimitive(kinEGamma);

	G4VPrimitiveScorer* kinEElectron = new G4PSIncidentKineticEnergy("kinEElectron", fCurrent_Out);
	kinEElectron->SetFilter(electronFilter);
    SourceScorer->RegisterPrimitive(kinEElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorAngle(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetDetectorAngle()", "G4PlasticDetector", 
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
    	G4Exception ("DetectorConstruction::SetSourceRadius()", "G4PlasticDetector", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		sourceRadius = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceInnerRadius(G4double val)
{
	sourceInnerRadius = val;
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
	// Return the source radius
	return sourceRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetSourceInnerRadius()
{
	// Return the source inner radius
	return sourceInnerRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/G4PlasticDetector/", "Geometry control");

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
