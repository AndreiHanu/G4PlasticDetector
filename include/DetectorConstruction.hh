#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
    public:
  	    // Constructor
        DetectorConstruction();
        // Destructor
        virtual ~DetectorConstruction();

        // Defines the detector geometry and returns a pointer to the physical World Volume
        virtual G4VPhysicalVolume* Construct();
    
        // Sensitive Detector 
	    virtual void ConstructSDandField();
    
        // Set Methods
        void SetDetectorAngle(G4double val);
        void SetSourceRadius(G4double val);
    
        // Get Methods
        G4double GetDetectorAngle();
        G4double GetSourceRadius();
    
    private:
        // Defines all the detector materials
        void DefineMaterials();
    
        // Define commands to change the geometry
        void DefineCommands();
    
        G4GenericMessenger* fMessenger;
        G4bool fCheckOverlaps;
    
        // Standard Materials
        G4Material* fMatWorld;
        G4Material* fMatDetHousing;
        G4Material* fMatEntranceWindow;
        G4Material* fMatEJ204;
        G4Material* fMatLightGuide;
        G4Material* fMatDetInterior;
        G4Material* fMatPMT;
        G4Material* fMatPMTInterior;
    
        // Logical Volumes
        G4LogicalVolume* WorldLogical;
        G4LogicalVolume* SourceLogical;
        G4LogicalVolume* DetHousingLogical;
        G4LogicalVolume* DetInteriorLogical;
        G4LogicalVolume* EntranceWindowLogical;
        G4LogicalVolume* ScintillatorLogical;
        G4LogicalVolume* LightGuideLogical;
        G4LogicalVolume* PMTLogical;
        G4LogicalVolume* PMTInteriorLogical;
    
        // Physical Volumes
        G4VPhysicalVolume* WorldPhysical;
        G4VPhysicalVolume* SourcePhysical;
        G4VPhysicalVolume* DetHousingPhysical;
        G4VPhysicalVolume* DetInteriorPhysical;
        G4VPhysicalVolume* EntranceWindowPhysical;
        G4VPhysicalVolume* ScintillatorPhysical;
        G4VPhysicalVolume* LightGuidePhysical;
        G4VPhysicalVolume* PMTPhysical;
        G4VPhysicalVolume* PMTInteriorPhysical;

	    // Rotation Angles
	    G4double rotX;
        G4double sourceRadius;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

