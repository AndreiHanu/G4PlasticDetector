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
        G4Material* fMatHousing;
        G4Material* fMatReflector;
        G4Material* fMatLaBr3;
        G4Material* fMatLightGuide;
    
        // Logical Volumes
        G4LogicalVolume* WorldLogical;
        G4LogicalVolume* SourceLogical;
        G4LogicalVolume* HousingLogical;
        G4LogicalVolume* ReflectorLogical;
        G4LogicalVolume* LaBr3Logical;
        G4LogicalVolume* LightGuideLogical;
    
        // Physical Volumes
        G4VPhysicalVolume* WorldPhysical;
        G4VPhysicalVolume* SourcePhysical;
        G4VPhysicalVolume* HousingPhysical;
        G4VPhysicalVolume* ReflectorPhysical;
        G4VPhysicalVolume* LaBr3Physical;
        G4VPhysicalVolume* LightGuidePhysical;
    
        // Geometry Parameters
        G4double fLaBr3Diameter;
        G4double fLaBr3Length;
        G4double fHousingThickness;
        G4double fReflectorThickness;
        G4double fLightGuideDiameter;
        G4double fLightGuideThickness;

	    // Rotation Angles
	    G4double rotX;
        G4double sourceRadius;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

