#ifndef ExG4DetectorConstruction01_h
#define ExG4DetectorConstruction01_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "ExG4DetectorSD.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

class ExG4DetectorConstruction01 : public G4VUserDetectorConstruction
{
  public:
    
    ExG4DetectorConstruction01();
    ~ExG4DetectorConstruction01();
    
    virtual G4VPhysicalVolume* Construct();
    
private:
    G4LogicalVolume *logicDet11;
    G4LogicalVolume *logicDetector;
    virtual void ConstructSDandField();
    G4OpticalSurface *mirrorSurface;

};


#endif
