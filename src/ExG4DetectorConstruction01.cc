#include "ExG4DetectorConstruction01.hh"
#include "ExG4DetectorSD.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>

ExG4DetectorConstruction01::ExG4DetectorConstruction01()
: G4VUserDetectorConstruction()
{ }


ExG4DetectorConstruction01::~ExG4DetectorConstruction01()
{ }



G4VPhysicalVolume* ExG4DetectorConstruction01::Construct()
{  

  G4NistManager* nist = G4NistManager::Instance();
    

//  G4double det_sizeXY = 25*cm, det_sizeZ = 0.15*cm;


  G4Material* det_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  det_mat->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double energy[9] = {1.239841939*eV/0.10, 1.239841939*eV/0.20, 1.239841939*eV/0.30, 1.239841939*eV/0.40, 1.239841939*eV/0.50, 1.239841939*eV/0.60, 1.239841939*eV/0.70, 1.239841939*eV/0.80, 1.239841939*eV/0.90};
  G4double rindexPolystyrene[9] = {1.59,1.59,1.59,1.59,1.59,1.59,1.59,1.59,1.59};

  G4double absorptionLen[9] = {3.50*m,3.51*m,3.52*m,3.53*m,3.54*m,3.55*m,3.56*m,3.57*m,3.58*m};
  
  G4MaterialPropertiesTable *mptPolystyrene = new G4MaterialPropertiesTable();
  mptPolystyrene->AddProperty("RINDEX", energy, rindexPolystyrene, 9);
  mptPolystyrene->AddProperty("ABSLENGTH", energy, absorptionLen, 9);

  G4double ScintFast[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0};
  
  mptPolystyrene->AddProperty("SCINTILLATIONCOMPONENT1", energy, ScintFast, 9, true);

  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD", 8./keV);
  mptPolystyrene->AddConstProperty("RESOLUTIONSCALE", 1.);
  mptPolystyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.8*ns, true);
  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD1", 1.);

  det_mat->SetMaterialPropertiesTable(mptPolystyrene);

  G4bool checkOverlaps = false;

  G4double innerRadius = 0.*cm;
  G4double R = 0.0500*cm;
  G4double D = 0.1000*cm;
  G4double hx = 10*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  
  G4RotationMatrix rotm  = G4RotationMatrix();
  G4RotationMatrix rotm1  = G4RotationMatrix();
  rotm.rotateY(90*deg);
  rotm1.rotateY(0*deg);

  G4double world_sizeXY = 2*m;
  G4double world_sizeZ  = 2*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4double rindexWorld[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0};
  G4double absorptionLenWorld[9] = {0.03660256,0.03660256,0.03660256,0.03660256,0.03660256,0.03660256,0.03660256,0.03660256,0.03660256}; //for 1  mm diameter
  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("ABSLENGTH", energy, absorptionLenWorld, 9);
    mptWorld->AddProperty("RINDEX", energy, rindexWorld, 9);

    world_mat->SetMaterialPropertiesTable(mptWorld);

  // Создание объема для мира, определяется просто сама форма объема, берем параллелепипед
  G4Box* solidWorld =
    new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
    
    

  // Логический объем, здесь подключается материал, из которого сделан объем
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");            
                                             
                                             
                                             

  //Физический объем, а теперь наш логический объем помещаем в "ральный" мир
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);        
                     
  ///////////////////////////////////////////////////////////////////////////////////////////

  G4Tubs* solidDet11 = new G4Tubs("solidDet11", innerRadius, R*0.94, hx, startAngle, spanningAngle);

  logicDet11 = new G4LogicalVolume(solidDet11, det_mat, "logicDet11");


    mirrorSurface = new G4OpticalSurface("mirrorSurface");

    mirrorSurface->SetType(dielectric_dielectric);
    mirrorSurface->SetFinish(ground);
    mirrorSurface->SetModel(unified);

    G4double reflectivity[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0};

    G4MaterialPropertiesTable *mptMirror = new G4MaterialPropertiesTable();
    mptMirror->AddProperty("REFLECTIVITY", energy, reflectivity, 9);

    mirrorSurface->SetMaterialPropertiesTable(mptMirror);

  G4LogicalSkinSurface *skin = new G4LogicalSkinSurface("skin", logicWorld, mirrorSurface);  

  G4Element *elC = nist->FindOrBuildElement("C");
  G4Element *elH = nist->FindOrBuildElement("H");
  G4Element *elO = nist->FindOrBuildElement("O");
  G4Element *elF = nist->FindOrBuildElement("F");

  G4Material *C5H8O2 = new G4Material("C5H8O2", 1.190*g/cm3, 3);
  C5H8O2->AddElement(elC, 5);
  C5H8O2->AddElement(elH, 8);
  C5H8O2->AddElement(elO, 2);
  C5H8O2->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double rindexPMMA[9] = {1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49, 1.49};
  
  G4MaterialPropertiesTable *mptPMMA = new G4MaterialPropertiesTable();
  mptPMMA->AddProperty("RINDEX", energy, rindexPMMA, 9);

  C5H8O2->SetMaterialPropertiesTable(mptPMMA);

  G4Tubs* PMMA = new G4Tubs("PMMA", R*0.94, R*0.97, hx, startAngle, spanningAngle);
  G4LogicalVolume *logicPMMA = new G4LogicalVolume(PMMA, C5H8O2, "logicPMMA");

  G4Material *C6H7F3O2 = new G4Material("C6H7F3O2", 1.43*g/cm3, 4);
  C6H7F3O2->AddElement(elC, 6);
  C6H7F3O2->AddElement(elH, 7);
  C6H7F3O2->AddElement(elF, 3);
  C6H7F3O2->AddElement(elO, 2);
  C6H7F3O2->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double rindexFP[9] = {1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42, 1.42};
  
  G4MaterialPropertiesTable *mptFP = new G4MaterialPropertiesTable();
  mptFP->AddProperty("RINDEX", energy, rindexFP, 9);

  C6H7F3O2->SetMaterialPropertiesTable(mptFP);

  G4Tubs* FP = new G4Tubs("FP", R*0.97, R, hx, startAngle, spanningAngle);
  G4LogicalVolume *logicFP = new G4LogicalVolume(FP, C6H7F3O2, "logicFP");
  
  G4double DX = 0.00575*2*2*cm;
  G4double DY = 0.00625*2*2*cm;
  G4double DZ = 0.00575*2*2*cm;
  
  int N = 200;
  double H_0 = -0.4*world_sizeXY;
  
    /*new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,0,(0))), logicDet11, "physDet", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,0,(0))), logicPMMA, "physPMMA", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,0,(0))), logicFP, "physFP", logicWorld, false, 0, checkOverlaps);*/
    
    double DELTA_Z = 75*D;

////////// first layer

  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0,-DELTA_Z+(D + i*D))), logicDet11, "physDet", logicWorld, false, i, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0,-DELTA_Z+(D + i*D))), logicPMMA, "physPMMA", logicWorld, false, i, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0,-DELTA_Z+(D + i*D))), logicFP, "physFP", logicWorld, false, i, checkOverlaps);
  }
  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicDet11, "physDet", logicWorld, false, i+N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicPMMA, "physPMMA", logicWorld, false, i+N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicFP, "physFP", logicWorld, false, i+N, checkOverlaps);

  }  

  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*2*D,-DELTA_Z+(D + i*D))), logicDet11, "physDet", logicWorld, false, i+2*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*2*D,-DELTA_Z+(D + i*D))), logicPMMA, "physPMMA", logicWorld, false, i+2*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*2*D,-DELTA_Z+(D + i*D))), logicFP, "physFP", logicWorld, false, i+2*N, checkOverlaps);
  }
  
  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*3*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicDet11, "physDet", logicWorld, false, i+3*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*3*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicPMMA, "physPMMA", logicWorld, false, i+3*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*3*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicFP, "physFP", logicWorld, false, i+3*N, checkOverlaps);
  }

  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*4*D,-DELTA_Z+(D + i*D))), logicDet11, "physDet", logicWorld, false, i+4*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*4*D,-DELTA_Z+(D + i*D))), logicPMMA, "physPMMA", logicWorld, false, i+4*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*4*D,-DELTA_Z+(D + i*D))), logicFP, "physFP", logicWorld, false, i+4*N, checkOverlaps);
  }

  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*5*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicDet11, "physDet", logicWorld, false, i+5*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*5*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicPMMA, "physPMMA", logicWorld, false, i+5*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H_0+0.866025404*5*D,-DELTA_Z+(D + i*D) + 0.5*D)), logicFP, "physFP", logicWorld, false, i+5*N, checkOverlaps);
  }  
  
////////// second layer  
  
  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D, H_0+29*DZ, 25*D)), logicDet11, "physDet", logicWorld, false, i+6*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+6*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ,25*D)), logicFP, "physFP", logicWorld, false, i+6*N, checkOverlaps);
  }
  
  for(G4int i = 0; i < N; i++){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+0.866025404*D,25*D)), logicDet11, "physDet", logicWorld, false, i+7*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+0.866025404*D,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+7*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+0.866025404*D,25*D)), logicFP, "physFP", logicWorld, false, i+7*N, checkOverlaps);
  }
  
  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D, H_0+29*DZ+2*0.866025404*D, 25*D)), logicDet11, "physDet", logicWorld, false, i+8*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ+2*0.866025404*D,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+8*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ+2*0.866025404*D,25*D)), logicFP, "physFP", logicWorld, false, i+8*N, checkOverlaps);
  }
  
  for(G4int i = 0; i < N; i++){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+3*0.866025404*D,25*D)), logicDet11, "physDet", logicWorld, false, i+9*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+3*0.866025404*D,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+9*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+3*0.866025404*D,25*D)), logicFP, "physFP", logicWorld, false, i+9*N, checkOverlaps);
  }

  for(G4int i = 0; i < N; i++)
  {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D, H_0+29*DZ+4*0.866025404*D, 25*D)), logicDet11, "physDet", logicWorld, false, i+10*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ+4*0.866025404*D,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+10*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.5*D,H_0+29*DZ+4*0.866025404*D,25*D)), logicFP, "physFP", logicWorld, false, i+10*N, checkOverlaps);
  }
  
  for(G4int i = 0; i < N; i++){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+5*0.866025404*D,25*D)), logicDet11, "physDet", logicWorld, false, i+11*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+5*0.866025404*D,25*D)), logicPMMA, "physPMMA", logicWorld, false, i+11*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector((D + i*D)-hx-4.0*D,H_0+29*DZ+5*0.866025404*D,25*D)), logicFP, "physFP", logicWorld, false, i+11*N, checkOverlaps);
  }
 
////////// third layer
    
    double H = 10 * cm;
  
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }

////////// fourth layer                                   (D + i*D)-hx-4.0*D
    double H_4 = H + 29*DZ;
    
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_4+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
  
////////// fifth layer

    double H5 = 10 * cm;
  
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
////////// sixth layer                                   (D + i*D)-hx-4.0*D
    double H_6 = H + 29*DZ+H5;
    
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_6+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
  
////////// seventh layer
    double H7 = 10 * cm;
  
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }

////////// eighth layer                                   (D + i*D)-hx-4.0*D
    double H_8 = H + 29*DZ+H5+H5;
    
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_8+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
  
////////// nineth layer
    double H9 = 10 * cm;
  
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
  
////////// tenth layer                                   (D + i*D)-hx-4.0*D
    double H_10 = H + 29*DZ+H5+H5+H5;
    
    for(G4int i = 0; i < 6*N; i++)
  {
    if(static_cast<int>(floor(i/N)) % 2 == 0){
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
    else
    {
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
  
    for(G4int j = 0; j<10; j++){
        if(j % 2 == 0){ //11th 13th 15th 17th 19th layer j = 0, 2, 4, 6, 8
                for(G4int i = 0; i < 6*N; i++)
                {
                    if(static_cast<int>(floor(i/N)) % 2 == 0){
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
                    else
                    {
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0,H9*(j+2)/2+H9+H7+H5+H+H_0+ 0.866025404*D*floor((i)/N),-DELTA_Z+0.5*D+(D + (i-N*floor((i+1)/N))*D))), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
        }
        else{                                   //// double H_10 = H + 29*DZ+H5+H5+H5  j = 1, 3, 5, 7, 9 12th 14th 16th 18th 20th
                for(G4int i = 0; i < 6*N; i++)
                {
                    if(static_cast<int>(floor(i/N)) % 2 == 0){
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
                    else
                    {
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicDet11, "physDet", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicPMMA, "physPMMA", logicWorld, false, i+12*N, checkOverlaps);
                        new G4PVPlacement(G4Transform3D(rotm1, G4ThreeVector(-4.0*D-hx+0.5*D+(D + (i-N*floor((i+1)/N))*D),H5*(j+1)/2+H_10+H_0+ 0.866025404*D*floor((i)/N),25*D)), logicFP, "physFP", logicWorld, false, i+12*N, checkOverlaps);
    }
    
  }
        }
    }
 /* double wolfRamWidth = 0.35*mm;
  G4Material* target_mat = nist->FindOrBuildMaterial("G4_W");
  G4Box *solidWolframTarget = new G4Box("solidWolframTarget", hx, wolfRamWidth*0.5, 50*D);
  G4LogicalVolume *logicWolframTarget = new G4LogicalVolume(solidWolframTarget, target_mat, "logicWolframTarget");
  new G4PVPlacement(0, G4ThreeVector(0, 0, 25*D), logicWolframTarget, "physTarget", logicWorld, false, 0, checkOverlaps);
 */ 

  G4Box *solidDetector = new G4Box("solidDetector", 0.5*DX, 0.5*DY, 0.5*DZ);

  logicDetector = new G4LogicalVolume(solidDetector, world_mat, "logicDetector");
//first layer
    for(G4int i=0; i<26; i++)
    {
        for(G4int j=0; j<900; j++)
        {
            new G4PVPlacement(0, G4ThreeVector(-0.5*DX - hx, H_0+0.5*DY + i*DY - R - 2*DY, -DELTA_Z+0.5*DZ + j*DZ), logicDetector, "physDetector", logicWorld, false, j + 225*i,  checkOverlaps);
        }    
    }
    
    //G4cout << H_0+0.5*DY  - R - 2*DY-DY*0.5 << " -- " << H_0+0.5*DY + 25*DY - R - 2*DY+DY*0.5 << " look here mf " << G4endl;
    //G4cout << H_0+29*DZ+0.5*DY - R - 2*DY-0.5*DY - (H_0+0.5*DY + 25*DY - R - 2*DY+DY*0.5) << " this dy mf " << G4endl;
    //G4cout << -DELTA_Z+0.5*DZ + 899*DZ << " HERE " << G4endl;
    
//second layer
    
    double DELTA_X = 80*D;
    
    for(G4int i=0; i<26; i++)
    {
        for(G4int j=0; j<900; j++)
        {
            new G4PVPlacement(0, G4ThreeVector(-DELTA_X+0.5*DZ + j*DZ-25*D, H_0+29*DZ+0.5*DY + i*DY - R - 2*DY, -0.5*DX - hx+25*D), logicDetector, "physDetector", logicWorld, false, j + 225*i+5850,  checkOverlaps);
        }    
    }
//third layer
    for(G4int i=0; i<26; i++)
    {
        for(G4int j=0; j<900; j++)
        {
            new G4PVPlacement(0, G4ThreeVector(-0.5*DX - hx, H+H_0+0.5*DY + i*DY - R - 2*DY, -DELTA_Z+0.5*DZ + j*DZ), logicDetector, "physDetector", logicWorld, false, j + 225*i,  checkOverlaps);
        }    
    }

//fourth layer
    for(G4int i=0; i<26; i++)
    {
        for(G4int j=0; j<900; j++)
        {
            new G4PVPlacement(0, G4ThreeVector(-DELTA_X+0.5*DZ + j*DZ-25*D, H_4+H_0+0.5*DY + i*DY - R - 2*DY, -0.5*DX - hx+25*D), logicDetector, "physDetector", logicWorld, false, j + 225*i+5850,  checkOverlaps);
        }    
    }
    
    
    
    for(G4int k=0; k<16;k++){
        if(k % 2 == 0){ //5th 7th ... 17th 19th layer k = 0, 2, 4, 6, 8, 10, 12, 14
                for(G4int i=0; i<26; i++)
                {
                    for(G4int j=0; j<900; j++)
                    {
                        new G4PVPlacement(0, G4ThreeVector(-0.5*DX - hx, H*(k+2)/2+H+H_0+0.5*DY + i*DY - R - 2*DY, -DELTA_Z+0.5*DZ + j*DZ), logicDetector, "physDetector", logicWorld, false, j + 225*i,  checkOverlaps);
                    }    
                }
        }
        
        else{//6th 8th ... 16th 18th layer k = 1, 3, 5, 7, 9, 11, 13, 15
                for(G4int i=0; i<26; i++)
                {
                    for(G4int j=0; j<900; j++)
                    {
                        new G4PVPlacement(0, G4ThreeVector(-DELTA_X+0.5*DZ + j*DZ-25*D, H*(k+1)/2+H_4+H_0+0.5*DY + i*DY - R - 2*DY, -0.5*DX - hx+25*D), logicDetector, "physDetector", logicWorld, false, j + 225*i+5850,  checkOverlaps);
                    }    
                }
        }
    }
    
  return physWorld;
}




void ExG4DetectorConstruction01::ConstructSDandField()
{
    ExG4DetectorSD *sensDet = new ExG4DetectorSD("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}
