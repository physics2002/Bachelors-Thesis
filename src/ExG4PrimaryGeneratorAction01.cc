#include "ExG4PrimaryGeneratorAction01.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>


int iRun = 0;


ExG4PrimaryGeneratorAction01::ExG4PrimaryGeneratorAction01()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{

  fParticleGun  = new G4ParticleGun(1);
    
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  G4ParticleDefinition* particle
    = particleTable->FindParticle("proton");
  
  fParticleGun->SetParticleDefinition(particle);

  
}

ExG4PrimaryGeneratorAction01::~ExG4PrimaryGeneratorAction01()
{
  delete fParticleGun;
}


void ExG4PrimaryGeneratorAction01::GeneratePrimaries(G4Event* anEvent)
{
  G4double x0 = 0 *cm;
  G4double y0 = 25 *cm;
  G4double z0 = 13.1885 *cm;
  //G4double y0 = 3.0 *cm;
  //G4double z0 = 4.0 *cm;
  //G4double pi = 3.141592653589793;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1., -0.15838444));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-cos(5*pi/180)/sin(5*pi/180),-1., -cos(5*pi/180)/sin(5*pi/180)));

  fParticleGun->SetParticleEnergy(100*GeV);

  fParticleGun->GeneratePrimaryVertex(anEvent);
  

    G4cout << iRun << G4endl;
    iRun++;
}
