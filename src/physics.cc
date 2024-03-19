#include "physics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

MyPhysicsList::MyPhysicsList()
{
  SetVerboseLevel(1);
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4OpticalPhysics());
}

MyPhysicsList::~MyPhysicsList()
{}

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

void MyPhysicsList::SetCuts()
{ 
 // fixe lower limit for cut
 //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

 // call base class method to set cuts which default value can be
 // modified via /run/setCut/* commands
 
 
 G4VUserPhysicsList::SetCuts();

 DumpCutValuesTable();
}
