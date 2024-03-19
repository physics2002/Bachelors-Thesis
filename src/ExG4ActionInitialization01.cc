#include "ExG4ActionInitialization01.hh"
#include "ExG4PrimaryGeneratorAction01.hh"




ExG4ActionInitialization01::ExG4ActionInitialization01()
 : G4VUserActionInitialization()
{}

ExG4ActionInitialization01::~ExG4ActionInitialization01()
{}

void ExG4ActionInitialization01::BuildForMaster() const
{
    MyRunAction *runAction = new MyRunAction();
    SetUserAction(runAction);
}

void ExG4ActionInitialization01::Build() const
{
SetUserAction(new ExG4PrimaryGeneratorAction01);
    
    MyRunAction *runAction = new MyRunAction();
    SetUserAction(runAction);
}
