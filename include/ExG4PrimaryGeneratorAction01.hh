#ifndef B1PrimaryGeneratorAction_h
#define B1PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4AnalysisManager.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

class ExG4PrimaryGeneratorAction01 : public G4VUserPrimaryGeneratorAction
{
  public:
    ExG4PrimaryGeneratorAction01();
    virtual ~ExG4PrimaryGeneratorAction01();

    
    virtual void GeneratePrimaries(G4Event*);
    
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    
  
  private:
    G4ParticleGun*  fParticleGun; 
};
#endif
