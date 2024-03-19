#ifndef ExG4ActionInitialization01_h
#define ExG4ActionInitialization01_h 1

#include "G4VUserActionInitialization.hh"
#include "run.hh"

class ExG4ActionInitialization01 : public G4VUserActionInitialization
{
  public:
    ExG4ActionInitialization01();
    virtual ~ExG4ActionInitialization01();
    virtual void Build() const;
    virtual void BuildForMaster() const;
};


#endif
