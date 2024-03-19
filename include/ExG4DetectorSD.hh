#include <G4VSensitiveDetector.hh>

#ifndef ExG4DetectorSD_h
#define ExG4DetectorSD_h 1
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"

   class ExG4DetectorSD: public G4VSensitiveDetector
     {
     public:
       
        ExG4DetectorSD(G4String SDname);
        ~ExG4DetectorSD();       

     private:
        virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);


        
     };
#endif /* SENSITIVEDETECTOR */
