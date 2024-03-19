#include<G4Step.hh>
#include<fstream>
#include<iostream>
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "ExG4DetectorSD.hh"
#include <vector>
#include "run.hh"


#include<G4UserRunAction.hh>
#include<G4Run.hh>

ExG4DetectorSD::ExG4DetectorSD(G4String name):G4VSensitiveDetector(name)
{
}

G4bool ExG4DetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4Track *track = aStep->GetTrack();
    track->SetTrackStatus(fStopAndKill);
    
    G4TouchableHandle touchable = aStep->GetPreStepPoint()->GetTouchableHandle();
    G4VPhysicalVolume* volume = touchable->GetVolume();
    G4ThreeVector posDetector = volume->GetTranslation();
   
//    G4cout << copyNo << " " << posDetector.z()/cm << " " << posDetector.y()/cm<<G4endl;

    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
   /* analysisManager->FillH1(0, posDetector.z()/cm);
    analysisManager->FillH1(1, posDetector.y()/cm);
    analysisManager->FillH1(2, evt);
    analysisManager->FillH2(0, posDetector.z()/cm, posDetector.y()/cm);*/
    
    double dy = 0.017*cm;
    double Y0 = -80.1*cm;
    double Y  = -79.45*cm;
    double dY = Y - Y0;
    double h  = 10*cm;
    
    //first layer
    if(posDetector.y() > Y0 && posDetector.y() < dY+Y0){
        analysisManager->FillNtupleDColumn(0, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(0, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(0, 2, evt);
        analysisManager->AddNtupleRow(0);
    }
    
    //perpendicular first layer
    if(posDetector.y() > dY+Y0+dy && posDetector.y() < dY+Y0+dy+dY){
        analysisManager->FillNtupleDColumn(1, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(1, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(1, 2, evt);
        analysisManager->AddNtupleRow(1);
    }
    
    //second layer
    if(posDetector.y() > Y0 +h && posDetector.y() < dY+Y0+h){
        analysisManager->FillNtupleDColumn(2, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(2, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(2, 2, evt);
        analysisManager->AddNtupleRow(2);
    }
    //perpendicular second layer
    if(posDetector.y() > dY+Y0+dy+h && posDetector.y() < dY+Y0+dy+dY+h){
        analysisManager->FillNtupleDColumn(3, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(3, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(3, 2, evt);
        analysisManager->AddNtupleRow(3);
    }
    
    //third layer
    if(posDetector.y() > Y0 +2*h && posDetector.y() < dY+Y0+2*h){
        analysisManager->FillNtupleDColumn(4, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(4, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(4, 2, evt);
        analysisManager->AddNtupleRow(4);
    }
    //perpendicular third layer
    if(posDetector.y() > dY+Y0+dy+2*h && posDetector.y() < dY+Y0+dy+dY+2*h){
        analysisManager->FillNtupleDColumn(5, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(5, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(5, 2, evt);
        analysisManager->AddNtupleRow(5);
    }

    //fourth layer
    if(posDetector.y() > Y0 +3*h && posDetector.y() < dY+Y0+3*h){
        analysisManager->FillNtupleDColumn(6, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(6, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(6, 2, evt);
        analysisManager->AddNtupleRow(6);
    }
    //perpendicular fourth layer
    if(posDetector.y() > dY+Y0+dy+3*h && posDetector.y() < dY+Y0+dy+dY+3*h){
        analysisManager->FillNtupleDColumn(7, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(7, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(7, 2, evt);
        analysisManager->AddNtupleRow(7);
    }
    
    //fifth layer
    if(posDetector.y() > Y0 +4*h && posDetector.y() < dY+Y0+4*h){
        analysisManager->FillNtupleDColumn(8, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(8, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(8, 2, evt);
        analysisManager->AddNtupleRow(8);
    }
    //perpendicular fifth layer
    if(posDetector.y() > dY+Y0+dy+4*h && posDetector.y() < dY+Y0+dy+dY+4*h){
        analysisManager->FillNtupleDColumn(9, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(9, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(9, 2, evt);
        analysisManager->AddNtupleRow(9);
    }
    
    //sixth layer
    if(posDetector.y() > Y0 +5*h && posDetector.y() < dY+Y0+5*h){
        analysisManager->FillNtupleDColumn(10, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(10, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(10, 2, evt);
        analysisManager->AddNtupleRow(10);
    }
    //perpendicular sixth layer
    if(posDetector.y() > dY+Y0+dy+5*h && posDetector.y() < dY+Y0+dy+dY+5*h){
        analysisManager->FillNtupleDColumn(11, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(11, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(11, 2, evt);
        analysisManager->AddNtupleRow(11);
    }
    
    //seventh layer
    if(posDetector.y() > Y0 +6*h && posDetector.y() < dY+Y0+6*h){
        analysisManager->FillNtupleDColumn(12, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(12, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(12, 2, evt);
        analysisManager->AddNtupleRow(12);
    }
    //perpendicular seventh layer
    if(posDetector.y() > dY+Y0+dy+6*h && posDetector.y() < dY+Y0+dy+dY+6*h){
        analysisManager->FillNtupleDColumn(13, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(13, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(13, 2, evt);
        analysisManager->AddNtupleRow(13);
    }

    //eighth layer
    if(posDetector.y() > Y0 +7*h && posDetector.y() < dY+Y0+7*h){
        analysisManager->FillNtupleDColumn(14, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(14, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(14, 2, evt);
        analysisManager->AddNtupleRow(14);
    }
    //perpendicular eight layer
    if(posDetector.y() > dY+Y0+dy+7*h && posDetector.y() < dY+Y0+dy+dY+7*h){
        analysisManager->FillNtupleDColumn(15, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(15, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(15, 2, evt);
        analysisManager->AddNtupleRow(15);
    }
    
    //nineth layer
    if(posDetector.y() > Y0 +8*h && posDetector.y() < dY+Y0+8*h){
        analysisManager->FillNtupleDColumn(16, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(16, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(16, 2, evt);
        analysisManager->AddNtupleRow(16);
    }
    //perpendicular nineth layer
    if(posDetector.y() > dY+Y0+dy+8*h && posDetector.y() < dY+Y0+dy+dY+8*h){
        analysisManager->FillNtupleDColumn(17, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(17, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(17, 2, evt);
        analysisManager->AddNtupleRow(17);
    }
    
    //tenth layer
    if(posDetector.y() > Y0 +9*h && posDetector.y() < dY+Y0+9*h){
        analysisManager->FillNtupleDColumn(18, 0, posDetector.z()/cm);
        analysisManager->FillNtupleDColumn(18, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(18, 2, evt);
        analysisManager->AddNtupleRow(18);
    }
    //perpendicular tenth layer
    if(posDetector.y() > dY+Y0+dy+9*h && posDetector.y() < dY+Y0+dy+dY+9*h){
        analysisManager->FillNtupleDColumn(19, 0, posDetector.x()/cm);
        analysisManager->FillNtupleDColumn(19, 1, posDetector.y()/cm);
        analysisManager->FillNtupleIColumn(19, 2, evt);
        analysisManager->AddNtupleRow(19);
    }
    /*for(G4int i=0; i<20; i++){
        if(posDetector.y() > Y0 + i*h && posDetector.y() < dY+Y0 + i*h){
            analysisManager->FillNtupleDColumn(i, 0, posDetector.z()/cm);
            analysisManager->FillNtupleDColumn(i, 1, posDetector.y()/cm);
            analysisManager->AddNtupleRow(i);
        }
        if(posDetector.y() > dY+Y0+dy+i*h && posDetector.y() < dY+Y0+dy+dY+i*h){
            analysisManager->FillNtupleDColumn(i, 1, posDetector.x()/cm);
            analysisManager->FillNtupleDColumn(i, 1, posDetector.y()/cm);
            analysisManager->AddNtupleRow(i);
        }
    }*/
    
    return true;
}


ExG4DetectorSD::~ExG4DetectorSD()
{}

