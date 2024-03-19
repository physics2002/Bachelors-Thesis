#include "run.hh"



MyRunAction::MyRunAction()
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    /*analysisManager->SetVerboseLevel(1);
    analysisManager->SetActivation(true);
    G4int ih = analysisManager->CreateH1("Z", "Z", 18, 0.9, 1.);
    analysisManager->SetH1Activation(ih, true);
    
    G4int ih1 = analysisManager->CreateH1("Y", "Y", 26, -4.021875, -3.865625);
    analysisManager->SetH1Activation(ih1, true);
    
    const G4Run *run = new G4Run();
    G4int runID = run->GetNumberOfEvent();
    G4cout << runID << G4endl;
    G4int ih3 = analysisManager->CreateH1("Nphotons", "Nphotons", runID, 0, 300);
    analysisManager->SetH1Activation(ih3, true);
    
    G4int ih2 = analysisManager->CreateH2("ZY", "ZY",18, 0.9, 1., 26, -4.021875, -3.865625);
    analysisManager->SetH2Activation(ih2, true); */
    for(G4int i=0; i<20; i++){
        std::stringstream layerID;
        layerID << i+1; 
        analysisManager->CreateNtuple("Position"+layerID.str(), "Position"+layerID.str());
        analysisManager->CreateNtupleDColumn(i, "fZX");
        analysisManager->CreateNtupleDColumn(i, "fY");
        analysisManager->CreateNtupleIColumn(i, "N");
        analysisManager->FinishNtuple(i);
    }

}

MyRunAction::~MyRunAction()
{}


void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    G4int runID = run->GetRunID();
    std::stringstream strRunID;
    strRunID << runID;
    analysisManager->OpenFile("output"+strRunID.str()+".root");    
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

}
