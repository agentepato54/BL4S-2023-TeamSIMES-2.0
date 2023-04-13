//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EventAction.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "HodoscopeHit.hh"
#include "DriftChamberHit.hh"
#include "EmCalorimeterHit.hh"
#include "HadCalorimeterHit.hh"
#include "Analysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(), 
  fHHC1ID(-1),
  fHHC2ID(-1),
  fDHC1ID(-1),
  fDHC2ID(-1),
  fECHCID(-1),
  fCC1CID(-1),
  fCC2CID(-1),
  fCC3CID(-1),
  fHCHCID(-1)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    if (fHHC1ID==-1) {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      fHHC1ID = sdManager->GetCollectionID("hodoscope1/EMcalorimeterColl");
      fHHC2ID = sdManager->GetCollectionID("hodoscope2/EMcalorimeterColl");
      fDHC1ID = sdManager->GetCollectionID("chamber1/driftChamberColl");
      fDHC2ID = sdManager->GetCollectionID("chamber2/driftChamberColl");
      fECHCID = sdManager->GetCollectionID("EMcalorimeter/EMcalorimeterColl");

      fCC1CID = sdManager->GetCollectionID("EMcalorimeterCol1/EMcalorimeterColl");
      fCC2CID = sdManager->GetCollectionID("EMcalorimeterCol2/EMcalorimeterColl");
      fCC3CID = sdManager->GetCollectionID("EMcalorimeterCol3/EMcalorimeterColl");

      fHCHCID = sdManager->GetCollectionID("HadCalorimeter/HadCalorimeterColl");
      
    }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) 
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n"; 
        G4Exception("EventAction::EndOfEventAction()",
                    "Code001", JustWarning, msg);
        return;
    }   


    // Get hits collections 
    EmCalorimeterHitsCollection* hHC1 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fHHC1ID));
      
    EmCalorimeterHitsCollection* hHC2 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fHHC2ID));
      
    DriftChamberHitsCollection* dHC1 
      = static_cast<DriftChamberHitsCollection*>(hce->GetHC(fDHC1ID));
      
    DriftChamberHitsCollection* dHC2 
      = static_cast<DriftChamberHitsCollection*>(hce->GetHC(fDHC2ID));
      
    EmCalorimeterHitsCollection* ecHC 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fECHCID));


    EmCalorimeterHitsCollection* CHC1 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fCC1CID));

    EmCalorimeterHitsCollection* CHC2 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fCC2CID));

    EmCalorimeterHitsCollection* CHC3 
      = static_cast<EmCalorimeterHitsCollection*>(hce->GetHC(fCC3CID));

      
    HadCalorimeterHitsCollection* hcHC 
      = static_cast<HadCalorimeterHitsCollection*>(hce->GetHC(fHCHCID));
      
    if ( (!hHC1) || (!hHC2) || (!dHC1) || (!dHC2) || (!ecHC) || (!hcHC) || (!CHC1) || (!CHC2) || (!CHC3)) 
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found.\n"; 
        G4Exception("EventAction::EndOfEventAction()",
                    "Code001", JustWarning, msg);
        return;
    }   
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    G4int n_hit;

    for (G4int i=0;i<n_hit;i++)
    {
       DriftChamberHit* hit = (*dHC1)[i];
       G4ThreeVector localPos = hit->GetLocalPos();
       analysisManager->FillH2(0, localPos.x(), localPos.y());
    }
        
    // Fill ntuple
    
    G4int totalEmHit = 0;
    G4double totalEmE = 0.;
    
    G4double totalEmSC1 = 0.;
    G4double totalEmC = 0.;
    G4double totalEmSC2 = 0.;

    G4double totalEmCol1_1 = 0.;
    G4double totalEmCol1_2 = 0.;
    G4double totalEmCol1_3 = 0.;

    G4double totalEmCol2_1 = 0.;
    G4double totalEmCol2_2 = 0.;
    G4double totalEmCol2_3 = 0.;

    G4double totalEmCol3_1 = 0.;
    G4double totalEmCol3_2 = 0.;
    G4double totalEmCol3_3 = 0.;

    EmCalorimeterHit* hitSC1 = (*hHC1)[0];
    G4double eDepSC1 = hitSC1->GetEdep();
    totalEmSC1 += eDepSC1;
    totalEmE += eDepSC1;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmSC1);
    totalEmHit++;

    EmCalorimeterHit* hit = (*ecHC)[0];
    G4double eDepC = hit->GetEdep();
    totalEmC += eDepC;
    totalEmE += eDepC;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmC);
    totalEmHit++;

    EmCalorimeterHit* hitSC2 = (*hHC2)[0];
    G4double eDepSC2 = hitSC2->GetEdep();
    totalEmSC2 += eDepSC2;
    totalEmE += eDepSC2;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmSC2);
    totalEmHit++;

    EmCalorimeterHit* hitc11 = (*CHC1)[0];
    G4double eDepC1_1 = hitc11->GetEdep();
    totalEmCol1_1 += eDepC1_1;
    totalEmE += eDepC1_1;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol1_1);
    totalEmHit++;

    EmCalorimeterHit* hitc12 = (*CHC1)[1];
    G4double eDepC1_2 = hitc12->GetEdep();
    totalEmCol1_2 += eDepC1_2;
    totalEmE += eDepC1_2;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol1_2);
    totalEmHit++;

    EmCalorimeterHit* hitc13 = (*CHC1)[2];
    G4double eDepC1_3 = hitc13->GetEdep();
    totalEmCol1_3 += eDepC1_3;
    totalEmE += eDepC1_3;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol1_3);
    totalEmHit++;

    EmCalorimeterHit* hitc21 = (*CHC2)[0];
    G4double eDepC2_1 = hitc21->GetEdep();
    totalEmCol2_1 += eDepC2_1;
    totalEmE += eDepC2_1;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol2_1);
    totalEmHit++;

    EmCalorimeterHit* hitc22 = (*CHC2)[1];
    G4double eDepC2_2 = hitc22->GetEdep();
    totalEmCol2_2 += eDepC2_2;
    totalEmE += eDepC2_2;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol2_2);
    totalEmHit++;

    EmCalorimeterHit* hitc23 = (*CHC2)[2];
    G4double eDepC2_3 = hitc23->GetEdep();
    totalEmCol2_3 += eDepC2_3;
    totalEmE += eDepC2_3;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol2_3);
    totalEmHit++;

    EmCalorimeterHit* hitc31 = (*CHC3)[0];
    G4double eDepC3_1 = hitc31->GetEdep();
    totalEmCol3_1 += eDepC3_1;
    totalEmE += eDepC3_1;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol3_1);
    totalEmHit++;

    EmCalorimeterHit* hitc32 = (*CHC3)[1];
    G4double eDepC3_2 = hitc32->GetEdep();
    totalEmCol3_2 += eDepC3_2;
    totalEmE += eDepC3_2;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol3_2);
    totalEmHit++;

    EmCalorimeterHit* hitc33 = (*CHC3)[2];
    G4double eDepC3_3 = hitc33->GetEdep();
    totalEmCol3_3 += eDepC3_3;
    totalEmE += eDepC3_3;
    analysisManager->FillNtupleDColumn(totalEmHit, totalEmCol3_3);
    totalEmHit++;

    analysisManager->FillNtupleDColumn(totalEmHit, totalEmE);
    totalEmHit++;

    G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout << G4endl
           << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
           << primary->GetG4code()->GetParticleName()
           << " " << primary->GetMomentum() << G4endl;
    G4int p;
    p = primary->GetMomentum()[2];
    analysisManager->FillNtupleDColumn(totalEmHit, p);
    
      
    analysisManager->AddNtupleRow();  
    
   // ------------------------------------------------------------------------------

    G4cout << "Scintillator 1. Local Edep is " << eDepSC1/MeV << " (MeV)" << G4endl;
    G4cout << "EM Calorimeter (Single). Local Edep is " << eDepC/MeV << " (MeV)" << G4endl;
    G4cout << "Scintillator 2. Local Edep is " << eDepSC2/MeV << " (MeV)" << G4endl;

    for (G4int i=0;i<3;i++)
    {
        EmCalorimeterHit* hit = (*CHC1)[i];
        G4double eDep = hit->GetEdep();
        G4cout << "EM Calorimeter " << i+1 << ". Local Edep is " << eDep/MeV << " (MeV)" << G4endl;
    }

    for (G4int i=0;i<3;i++)
    {
        EmCalorimeterHit* hit = (*CHC2)[i];
        G4double eDep = hit->GetEdep();
        G4cout << "EM Calorimeter " << i+1 << ". Local Edep is " << eDep/MeV << " (MeV)" << G4endl;
    }

    for (G4int i=0;i<3;i++)
    {
        EmCalorimeterHit* hit = (*CHC3)[i];
        G4double eDep = hit->GetEdep();
        G4cout << "EM Calorimeter " << i+1 << ". Local Edep is " << eDep/MeV << " (MeV)" << G4endl;
    } 

    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
    

    // EM calorimeter
    G4cout << "EM Calorimeter has " << totalEmHit << " hits. Total Edep is "
    << totalEmE/MeV << " (MeV)" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
