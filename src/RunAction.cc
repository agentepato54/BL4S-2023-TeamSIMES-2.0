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
// $Id: RunAction.cc 74204 2013-10-01 07:04:43Z ihrivnac $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"
#include "Run.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{
  
  //=================================
  // Exercise 3 Step 1:
  // Create an output file containing
  // the histograms and the ntuple
  // Histograms:
  // 1D -> Drift Chamber 1 containing number of hits (50 bins from 0 to 50)
  // 1D -> Drift Chamber 2 containing number of hits (50 bins from 0 to 50)
  // 2D -> Drift Chamber 1 X vs Y binning: 50[-1000,1000] x 50[-300,300]
  // 2D -> Drift Chamber 2 X vs Y binning: 50[-1000,1000] x 50[-300,300]
  // Ntuple:
  // Integer Column: Number Hits Drift Chamber 1
  // Integer Column: Number Hits Drift Chamber 2
  // Double Column : Energy in EM Calo
  // Double Column : Energy in HAD Calo
  // Double Column : Time on hodoscope 1
  // Double Column : Time on hodoscope 2
    
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("Tutorial");

  // Creating 2D histograms
  analysisManager                                                
    ->CreateH2("CHC1_XY","Em Cal Colum 1 X vs Y",           // h2 Id = 0
               50, -1000., 1000, 50, -300., 300.); 
  analysisManager                                                
    ->CreateH2("CHC2_XY","Em Cal Colum 2 X vs Y",           // h2 Id = 1
               50, -1500., 1500, 50, -300., 300.);


  // Creating ntuple
  //
  analysisManager->CreateNtuple("Tutorial", "Hits");
  analysisManager->CreateNtupleDColumn("SC1Energy"); // column Id = 2
  analysisManager->CreateNtupleDColumn("EmEnergy"); // column Id = 3
  analysisManager->CreateNtupleDColumn("SC2Energy");    // column Id = 4
  analysisManager->CreateNtupleDColumn("EmC11Energy");    // column Id = 4
  analysisManager->CreateNtupleDColumn("EmC12Energy"); 
  analysisManager->CreateNtupleDColumn("EmC13Energy");    // column Id = 4
  analysisManager->CreateNtupleDColumn("EmC21Energy"); 
  analysisManager->CreateNtupleDColumn("EmC22Energy"); 
  analysisManager->CreateNtupleDColumn("EmC23Energy"); 
  analysisManager->CreateNtupleDColumn("EmC31Energy"); 
  analysisManager->CreateNtupleDColumn("EmC32Energy"); 
  analysisManager->CreateNtupleDColumn("EmC33Energy"); 
  analysisManager->CreateNtupleDColumn("EmTotalEnergy");
  analysisManager->CreateNtupleDColumn("RealTotalEnergy");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{

    
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //=================================
  // Exercise 3 Step 2:
  // Open output file at each new run

    // Open an output file
  // The default file name is set in RunAction::RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 2:
// Collect the energy accumulated in the local Run
// And dump on screen the needed physics quantities
// for this particular run.
void RunAction::EndOfRunAction(const G4Run* run)
{
    const Run* myrun = dynamic_cast<const Run*>(run);
    if ( myrun )
    {
        G4int nEvents = myrun->GetNumberOfEvent();
        if ( nEvents < 1 )
        {
            G4ExceptionDescription msg;
            msg << "Run consists of 0 events";
            G4Exception("RunAction::EndOfRunAction()",
                        "Code001", JustWarning, msg);
        }
        G4double em_ene = myrun->GetEmEnergy();
        G4double had_ene = myrun->GetHadEnergy();
        G4double shower_shape = myrun->GetShowerShape();
	G4int safety = ( nEvents > 0 ? nEvents : 1);//To avoid divisions by zero
    } else {
        G4ExceptionDescription msg;
        msg << "Run is not of correct type, skypping analysis via RunAction";
        G4Exception("RunAction::EndOfRunAction()",
                    "Code001", JustWarning, msg);
    }
    
  //=================================
  // Exercise 3 Step 3:
  // Write and close output file
    
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 1 Step 4:
// Implement this method,
// instead of creating a generic G4Run
// create a user-defined Run
G4Run* RunAction::GenerateRun() {
    return new Run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
