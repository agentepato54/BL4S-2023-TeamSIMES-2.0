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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{

  //=================================
  // Exercise 3 Step 2:
  // Open output file at each new run

    // Open an output file
  // The default file name is set in RunAction::RunAction(),
  // it can be overwritten in a macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 2:
// Collect the energy accumulated in the local Run
// And dump on screen the physics quantities
// for this particular run.
void RunAction::EndOfRunAction(const G4Run* run)
{
    
  //=================================
  // Exercise 3 Step 3:
  // Write and close output file
    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=================================
// Exercise 1 Step 4:
// Implement this method,
// instead of creating a generic G4Run
// create a user-defined Run
G4Run* RunAction::GenerateRun() {
    return new G4Run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
