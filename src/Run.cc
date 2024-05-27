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

#include "Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "HadCalorimeterHit.hh"
#include "EmCalorimeterHit.hh"


Run::Run() :
    G4Run()
{ }


//=================================
// Exercise 1 Step 2:
// Get the hits collection from the event
// For the EM and HAD calorimeter sensitive detectors
// Accumulate on them
// Do not forget to forward the call to the base class RecordEvent
// method

void Run::RecordEvent(const G4Event* evt)
{
}

//=================================
// Exercise 1 Step 3:
// Implement a Merge method.
// The meaning of this will become clear
// When Multi-threading is discussed.
// This method of a MASTER instance is called
// By all worker threads passing WORKER run as
// parameter. Get the energy from WORKER run and
// Accumulate in the MASTER
// Do not forget to forward the call to the base class
// method to accumulate the number of processed events!
void Run::Merge(const G4Run* aRun)
{
}