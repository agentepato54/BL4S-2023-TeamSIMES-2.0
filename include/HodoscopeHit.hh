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
// $Id: HodoscopeHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file HodoscopeHit.hh
/// \brief Definition of the HodoscopeHit class

#ifndef HodoscopeHit_h
#define HodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hodoscope hit
///
/// It records:
/// - the strip ID
/// - the particle time
/// - the strip logical volume, its position and rotation

class HodoscopeHit : public G4VHit
{
public:
    HodoscopeHit(G4int i,G4double t);
    HodoscopeHit(const HodoscopeHit &right);
    virtual ~HodoscopeHit();

    const HodoscopeHit& operator=(const HodoscopeHit &right);
    int operator==(const HodoscopeHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    
    G4int GetID() const { return fId; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }
    
private:
    G4int fId;
    G4double fTime;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV;
};

typedef G4THitsCollection<HodoscopeHit> HodoscopeHitsCollection;

extern G4ThreadLocal G4Allocator<HodoscopeHit>* HodoscopeHitAllocator;

inline void* HodoscopeHit::operator new(size_t)
{
    if (!HodoscopeHitAllocator)
        HodoscopeHitAllocator = new G4Allocator<HodoscopeHit>;
    return (void*)HodoscopeHitAllocator->MallocSingle();
}

inline void HodoscopeHit::operator delete(void*aHit)
{
    HodoscopeHitAllocator->FreeSingle((HodoscopeHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
