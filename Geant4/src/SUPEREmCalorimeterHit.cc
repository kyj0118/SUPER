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
//
/// \file SUPEREmCalorimeterHit.cc
/// \brief Implementation of the SUPEREmCalorimeterHit class

#include "SUPEREmCalorimeterHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<SUPEREmCalorimeterHit>* SUPEREmCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREmCalorimeterHit::SUPEREmCalorimeterHit()
  : G4VHit(), 
    fCellID(-1),fLayerID(-1), fSegmentID(-1), fEdep(0.), fPos(0.), fTime(0.), fPLogV(nullptr),
    fParticlePx(0), fParticlePy(0), fParticlePz(0), fParticleTrackID(0),fParticleParentID(0),
    fParticleCharge(0), fParticleMass(0), fParticlePDGID(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREmCalorimeterHit::SUPEREmCalorimeterHit(G4int cellID)
  : G4VHit(), 
    fCellID(cellID),fLayerID(-1), fSegmentID(-1), fEdep(0.), fPos(0.), fTime(0.), fPLogV(nullptr),
    fParticlePx(0), fParticlePy(0), fParticlePz(0), fParticleTrackID(0),fParticleParentID(0),
    fParticleCharge(0), fParticleMass(0), fParticlePDGID(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREmCalorimeterHit::~SUPEREmCalorimeterHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......







