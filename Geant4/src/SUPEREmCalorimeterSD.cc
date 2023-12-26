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
/// \file SUPEREmCalorimeterSD.cc
/// \brief Implementation of the SUPEREmCalorimeterSD class

#include "SUPEREmCalorimeterSD.hh"
#include "SUPEREmCalorimeterHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREmCalorimeterSD::SUPEREmCalorimeterSD(G4String name, G4int ctype)
  : G4VSensitiveDetector(name), fNameSD(name), fHitsCollection(nullptr), fHCID(-1), fCrystalType(ctype)
{
  TotalPhi = 48;
  fEdep.clear();
  fEweightedx.clear();
  fEweightedy.clear();
  fEweightedz.clear();
  fEweightedt.clear();
  
  fEdep.resize(TotalPhi,0);
  fEweightedx.resize(TotalPhi,0);
  fEweightedy.resize(TotalPhi,0);
  fEweightedz.resize(TotalPhi,0);
  fEweightedt.resize(TotalPhi,0);
  
  collectionName.insert("EMCalHitCollection"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREmCalorimeterSD::~SUPEREmCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPEREmCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new SUPEREmCalorimeterHitsCollection(fNameSD,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SUPEREmCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep != 0){
    auto prepoint = step -> GetPreStepPoint();
    auto postpoint = step -> GetPostStepPoint();
    auto physVolume = prepoint -> GetPhysicalVolume();

    auto touchable = prepoint -> GetTouchable();
    auto physVolumeMother = touchable -> GetVolume(1);
    
    //G4int cpid = physVolume -> GetCopyNo();
    G4int cpid_me = touchable -> GetCopyNumber(0);
    G4int cpid_mother = physVolumeMother -> GetCopyNo();
    G4int cpid = cpid_me;
    //G4int cpid = cpid_me+cpid_mother*25;
    //G4cout << physVolume -> GetName() << G4endl;
    //std::cout << cpid << std::endl;
    //G4cout << "me : " << cpid_me << ", mother : " <<  cpid_mother << G4endl;
    
    //G4cout << "vol1 : " << physVolume0 -> GetName() << G4endl;
    //G4cout << "vol2 : " << physVolume1 -> GetName() << G4endl;
    
    G4double prex = (prepoint -> GetPosition()).x();
    G4double prey = (prepoint -> GetPosition()).y();
    G4double prez = (prepoint -> GetPosition()).z();
    G4double pret = prepoint -> GetGlobalTime();
    
    G4double postx = (postpoint -> GetPosition()).x();
    G4double posty = (postpoint -> GetPosition()).y();
    G4double postz = (postpoint -> GetPosition()).z();
    G4double postt = prepoint -> GetGlobalTime();
    
    G4double x = (prex + postx)/2.0;
    G4double y = (prey + posty)/2.0;
    G4double z = (prez + postz)/2.0;
    G4double t = (pret + postt)/2.0;
    
    fEdep[cpid] += edep;
    fEweightedx[cpid] += x * edep;
    fEweightedy[cpid] += y * edep;
    fEweightedz[cpid] += z * edep;
    fEweightedt[cpid] += t * edep; 
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPEREmCalorimeterSD::EndOfEvent(G4HCofThisEvent* hce){
  for (int i = 0; i < TotalPhi; i++){
    if (fEdep[i] == 0) continue;
    
    fEweightedx[i]/=fEdep[i];
    fEweightedy[i]/=fEdep[i];
    fEweightedz[i]/=fEdep[i];
    fEweightedt[i]/=fEdep[i];
    
    fHitsCollection->insert(new SUPEREmCalorimeterHit(fHCID));
    G4int CurrentHitID = fHitsCollection->GetSize()-1;
    auto hit = (SUPEREmCalorimeterHit*) ((hce -> GetHC(fHCID)) -> GetHit(CurrentHitID));
    hit -> SetXYZTE(fEweightedx[i], fEweightedy[i], fEweightedz[i], fEweightedt[i], fEdep[i]);
    hit -> SetLayerID(fCrystalType);
    hit -> SetSegmentID(i);
    hit -> SetCellID(0);
    
    fEdep[i] = 0.0; fEweightedx[i] = 0.0; fEweightedy[i] = 0.0; fEweightedz[i] = 0.0; fEweightedt[i] = 0.0;
  }
}
