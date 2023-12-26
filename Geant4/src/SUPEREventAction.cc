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
/// \file SUPEREventAction.cc
/// \brief Implementation of the SUPEREventAction class

#include "SUPEREventAction.hh"
#include "SUPERRunAction.hh"
#include "SUPEREmCalorimeterHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "SUPERPrimaryGeneratorAction.hh"
#include "TClonesArray.h"
#include "TTree.h"
#include "TObject.h"
#include <Riostream.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4ThreeVector gPrimaryParticlePosition;
extern G4ThreeVector gPrimaryParticleMomentumDirection;
extern int gPrimaryParticlePDG;
extern double gPrimaryParticleEnergy;
extern double gPrimaryParticleMass;

extern G4ThreeVector gPrimaryParticlePosition2;
extern G4ThreeVector gPrimaryParticleMomentumDirection2;
extern int gPrimaryParticlePDG2;
extern double gPrimaryParticleEnergy2;
extern double gPrimaryParticleMass2;


extern vector<G4int> gTrackIDvector;
extern vector<G4int> gPDGEncodingvector;
extern vector<G4ThreeVector> gPositionvector;
extern vector<G4ThreeVector> gMomentumvector;
extern vector<G4double> gEtotvector;
extern vector<G4double> gMassvector;

extern int gGenGamma;

SUPEREventAction::SUPEREventAction(SUPERRunAction *runAction, TTree *tr)
  : G4UserEventAction(), fRunAction(runAction), fTree(tr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPEREventAction::~SUPEREventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPEREventAction::BeginOfEventAction(const G4Event*)
{
  SetRunID(fRunAction -> GetRunID());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPEREventAction::EndOfEventAction(const G4Event* event)
{
  

  EventInfo.eventID = event -> GetEventID();
  if (EventInfo.eventID %1000 == 0){
    G4cout << "event ID : " << EventInfo.eventID << G4endl;
  }
  auto hce = event -> GetHCofThisEvent();
  
  int iarrayEMHit = 0;
  for (int i = 0; i < hce -> GetCapacity(); i++){
    if (hce -> GetHC(i) -> GetSize() == 0) continue;
    G4String iHCName = hce -> GetHC(i) -> GetName();
    // EM Hit
    if (iHCName == "EMCalHitCollection"){
      int nEMHitsInLayer = hce -> GetHC(i) -> GetSize();
      for (int ih = 0; ih < nEMHitsInLayer; ih++){
        auto hit = (SUPEREmCalorimeterHit*) (hce -> GetHC(i) -> GetHit(ih));
	
        EMHit.cid[iarrayEMHit] = hit -> GetCellID();
        EMHit.lid[iarrayEMHit] = hit -> GetLayerID();
        EMHit.segid[iarrayEMHit] = hit -> GetSegmentID();
        double xx,yy,zz,tt,ee;
        hit -> GetXYZTE(xx,yy,zz,tt,ee);
	
        EMHit.one[iarrayEMHit] = 1;
        EMHit.x[iarrayEMHit] = xx;
        EMHit.y[iarrayEMHit] = yy;
        EMHit.z[iarrayEMHit] = zz;
        EMHit.t[iarrayEMHit] = tt;
        EMHit.e[iarrayEMHit] = ee;
	
        iarrayEMHit++;
      }
    }
  }
  
  EMHit.nhit = iarrayEMHit;
  
  PrimaryParticle.x = gPrimaryParticlePosition.getX();
  PrimaryParticle.y = gPrimaryParticlePosition.getY();
  PrimaryParticle.z = gPrimaryParticlePosition.getZ();
  
  PrimaryParticle.e = gPrimaryParticleEnergy + gPrimaryParticleMass;
  PrimaryParticle.m = gPrimaryParticleMass;
  PrimaryParticle.p = sqrt(PrimaryParticle.e*PrimaryParticle.e - PrimaryParticle.m*PrimaryParticle.m); // magnitude of momentum
  
  PrimaryParticle.px = gPrimaryParticleMomentumDirection.getX() * PrimaryParticle.p;
  PrimaryParticle.py = gPrimaryParticleMomentumDirection.getY() * PrimaryParticle.p;
  PrimaryParticle.pz = gPrimaryParticleMomentumDirection.getZ() * PrimaryParticle.p;

  // Second primary particle for two hard photons production
  if (gGenGamma == 2){
    PrimaryParticle2.x = gPrimaryParticlePosition2.getX();
    PrimaryParticle2.y = gPrimaryParticlePosition2.getY();
    PrimaryParticle2.z = gPrimaryParticlePosition2.getZ();
  
    PrimaryParticle2.e = gPrimaryParticleEnergy2 + gPrimaryParticleMass2;
    PrimaryParticle2.m = gPrimaryParticleMass2;
    PrimaryParticle2.p = sqrt(PrimaryParticle2.e*PrimaryParticle2.e - PrimaryParticle2.m*PrimaryParticle2.m); // magnitude of momentum
  
    PrimaryParticle2.px = gPrimaryParticleMomentumDirection2.getX() * PrimaryParticle2.p;
    PrimaryParticle2.py = gPrimaryParticleMomentumDirection2.getY() * PrimaryParticle2.p;
    PrimaryParticle2.pz = gPrimaryParticleMomentumDirection2.getZ() * PrimaryParticle2.p;
  }
  
  SecondaryParticle.nparticle = gTrackIDvector.size();
  for (int i = 0 ; i < gTrackIDvector.size(); i++){
    SecondaryParticle.id[i] = gTrackIDvector[i];
    SecondaryParticle.x[i] = gPositionvector[i].getX();
    SecondaryParticle.y[i] = gPositionvector[i].getY();
    SecondaryParticle.z[i] = gPositionvector[i].getZ();

    SecondaryParticle.e[i] = gEtotvector[i];
    SecondaryParticle.m[i] = gMassvector[i];
    SecondaryParticle.p[i] = sqrt(gEtotvector[i]*gEtotvector[i] - gMassvector[i]*gMassvector[i]);

    SecondaryParticle.px[i] = gMomentumvector[i].getX();
    SecondaryParticle.py[i] = gMomentumvector[i].getY();
    SecondaryParticle.pz[i] = gMomentumvector[i].getZ();
    SecondaryParticle.PDG[i] = gPDGEncodingvector[i];
  }

  fTree -> Fill();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPEREventAction::SetBranch(){
  fTree -> Branch("eventID",&EventInfo.eventID,"eventID/I");
  fTree -> Branch("runID",&EventInfo.runID,"runID/I");
  fTree -> Branch("randomSeed",&EventInfo.randomSeed,"randomSeed/L");
  
  fTree -> Branch("PrimaryParticle.x",&PrimaryParticle.x,"PrimaryParticle.x/D");
  fTree -> Branch("PrimaryParticle.y",&PrimaryParticle.y,"PrimaryParticle.y/D");
  fTree -> Branch("PrimaryParticle.z",&PrimaryParticle.z,"PrimaryParticle.z/D");

  fTree -> Branch("PrimaryParticle.px",&PrimaryParticle.px,"PrimaryParticle.px/D");
  fTree -> Branch("PrimaryParticle.py",&PrimaryParticle.py,"PrimaryParticle.py/D");
  fTree -> Branch("PrimaryParticle.pz",&PrimaryParticle.pz,"PrimaryParticle.pz/D");
  fTree -> Branch("PrimaryParticle.p",&PrimaryParticle.p,"PrimaryParticle.p/D");
  fTree -> Branch("PrimaryParticle.m",&PrimaryParticle.m,"PrimaryParticle.m/D");
  fTree -> Branch("PrimaryParticle.e",&PrimaryParticle.e,"PrimaryParticle.e/D");
  fTree -> Branch("PrimaryParticle.PDG",&PrimaryParticle.PDG,"PrimaryParticle.PDG/I");

  if (gGenGamma == 2){
    fTree -> Branch("PrimaryParticle2.x",&PrimaryParticle2.x,"PrimaryParticle2.x/D");
    fTree -> Branch("PrimaryParticle2.y",&PrimaryParticle2.y,"PrimaryParticle2.y/D");
    fTree -> Branch("PrimaryParticle2.z",&PrimaryParticle2.z,"PrimaryParticle2.z/D");

    fTree -> Branch("PrimaryParticle2.px",&PrimaryParticle2.px,"PrimaryParticle2.px/D");
    fTree -> Branch("PrimaryParticle2.py",&PrimaryParticle2.py,"PrimaryParticle2.py/D");
    fTree -> Branch("PrimaryParticle2.pz",&PrimaryParticle2.pz,"PrimaryParticle2.pz/D");
    fTree -> Branch("PrimaryParticle2.p",&PrimaryParticle2.p,"PrimaryParticle2.p/D");
    fTree -> Branch("PrimaryParticle2.m",&PrimaryParticle2.m,"PrimaryParticle2.m/D");
    fTree -> Branch("PrimaryParticle2.e",&PrimaryParticle2.e,"PrimaryParticle2.e/D");
    fTree -> Branch("PrimaryParticle2.PDG",&PrimaryParticle2.PDG,"PrimaryParticle2.PDG/I");
  }
  
  fTree -> Branch("nSecondary",&SecondaryParticle.nparticle,"nSecondary/I");
  fTree -> Branch("SecondaryParticle.px",SecondaryParticle.px,"SecondaryParticle.px[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.py",SecondaryParticle.py,"SecondaryParticle.py[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.pz",SecondaryParticle.pz,"SecondaryParticle.pz[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.p",SecondaryParticle.p,"SecondaryParticle.p[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.m",SecondaryParticle.m,"SecondaryParticle.m[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.e",SecondaryParticle.e,"SecondaryParticle.e[nSecondary]/D");
  fTree -> Branch("SecondaryParticle.PDG",SecondaryParticle.PDG,"SecondaryParticle.PDG[nSecondary]/I");
  fTree -> Branch("SecondaryParticle.id",SecondaryParticle.id,"SecondaryParticle.id[nSecondary]/I");

  fTree -> Branch("nEMHit",&EMHit.nhit,"nEMHit/I");
  fTree -> Branch("EMHit.one",EMHit.one,"EMHit.one[nEMHit]/I");
  fTree -> Branch("EMHit.CellID",EMHit.cid,"EMHit.CellID[nEMHit]/I");
  fTree -> Branch("EMHit.LayerID",EMHit.lid,"EMHit.LayerID[nEMHit]/I");
  fTree -> Branch("EMHit.SegmentID",EMHit.segid,"EMHit.SegmentID[nEMHit]/I");
  fTree -> Branch("EMHit.x",EMHit.x,"EMHit.x[nEMHit]/D");
  fTree -> Branch("EMHit.y",EMHit.y,"EMHit.y[nEMHit]/D");
  fTree -> Branch("EMHit.z",EMHit.z,"EMHit.z[nEMHit]/D");
  fTree -> Branch("EMHit.t",EMHit.t,"EMHit.t[nEMHit]/D");
  fTree -> Branch("EMHit.e",EMHit.e,"EMHit.e[nEMHit]/D");

}
void SUPEREventAction::SetRunID(G4int RunID){
  EventInfo.runID = RunID;
}
void SUPEREventAction::SetRandomSeed(long seed){
  EventInfo.randomSeed = seed;
}
