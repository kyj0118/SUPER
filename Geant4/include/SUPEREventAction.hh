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
/// \file SUPEREventAction.hh
/// \brief Definition of the SUPEREventAction class

#ifndef SUPEREventAction_h
#define SUPEREventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>
#include <array>
#include "SUPERRunAction.hh"

#include "TTree.h"
#include "TInterpreter.h"
#include "TSystem.h"

using namespace std;
/// Event action
const int kMaxScintillator = 1000;

struct EMHitStruct{
  int nhit;
  int nMaximumHits = kMaxScintillator;
  int one[kMaxScintillator];
  int cid[kMaxScintillator];
  int lid[kMaxScintillator];
  int segid[kMaxScintillator];
  double x[kMaxScintillator];
  double y[kMaxScintillator];
  double z[kMaxScintillator];
  double t[kMaxScintillator];
  double e[kMaxScintillator];
};

struct EventInfoStruct{
  int eventID;
  int runID;
  long randomSeed;
};

struct PrimaryParticleInfoStruct{
  double x,y,z;
  double px,py,pz;
  double p,m,e;
  int PDG;
};

struct SecondaryInfoStruct{
  double x[1000];
  double y[1000];
  double z[1000];
  double px[1000];
  double py[1000];
  double pz[1000];
  double p[1000];
  double m[1000];
  double e[1000];
  int id[1000];
  int PDG[1000];
  int nparticle;
};



class SUPEREventAction : public G4UserEventAction
{
public:
  SUPEREventAction(SUPERRunAction *runAction, TTree *tr);
  virtual ~SUPEREventAction();
  void SetBranch();
  void SetRunID(G4int RunID);
  void SetRandomSeed(long seed);
  void SetSaveStepLevel(bool flag){fSaveStepLevel = flag;}
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  
  EMHitStruct EMHit;
  EventInfoStruct EventInfo;
  PrimaryParticleInfoStruct PrimaryParticle;
  PrimaryParticleInfoStruct PrimaryParticle2; // for two hard photon generation
  SecondaryInfoStruct SecondaryParticle;
private:
  SUPERRunAction* fRunAction;  
  TTree *fTree;
  
  bool fSaveStepLevel;



  vector<vector<G4double>> fEMStepEdep;
  vector<vector<G4double>> fEMPreStepx;
  vector<vector<G4double>> fEMPreStepy;
  vector<vector<G4double>> fEMPreStepz;
  vector<vector<G4double>> fEMPreStept;
  vector<vector<G4double>> fEMPostStepx;
  vector<vector<G4double>> fEMPostStepy;
  vector<vector<G4double>> fEMPostStepz;
  vector<vector<G4double>> fEMPostStept;

  vector<vector<G4double>> fEMParticlePx;
  vector<vector<G4double>> fEMParticlePy;
  vector<vector<G4double>> fEMParticlePz;
  vector<vector<G4int> > fEMParticleTrackID;
  vector<vector<G4int> > fEMParticleParentID;
  vector<vector<G4double>> fEMParticleCharge;
  vector<vector<G4double>> fEMParticleMass;
  vector<vector<G4int> > fEMParticlePDGID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
