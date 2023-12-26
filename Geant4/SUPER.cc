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
/// \file exampleSUPER.cc
/// \brief Main program of the analysis/SUPER example

#include "globals.hh"
#include "Randomize.hh"
#include "time.h"

// User Defined Detector
#include "SUPERDetectorConstruction.hh"
#include "SUPERPrimaryGeneratorAction.hh"
#include "SUPERDetectorConstruction.hh"
#include "SUPERActionInitialization.hh"

// Geant4
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4UImanager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// Physics list package
#include "FTFP_BERT.hh"

// Root
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TRandom3.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Global variables set by user
bool gSaveStepLevel = false;    // 
long gSeed = 0;                 // Random seed number. 0 for time seed 
bool gUseGPS = true;            //
bool gGenerateStepTheta;  //
bool gGenerateUniformPhi = true;  //
bool gGenearteUniformMomentum = false;
bool gGenearteUniformPosition = false;
int gGenGamma;
G4double gBeamMomentumMax;
G4double gBeamMomentumMin;
G4double gThetaLimitMin;  
G4double gThetaLimitMax;
G4double gGeneratePhi;
G4double gBeamMomentum;
G4String gParticle;
G4ThreeVector gPrimaryPosition;
G4double gPrimaryParticlePositionXmin;
G4double gPrimaryParticlePositionXmax;
G4double gPrimaryParticlePositionYmin;
G4double gPrimaryParticlePositionYmax;
G4double gPrimaryParticlePositionZ;

G4double gNsteps;
G4double gTheta_step;

int main(int argc,char** argv)
{
  if (argc != 1 && argc != 4){
    std::cout << "./exampleSUPER [1. Macro file name] [2. Output file name] [3. Random seed number]" << std::endl;
    std::cout << "ex) ./exampleSUPER run.mac example 1" << std::endl;
    return 0;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                      User Defined Parameters                                                //
  
  
  gSaveStepLevel = false;  // whether save all the step information or not
  gUseGPS = false;  // [true] : use General Particle Source described in your.mac file (/gps/position ...)   [false] : defined at bellow
  gGenerateStepTheta = false;
  gPrimaryPosition = G4ThreeVector(0,0,0); // Primary particle position  

  gParticle = "pi0"; // Particle name
  
  gGenGamma = 2; // 0: pi0,  1: photon for training,  2: two hard photons(brems) 
  
  //gGenGamma = true;

  if (argc == 1) gUseGPS = true; // use vis.mac
  if (argc == 4){
    //Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    gSeed = (long) atol(argv[3]);
    CLHEP::HepRandom::setTheSeed(gSeed);
    gRandom -> SetSeed(gSeed);
  }
  
  G4cout << "Save Step Level : " << gSaveStepLevel << G4endl;
  if (gUseGPS)
    G4cout << "Use GPS in .mac file" << G4endl;
  

  TString str_fname = argv[2];
  if (!str_fname.EndsWith(".root")){
    str_fname += ".root";
  }
  else if (str_fname == ""){
    str_fname = "VisMac.root";
  }

  
  TFile *tf = new TFile(str_fname,"RECREATE");
  auto tr = new TTree("tree","Geant4 output");
  tr -> SetAutoSave();
  auto physicsList = new FTFP_BERT();
  G4RunManager* runManager = new G4RunManager;
  runManager -> SetUserInitialization(physicsList);
  runManager -> SetUserInitialization(new SUPERActionInitialization(tr));
  runManager -> SetUserInitialization(new SUPERDetectorConstruction());
  //runManager -> SetUserAction(new SUPERPrimaryGeneratorAction());
  runManager -> Initialize();
  
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc != 1) {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager -> ApplyCommand(command+fileName);
  }
  else 
  {
    gUseGPS = true;
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager -> ApplyCommand("/control/execute vis.mac"); 
    ui -> SessionStart();
    
    delete ui;
  }

  tf -> cd();
  tr -> Write();
  tf -> Close();

  delete visManager;

  return 0;
}
