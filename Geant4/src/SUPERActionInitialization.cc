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
/// \file SUPERActionInitialization.cc
/// \brief Implementation of the SUPERActionInitialization class

#include "SUPERActionInitialization.hh"
#include "SUPERPrimaryGeneratorAction.hh"
#include "SUPERRunAction.hh"
#include "SUPERStackingAction.hh"
#include "SUPEREventAction.hh"


#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern bool gSaveStepLevel;

SUPERActionInitialization::SUPERActionInitialization(TTree* tr)
  : G4VUserActionInitialization(), fTree(tr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERActionInitialization::~SUPERActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERActionInitialization::BuildForMaster() const
{
  SetUserAction(new SUPERRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERActionInitialization::Build() const
{
  SetUserAction(new SUPERPrimaryGeneratorAction);
  auto runAction = new SUPERRunAction();
  auto eventAction = new SUPEREventAction(runAction,fTree);
  auto stackingAction = new SUPERStackingAction;
  SetUserAction(eventAction);
  SetUserAction(runAction);
  SetUserAction(stackingAction);
  eventAction -> SetRandomSeed(CLHEP::HepRandom::getTheSeed());
  eventAction -> SetSaveStepLevel(gSaveStepLevel);
  eventAction -> SetBranch();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
