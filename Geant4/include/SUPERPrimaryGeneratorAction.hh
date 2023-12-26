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
/// \file SUPERPrimaryGeneratorAction.hh
/// \brief Definition of the SUPERPrimaryGeneratorAction class

#ifndef SUPERPrimaryGeneratorAction_h
#define SUPERPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <G4GeneralParticleSource.hh>
#include "TF1.h"
#include "TRandom3.h"


class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

/// Primary generator
///
/// A single particle is generated.
/// User can select 
/// - the initial momentum and angle
/// - the momentum and angle spreads
/// - random selection of a particle type from proton, kaon+, pi+, muon+, e+ 


class SUPERPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  SUPERPrimaryGeneratorAction();
  virtual ~SUPERPrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event*);

  
private:
  G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* fGeneralParticleSource;
  G4ParticleDefinition* fParticle;
  G4ParticleDefinition* fAl27;
  G4ParticleDefinition* fN14;
  G4double fMomentum;
  double fpi0mass;
  double N14mass;
  double Al27mass;
  //G4double pi0pdist(G4double *x,G4double *p);
  //G4double pi0Angulardist(G4double *x,G4double *p);
  TF1 *fpi0pdist = NULL;
  TF1 *fpi0Angulardist = NULL;;

  TF1 *fgammaEdist = NULL;
  TF1 *fgammaAngulardist = NULL;;

  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
