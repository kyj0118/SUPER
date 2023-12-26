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
/// \file SUPERPrimaryGeneratorAction.cc
/// \brief Implementation of the SUPERPrimaryGeneratorAction class

#include "SUPERPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "TLorentzVector.h"
#include "Randomize.hh"
#include "TRandom3.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector gPrimaryParticlePosition;
G4ThreeVector gPrimaryParticleMomentumDirection;
int gPrimaryParticlePDG;
double gPrimaryParticleEnergy;
double gPrimaryParticleMass;

G4ThreeVector gPrimaryParticlePosition2;
G4ThreeVector gPrimaryParticleMomentumDirection2;
int gPrimaryParticlePDG2;
double gPrimaryParticleEnergy2;
double gPrimaryParticleMass2;


extern bool gUseGPS;
extern bool gGenerateStepTheta;
extern bool gGenerateUniformPhi;

extern bool gGenearteUniformPosition;
extern G4ThreeVector gPrimaryPosition;
extern G4double gPrimaryParticlePositionXmin;
extern G4double gPrimaryParticlePositionXmax;
extern G4double gPrimaryParticlePositionYmin;
extern G4double gPrimaryParticlePositionYmax;
extern G4double gPrimaryParticlePositionZ;

extern bool gGenearteUniformMomentum;
extern G4double gBeamMomentumMax;
extern G4double gBeamMomentumMin; 

extern int gGenGamma;

extern G4double gNsteps;
extern G4double gTheta_step;
extern G4double gThetaLimitMin;  
extern G4double gThetaLimitMax;
extern G4double gGeneratePhi;
extern G4double gBeamMomentum;

extern G4String gParticle;
using namespace std;
G4double gpi0mass = 0;

G4double pi0pdist(G4double *x,G4double *p){
  double Pt = x[0];
  double T0 = p[0];
  double pi0mass = p[1];
  double Et = sqrt(Pt*Pt + pi0mass*pi0mass);
  return Pt*sqrt(Et)*exp(-Et/T0);
}

G4double pi0Edist(G4double *x,G4double *p){
  G4double Ek = x[0];
  G4double pi0mass = p[0];
  double T0 = 10.0 * MeV;
  double Et = sqrt(Ek*Ek + pi0mass*pi0mass);
  return Et*sqrt(Et)*exp(-Et/T0);
}

G4double PhotonEdist(G4double *x,G4double *p){
  G4double Ek = x[0];
  double E0 = 30.0 * MeV;
  //double Et = sqrt(Ek*Ek + pi0mass*pi0mass);
  //return Et*sqrt(Et)*exp(-Et/T0);
  return exp(-Ek/E0);
}

G4double pi0Angulardist(G4double *x,G4double *p){
  double costh = x[0];
  double A2 = 0.31;
  double P2 = (3.0*costh*costh - 1.0)/2.0;
  return 1.0 + A2*P2;
}

SUPERPrimaryGeneratorAction::SUPERPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(nullptr),
  fGeneralParticleSource(nullptr), 
  fParticle(nullptr),
  fMomentum(1000.*MeV)
{
  fMomentum = gBeamMomentum;

  double au =  0.93149410242 * GeV;
  N14mass = 14.003074004251 * au;
  Al27mass = 26.98153841 * au;
  
  if(gUseGPS){
    fGeneralParticleSource  = new G4GeneralParticleSource();
  }
  else {
    if (gGenGamma == 0){
    auto particleTable = G4ParticleTable::GetParticleTable();
    auto ionTable = particleTable -> GetIonTable();
    fParticle = particleTable->FindParticle(gParticle);
    fParticleGun  = new G4ParticleGun(fParticle);
    fParticleGun->SetParticleMomentum(fMomentum);

    //fAl27 = particleTable->FindParticle("Al27");
    //fN14 = particleTable->FindParticle("N14");

    //fAl27 = ionTable->FindIon(13,27,0);
    //fN14 = ionTable->FindIon(7,14,0);

    //fAl27 = particleTable->GetIon(13,27,0);
    //fN14 = particleTable->GetIon(7,14,0);

    fpi0mass = particleTable->FindParticle("pi0")->GetPDGMass();
    
    fpi0pdist = new TF1("fpdist",pi0pdist,0,300,2);
    fpi0pdist -> SetParameter(0,10*MeV);
    fpi0pdist -> SetParameter(1,fpi0mass);

    fpi0Angulardist = new TF1("fangdist",pi0Angulardist,-1,1,0);
    }
    else {
      auto particleTable = G4ParticleTable::GetParticleTable();
      fParticle = particleTable->FindParticle("gamma");
      fParticleGun  = new G4ParticleGun(fParticle);
      if (gGenGamma == 2){
	fgammaEdist = new TF1("fgammaEdist",PhotonEdist,10,200,0);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERPrimaryGeneratorAction::~SUPERPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGeneralParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if(gUseGPS){
    gPrimaryParticlePosition = fGeneralParticleSource -> GetParticlePosition();
    gPrimaryParticleEnergy = fGeneralParticleSource -> GetParticleEnergy();
    gPrimaryParticleMomentumDirection = fGeneralParticleSource -> GetParticleMomentumDirection();
    gPrimaryParticlePDG = fGeneralParticleSource -> GetParticleDefinition() -> GetPDGEncoding();
    gPrimaryParticleMass = fGeneralParticleSource -> GetParticleDefinition() -> GetPDGMass();
    fGeneralParticleSource -> GeneratePrimaryVertex(event);
  }
  else{
    if (gGenGamma == 0){
      fParticleGun->SetParticleDefinition(fParticle);
      fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
      fParticleGun->SetParticleTime(0.0*ns);
      // random generation

      double pi0MomentumCM = fpi0pdist -> GetRandom();
      double CosThetaCM = fpi0Angulardist -> GetRandom();
      double phi = gRandom -> Uniform(0,360)*deg; // uniform phi
      
      double dx,dy,dz;
      dz = CosThetaCM;
      double sin_theta= sqrt(1.0-dz*dz);
      dx = sin_theta * cos(phi);
      dy = sin_theta * sin(phi);
      TLorentzVector pi0lv(dx*pi0MomentumCM,dy*pi0MomentumCM,dz*pi0MomentumCM,sqrt(pi0MomentumCM*pi0MomentumCM+fpi0mass*fpi0mass));
      double N14KineticEnergy = 14*40*MeV;
      
      double N14Momentum = sqrt(N14KineticEnergy*N14KineticEnergy + 2.0*N14KineticEnergy*N14mass);
      TLorentzVector NAlSystemCM(0,0,N14Momentum,N14KineticEnergy + N14mass + Al27mass);
      auto boostVec = NAlSystemCM.BoostVector();
      pi0lv.Boost(boostVec);
      G4double px = pi0lv.Px()/pi0lv.P();
      G4double py = pi0lv.Py()/pi0lv.P();
      G4double pz = pi0lv.Pz()/pi0lv.P();
      G4double Ek = pi0lv.E() - fpi0mass;
      
      //G4ParticleMomentum momentum(pi0lv.Px(),pi0lv.Py(),pi0lv.Pz());
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
      fParticleGun->SetParticleEnergy(Ek);

      gPrimaryParticlePosition = fParticleGun -> GetParticlePosition();
      gPrimaryParticleEnergy = fParticleGun -> GetParticleEnergy();
      gPrimaryParticleMomentumDirection = fParticleGun -> GetParticleMomentumDirection();
      gPrimaryParticlePDG = fParticleGun -> GetParticleDefinition() -> GetPDGEncoding();
      gPrimaryParticleMass = fParticleGun -> GetParticleDefinition() -> GetPDGMass();
      fParticleGun->GeneratePrimaryVertex(event);
      
    }
    else if (gGenGamma == 1){
      G4double dx,dy,dz;
      G4double Ek = gRandom -> Uniform(0,200) * MeV;
      gRandom -> Sphere(dx,dy,dz,1);
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx,dy,dz));
      fParticleGun->SetParticleEnergy(Ek);

      gPrimaryParticlePosition = fParticleGun -> GetParticlePosition();
      gPrimaryParticleEnergy = fParticleGun -> GetParticleEnergy();
      gPrimaryParticleMomentumDirection = fParticleGun -> GetParticleMomentumDirection();
      gPrimaryParticlePDG = fParticleGun -> GetParticleDefinition() -> GetPDGEncoding();
      gPrimaryParticleMass = fParticleGun -> GetParticleDefinition() -> GetPDGMass();
      fParticleGun->GeneratePrimaryVertex(event);
    }
    else if (gGenGamma == 2){
      fParticleGun -> SetNumberOfParticles(2);
      double gamma1E = fgammaEdist -> GetRandom();
      double gamma2E = fgammaEdist -> GetRandom();
      //double gamma_cos1 = fpi0Angulardist -> GetRandom();
      //double gamma_cos2 = fpi0Angulardist -> GetRandom();
      
      //double phi1 = gRandom -> Uniform(0,360)*deg; // uniform phi
      //double phi2 = gRandom -> Uniform(0,360)*deg; // uniform phi

      G4double dx1,dy1,dz1;
      G4double dx2,dy2,dz2;

      gRandom -> Sphere(dx1,dy1,dz1,1);
      gRandom -> Sphere(dx2,dy2,dz2,1);
      
      //dz1 = gamma_cos1;
      //dz2 = gamma_cos2;
      
      //dx1 = sqrt(1.0-dz1*dz1) * cos(phi1);
      //dx2 = sqrt(1.0-dz2*dz2) * cos(phi2);
      
      TLorentzVector gamma1(dx1*gamma1E,dy1*gamma1E,dz1*gamma1E,gamma1E);
      TLorentzVector gamma2(dx2*gamma2E,dy2*gamma2E,dz2*gamma2E,gamma2E);
      
      double N14KineticEnergy = 14*40*MeV;
      double N14Momentum = sqrt(N14KineticEnergy*N14KineticEnergy + 2.0*N14KineticEnergy*N14mass);
      TLorentzVector NAlSystemCM(0,0,N14Momentum,N14KineticEnergy + N14mass + Al27mass);
      auto boostVec = NAlSystemCM.BoostVector();

      gamma1.Boost(boostVec);
      gamma2.Boost(boostVec);

      G4double px1 = gamma1.Px();
      G4double py1 = gamma1.Py();
      G4double pz1 = gamma1.Pz();
      G4double E1 = gamma1.E();

      G4double px2 = gamma2.Px();
      G4double py2 = gamma2.Py();
      G4double pz2 = gamma2.Pz();
      G4double E2 = gamma2.E();

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx1,dy1,dz1));
      fParticleGun->SetParticleEnergy(E1);

      gPrimaryParticlePosition = fParticleGun -> GetParticlePosition();
      gPrimaryParticleEnergy = fParticleGun -> GetParticleEnergy();
      gPrimaryParticleMomentumDirection = fParticleGun -> GetParticleMomentumDirection();
      gPrimaryParticlePDG = fParticleGun -> GetParticleDefinition() -> GetPDGEncoding();
      gPrimaryParticleMass = fParticleGun -> GetParticleDefinition() -> GetPDGMass();
      fParticleGun->GeneratePrimaryVertex(event);

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx2,dy2,dz2));
      fParticleGun->SetParticleEnergy(E2);
      
      gPrimaryParticlePosition2 = fParticleGun -> GetParticlePosition();
      gPrimaryParticleEnergy2 = fParticleGun -> GetParticleEnergy();
      gPrimaryParticleMomentumDirection2 = fParticleGun -> GetParticleMomentumDirection();
      gPrimaryParticlePDG2 = fParticleGun -> GetParticleDefinition() -> GetPDGEncoding();
      gPrimaryParticleMass2 = fParticleGun -> GetParticleDefinition() -> GetPDGMass();
      fParticleGun->GeneratePrimaryVertex(event);

      
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
