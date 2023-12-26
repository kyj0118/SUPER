#include "SUPERStackingAction.hh"
#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4int> gTrackIDvector;
std::vector<G4int> gPDGEncodingvector;
std::vector<G4ThreeVector> gPositionvector;
std::vector<G4ThreeVector> gMomentumvector;
std::vector<G4double> gEtotvector;
std::vector<G4double> gMassvector;
//std::vector<G4ThreeVector> gPositionvector;
extern int gGenGamma;
SUPERStackingAction::SUPERStackingAction()
  : G4UserStackingAction(),
    fTrackCount(0), fGammaCount(0), felectronCount(0), fpositronCount(0), fOtherCount(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERStackingAction::~SUPERStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
SUPERStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if (gGenGamma == 0){
    G4int trackid = aTrack -> GetTrackID();
    G4int pdgnum = aTrack -> GetParticleDefinition() -> GetPDGEncoding();
    G4int parentid = 	aTrack -> GetParentID();
    if (parentid == 1){
      //auto pos = aTrack -> GetVertexPosition();
      auto pos = aTrack -> GetPosition();
      //auto mom = aTrack -> GetMomentumDirection();
      auto mom = aTrack -> GetMomentum();
      auto Ek = aTrack -> GetKineticEnergy();
      //auto Ek = aTrack -> GetVertexKineticEnergy();
      //std::cout << "pos : (" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << std::endl;
      //std::cout << "(" << mom.x() << ", " << mom.y() << ", " << mom.z() << ")" << std::endl;

      auto mass = aTrack -> GetParticleDefinition() -> GetPDGMass();
      G4double Etot = Ek+mass;
      auto p = sqrt(Etot*Etot - mass*mass);
      //mom = mom * p;

      //std::cout << "pdg : " << pdgnum <<" : " <<  Ek << ", " << p << std::endl;
      //std::cout << "(" << mom.x() << ", " << mom.y() << ", " << mom.z() << ")" << std::endl;
      gTrackIDvector.push_back(trackid);
      gPDGEncodingvector.push_back(pdgnum);
      gPositionvector.push_back(pos);
      gMomentumvector.push_back(mom);
      gMassvector.push_back(mass);
      gEtotvector.push_back(Etot);

    }
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERStackingAction::NewStage()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERStackingAction::PrepareNewEvent()
{
  gTrackIDvector.clear();
  gPDGEncodingvector.clear();
  gPositionvector.clear();
  gMomentumvector.clear();
  gMassvector.clear();
  gEtotvector.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
