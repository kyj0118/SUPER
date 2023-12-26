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
/// \file SUPEREmCalorimeterSD.hh
/// \brief Definition of the SUPEREmCalorimeterSD class

#ifndef SUPEREmCalorimeterSD_h
#define SUPEREmCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SUPEREmCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// EM calorimeter sensitive detector
using namespace std;
class SUPEREmCalorimeterSD : public G4VSensitiveDetector
{   
public:
  SUPEREmCalorimeterSD(G4String name, G4int ctype);
  virtual ~SUPEREmCalorimeterSD();
  
  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*);
private:
  G4String fNameSD;

  vector<G4double> fEdep;
  vector<G4double> fEweightedx;
  vector<G4double> fEweightedy;
  vector<G4double> fEweightedz;
  vector<G4double> fEweightedt;
  
  SUPEREmCalorimeterHitsCollection* fHitsCollection;
  G4int TotalPhi;
  G4int fHCID;
  G4int fPhiID;
  G4int fCrystalType;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
