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
/// \file SUPERDetectorConstruction.cc
/// \brief Implementation of the SUPERDetectorConstruction class

#include "SUPERDetectorConstruction.hh"
#include "SUPEREmCalorimeterSD.hh"

#include "G4TransportationManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"

#include "G4Box.hh"
//#include "G4Trap.hh"
#include "G4GenericTrap.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4tgbRotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

// Root classes
//#include "TString.h"
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int gNLayers;
G4int gNStrips;
G4int gTotalNScintillators;

using namespace std;
SUPERDetectorConstruction::SUPERDetectorConstruction()
  : G4VUserDetectorConstruction()
{
  gTotalNScintillators = 1000;
  fScintLength = 25*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERDetectorConstruction::~SUPERDetectorConstruction()
{
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SUPERDetectorConstruction::Construct(){
  G4NistManager* nist = G4NistManager::Instance();

  // -----------------------------------------------------
  // World
  
  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");
  G4double world_size = 5.*m;

  G4Box* solidWorld =    
    new G4Box("World",                       // its name
              0.5*world_size,                // half x
              0.5*world_size,                // half y
              0.5*world_size);               // half z
  
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      false);                 //overlaps checking


  // -----------------------------------------------------

  // Vis attributer
  auto visScint = new G4VisAttributes;
  //visScint -> SetColor(0.0,0.0,0.0,0.5);
  visScint -> SetColor(1.0,1.0,1.0);
  visScint -> SetLineWidth(2);
  
  // Detector
  G4Material* Material_CsI = nist -> FindOrBuildMaterial("G4_CESIUM_IODIDE");
  int icopynum = 0;
  for (int itype = 0; itype < 10; itype++){
    //for (int itype = 0; itype < 1; itype++){
    int iitype = itype+1;

    vector<G4TwoVector> vertices =   getGenericTrapVertices(iitype);
    vector<G4TwoVector> vertices_r = getGenericTrapVertices(-iitype);

    G4double z = fScintLength;
    auto crystal   = new G4GenericTrap("crystal",z/2.0,vertices);
    auto crystal_r = new G4GenericTrap("crystal_r",z/2.0,vertices_r);

    logicCsI1[itype] = new G4LogicalVolume(crystal,Material_CsI,"logic_crystal");
    logicCsI2[itype] = new G4LogicalVolume(crystal_r,Material_CsI,"logic_crystal_r");

    logicCsI1[itype] -> SetVisAttributes(visScint);
    logicCsI2[itype] -> SetVisAttributes(visScint);
    
    //for (int iphi = 0; iphi < 1; iphi++){
    //for (int iphi = 0; iphi < 3; iphi++){
    for (int iphi = 0; iphi < 48; iphi++){
    //for (int iphi = 0; iphi < 12; iphi++){
      if (itype >= 0 && itype <= 2 && ((int) (iphi/2)) % 2 == 0) continue;
      G4ThreeVector Zaxis1(0,0,1);
      G4ThreeVector Zaxis2(0,0,1);
      G4RotationMatrix *rotM1 = new G4RotationMatrix();
      G4RotationMatrix *rotM2 = new G4RotationMatrix();

      rotM1 -> rotateY(90*deg);
      rotM2 -> rotateY(90*deg);

      G4double InnerR = 20*cm;
      G4double InnerHolderR = 13.25*cm;
      G4double dTheta = 7.5*deg;
      G4double dPhi = 7.5*deg;
      G4double theta_i = (((G4double) itype) + 0.5) * dTheta;
      G4double phi_i = dPhi * ((G4double) iphi);
      G4double L_i = InnerR / ( (1.0 - tan(dTheta/2.0)*tan(theta_i))  * cos(theta_i));

      if (itype >= 6){
        double theta_i2 = (((G4double) itype) + 1.0) * dTheta;
        double holderSlope = (InnerR - InnerHolderR) / (InnerHolderR/tan(dTheta*3.0) - InnerR);
        double dx = (InnerR * tan(theta_i2) - InnerR) / (1.0 + tan(theta_i2)*holderSlope);
        InnerR = (InnerR - dx * holderSlope);
        //std::cout << InnerR << std::endl;
        L_i = InnerR / ( (1.0 - tan(dTheta/2.0)*tan(theta_i))  * cos(theta_i));
      }


      G4double xpos = L_i + fScintLength/2.0;
      G4double ypos = 0;
      G4double zpos = 0;

      G4double xrot = cos(theta_i) * xpos - sin(theta_i) * zpos;
      G4double yrot = 0;
      G4double zrot = sin(theta_i) * xpos + cos(theta_i) * zpos;

      G4double xrrot = cos(phi_i) * xrot - sin(phi_i) * yrot;
      G4double yrrot = sin(phi_i) * xrot + cos(phi_i) * yrot;
      G4double zrrot = zrot;

      G4double x2rot = cos(-theta_i) * xpos - sin(-theta_i) * zpos;
      G4double y2rot = 0;
      G4double z2rot = sin(-theta_i) * xpos + cos(-theta_i) * zpos;

      G4double x2rrot = cos(phi_i) * x2rot - sin(phi_i) * y2rot;
      G4double y2rrot = sin(phi_i) * x2rot + cos(phi_i) * y2rot;
      G4double z2rrot = z2rot;
      //G4double x2rrot = cos(-phi_i) * x2rot - sin(-phi_i) * y2rot;
      //G4double y2rrot = sin(-phi_i) * x2rot + cos(-phi_i) * y2rot;

      G4ThreeVector pos1(xrrot,yrrot,zrrot);
      G4ThreeVector v1 = Zaxis1.cross(-pos1);
      G4RotationMatrix *rotMtest = new G4RotationMatrix();
      //rotMtest->rotate(pos1.angle(Zaxis1)+180*deg,v1);
      rotMtest->rotate(pos1.angle(Zaxis1),v1);

      G4ThreeVector pos2(x2rrot,y2rrot,z2rrot);
      G4ThreeVector v2 = Zaxis2.cross(-pos2);
      G4RotationMatrix *rotMtest2 = new G4RotationMatrix();
      rotMtest2->rotate(pos2.angle(Zaxis2)+180*deg,v2);


      rotM1 -> rotateY(theta_i);
      Zaxis1 *= (*rotM1);
      rotM1 -> rotate(-phi_i,Zaxis1);

      rotM2 -> rotateY(-theta_i);
      Zaxis2 *= (*rotM2);
      rotM2 -> rotate(-phi_i,Zaxis2);

      //cout << Zaxis1.x() <<", " << Zaxis1.y() << ", " << Zaxis1.z() << endl;
      //cout << Zaxis2.x() <<", " << Zaxis2.y() << ", " << Zaxis2.z() << endl;
      //Zaxis2 *= (rotM2->inverse());
      //rotM2 -> rotate(phi_i,Zaxis2);

      /*
      new G4PVPlacement(rotMtest,
                        pos1,
                        logicCsI,
                        "physicalCsI",
                        logicWorld,
                        false,
                        0,
                        false);
      */
      /*
      new G4PVPlacement(rotMtest2,
                        pos2,
                        logicCsI,
                        "physicalCsI",
                        logicWorld,
                        false,
                        0,
                        false);
      */


    
      new G4PVPlacement(rotM1,
                        G4ThreeVector(xrrot,yrrot,zrrot),
                        logicCsI1[itype],
                        "physicalCsI1",
                        logicWorld,
                        false,
                        iphi,
                        false);


      new G4PVPlacement(rotM2,
                        G4ThreeVector(x2rrot,y2rrot,z2rrot),
                        logicCsI2[itype],
                        "physicalCsI2",
                        logicWorld,
                        false,
                        iphi,
                        false);

      icopynum++;
    }
  }
  return physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERDetectorConstruction::ConstructSDandField()
{

  for (int itype = 0; itype < 10; itype++){
    int iitype = itype+1;
    G4String name1 = "ScintSD1" + std::to_string(iitype);
    G4String name2 = "ScintSD2" + std::to_string(iitype);
    SUPEREmCalorimeterSD* ScintillatorSD1 = new SUPEREmCalorimeterSD(name1,iitype);
    SUPEREmCalorimeterSD* ScintillatorSD2 = new SUPEREmCalorimeterSD(name2,-iitype);
    G4SDManager::GetSDMpointer() -> AddNewDetector(ScintillatorSD1);
    G4SDManager::GetSDMpointer() -> AddNewDetector(ScintillatorSD2);

    SetSensitiveDetector(logicCsI1[itype], ScintillatorSD1);
    SetSensitiveDetector(logicCsI2[itype], ScintillatorSD2);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<G4TwoVector> SUPERDetectorConstruction::getGenericTrapVertices(int iitype)
{

  bool is_reversed = false;
  if (iitype < 0) is_reversed = true;
  iitype = abs(iitype);
  int itype = iitype - 1;
  /*
  double A[10] = {61.14, 59.76, 57.42, 53.30, 50.92, 48.10, 42.26, 35.77, 28.52, 42.15};
  double B[10] = {61.98, 62.31, 61.71, 56.87, 55.61, 54.01, 49.16, 43.69, 37.51, 62.21};
  double H[10] = {57.76, 58.46, 59.69, 61.58, 64.35, 68.37, 70.05, 72.64, 76.42, 82.42};
  double a[10] = {5.63,  5.64,  5.66,  24.86, 24.87, 24.88, 22.27, 19.35, 15.96, 25.08};
  double b[10] = {6.00,  6.77,  7.60,  26.53, 27.17, 27.95, 25.94, 23.70, 21.09, 37.12};
  double h[10] = {24.99, 25.68, 26.92, 28.81, 31.58, 35.60, 37.28, 39.87, 43.65, 49.46};

  G4double xb1 = H[itype]*mm;
  G4double xb2 = H[itype]*mm;
  G4double yb1 = A[itype]*mm;
  G4double yb2 = B[itype]*mm;

  G4double xt1 = h[itype]*mm;
  G4double xt2 = h[itype]*mm;
  G4double yt1 = a[itype]*mm;
  G4double yt2 = b[itype]*mm;
  */

  G4double InnerR = 20*cm;
  G4double InnerHolderR = 13.25*cm;
  G4double dTheta = 7.5*deg;
  G4double dPhi = 7.5*deg;
  G4double theta_i = (((G4double) itype) + 0.5) * dTheta;
  G4double L_i = InnerR / ( (1.0 - tan(dTheta/2.0)*tan(theta_i))  * cos(theta_i));

  if (itype >= 6){
    double theta_i2 = (((G4double) itype) + 1.0) * dTheta;
    double holderSlope = (InnerR - InnerHolderR) / (InnerHolderR/tan(dTheta*3.0) - InnerR);
    //std::cout << holderSlope << std::endl;
    double dx = (InnerR * tan(theta_i2) - InnerR) / (1.0 + tan(theta_i2)*holderSlope);
    InnerR = (InnerR - dx * holderSlope);
    L_i = InnerR / ( (1.0 - tan(dTheta/2.0)*tan(theta_i))  * cos(theta_i));
  }

//G4double h = 2.0 * L_i * tan(dTheta/2.0);
  //G4double H = 2.0 * (L_i + 25.0*cm) * tan(dTheta/2.0);

  G4double h = 2.0 * L_i * tan(dTheta/2.0);
  G4double H = 2.0 * (L_i + fScintLength) * tan(dTheta/2.0);

  G4double a = 2.0 * (InnerR + h*sin(theta_i)) * tan(dPhi/2.0);
  G4double b = 2.0 * InnerR * tan(dPhi/2.0);

  //G4double A = 2.0 * (L_i + 25.0*cm) * tan(dTheta/2.0);
  //G4double B = 2.0 * (L_i + 25.0*cm) * tan(dTheta/2.0);

  G4double A = 2.0 * ((L_i + fScintLength) * cos(theta_i) + H*sin(theta_i)/2.0 ) * tan(dPhi/2.0);
  G4double B = 2.0 * ((L_i + fScintLength) * cos(theta_i) - H*sin(theta_i)/2.0 ) * tan(dPhi/2.0);


  //G4double a = 2.0 * L_i * tan(dTheta/2.0);
  //G4double a = 1.0 * InnerR * tan(dTheta/2.0);
  //G4double b = 2.0 * L_i * tan(dTheta/2.0);
  //G4double b = 2.0 * (InnerR + h*sin(theta_i)) * tan(dTheta/2.0);



  //G4double a =  L_i * tan(dTheta/2.0);
  //G4double A =  (L_i + 25.0*cm) * tan(dTheta/2.0);

  //G4double b =  L_i * tan(dTheta/2.0);
  //G4double B =  (L_i + 25.0*cm) * tan(dTheta/2.0);

  /*
  G4double xb1 = A;
  G4double xb2 = B;
  G4double yb1 = H;
  G4double yb2 = H;

  G4double xt1 = a;
  G4double xt2 = b;
  G4double yt1 = h;
  G4double yt2 = h;
  */

  G4double xb1 = H;
  G4double xb2 = H;
  G4double yb1 = A;
  G4double yb2 = B;

  G4double xt1 = h;
  G4double xt2 = h;
  G4double yt1 = a;
  G4double yt2 = b;
  if (is_reversed){
    yb1 = B;
    yb2 = A;
    yt1 = b;
    yt2 = a;
  }

  vector<G4TwoVector> vertices;
  vertices.clear();
  /*
  vertices.push_back(G4TwoVector(xb1/2.0,-yb1/2.0));
  vertices.push_back(G4TwoVector(-xb1/2.0,-yb1/2.0));
  vertices.push_back(G4TwoVector(-xb2/2.0,yb2/2.0));
  vertices.push_back(G4TwoVector(xb2/2.0,yb2/2.0));

  vertices.push_back(G4TwoVector(xt1/2.0,-yt1/2.0));
  vertices.push_back(G4TwoVector(-xt1/2.0,-yt1/2.0));
  vertices.push_back(G4TwoVector(-xt2/2.0,yt2/2.0));
  vertices.push_back(G4TwoVector(xt2/2.0,yt2/2.0));
  */


  vertices.push_back(G4TwoVector(-xb1/2.0,-yb1/2.0));
  vertices.push_back(G4TwoVector(-xb2/2.0,yb1/2.0));
  vertices.push_back(G4TwoVector(xb2/2.0,yb2/2.0));
  vertices.push_back(G4TwoVector(xb1/2.0,-yb2/2.0));

  vertices.push_back(G4TwoVector(-xt1/2.0,-yt1/2.0));
  vertices.push_back(G4TwoVector(-xt2/2.0,yt1/2.0));
  vertices.push_back(G4TwoVector(xt2/2.0,yt2/2.0));
  vertices.push_back(G4TwoVector(xt1/2.0,-yt2/2.0));



  //vertices.push_back(G4TwoVector(-xt1/2.0,-yt1/2.0));
  //vertices.push_back(G4TwoVector(xt1/2.0,-yt2/2.0));
  //vertices.push_back(G4TwoVector(xt2/2.0,yt2/2.0));
  //vertices.push_back(G4TwoVector(-xt2/2.0,yt1/2.0));

  return vertices;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
