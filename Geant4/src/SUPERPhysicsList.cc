#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "SUPERPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4EmStandardPhysics.hh"
// For setting smallest step length
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4Decay.hh"
#include "G4Threading.hh"

#include "G4PhotoNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4MuonNuclearProcess.hh"

G4ThreadLocal G4int SUPERPhysicsList::fVerboseLevel = 1;
G4ThreadLocal G4int SUPERPhysicsList::fMaxNumPhotonStep = 20;
G4ThreadLocal G4Cerenkov* SUPERPhysicsList::fCerenkovProcess = 0;
G4ThreadLocal G4Scintillation* SUPERPhysicsList::fScintillationProcess = 0;
G4ThreadLocal G4OpAbsorption* SUPERPhysicsList::fAbsorptionProcess = 0;
G4ThreadLocal G4OpRayleigh* SUPERPhysicsList::fRayleighScatteringProcess = 0;
G4ThreadLocal G4OpMieHG* SUPERPhysicsList::fMieHGScatteringProcess = 0;
G4ThreadLocal G4OpBoundaryProcess* SUPERPhysicsList::fBoundaryProcess = 0;

G4ThreadLocal G4OpWLS* SUPERPhysicsList::fWLSProcess = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERPhysicsList::SUPERPhysicsList() 
 : G4VUserPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SUPERPhysicsList::~SUPERPhysicsList() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4BosonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  G4BaryonConstructor rConstructor;
  rConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructDecay();
  G4VPhysicsConstructor* emList = new G4EmStandardPhysics();
  emList->ConstructProcess(); 
  //ConstructNuclearProcess();  
  //ConstructEM();
  //ConstructOp();
}

void SUPERPhysicsList::ConstructNuclearProcess(){
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      pmanager->AddDiscreteProcess(new G4PhotoNuclearProcess());
    } else if (particleName == "e-") {
      pmanager->AddDiscreteProcess(new G4ElectronNuclearProcess());
    } else if (particleName == "e+") {
      pmanager->AddDiscreteProcess(new G4PositronNuclearProcess());
    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      pmanager->AddDiscreteProcess(new G4MuonNuclearProcess()); 
    }
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
    } else if (particleName == "e-") {
      //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);

    } else if (particleName == "e+") {
      //positron
      // Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1, 4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      //muon
      // Construct processes for muon
      pmanager->AddProcess(new G4MuMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation(),      -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung(),  -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction(),  -1, 4, 4);

    } else {
      if ((particle->GetPDGCharge() != 0.0) &&
          (particle->GetParticleName() != "chargedgeantino") &&
          !particle->IsShortLived()) {
	// all others charged particles except geantino
	pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
	pmanager->AddProcess(new G4hIonisation(),       -1,2,2);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void SUPERPhysicsList::ConstructOp()
{
  fCerenkovProcess = new G4Cerenkov("Cerenkov");
  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  fCerenkovProcess->SetTrackSecondariesFirst(true);
  
  fScintillationProcess = new G4Scintillation("Scintillation");
  fScintillationProcess->SetScintillationYieldFactor(1.);
  //  fScintillationProcess->SetScintillationExcitationRatio(0.0);
  // ScintillationExcitation ratio is used when there are a lot of excitation level, so you need to set where to go.
  fScintillationProcess->SetTrackSecondariesFirst(true);

  fAbsorptionProcess = new G4OpAbsorption();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  fBoundaryProcess = new G4OpBoundaryProcess();

  fWLSProcess = new G4OpWLS();
  //fWLSProcess->UseTimeProfile("delta");
  fWLSProcess->UseTimeProfile("exponential");

  //  fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  //  fScintillationProcess->SetVerboseLevel(fVerboseLevel);
  //  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  // fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  // fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  // fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
  // fWLSProcess->SetVerboseLevel(fVerboseLevel);

  // Use Birks Correction in the Scintillation process
  if(G4Threading::IsMasterThread())
    {
      G4EmSaturation* emSaturation =
	G4LossTableManager::Instance()->EmSaturation();
      fScintillationProcess->AddSaturation(emSaturation);
    }

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (fCerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(fCerenkovProcess);
      pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    if (fScintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(fScintillationProcess);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
    }
    /*
      if (fWLSProcess->IsApplicable(*particle))
      {
      pmanager->AddProcess(fWLSProcess);
      pmanager->SetProcessOrderingToLast(fWLSProcess, idxPostStep);
      }
    */
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(fAbsorptionProcess);
      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
      pmanager->AddDiscreteProcess(fBoundaryProcess);
      pmanager->AddDiscreteProcess(fWLSProcess);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::SetVerbose(G4int verbose)
{/*
   fVerboseLevel = verbose;

   fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
   fScintillationProcess->SetVerboseLevel(fVerboseLevel);
   fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
   fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
   fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
   fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
 */}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
  fMaxNumPhotonStep = MaxNumber;

  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SUPERPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  SetCutsWithDefault();
  /*
    defaultCutValue = 0.1*mm;
    
    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");
  */
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
