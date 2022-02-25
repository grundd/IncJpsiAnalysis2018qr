// ReadKinematicsTree.C
// dgrund, 01-02-2022 (Jan 2, 2022)

// cpp headers
#include "fstream"
#include "stdio.h"
#include "iomanip" // std::setprecision()
// aliroot headers
#include "AliRunLoader.h"
#include "AliStack.h"
// root headers
#include "TString.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TMath.h"

// ***************
// Options to set:
Int_t opt = 7;
Bool_t write_to_txt = kTRUE;
// ***************

void ReadKinematicsTree() 
{
    Bool_t FD = kFALSE;

    TString dir = "";
    if(opt == 0){ dir = "kCohPsi2sToMuPi_12-30-2021"; FD = kTRUE; }
    if(opt == 1){ dir = "kCohPsi2sToMuPi"; FD = kTRUE; }
    if(opt == 2){ dir = "kCohJpsiToMu"; }
    if(opt == 3){ dir = "officialALICE_kCohPsi2sToMuPi"; FD = kTRUE; }
    if(opt == 4){ dir = "officialALICE_kCohJpsiToMu"; }
    if(opt == 5){ dir = "officialALICE_kTwoGammaToMuMedium"; }
    if(opt == 6){ dir = "officialALICE_kCohJpsiToElRad"; }
    if(opt == 7){ dir = "officialALICE_kCohJpsiToEl"; }
    // Get the directory
    TString gAliceFile = dir.Data();
    gAliceFile += "/";
    gAliceFile += "galice.root";

    // Open the working directory
    AliRunLoader* runLoader = AliRunLoader::Open(gAliceFile);
    runLoader->LoadKinematics();

    // Loop over events
    Int_t nEvents = runLoader->GetNumberOfEvents();
    Printf("TOTAL NUMBER OF EVENTS: %i", nEvents);

    TLorentzVector vGen;        // current particle from Kinematics.root
    TLorentzVector vJpsi_orig;  // original J/psi in Kinematics.root
    TLorentzVector vJpsi_muon;  // J/psi reconstructed from mu and antimu 
    TLorentzVector vJpsi_muga;  // J/psi reconstructed from muons and photons (with indices larger than those of muons)
    TLorentzVector vPsi2s;      // Psi(2s) reconstructed from all decay products
    // Create output text file
    ofstream out_txt;
    out_txt.open(dir + "/out.txt");
    // Analysis of photons
    ofstream out_txt_ph;
    if(FD) out_txt_ph.open(dir + "/out_photons.txt");
    Int_t nPh0(0), nPh1(0), nPh2(0);
    // Create output root file
    TFile out_file(dir + "/tPtGen.root","RECREATE");
    // Create output tree
    // To this tree, we will store the transverse momentum of vJpsi_muon 
    Double_t fPtGen = 0;
    TTree *tPtGen = new TTree("tPtGen", "tPtGen");
    tPtGen->Branch("fPtGen", &fPtGen, "fPtGen/D");

    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
        vJpsi_muon.SetXYZM(0.,0.,0.,0.);
        vJpsi_muga.SetXYZM(0.,0.,0.,0.);
        vPsi2s.SetXYZM(0.,0.,0.,0.);

        if(write_to_txt) out_txt << "++++ EVENT no. " << iEvent+1 << " ++++" << endl;

        if((iEvent+1 % 1000) == 0) Printf("%i events analysed", iEvent+1);

        runLoader->GetEvent(iEvent);	    
        AliStack* stack = runLoader->Stack();
        Int_t nParticles = stack->GetNtrack();
        Int_t nGammaPsi2s = 0; // number of photons from a decay of Psi(2s)
        Bool_t alreadyMuon = kFALSE; // if a muon already seen

        // loop over particles in the event
        for(Int_t iPart = 0; iPart < nParticles; iPart++){
            // Get current particle
            TParticle *part = stack->Particle(iPart);
            vGen.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
            // muons
            if(TMath::Abs(part->GetPdgCode()) == 13){
                alreadyMuon = kTRUE;
                // get mother
                // if not mother (kCohJpsiToMu, kIncohJpsiToMu, kTwoGammaToMuMedium)
                if(part->GetFirstMother() == -1){
                    vJpsi_muon += vGen;
                } else {
                    TParticle *mother = stack->Particle(part->GetFirstMother());
                    // if mother is Psi(2s) (kCohPsi2sToMuPi, kIncohPsi2sToMuPi)
                    if(TMath::Abs(mother->GetPdgCode()) == 100443){
                        vJpsi_muon += vGen;
                        vJpsi_muga += vGen;
                        vPsi2s += vGen;
                    }
                }
            }
            if(FD){
                // photons
                if(part->GetPdgCode() == 22){
                    // get mother
                    TParticle *mother = stack->Particle(part->GetFirstMother());
                    // if mother is Psi(2s)
                    if(TMath::Abs(mother->GetPdgCode()) == 100443){
                        if(alreadyMuon) vJpsi_muga += vGen;
                        else{
                            Printf("Ev %i has a photon from the decay of Psi(2S).", iEvent+1);
                            nGammaPsi2s++; 
                        }
                        vPsi2s += vGen;  
                    }
                } 
                // pions
                if(TMath::Abs(part->GetPdgCode()) == 211){
                    // get mother
                    TParticle *mother = stack->Particle(part->GetFirstMother());
                    // if mother is Psi(2s)
                    if(TMath::Abs(mother->GetPdgCode()) == 100443){
                        vPsi2s += vGen;                    
                    }                
                }
                // J/psi
                if(TMath::Abs(part->GetPdgCode()) == 443) vJpsi_orig.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
            }
            // print to text file
            if(write_to_txt){
                out_txt << " i = " << iPart 
                        << "\tmother = " << part->GetMother(0) 
                        //<< "\tstatus = " << part->GetStatusCode() 
                        << "\tpdg = " << part->GetPdgCode()
                        << "\tname = " << part->GetName() 
                        << "\tdaughter = " << part->GetDaughter(0)
                        << std::fixed << std::setprecision(5)
                        << "\teta = " << part->Eta() 
                        << "\tp_x = " << part->Px()
                        << "\tp_y = " << part->Py() 
                        << "\tp_z = " << part->Pz()
                        << "\tm = " << part->GetCalcMass() << endl;
            }
        } // particle loop
        // count photons from the decay of Psi(2s)
        if(FD){
            if(nGammaPsi2s > 0){
                out_txt_ph << "Ev " << iEvent+1 << " has " << nGammaPsi2s << " gamma." << endl;
            }
            if(nGammaPsi2s == 0) nPh0++;
            if(nGammaPsi2s == 1) nPh1++;
            if(nGammaPsi2s == 2) nPh2++;
        }
        // save p_T of vJpsi_muon
        fPtGen = vJpsi_muon.Pt();
        tPtGen->Fill();
        // print info to text files
        if(FD && write_to_txt){
            out_txt << "J/psi (orig): " << std::fixed << std::setprecision(5)
                    << "\teta = " << vJpsi_orig.Eta() 
                    << "\tp_x = " << vJpsi_orig.Px()
                    << "\tp_y = " << vJpsi_orig.Py() 
                    << "\tp_z = " << vJpsi_orig.Pz()
                    << "\tp_T = " << vJpsi_orig.Pt()
                    << "\tm = " << vJpsi_orig.M() << endl;
            out_txt << "J/psi (muon): " << std::fixed << std::setprecision(5)
                    << "\teta = " << vJpsi_muon.Eta() 
                    << "\tp_x = " << vJpsi_muon.Px()
                    << "\tp_y = " << vJpsi_muon.Py() 
                    << "\tp_z = " << vJpsi_muon.Pz()
                    << "\tp_T = " << vJpsi_muon.Pt()
                    << "\tm = " << vJpsi_muon.M() << endl;
            out_txt << "J/psi (mu+g): " << std::fixed << std::setprecision(5)
                    << "\teta = " << vJpsi_muga.Eta() 
                    << "\tp_x = " << vJpsi_muga.Px()
                    << "\tp_y = " << vJpsi_muga.Py() 
                    << "\tp_z = " << vJpsi_muga.Pz()
                    << "\tp_T = " << vJpsi_muga.Pt()
                    << "\tm = " << vJpsi_muga.M() << endl;
            out_txt << "Psi(2s) (all):" << std::fixed << std::setprecision(5)
                    << "\teta = " << vPsi2s.Eta() 
                    << "\tp_x = " << vPsi2s.Px()
                    << "\tp_y = " << vPsi2s.Py() 
                    << "\tp_z = " << vPsi2s.Pz()
                    << "\tp_T = " << vPsi2s.Pt()
                    << "\tm = " << vPsi2s.M() << endl;
        }
    } // event loop

    if(FD) out_txt_ph << endl << "nPh0 = " << nPh0 << endl << "nPh1 = " << nPh1 << endl << "nPh2 = " << nPh2 << endl;

    runLoader->UnloadAll();

    out_file.Write("",TObject::kWriteDelete);
    out_txt.close();
    if(FD) out_txt_ph.close();

    return;
}

