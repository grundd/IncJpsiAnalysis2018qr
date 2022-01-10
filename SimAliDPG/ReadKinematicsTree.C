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
Int_t opt = 0;
// = 0 => kCohPsi2sToMuPi, 12-30-2021 (1000 events, std RA = 6.602 fm)
// = 1 => kCohPsi2sToMuPi, 01-03-2022 (1e4 events, std RA = 6.602 fm)
Bool_t write_to_txt_file = kTRUE;
// ***************

void ReadKinematicsTree() 
{
    TString dir = "";
    if(opt == 0) dir = "kCohPsi2sToMuPi_12-30-2021";
    if(opt == 1) dir = "kCohPsi2sToMuPi";
    if(opt == 2) dir = "kCohJpsiToMu";
    TString gAliceFile = dir.Data();
    gAliceFile += "/";
    gAliceFile += "galice.root";

    // Open the working directory
    AliRunLoader* runLoader = AliRunLoader::Open(gAliceFile);
    runLoader->LoadKinematics();

    // Loop over events
    Int_t nEvents = runLoader->GetNumberOfEvents();
    Printf("TOTAL NUMBER OF EVENTS: %i", nEvents);

    TLorentzVector vGen;
    TLorentzVector vJpsiMuons;// J/psi reconstructed from mu and antimu 
    TLorentzVector vJpsiMuPho;// J/psi reconstructed from muons and photons (with indices larger than those of muons)
    TLorentzVector vPsi2s;    // Psi(2s) reconstructed from all decay products
    Double_t pT_Jpsi = 0;
    // Create output text file
    ofstream outfile;
    outfile.open(dir + "/out.txt");
    // Create output root file
    TFile fFileOut(dir + "/tPtGen.root","RECREATE");
    // Create output tree
    // To this tree, we will store the transverse momentum of the 
    // J/psi reconstructed from muon pairs (i.e., p_T of vJpsiMuons)
    Double_t fPtGen = 0;
    TTree *tPtGen = new TTree("tPtGen", "tPtGen");
    tPtGen->Branch("fPtGen", &fPtGen, "fPtGen/D");

    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
        vJpsiMuons.SetXYZM(0.,0.,0.,0.);
        vJpsiMuPho.SetXYZM(0.,0.,0.,0.);
        vPsi2s.SetXYZM(0.,0.,0.,0.);

        if(write_to_txt_file) outfile << "++++ EVENT no. " << iEvent+1 << " ++++" << endl;

        if((iEvent+1 % 1000) == 0) Printf("%i events analysed", iEvent+1);

        runLoader->GetEvent(iEvent);	    
        AliStack* stack = runLoader->Stack();
        Int_t nParticles = stack->GetNtrack();

        Int_t iMuons = 0;

        // loop over particles in the event
        for(Int_t iPart = 0; iPart < nParticles; iPart++){
            TParticle *part = stack->Particle(iPart);
            vGen.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
            // muons
            if(TMath::Abs(part->GetPdgCode()) == 13){
                vJpsiMuons += vGen;
                vJpsiMuPho += vGen;
                vPsi2s += vGen;
                iMuons++;
            }
            // photons
            if(part->GetPdgCode() == 22){
                vPsi2s += vGen;
                if(iMuons == 2) vJpsiMuPho += vGen;
            } 
            // pions
            if(TMath::Abs(part->GetPdgCode()) == 211) vPsi2s += vGen;
            // print to text file
            if(write_to_txt_file){
                outfile << " i = " << iPart 
                        << "\tmother = " << part->GetMother(0) 
                        << "\tstatus = " << part->GetStatusCode() 
                        << "\tpdg = " << part->GetPdgCode()
                        << "\tname = " << part->GetName() 
                        << "\tdaughter = " << part->GetDaughter(0)
                        << std::fixed << std::setprecision(5)
                        << "\teta = " << part->Eta() 
                        << "\tp_x = " << part->Px()
                        << "\tp_y = " << part->Py() 
                        << "\tp_z = " << part->Pz()
                        << std::fixed << std::setprecision(3)
                        << "\tm = " << part->GetCalcMass() << endl;
            }
            // save J/psi transverse momentum
            if(part->GetPdgCode() == 443) pT_Jpsi = part->Pt();
        }
        // save p_T of reconstructed J/psi to the tree
        fPtGen = vJpsiMuons.Pt();
        tPtGen->Fill();
        // print info to text files
        if(write_to_txt_file){
            ///*
            outfile << "J/psi (muons):" << std::fixed << std::setprecision(5)
                    << "\teta = " << vJpsiMuons.Eta() 
                    << "\tp_x = " << vJpsiMuons.Px()
                    << "\tp_y = " << vJpsiMuons.Py() 
                    << "\tp_z = " << vJpsiMuons.Pz()
                    << std::fixed << std::setprecision(3)
                    << "\tm = " << vJpsiMuons.M() << endl;
            outfile << "J/psi (mu+ph):" << std::fixed << std::setprecision(5)
                    << "\teta = " << vJpsiMuPho.Eta() 
                    << "\tp_x = " << vJpsiMuPho.Px()
                    << "\tp_y = " << vJpsiMuPho.Py() 
                    << "\tp_z = " << vJpsiMuPho.Pz()
                    << std::fixed << std::setprecision(3)
                    << "\tm = " << vJpsiMuPho.M() << endl;
            outfile << "Psi(2s) (all):" << std::fixed << std::setprecision(5)
                    << "\teta = " << vPsi2s.Eta() 
                    << "\tp_x = " << vPsi2s.Px()
                    << "\tp_y = " << vPsi2s.Py() 
                    << "\tp_z = " << vPsi2s.Pz()
                    << std::fixed << std::setprecision(3)
                    << "\tm = " << vPsi2s.M() << endl;
            //*/
            outfile << std::fixed << std::setprecision(5)
                    << "J/psi (orig): p_T = " << pT_Jpsi << endl
                    << "J/psi (muon): p_T = " << vJpsiMuons.Pt() << endl
                    << "diff = " << TMath::Abs(pT_Jpsi - vJpsiMuons.Pt()) << endl;
        }
    } //event loop

    runLoader->UnloadAll();

    fFileOut.Write("",TObject::kWriteDelete);
    
    return;
}

