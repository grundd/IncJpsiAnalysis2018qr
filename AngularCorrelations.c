// AngularCorrelations.c
// David Grund, Sep 9, 2021
// To prepare a new tree with px, py and pz for both muons to investigate angular correlations 

#include "TFile.h"
#include "TLorentzVector.h"

#include "fTreeJPsiManager.h"

Double_t fPx1, fPy1, fPz1, fPx2, fPy2, fPz2;
Double_t fMuonMass = 0.105658; // GeV/c^2

Bool_t EventPassedAngCorr();

void AngularCorrelations(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    TFile fFileOut("Trees/AngularCorrelations/AngularCorrelations.root","RECREATE");

    TTree *tAngCorr = new TTree("tAngCorr", "tAngCorr");
    tAngCorr->Branch("fPx1", &fPx1, "fPx1/D");
    tAngCorr->Branch("fPy1", &fPy1, "fPy1/D");
    tAngCorr->Branch("fPz1", &fPz1, "fPz1/D");
    tAngCorr->Branch("fPx2", &fPx2, "fPx2/D");
    tAngCorr->Branch("fPy2", &fPy2, "fPy2/D");
    tAngCorr->Branch("fPz2", &fPz2, "fPz2/D");
    tAngCorr->Branch("fPt", &fPt, "fPt/D");
    tAngCorr->Branch("fM", &fM, "fM/D");

    TLorentzVector fMuon1;
    TLorentzVector fMuon2;

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(-1,-1)){
            fMuon1.SetPtEtaPhiM(fPt1, fEta1, fPhi1, fMuonMass);
            fMuon2.SetPtEtaPhiM(fPt2, fEta2, fPhi2, fMuonMass);
            // Get the momentum coordinates
            fPx1 = fMuon1.Px();
            fPy1 = fMuon1.Py();
            fPz1 = fMuon1.Pz();
            fPx2 = fMuon2.Px();
            fPy2 = fMuon2.Py();
            fPz2 = fMuon2.Pz();
            // Fill the tree
            tAngCorr->Fill();
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}