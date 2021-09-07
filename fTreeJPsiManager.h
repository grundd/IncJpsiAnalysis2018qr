// TreesManager.h
// David Grund, Sep 7, 2021
// Definition of the function that can be used to connect variables from the fTreeJPsi tree
// And of the function that enables to perform selections of events

#include <stdio.h> // printf
#include "TTree.h"

Int_t fRunNumber;
TString *fTriggerName = NULL;
Double_t fTrk1SigIfMu, fTrk1SigIfEl, fTrk2SigIfMu, fTrk2SigIfEl;
Double_t fPt, fM, fY, fPhi;
Double_t fPt1, fPt2, fEta1, fEta2, fPhi1, fPhi2, fQ1, fQ2;
Double_t fZNA_energy, fZNC_energy;
Double_t fZNA_time[4], fZNC_time[4];
Double_t fV0A_time, fV0C_time, fADA_time, fADC_time;
Int_t fV0A_dec, fV0C_dec, fADA_dec, fADC_dec;
Bool_t fMatchingSPD;

void ConnectTreeVariables(TTree *t){
    // Set Branch Addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    t->SetBranchAddress("fTriggerName", &fTriggerName);
    // PID, sigmas:
    t->SetBranchAddress("fTrk1SigIfMu", &fTrk1SigIfMu);
    t->SetBranchAddress("fTrk1SigIfEl", &fTrk1SigIfEl);
    t->SetBranchAddress("fTrk2SigIfMu", &fTrk2SigIfMu);
    t->SetBranchAddress("fTrk2SigIfEl", &fTrk2SigIfEl);
    // Kinematics:
    t->SetBranchAddress("fPt", &fPt);
    t->SetBranchAddress("fPhi", &fPhi);
    t->SetBranchAddress("fY", &fY);
    t->SetBranchAddress("fM", &fM);
    // Two tracks:
    t->SetBranchAddress("fPt1", &fPt1);
    t->SetBranchAddress("fPt2", &fPt2);
    t->SetBranchAddress("fEta1", &fEta1);
    t->SetBranchAddress("fEta2", &fEta2);
    t->SetBranchAddress("fPhi1", &fPhi1);
    t->SetBranchAddress("fPhi2", &fPhi2);
    t->SetBranchAddress("fQ1", &fQ1);
    t->SetBranchAddress("fQ2", &fQ2);
    // ZDC:
    t->SetBranchAddress("fZNA_energy", &fZNA_energy);
    t->SetBranchAddress("fZNC_energy", &fZNC_energy);
    t->SetBranchAddress("fZNA_time", &fZNA_time);
    t->SetBranchAddress("fZNC_time", &fZNC_time);
    // V0:
    t->SetBranchAddress("fV0A_dec", &fV0A_dec);
    t->SetBranchAddress("fV0C_dec", &fV0C_dec);
    t->SetBranchAddress("fV0A_time", &fV0A_time);
    t->SetBranchAddress("fV0C_time", &fV0C_time);
    // AD:
    t->SetBranchAddress("fADA_dec", &fADA_dec);
    t->SetBranchAddress("fADC_dec", &fADC_dec);
    t->SetBranchAddress("fADA_time", &fADA_time);
    t->SetBranchAddress("fADC_time", &fADC_time);
    // Matching SPD clusters with FOhits:
    t->SetBranchAddress("fMatchingSPD", &fMatchingSPD);

    Printf("Variables from fTreeJPsi connected.");
    return;
}

Bool_t EventPassed(){

    return kTRUE;
}