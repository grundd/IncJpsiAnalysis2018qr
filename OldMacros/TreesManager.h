// TreesManager.h
// David Grund, Sep 7, 2021
// Definition of the functions that can be used to connect variables from the tree:
// fTreeJPsi, fTreeJPsiMCRec and fTreeJPsiMCGen
// And of the functions that enables to perform selections of events

#include <stdio.h> // printf
#include "TTree.h"

Int_t fRunNumber;
TString *fTriggerName = NULL;
Bool_t fTriggerInputsMC[11];
Double_t fTrk1SigIfMu, fTrk1SigIfEl, fTrk2SigIfMu, fTrk2SigIfEl;
Double_t fPt, fM, fY, fPhi;
Double_t fPt1, fPt2, fEta1, fEta2, fPhi1, fPhi2, fQ1, fQ2;
Double_t fZNA_energy, fZNC_energy;
Double_t fZNA_time[4], fZNC_time[4];
Double_t fV0A_time, fV0C_time, fADA_time, fADC_time;
Int_t fV0A_dec, fV0C_dec, fADA_dec, fADC_dec;
Bool_t fMatchingSPD;
Double_t fPtGen, fYGen, fMGen, fPhiGen;

void ConnectTreeVariables(TTree *t){
    // Set branch addresses
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

    Printf("Variables from %s connected.", t->GetName());
    return;
}

void ConnectTreeVariablesMCRec(TTree *t){
    // Set branch addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    t->SetBranchAddress("fTriggerInputsMC", &fTriggerInputsMC);
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
    // MC kinematics on generated level
    t->SetBranchAddress("fPtGen", &fPtGen);
    t->SetBranchAddress("fPhiGen", &fPhiGen);
    t->SetBranchAddress("fYGen", &fYGen);
    t->SetBranchAddress("fMGen", &fMGen);

    Printf("Variables from %s connected.", t->GetName());
    return;
}

void ConnectTreeVariablesMCGen(TTree *t){
    // Set branch addresses
    // Basic things:
    t->SetBranchAddress("fRunNumber", &fRunNumber);
    // MC kinematics on generated level
    t->SetBranchAddress("fPtGen", &fPtGen);
    t->SetBranchAddress("fPhiGen", &fPhiGen);
    t->SetBranchAddress("fYGen", &fYGen);
    t->SetBranchAddress("fMGen", &fMGen);

    Printf("Variables from %s connected.", t->GetName());
    return;
}

Bool_t EventPassed(Int_t iMassCut = 0, Int_t iPtCut = 0){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD

    // 5) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;

    // 8) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;

    // 10) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;

    // 11) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case -1: // No inv mass cut
            bMassCut = kTRUE;
            break;
        case 0:
            if(fM > 2.2 && fM < 4.5) bMassCut = kTRUE;
            break;
        case 1:
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
    }
    if(!bMassCut) return kFALSE;

    // 13) Transverse momentum cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // No pt cut
            bPtCut = kTRUE;
            break;
        case 0: // Incoherent-enriched sample
            if(fPt > 0.20) bPtCut = kTRUE;
            break;
        case 1: // Coherent-enriched sample
            if(fPt < 0.11) bPtCut = kTRUE;
            break;
        case 2: // Total sample (pt < 2.0 GeV/c)
            if(fPt < 2.00) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCRec(Int_t iMassCut = -1, Int_t iPtCut = -1, Int_t iPtBin = -1){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    Bool_t CCUP31 = kFALSE;
    if(
        !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
        !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
        !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
        !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
        fTriggerInputsMC[10] &&  //  0STG (SPD topological)
        fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
    ) CCUP31 = kTRUE;
    if(!CCUP31) return kFALSE;

    // 5) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;

    // 8) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;

    // 10) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;

    // 11) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case -1: // No inv mass cut
            bMassCut = kTRUE;
            break;
        case 0:
            if(fM > 2.2 && fM < 4.5) bMassCut = kTRUE;
            break;
        case 1:
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
    }
    if(!bMassCut) return kFALSE;

    // 13) Transverse momentum cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // No pt cut
            bPtCut = kTRUE;
            break;
        case 0: // Incoherent-enriched sample
            if(fPt > 0.20) bPtCut = kTRUE;
            break;
        case 1: // Coherent-enriched sample
            if(fPt < 0.11) bPtCut = kTRUE;
            break;
        case 2: // Total sample (pt < 2.0 GeV/c)
            if(fPt < 2.00) bPtCut = kTRUE;
            break;
        case 3: // Sample with pt from 0.2 to 1 GeV/c 
            if(fPt > 0.20 && fPt < 1.00) bPtCut = kTRUE;
            break;
        case 4: // Pt bins
            //if(fPt > ptBoundaries[iPtBin-1] && fPt < ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCGen(Int_t iPtCut = -1, Int_t iPtBin = -1){
    // 1) Dilepton rapidity |y| < 0.8
    if(!(abs(fYGen) < 0.8)) return kFALSE;    

    // 2) Transverse momentum cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case -1: // No pt cut
            bPtCut = kTRUE;
            break;
        case 0: // Pt bins
            //if(fPt > ptBoundaries[iPtBin-1] && fPt < ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}