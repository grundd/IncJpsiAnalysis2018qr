// RaphaelleComments.C
// Macros to answer Raphaelle's questions

// root headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"

void R02_AcceptanceDimuons(Bool_t etaCut);
void R15_DeltaPhiVsPt(Bool_t fromMC);

void RaphaelleComments(){

    //R02_AcceptanceDimuons(kTRUE);
    //R02_AcceptanceDimuons(kFALSE);

    //R15_DeltaPhiVsPt(kTRUE);
    //R15_DeltaPhiVsPt(kFALSE);

    return;
}

void R02_AcceptanceDimuons(Bool_t etaCut){

    TFile *f = NULL;
    TTree *t = NULL;    
    f = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(f) Printf("Input data loaded.");
    t = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsi"));
    if(t) Printf("Input tree loaded.");
    ConnectTreeVariables(t); 

    TH1D *hAccDimuons = new TH1D("hAccDimuons","hAccDimuons",200,-1.,1.);

    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++){
        t->GetEntry(iEntry);
        
        // 5) SPD cluster matches FOhits
        if(!(fMatchingSPD == kTRUE)) continue;
        // 6a) ADA offline veto (no effect on MC)
        else if(!(fADA_dec == 0)) continue;
        // 6b) ADC offline veto (no effect on MC)
        else if(!(fADC_dec == 0)) continue;
        // 7a) V0A offline veto (no effect on MC)
        else if(!(fV0A_dec == 0)) continue;
        // 7b) V0C offline veto (no effect on MC)
        else if(!(fV0C_dec == 0)) continue;
        // 8) dilepton rapidity |y| < 0.8
        // (skipped)
        // 9) pseudorapidity of both tracks |eta| < 0.8
        else if(etaCut && !(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // 10) tracks have opposite charges
        else if(!(fQ1 * fQ2 < 0)) continue;
        // 11) muon pairs only
        else if(!(TMath::Power(fTrk1SigIfMu,2) + TMath::Power(fTrk2SigIfMu,2) < TMath::Power(fTrk1SigIfEl,2) + TMath::Power(fTrk2SigIfEl,2))) continue;
        // 12) invariant mass between 2.2 and 4.5 GeV/c^2
        else if(!(fM > 2.2 && fM < 4.5)) continue;
        // 13) transverse momentum cut
        // (skipped)
        // event passed all the selections =>
        hAccDimuons->Fill(fY);
    }        

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    // horizontal axis
    hAccDimuons->GetXaxis()->SetTitle("#it{y} (-)");
    hAccDimuons->GetXaxis()->SetTitleSize(0.05);
    hAccDimuons->GetXaxis()->SetTitleOffset(1.18);
    hAccDimuons->GetXaxis()->SetLabelSize(0.05);
    hAccDimuons->GetXaxis()->SetLabelOffset(0.015);
    hAccDimuons->GetXaxis()->SetDecimals(1);
    // vertical axis
    hAccDimuons->GetYaxis()->SetTitle("Counts per bin");
    hAccDimuons->GetYaxis()->SetTitleSize(0.05);
    hAccDimuons->GetYaxis()->SetTitleOffset(0.99);
    hAccDimuons->GetYaxis()->SetLabelSize(0.05);
    //hAccDimuons->GetYaxis()->SetDecimals(1);

    hAccDimuons->Draw();

    if(etaCut){
        c->Print("Results/RaphaelleComments/AccDimuons.pdf");
        c->Print("Results/RaphaelleComments/AccDimuons.png");
    } else {
        c->Print("Results/RaphaelleComments/AccDimuons_noEtaCut.pdf");
        c->Print("Results/RaphaelleComments/AccDimuons_noEtaCut.png");
    }


    return;    
}

void R15_DeltaPhiVsPt(Bool_t fromMC)
{
    TFile *f = NULL;
    TTree *t = NULL;
    if(fromMC){
        f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
        if(f) Printf("MC rec file loaded.");
        t = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsiMCRec"));
        if(t) Printf("MC rec tree loaded.");
        ConnectTreeVariablesMCRec(t);
    } else {
        f = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
        if(f) Printf("Input data loaded.");
        t = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsi"));
        if(t) Printf("Input tree loaded.");
        ConnectTreeVariables(t);     
    }

    // horizontal axis = pT, vertical axis = DeltaPhi
    TH2D *hDeltaPhiVsPt = new TH2D("hDeltaPhiVsPt","hDeltaPhiVsPt",160,0.,1.6,80,0.,3.2);

    Double_t DeltaPhi;
    ///*
    for(Int_t iEntry = 0; iEntry < t->GetEntries(); iEntry++){
        t->GetEntry(iEntry);

        // 10) tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;
        // 11) muon pairs only
        else if(!(TMath::Power(fTrk1SigIfMu,2) + TMath::Power(fTrk2SigIfMu,2) < TMath::Power(fTrk1SigIfEl,2) + TMath::Power(fTrk2SigIfEl,2))) continue;
        // 12) invariant mass between 2.2 and 4.5 GeV/c^2
        else if(!(fM > 2.2 && fM < 4.5)) continue;

        DeltaPhi = TMath::Abs(fPhi1 - fPhi2);
        if(DeltaPhi > TMath::Pi()) DeltaPhi = 2*TMath::Pi() - DeltaPhi;
        hDeltaPhiVsPt->Fill(fPt,DeltaPhi);
    }
    //*/

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.11);
    c->SetLeftMargin(0.95);
    c->SetLogz();

    // horizontal axis
    hDeltaPhiVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDeltaPhiVsPt->GetXaxis()->SetTitleSize(0.05);
    hDeltaPhiVsPt->GetXaxis()->SetTitleOffset(1.18);
    hDeltaPhiVsPt->GetXaxis()->SetLabelSize(0.05);
    hDeltaPhiVsPt->GetXaxis()->SetLabelOffset(0.015);
    hDeltaPhiVsPt->GetXaxis()->SetDecimals(1);
    // vertical axis
    hDeltaPhiVsPt->GetYaxis()->SetTitle("#Delta#phi (-)");
    hDeltaPhiVsPt->GetYaxis()->SetTitleSize(0.05);
    hDeltaPhiVsPt->GetYaxis()->SetTitleOffset(0.9);
    hDeltaPhiVsPt->GetYaxis()->SetLabelSize(0.05);
    hDeltaPhiVsPt->GetYaxis()->SetDecimals(1);

    hDeltaPhiVsPt->Draw("COLZ");
    if(fromMC){
        c->Print("Results/RaphaelleComments/DeltaPhiVsPt_MC.pdf");
        c->Print("Results/RaphaelleComments/DeltaPhiVsPt_MC.png");
    } else {
        c->Print("Results/RaphaelleComments/DeltaPhiVsPt.pdf");
        c->Print("Results/RaphaelleComments/DeltaPhiVsPt.png");
    }

    return;    
}

