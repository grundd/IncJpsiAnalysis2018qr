// ElectronsMuonsPID.c
// David Grund, Dec 2, 2021

// root headers
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"

TString path_out = "Results/ElectronsMuonsPID/";

Int_t nBins = 100;
Double_t sigma_low = 0.0;
Double_t sigma_upp = 20.;
TH2D *hSigmasTPC = new TH2D("hSigmasTPC","hSigmasTPC",nBins,sigma_low,sigma_upp,nBins,sigma_low,sigma_upp);
// horizontal axis = sigma if electron
// vertical axis = sigma if muon

void Plot2DHistTPCSigmas();
Bool_t EventPassedLocal();

void ElectronsMuonsPID(){

    Plot2DHistTPCSigmas();

    return;
}

void Plot2DHistTPCSigmas(){

    // Load data
    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    ///*
    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        Double_t SigmaIfEls = TMath::Sqrt(fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl);
        Double_t SigmaIfMus = TMath::Sqrt(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu);
        if(EventPassedLocal()) hSigmasTPC->Fill(SigmaIfEls,SigmaIfMus);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    //*/

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(1);

    TCanvas *c = new TCanvas("c","c",900,800);
    c->SetGrid();
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.145);
    c->SetRightMargin(0.11);
    c->SetLeftMargin(0.13);

    // horizontal axis
    hSigmasTPC->GetXaxis()->SetTitle("#sqrt{#it{N}#sigma_{e}^{2}(+) + #it{N}#sigma_{e}^{2}(-)}");
    hSigmasTPC->GetXaxis()->SetTitleSize(0.05);
    hSigmasTPC->GetXaxis()->SetTitleOffset(1.15);
    hSigmasTPC->GetXaxis()->SetLabelSize(0.05);
    hSigmasTPC->GetXaxis()->SetDecimals(0);
    // vertical axis
    hSigmasTPC->GetYaxis()->SetTitle("#sqrt{#it{N}#sigma_{#mu}^{2}(+) + #it{N}#sigma_{#mu}^{2}(-)}");
    hSigmasTPC->GetYaxis()->SetTitleSize(0.05);
    hSigmasTPC->GetYaxis()->SetTitleOffset(1.15);
    hSigmasTPC->GetYaxis()->SetLabelSize(0.05);
    hSigmasTPC->GetYaxis()->SetDecimals(0);
    // Z-axis
    //hSigmasTPC->GetZaxis()->SetLabelSize(0.05);
    // Set ranges and draw
    hSigmasTPC->GetXaxis()->SetRangeUser(0.0,20.0);
    hSigmasTPC->GetYaxis()->SetRangeUser(0.0,20.0);
    hSigmasTPC->Draw("COLZ");

    // Draw dashed line y = x
    TLine *line = new TLine(0.0,0.0,20.0,20.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(9);
    line->Draw("SAME");

    // Legends
    Double_t x_el = 0.55;
    Double_t y_el = 0.85;
    TLegend *leg_el = new TLegend(x_el,y_el,x_el+0.18,y_el+0.06);
    leg_el->AddEntry((TObject*)0,Form("electrons"),""); 
    leg_el->SetTextSize(0.05);
    leg_el->SetBorderSize(0);
    leg_el->SetMargin(0.05);
    //leg_el->SetFillColor(kBlue);
    leg_el->Draw();
    Double_t x_mu = 0.72;
    Double_t y_mu = 0.7;
    TLegend *leg_mu = new TLegend(x_mu,y_mu,x_mu+0.135,y_mu+0.06);
    leg_mu->AddEntry((TObject*)0,Form("muons"),""); 
    leg_mu->SetTextSize(0.05);
    leg_mu->SetBorderSize(0);
    leg_mu->SetMargin(0.05);
    //leg_mu->SetFillColor(kBlue);
    leg_mu->Draw();

    c->Print((path_out + "Hist2DSigmasTPC.pdf").Data());
    c->Print((path_out + "Hist2DSigmasTPC.png").Data());

    return;
}

Bool_t EventPassedLocal(){

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
    // (no cut)

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;

    // 13) Transverse momentum cut
    if(!(fPt > 0.2 && fPt < 1.0)) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}