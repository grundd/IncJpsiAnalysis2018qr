// ZNPlots_Guillermo.c
// David Grund, Feb 25, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
// my headers
#include "AnalysisManager.h"

Double_t fZNA_n, fZNC_n; // number of neutrons (ZN energy divided by 2510 GeV)

void PlotHistogram(Int_t iPtRange, Int_t iMassRange);
// iPtRange = 0 => pT in (0, 0.11) GeV/c
//          = 1 => pT in (0.2,0.4) GeV/c
//          = 2 => pT in (0.4,1.0) GeV/c
// iMassRange = 0 => M in (3.0,3.2) GeV/c^2
//            = 1 => M in (2.2,2.9) GeV/c^2
void PrepareData();

void ZNPlots_Guillermo(){

    //PrepareData();

    for(Int_t iM = 0; iM < 2; iM++){
        for(Int_t iPt = 0; iPt < 3; iPt++){
            PlotHistogram(iPt,iM);
        }
    }

    return;
}

void PlotHistogram(Int_t iPtRange, Int_t iMassRange){

    TFile *f_in = TFile::Open("Trees/ZNPlots_Guillermo/ZNPlots_Guillermo.root", "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("tZN_n"));
    if(t_in) Printf("Input tree loaded.");    

    t_in->SetBranchAddress("fZNA_time", &fZNA_time);
    t_in->SetBranchAddress("fZNC_time", &fZNC_time);
    t_in->SetBranchAddress("fZNA_n", &fZNA_n);
    t_in->SetBranchAddress("fZNC_n", &fZNC_n);
    t_in->SetBranchAddress("fPt", &fPt);
    t_in->SetBranchAddress("fM", &fM);

    TH1F *hZNA = new TH1F("hZNA","hZNA",50,0.,10.);
    TH1F *hZNC = new TH1F("hZNC","hZNC",50,0.,10.);

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        if(iMassRange == 0 && !(fM > 3.0 && fM < 3.2)) continue;
        if(iMassRange == 1 && !(fM > 2.2 && fM < 2.9)) continue;
        if(iPtRange == 0 && !(fPt > 0.00 && fPt < 0.11)) continue;
        if(iPtRange == 1 && !(fPt > 0.20 && fPt < 0.40)) continue;
        if(iPtRange == 2 && !(fPt > 0.40 && fPt < 1.00)) continue;

        Bool_t fZNA_hit = kFALSE;
        Bool_t fZNC_hit = kFALSE;
        for(Int_t i = 0; i < 4; i++){
            // hit in ZNA but not in ZNC
            if(TMath::Abs(fZNA_time[i]) < 2) fZNA_hit = kTRUE;
            // hit in ZNC but not in ZNA
            if(TMath::Abs(fZNC_time[i]) < 2) fZNC_hit = kTRUE;
        }
        if(fZNA_hit && !fZNC_hit) hZNA->Fill(fZNA_n);
        if(fZNC_hit && !fZNA_hit) hZNC->Fill(fZNC_n);
    }    

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.11);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.04);
    //c->SetLogy();
    // X-axis
    hZNA->GetXaxis()->SetTitle("# of neutrons (ZN energy/2510 GeV)");
    hZNA->GetXaxis()->SetTitleSize(0.05);
    hZNA->GetXaxis()->SetLabelSize(0.05);
    // Y-axis
    hZNA->GetYaxis()->SetTitle("Counts per 0.2");
    hZNA->GetYaxis()->SetTitleSize(0.05);
    hZNA->GetYaxis()->SetLabelSize(0.05);
    hZNA->GetYaxis()->SetTitleOffset(1.2);
    // Style hist ZNA
    hZNA->SetLineColor(kRed);
    hZNA->SetLineWidth(1);
    hZNA->SetMarkerStyle(21);
    hZNA->SetMarkerColor(kRed);
    hZNA->SetMarkerSize(0.5);
    // Style hist ZNC
    hZNC->SetLineColor(kBlue);
    hZNC->SetLineWidth(1);
    hZNC->SetMarkerStyle(21);
    hZNC->SetMarkerColor(kBlue);
    hZNC->SetMarkerSize(0.5);
    // Draw
    hZNA->Draw("E0");
    hZNC->Draw("SAME E0");
    // Legend 1
    TLegend *l1 = new TLegend(0.215,0.75,0.99,0.96);
    l1->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    if(iMassRange == 0) l1->AddEntry((TObject*)0,"3.0 < #it{m}_{#mu#mu} < 3.2 GeV/#it{c}^{2}","");
    if(iMassRange == 1) l1->AddEntry((TObject*)0,"2.2 < #it{m}_{#mu#mu} < 2.9 GeV/#it{c}^{2}","");
    if(iPtRange == 0) l1->AddEntry((TObject*)0,"#it{p}_{T} < 0.11 GeV/#it{c}","");
    if(iPtRange == 1) l1->AddEntry((TObject*)0,"0.2 < #it{p}_{T} < 0.4 GeV/#it{c}","");
    if(iPtRange == 2) l1->AddEntry((TObject*)0,"0.4 < #it{p}_{T} < 1.0 GeV/#it{c}","");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.4,0.60,0.82,0.72);
    l2->AddEntry(hZNA,Form("ZNA (total: %.0f events)", hZNA->Integral()),"PL");
    l2->AddEntry(hZNC,Form("ZNC (total: %.0f events)", hZNC->Integral()),"PL");
    l2->SetTextSize(0.05);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    TString str_out = Form("Results/ZNPlots_Guillermo/Mass%i_pT%i", iMassRange, iPtRange);
    c->Print((str_out + ".pdf").Data());
    c->Print((str_out + ".png").Data());

    return;
}

void PrepareData(){

    TString str_file = "Trees/AnalysisData_pass3/AnalysisResults.root";
    TFile *f_in = TFile::Open(str_file.Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get("AnalysisOutput/fTreeJpsi"));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariables(t_in, kTRUE);

    // Create new data tree with applied cuts
    TFile f_out("Trees/ZNPlots_Guillermo/ZNPlots_Guillermo.root","RECREATE");

    // pT in (0, 0.11) GeV/c, M in (3.0,3.2) GeV/c^2
    TTree *tZN_n = new TTree("tZN_n", "tZN_n");
    tZN_n->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    tZN_n->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    tZN_n->Branch("fZNA_n", &fZNA_n, "fZNA_n/D");
    tZN_n->Branch("fZNC_n", &fZNC_n, "fZNC_n/D");
    tZN_n->Branch("fPt", &fPt, "fPt/D");
    tZN_n->Branch("fM",  &fM,  "fM/D");

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        // pass3:
        // 0) fEvent non-empty
        // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
        // 2) Central UPC trigger CCUP31:
        // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
        // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD   
        // 3) At least two tracks associated with the vertex
        if(fVertexContrib < 2) continue;
        // 4) Distance from the IP lower than 15 cm
        if(fVertexZ > 15) continue;
        // 5a) ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;
        // 5b) ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;
        // 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;
        // 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;
        // 7) SPD cluster matches FOhits
        if(!(fMatchingSPD == kTRUE)) continue;
        // 8) Muon pairs only
        if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) continue;
        // 9) Dilepton rapidity |y| < 0.8
        if(!(abs(fY) < 0.8)) continue;
        // 10) Pseudorapidity of both tracks |eta| < 0.8
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // 11) Tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;
        // 12) Invariant mass cut
        if(!(fM > 2.0 && fM < 4.0)) continue;
        // 13) Transverse momentum cut
        if(!(fPt < 1.0)) continue;
        
        // Save event to the tree
        fZNA_n = fZNA_energy / 2510.;
        fZNC_n = fZNC_energy / 2510.;
        tZN_n->Fill();
    }    

    f_out.Write("",TObject::kWriteDelete);

    return;
}