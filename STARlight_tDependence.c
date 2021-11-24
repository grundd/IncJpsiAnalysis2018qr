// STARlight_tDependence.c
// David Grund, Nov 21, 2021

// root headers
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
// my headers
#include "STARlight_Utilities.h"

TString str_folder = "";
Double_t sig_gPb, rap_cut;
Double_t pT_low, pT_upp;
Double_t nEv_tot = 0;
const Int_t t_nBins = 50;
Double_t sigma[t_nBins] = { 0 };
Double_t t_abs[t_nBins] = { 0 };

Int_t iProd = 1;

void CalculateDependence();

void STARlight_tDependence(){

    if(iProd == 1){
        str_folder = "Trees/STARlight/inc_tDep/";
        sig_gPb = 0.015499; // see output.txt, line y = +0.000
        rap_cut = 0.01;
        pT_low = 0.15;
        pT_upp = 1.02;
    } 

    CalculateDependence();

    return;
}

void CalculateDependence(){

    Double_t t_low = TMath::Power(pT_low, 2);
    Double_t t_upp = TMath::Power(pT_upp, 2);

    TFile *fSL = TFile::Open((str_folder + "trees_starlight.root").Data(), "read");
    if(fSL) Printf("Input file SL loaded.");

    TTree *tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
    if(tSL) Printf("Tree %s loaded.", tSL->GetName());
    ConnectTreeVariables_tSL(tSL);

    TFile *fPt = TFile::Open((str_folder + "trees_tPtGammaVMPom.root").Data(), "read");
    if(fPt) Printf("Input file Pt loaded.");

    TList *lPt = (TList*) fPt->Get("TreeList");
    if(lPt) Printf("List %s loaded.", lPt->GetName()); 

    TTree *tPtGammaVMPom = (TTree*)lPt->FindObject("tPtGammaVMPom");
    if(tPtGammaVMPom) Printf("Tree %s loaded.", tPtGammaVMPom->GetName());
    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);

    Printf("STARlight tree contains %lli entries.", tSL->GetEntries());
    Printf("tPtGammaVMPom tree contains %lli entries.", tPtGammaVMPom->GetEntries());

    TH1D *hSigmaPhotoNuc = new TH1D("hSigmaPhotoNuc","hSigmaPhotoNuc", t_nBins, t_low, t_upp);

    for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++){
        tSL->GetEntry(iEntry);
        tPtGammaVMPom->GetEntry(iEntry);

        // if the values differ by more than 10% => something wrong
        if(TMath::Abs(fPtVM - parent->Pt())/fPtVM > 0.05){
            Printf("Entry no. %i: fPtVM = %.6f, from SL = %.6f", iEntry+1, fPtVM, parent->Pt());
        }

        if(TMath::Abs(parent->Rapidity()) < rap_cut){
            nEv_tot++;
            if(parent->Pt() > pT_low && parent->Pt() < pT_upp){
                hSigmaPhotoNuc->Fill(fPtPm*fPtPm);
            }
        }
    }

    Printf("No. of events with |y| < %.2f: %.0f", rap_cut, nEv_tot);
    Printf("No. of entries in histogram: %.0f", hSigmaPhotoNuc->GetEntries());
    // Normalize the histogram
    hSigmaPhotoNuc->Scale(1.0/nEv_tot);
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));
    hSigmaPhotoNuc->Scale(sig_gPb);
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));
    hSigmaPhotoNuc->Scale(1.0, "width");
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral());
    Printf("Integral is %.5f", hSigmaPhotoNuc->Integral("width"));
    // Check the value of the integral
    Double_t integral = 0;
    for(Int_t iBin = 1; iBin <= t_nBins; iBin++){
        integral += hSigmaPhotoNuc->GetBinContent(iBin) * (hSigmaPhotoNuc->GetBinLowEdge(iBin+1) - hSigmaPhotoNuc->GetBinLowEdge(iBin));
    }
    Printf("Manually calculated integral is %.5f", integral);

    // Fill arrays
    for(Int_t iBin = 1; iBin <= t_nBins; iBin++){
        t_abs[iBin-1] = hSigmaPhotoNuc->GetBinCenter(iBin);
        sigma[iBin-1] = hSigmaPhotoNuc->GetBinContent(iBin);
    }
    // Define TGraph
    TGraph *grSL = new TGraph(t_nBins, t_abs, sigma);
    grSL->SetLineStyle(1);
    grSL->SetLineColor(215);
    grSL->SetLineWidth(2);
    grSL->GetXaxis()->SetRangeUser(t_low,t_upp);

    grSL->Print();
    
    TCanvas *c1 = new TCanvas("c1","c1",900,600);
    c1->SetLogy();
    grSL->Draw("AL");

    nEv_tot = 0;

    // Print the results to text file
    TString str_out = "";
    if(iProd == 1) str_out = "PhenoPredictions/STARlight/inc_tDep.txt";
    ofstream outfile(str_out.Data());
    for(Int_t iBin = 0; iBin < t_nBins; iBin++){
        outfile << Form("%.4f \t%.6f \n", t_abs[iBin], sigma[iBin]);
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str_out.Data());
    
    return;
}