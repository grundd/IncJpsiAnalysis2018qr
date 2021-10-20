// Plot2DHistogramPtAndT.c
// David Grund, Oct 19, 2021

// cpp headers
#include <fstream>
#include <stdio.h> // printf
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

Double_t fPt2Gm, fPt2VM, fPt2Pm;
Double_t fPtGm, fPtVM, fPtPm;
Int_t nGenEv = 500000;
Int_t nBins = 1000;
Double_t fPt2Low, fPt2Upp;

TString TxtPtGamma = "SL_simulations_10-19-2021/PtGamma.txt";
TString TxtPtVMPom = "SL_simulations_10-19-2021/PtVMpomeron.txt";

void PlotResults(Int_t opt);
void PrepareTree();

void Plot2DHistogramPtAndT(){

    //PrepareTree();

    PlotResults(1);

    return;
}

void PlotResults(Int_t opt){

    if(opt == 0){
        fPt2Low = 0.00; // GeV^2
        fPt2Upp = 2.56; // GeV^2
    } else if(opt == 1){
        fPt2Low = 0.04; // GeV^2
        fPt2Upp = 1.00; // GeV^2
    }

    TFile *f = TFile::Open("SL_simulations_10-19-2021/tree.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());

    TList *l = (TList*) f->Get("TreeList");
    if(l) Printf("List %s loaded.", l->GetName()); 

    TTree *tPtGamma = (TTree*)l->FindObject("tPtGamma");
    if(tPtGamma) Printf("Tree %s loaded.", tPtGamma->GetName());

    TTree *tPtVMPom = (TTree*)l->FindObject("tPtVMPom");
    if(tPtVMPom) Printf("Tree %s loaded.", tPtVMPom->GetName());

    tPtGamma->SetBranchAddress("fPt2Gm", &fPt2Gm);
    tPtGamma->SetBranchAddress("fPtGm", &fPtGm);

    tPtVMPom->SetBranchAddress("fPt2VM", &fPt2VM);
    tPtVMPom->SetBranchAddress("fPt2Pm", &fPt2Pm);
    tPtVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtVMPom->SetBranchAddress("fPtPm", &fPtPm);

    TH2D *Hist = new TH2D("Hist", "#it{t} vs #it{p}_{T}^{2} of J/#psi", nBins, fPt2Low, fPt2Upp, nBins, fPt2Low, fPt2Upp);
    // on a horizontal axis: J/psi transverse momentum squared, i.e. fPt2VM
    // on a vertical axis: Mandelstam t, i.e. fPt2Pm

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++){
        tPtGamma->GetEntry(iEntry);
        tPtVMPom->GetEntry(iEntry);
        //Printf("fPt2VM: %.4f, fPt2Gm+fPt2Pm: %.4f, fPt2Gm: %.4f, fPt2Pm: %.4f", fPt2VM, fPt2Gm+fPt2Pm, fPt2Gm, fPt2Pm);

        Hist->Fill(fPt2VM, fPt2Pm);
    }

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetLogz();
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.145);
    c->SetRightMargin(0.11);
    c->SetLeftMargin(0.1);

    // a vertical axis
    Hist->GetYaxis()->SetTitle("#it{p}_{T, pom}^{2} (GeV/#it{c})");
    Hist->GetYaxis()->SetTitleSize(0.05);
    Hist->GetYaxis()->SetLabelSize(0.05);
    Hist->GetYaxis()->SetTitleOffset(0.915);
    Hist->GetYaxis()->SetDecimals(1);
    // a horizontal axis
    Hist->GetXaxis()->SetTitle("#it{p}_{T, J/#psi}^{2} (GeV/#it{c})");
    Hist->GetXaxis()->SetTitleSize(0.05);
    Hist->GetXaxis()->SetTitleOffset(1.3);
    Hist->GetXaxis()->SetLabelSize(0.05);
    Hist->GetXaxis()->SetLabelOffset(0.02);
    Hist->GetXaxis()->SetDecimals(1);
    // draw the histogram
    Hist->GetZaxis()->SetLabelSize(0.05);
    Hist->Draw("COLZ");

    return;
}

void PrepareTree(){

    TTree *tPtGamma = new TTree("tPtGamma", "tPtGamma");
    tPtGamma->Branch("fPt2Gm", &fPt2Gm, "fPt2Gm/D");
    tPtGamma->Branch("fPtGm", &fPtGm, "fPtGm/D");

    TTree *tPtVMPom = new TTree("tPtVMPom", "tPtVMPom");
    tPtVMPom->Branch("fPt2VM", &fPt2VM, "fPt2VM/D");
    tPtVMPom->Branch("fPt2Pm", &fPt2Pm, "fPt2Pm/D");
    tPtVMPom->Branch("fPtVM", &fPtVM, "fPtVM/D");
    tPtVMPom->Branch("fPtPm", &fPtPm, "fPtPm/D");

    ifstream ifs;
    ifs.open(TxtPtGamma.Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < nGenEv; i++){
            ifs >> fPt2Gm;
            fPtGm = TMath::Sqrt(fPt2Gm);
            tPtGamma->Fill();
        }
        ifs.close();
    } else {
        Printf("File %s missing. Terminating.", TxtPtGamma.Data());
        return;
    }

    ifs.open(TxtPtVMPom.Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < nGenEv; i++){
            // see the src/eventfilewriter.cpp in the STARlight source code
            ifs >> fPt2VM >> fPt2Pm;
            fPtVM = TMath::Sqrt(fPt2VM);
            fPtPm = TMath::Sqrt(fPt2Pm);
            tPtVMPom->Fill();
        }
        ifs.close();
    } else {
        Printf("File %s missing. Terminating.", TxtPtVMPom.Data());
        return;
    }

    TList *l = new TList();
    l->Add(tPtGamma);
    l->Add(tPtVMPom);

    TFile *f = new TFile("SL_simulations_10-19-2021/tree.root","RECREATE");
    l->Write("TreeList", TObject::kSingleKey);
    f->ls();

    return;
}