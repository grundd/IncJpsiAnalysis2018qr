// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"

TString TxtPtGamma = "DependenceOnT/SL_simulations_11-15-2021_coherent/PtGamma.txt";
TString TxtPtVMPom = "DependenceOnT/SL_simulations_11-15-2021_coherent/PtVMpomeron.txt";

Double_t fPt2Gm, fPt2VM, fPt2Pm;
Double_t fPtGm, fPtVM, fPtPm;
Int_t nGenEv = 500000;

Double_t SigmaCohSL = 10.089;
const Int_t nBins = 1000;

void PrepareTree();
void CalculateSTARlightPredictionsCoherent();

void PhenoPredictions_SL(){

    PrepareTree();

    CalculateSTARlightPredictionsCoherent();

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

    TFile *f = new TFile("DependenceOnT/SL_simulations_11-15-2021_coherent/tree.root","RECREATE");
    l->Write("TreeList", TObject::kSingleKey);
    f->ls();

    return;
}

void CalculateSTARlightPredictionsCoherent(){

    TH1D *hSL = new TH1D("hSL","hSL",nBins,0.0,0.8);

    TFile *f = TFile::Open("DependenceOnT/SL_simulations_11-15-2021_coherent/tree.root","read");
    if(f) Printf("File %s loaded.", f->GetName());

    TList *l = (TList*) f->Get("TreeList");
    if(l) Printf("List %s loaded.", l->GetName()); 

    TTree *tPtVMPom = (TTree*)l->FindObject("tPtVMPom");
    if(tPtVMPom) Printf("Tree %s loaded.", tPtVMPom->GetName());

    tPtVMPom->SetBranchAddress("fPt2Pm", &fPt2Pm);
    tPtVMPom->SetBranchAddress("fPt2VM", &fPt2VM);

    Int_t nGenEv = tPtVMPom->GetEntries();
    Printf("Tree contains %i entries.", nGenEv);

    for(Int_t iEntry = 0; iEntry < nGenEv; iEntry++){
        tPtVMPom->GetEntry(iEntry);
        hSL->Fill(fPt2Pm);
    }

    Int_t opt = 1;

    if(opt == 0){
        // option 0 ~ Guillermo
        Double_t L = nGenEv / SigmaCohSL;
        for(Int_t i = 1; i <= nBins; i++){
            hSL->SetBinContent(i, hSL->GetBinContent(i) / L);
        }
        // Check the value of the integral
        Double_t Integral = 0;
        for(Int_t iBin = 1; iBin <= nBins; iBin++){
            Integral += hSL->GetBinContent(iBin) * (hSL->GetBinLowEdge(iBin+1) - hSL->GetBinLowEdge(iBin));
        }
        Printf("Integral is %.5f (%.5f).\n", Integral, hSL->Integral("width"));

    } else if(opt == 1){
        // option 1 ~ me
        hSL->Scale(SigmaCohSL/hSL->Integral(), "width");
        // Check the value of the integral
        Double_t Integral = 0;
        for(Int_t iBin = 1; iBin <= nBins; iBin++){
            Integral += hSL->GetBinContent(iBin) * (hSL->GetBinLowEdge(iBin+1) - hSL->GetBinLowEdge(iBin));
        }
        Printf("Integral is %.5f (%.5f).\n", Integral, hSL->Integral("width"));

    }

    Bool_t plot = kTRUE;
    if(plot){
        TCanvas *cSL = new TCanvas("cSL","cSL",900,600);
        cSL->SetLogy();
        hSL->Draw();
        hSL->GetXaxis()->SetRangeUser(0.0,0.012);
    }

    return;
}