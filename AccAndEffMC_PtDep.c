// AccAndEffMC_PtDep.c
// David Grund, 20-09-2021
// To investigate pt dependence of AxE

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h" // gStyle
#include "TLegend.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

Int_t nBins = 0;
// From "AxE/PtBinEdges.txt":
Double_t edges[25] = {0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.68, 0.76, 0.84, 0.92, 1.00, 1.15, 1.30, 1.45, 1.60};
TH1D *hNRec = NULL;
TH1D *hNGen = NULL;
TH1D* hAxE = NULL;

// Global folder for saving files
TString GlobPath = "AxE_PtDep/binning1/";

Int_t nCuts = 16;
Bool_t cuts[16] = {
    0,  // 0) pt cut rec (kFALSE) or gen (kTRUE)
    0,  // 1) !0VBA (no signal in the V0A) 
    0,  // 2) !0VBC (no signal in the V0C)
    0,  // 3) !0UBA (no signal in the ADA)
    0,  // 4) !0UBC (no signal in the ADC)
    1,  // 5) 0STG (SPD topological)
    1,  // 6) 0OMU (TOF two hits topology)
    0,  // 7) SPD cluster matches FOhits
    0,  // 8) AD offline veto
    0,  // 9) V0 offline veto
    0,  // 10) Rapidity
    0,  // 11) Pseudorapidity
    0,  // 12) Opposite charges
    0,  // 13) Muons only
    0,  // 14) Inv mass
    1   // 15) Transverse momentum
};

void FillHistNRec();
void FillHistNGen();
Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut, Int_t iPtBin);
Bool_t EventPassedMCGen_AxEPtDep(Int_t iPtBin);
TString ConvertCutsToString(Bool_t *cuts);
void SaveToFile(TH1D* hist, TString name);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

void AccAndEffMC_PtDep(){

    // Get number of bins
    nBins = sizeof(edges) / sizeof(edges[0]) - 1;
    Printf("%i pt bins defined.", nBins);
    // Define the histograms
    hNRec = new TH1D("hNRec","N rec per bin",nBins,edges);
    hNGen = new TH1D("hNGen","N gen per bin",nBins,edges);

    FillHistNRec();
    FillHistNGen();

    hAxE = (TH1D*)hNRec->Clone("hAxE");
    hAxE->SetTitle("AxE per bin");
    hAxE->Sumw2();
    hAxE->Divide(hNGen);

    // Draw the histogram:
    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.145);
    // gStyle
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // Marker and line
    hAxE->SetMarkerStyle(21);
    hAxE->SetMarkerColor(kBlue);
    hAxE->SetMarkerSize(1.0);
    hAxE->SetLineColor(kBlue);
    hAxE->SetLineWidth(1.0);
    // Vertical axis
    hAxE->GetYaxis()->SetTitle("#it{N}_{rec}/#it{N}_{gen}");
    hAxE->GetYaxis()->SetTitleSize(0.056);
    hAxE->GetYaxis()->SetTitleOffset(1.3);
    hAxE->GetYaxis()->SetLabelSize(0.056);
    hAxE->GetYaxis()->SetDecimals(3);
    hAxE->GetYaxis()->SetRangeUser(0.0,hAxE->GetBinContent(1)*1.1);
    // Horizontal axis
    hAxE->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAxE->GetXaxis()->SetTitleSize(0.056);
    hAxE->GetXaxis()->SetTitleOffset(1.2);
    hAxE->GetXaxis()->SetLabelSize(0.056);
    hAxE->GetXaxis()->SetLabelOffset(0.015);
    hAxE->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    hAxE->Draw("P E1");
    // Legend
    TLegend *l = new TLegend(0.52,0.77,0.85,0.97);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->SetTextSize(0.056);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    // Save the figures and print the results to txt file
    TString CutConfiguration = ConvertCutsToString(cuts);
    TString path((GlobPath + "fig/" + CutConfiguration).Data());
    c->Print((path + ".pdf").Data());
    c->Print((path + ".png").Data());
    ofstream outfile((path + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << "Bin \tPtLow \tPtUpp \tAxE [%%] \tAxE_err [%%] \n";
    for(Int_t i = 1; i <= nBins; i++){
        outfile << i << "\t" << edges[i-1] << "\t" << edges[i] << "\t" << hAxE->GetBinContent(i)*100 << "\t\t" << hAxE->GetBinError(i)*100 << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", (path + ".txt").Data());

    // Compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Double_t ErrRoot = 0;
        Double_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++){
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hNRec->GetBinContent(i),hNGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // Cross-check: calculate the total value of AxE
    Double_t NRecTot = 0;
    Double_t NGenTot = 0;
    for(Int_t i = 1; i <= nPtBins; i++){
        NRecTot += hNRec->GetBinContent(i);
        NGenTot += hNGen->GetBinContent(i);
    }
    Double_t AxETot = NRecTot / NGenTot;
    Double_t AxETot_err = CalculateErrorBayes(NRecTot, NGenTot);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxETot*100, AxETot_err*100);

    return;
}

void FillHistNRec(){
    // Check if the corresponding text file already exists
    TString CutConfiguration = ConvertCutsToString(cuts);
    TString file((GlobPath + CutConfiguration + ".txt").Data());

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNRec with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNRec->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
        if(tRec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCRec(tRec);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nBins; iPtBin++){
            Int_t NRec = 0;
            for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
                tRec->GetEntry(iEntry);
                if(EventPassedMCRec_AxEPtDep(0, iPtBin)) NRec++;
            }
            hNRec->SetBinContent(iPtBin, NRec);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNRec, file);

        return;
    }
}

void FillHistNGen(){
    // Check if the corresponding text file already exists
    TString file((GlobPath + "NGen.txt").Data());

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNGen with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNGen->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
        if(tGen) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nBins; iPtBin++){
            Int_t NGen = 0;
            for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
                tGen->GetEntry(iEntry);
                if(EventPassedMCGen_AxEPtDep(iPtBin)) NGen++;
            }
            hNGen->SetBinContent(iPtBin, NGen);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNGen, file);

        return;
    }
}

Bool_t EventPassedMCRec_AxEPtDep(Int_t iMassCut = 0, Int_t iPtBin = 0){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    if(cuts[1]){ // !0VBA (no signal in the V0A)
        if(fTriggerInputsMC[0]) return kFALSE;
    }
    if(cuts[2]){ // !0VBC (no signal in the V0C)
        if(fTriggerInputsMC[1]) return kFALSE;
    }
    if(cuts[3]){ // !0UBA (no signal in the ADA)
        if(fTriggerInputsMC[2]) return kFALSE;
    }
    if(cuts[4]){ // !0UBC (no signal in the ADC)
        if(fTriggerInputsMC[3]) return kFALSE;
    }
    if(cuts[5]){ // 0STG (SPD topological)
        if(!fTriggerInputsMC[10]) return kFALSE;
    }
    if(cuts[6]){ // 0OMU (TOF two hits topology)
        if(!fTriggerInputsMC[4]) return kFALSE;
    }

    // 5) SPD cluster matches FOhits
    if(cuts[7]){
        if(!(fMatchingSPD == kTRUE)) return kFALSE;
    }

    // 6) AD offline veto (negligible effect on MC)
    if(cuts[8]){
        if(!(fADA_dec == 0 && fADC_dec == 0)) return kFALSE;
    }

    // 7) V0 offline veto (negligible effect on MC)
    if(cuts[9]){
        if(!(fV0A_dec == 0 && fV0C_dec == 0)) return kFALSE;
    }

    // 8) Dilepton rapidity |y| < 0.8
    if(cuts[10]){
        if(!(abs(fY) < 0.8)) return kFALSE;
    }

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(cuts[11]){
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    }

    // 10) Tracks have opposite charges
    if(cuts[12]){
        if(!(fQ1 * fQ2 < 0)) return kFALSE;
    }
    
    // 11) Muon pairs only
    if(cuts[13]){
        if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    }

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(cuts[14]){
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
    }

    // 13) Transverse momentum cut
    if(cuts[15]){
        if(cuts[0] == kFALSE){ // vs PtRec
            if(!(fPt > edges[iPtBin-1] && fPt <= edges[iPtBin])) return kFALSE;
        } else if(cuts[0] == kTRUE){ // vs PtGen
            if(!(fPtGen > edges[iPtBin-1] && fPtGen < edges[iPtBin])) return kFALSE;
        }
    }

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCGen_AxEPtDep(Int_t iPtBin = 0){
    // 1) Dilepton rapidity |y| < 0.8
    if(!(abs(fYGen) < 0.8)) return kFALSE;

    // 2) Transverse momentum cut
    if(!(fPtGen > edges[iPtBin-1] && fPtGen < edges[iPtBin])) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

void SaveToFile(TH1D* hist, TString name){
    ofstream outfile (name.Data());
    for(Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
        outfile << iBin << "\t" << hist->GetBinContent(iBin) << "\n";
    }
    outfile.close();
    Printf("*** Saved to %s.***", name.Data());
}

TString ConvertCutsToString(Bool_t *cuts){
    TString s("cuts_");
    for(Int_t iCut = 0; iCut < nCuts; iCut++){
        if(cuts[iCut] == kTRUE) s.Append("1");
        else s.Append("0");
    }
    //cout << s << endl;
    return s;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}