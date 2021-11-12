// ZNClasses.c
// David Grund, Nov 3, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

// Main functions
void PrepareData();
void InvMassFitsInClasses();
void CalculateClasses(Int_t iMassCut);
void PlotEnergyDistribution();
// Supporting functions
Bool_t EventPassedZN(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin);
Double_t CalculateErrorBayes(Double_t k, Double_t n);
void DoInvMassFitMain(Int_t ZNClass, Double_t fPtCutLow, Double_t fPtCutUpp, Int_t iPtBin);
void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit);
void SetCanvas(TCanvas *c, Bool_t bLogScale);

// All events (J/psi and background mixed)
// No. of events
Double_t nEv_all_total[nPtBins] = { 0 };    // total number of events per bin
Double_t nEv_all_fired[nPtBins] = { 0 };    // number of events per bin with at least one of the ZN modules fired
Double_t nEv_all_ZNA[nPtBins] = { 0 };      // number of events per bin with fired ZNA
Double_t nEv_all_ZNC[nPtBins] = { 0 };      // number of events per bin with fired ZNC
Double_t nEv_all_Xn0n[nPtBins] = { 0 };     // number of events per bin with fired ZNA only
Double_t nEv_all_0nXn[nPtBins] = { 0 };     // number of events per bin with fired ZNC only
Double_t nEv_all_XnXn[nPtBins] = { 0 };     // number of events per bin with fired both ZN modules
Double_t nEv_all_0n0n[nPtBins] = { 0 };     // number of events per bin with no fired ZN modules
Double_t nEv_all_sums[8] = { 0 };           // sum of the channels above over bins
// Percentages (relative to the total number of events in a bin)
Double_t Perc_all_XnXn_val[nPtBins] = { 0 }; 
Double_t Perc_all_XnXn_err[nPtBins] = { 0 }; 
Double_t Perc_all_0n0n_val[nPtBins] = { 0 }; 
Double_t Perc_all_0n0n_err[nPtBins] = { 0 }; 
Double_t Perc_all_Xn0n_val[nPtBins] = { 0 };
Double_t Perc_all_Xn0n_err[nPtBins] = { 0 }; 
Double_t Perc_all_0nXn_val[nPtBins] = { 0 }; 
Double_t Perc_all_0nXn_err[nPtBins] = { 0 }; 
// Jpsi events
// No. of events
Double_t nEv_J_total_val[nPtBins] = { 0 }; 
Double_t nEv_J_total_err[nPtBins] = { 0 }; 
Double_t nEv_J_Xn0n_val[nPtBins] = { 0 }; 
Double_t nEv_J_Xn0n_err[nPtBins] = { 0 };     
Double_t nEv_J_0nXn_val[nPtBins] = { 0 }; 
Double_t nEv_J_0nXn_err[nPtBins] = { 0 }; 
Double_t nEv_J_XnXn_val[nPtBins] = { 0 }; 
Double_t nEv_J_XnXn_err[nPtBins] = { 0 };     
Double_t nEv_J_0n0n_val[nPtBins] = { 0 };
Double_t nEv_J_0n0n_err[nPtBins] = { 0 };          
// Percentages: J/psi
Double_t Perc_J_XnXn_val[nPtBins] = { 0 }; 
Double_t Perc_J_XnXn_err[nPtBins] = { 0 }; 
Double_t Perc_J_0n0n_val[nPtBins] = { 0 }; 
Double_t Perc_J_0n0n_err[nPtBins] = { 0 }; 
Double_t Perc_J_Xn0n_val[nPtBins] = { 0 };
Double_t Perc_J_Xn0n_err[nPtBins] = { 0 }; 
Double_t Perc_J_0nXn_val[nPtBins] = { 0 }; 
Double_t Perc_J_0nXn_err[nPtBins] = { 0 }; 
// Background events
// No. of events
Double_t nEv_bkg_total_val[nPtBins] = { 0 };
Double_t nEv_bkg_total_err[nPtBins] = { 0 }; 
Double_t nEv_bkg_Xn0n_val[nPtBins] = { 0 }; 
Double_t nEv_bkg_Xn0n_err[nPtBins] = { 0 };     
Double_t nEv_bkg_0nXn_val[nPtBins] = { 0 };
Double_t nEv_bkg_0nXn_err[nPtBins] = { 0 };     
Double_t nEv_bkg_XnXn_val[nPtBins] = { 0 };   
Double_t nEv_bkg_XnXn_err[nPtBins] = { 0 };     
Double_t nEv_bkg_0n0n_val[nPtBins] = { 0 };  
Double_t nEv_bkg_0n0n_err[nPtBins] = { 0 };
// Percentages: background
Double_t Perc_bkg_XnXn_val[nPtBins] = { 0 }; 
Double_t Perc_bkg_XnXn_err[nPtBins] = { 0 }; 
Double_t Perc_bkg_0n0n_val[nPtBins] = { 0 }; 
Double_t Perc_bkg_0n0n_err[nPtBins] = { 0 }; 
Double_t Perc_bkg_Xn0n_val[nPtBins] = { 0 };
Double_t Perc_bkg_Xn0n_err[nPtBins] = { 0 }; 
Double_t Perc_bkg_0nXn_val[nPtBins] = { 0 }; 
Double_t Perc_bkg_0nXn_err[nPtBins] = { 0 }; 

TString classes[4] = {"0n0n", "Xn0n", "0nXn", "XnXn"};

// #############################################
// Options and variables for the inv. mass fits:
Int_t UpToBin = nPtBins;
Int_t UpToClass = 4;
Bool_t SigmaFixed = kTRUE;
Double_t N_Jpsi_out[2];
Double_t N_bkgr_out[2];
// #############################################

void ZNClasses(){

    //PrepareData();

    //InvMassFitsInClasses();

    CalculateClasses(0); // 3.0 < m < 3.2 GeV
    CalculateClasses(1); // 2.5 < m < 3.0 GeV && 3.2 < m < 4.0 GeV

    PlotEnergyDistribution(); // 3.0 < m < 3.2 GeV

    return;
}

void CalculateClasses(Int_t iMassCut){

    TFile *fFileIn = TFile::Open("Trees/ZNClasses/ZNClasses.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tZNClasses"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("fPt", &fPt);
    fTreeIn->SetBranchAddress("fM", &fM);
    fTreeIn->SetBranchAddress("fZNA_energy", &fZNA_energy);
    fTreeIn->SetBranchAddress("fZNC_energy", &fZNC_energy);
    fTreeIn->SetBranchAddress("fZNA_time", &fZNA_time);
    fTreeIn->SetBranchAddress("fZNC_time", &fZNC_time); 

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;     

    SetPtBinning();

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            if(EventPassedZN(iMassCut, 1, iBin+1)){
                Bool_t ZN_fired = kFALSE;
                Bool_t ZNA_fired = kFALSE;
                Bool_t ZNC_fired = kFALSE;
                for(Int_t i = 0; i < 4; i++){
                    // Check ZN times: at least one time has to fall within pm 2 ns
                    if(TMath::Abs(fZNA_time[i]) < 2 || TMath::Abs(fZNC_time[i]) < 2) ZN_fired = kTRUE;
                    if(TMath::Abs(fZNA_time[i]) < 2) ZNA_fired = kTRUE;
                    if(TMath::Abs(fZNC_time[i]) < 2) ZNC_fired = kTRUE;
                }
                nEv_all_total[iBin]++;
                if(ZN_fired) nEv_all_fired[iBin]++;
                if(ZNA_fired && ZNC_fired) nEv_all_XnXn[iBin]++;
                if(ZNA_fired) nEv_all_ZNA[iBin]++;
                if(ZNC_fired) nEv_all_ZNC[iBin]++;
                if(ZNA_fired && !ZNC_fired) nEv_all_Xn0n[iBin]++;
                if(ZNC_fired && !ZNA_fired) nEv_all_0nXn[iBin]++;
                if(!ZNA_fired && !ZNC_fired) nEv_all_0n0n[iBin]++;
            }
        }
        if((iEntry+1) % 500 == 0){
            nEntriesAnalysed += 500;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Calculate percentages (all events)
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Perc_all_0n0n_val[iBin] = nEv_all_0n0n[iBin] / nEv_all_total[iBin] * 100;
        Perc_all_0n0n_err[iBin] = CalculateErrorBayes(nEv_all_0n0n[iBin], nEv_all_total[iBin]) * 100;
        Perc_all_Xn0n_val[iBin] = nEv_all_Xn0n[iBin] / nEv_all_total[iBin] * 100;
        Perc_all_Xn0n_err[iBin] = CalculateErrorBayes(nEv_all_Xn0n[iBin], nEv_all_total[iBin]) * 100;
        Perc_all_0nXn_val[iBin] = nEv_all_0nXn[iBin] / nEv_all_total[iBin] * 100;
        Perc_all_0nXn_err[iBin] = CalculateErrorBayes(nEv_all_0nXn[iBin], nEv_all_total[iBin]) * 100;
        Perc_all_XnXn_val[iBin] = nEv_all_XnXn[iBin] / nEv_all_total[iBin] * 100;
        Perc_all_XnXn_err[iBin] = CalculateErrorBayes(nEv_all_XnXn[iBin], nEv_all_total[iBin]) * 100;
    }

    // Load the values from the invariant mass fits
    TString str_in_J = Form("Results/ZNClasses/%ibins_J.txt", nPtBins);
    TString str_in_bkg = Form("Results/ZNClasses/%ibins_bkg.txt", nPtBins);
    ifstream ifs;
    ifs.open(str_in_J.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        Int_t iBin;
        ifs >> iBin 
            >> nEv_J_total_val[i]
            >> nEv_J_total_err[i]
            >> nEv_J_0n0n_val[i]
            >> nEv_J_0n0n_err[i]
            >> nEv_J_Xn0n_val[i]
            >> nEv_J_Xn0n_err[i]
            >> nEv_J_0nXn_val[i]
            >> nEv_J_0nXn_err[i] 
            >> nEv_J_XnXn_val[i]
            >> nEv_J_XnXn_err[i];
    }
    ifs.close();
    Printf("Values from %s loaded.", str_in_J.Data());
    ifs.open(str_in_bkg.Data());
    for(Int_t i = 0; i < nPtBins; i++){
        Int_t iBin;
        ifs >> iBin 
            >> nEv_bkg_total_val[i]
            >> nEv_bkg_total_err[i]
            >> nEv_bkg_0n0n_val[i]
            >> nEv_bkg_0n0n_err[i]
            >> nEv_bkg_Xn0n_val[i]
            >> nEv_bkg_Xn0n_err[i]
            >> nEv_bkg_0nXn_val[i]
            >> nEv_bkg_0nXn_err[i] 
            >> nEv_bkg_XnXn_val[i]
            >> nEv_bkg_XnXn_err[i];
    }
    ifs.close();
    Printf("Values from %s loaded.", str_in_bkg.Data());

    // Calculate percentages (J/psi and background)
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Perc_J_0n0n_val[iBin] = nEv_J_0n0n_val[iBin] / nEv_J_total_val[iBin] * 100;
        Perc_J_0n0n_err[iBin] = Perc_J_0n0n_val[iBin] * TMath::Sqrt(TMath::Power(nEv_J_0n0n_err[iBin]/nEv_J_0n0n_val[iBin], 2) + TMath::Power(nEv_J_total_err[iBin]/nEv_J_total_val[iBin], 2));
        Perc_J_Xn0n_val[iBin] = nEv_J_Xn0n_val[iBin] / nEv_J_total_val[iBin] * 100;
        Perc_J_Xn0n_err[iBin] = Perc_J_Xn0n_val[iBin] * TMath::Sqrt(TMath::Power(nEv_J_Xn0n_err[iBin]/nEv_J_Xn0n_val[iBin], 2) + TMath::Power(nEv_J_total_err[iBin]/nEv_J_total_val[iBin], 2));
        Perc_J_0nXn_val[iBin] = nEv_J_0nXn_val[iBin] / nEv_J_total_val[iBin] * 100;
        Perc_J_0nXn_err[iBin] = Perc_J_0nXn_val[iBin] * TMath::Sqrt(TMath::Power(nEv_J_0nXn_err[iBin]/nEv_J_0nXn_val[iBin], 2) + TMath::Power(nEv_J_total_err[iBin]/nEv_J_total_val[iBin], 2));
        Perc_J_XnXn_val[iBin] = nEv_J_XnXn_val[iBin] / nEv_J_total_val[iBin] * 100;
        Perc_J_XnXn_err[iBin] = Perc_J_XnXn_val[iBin] * TMath::Sqrt(TMath::Power(nEv_J_XnXn_err[iBin]/nEv_J_XnXn_val[iBin], 2) + TMath::Power(nEv_J_total_err[iBin]/nEv_J_total_val[iBin], 2));
    }    
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Perc_bkg_0n0n_val[iBin] = nEv_bkg_0n0n_val[iBin] / nEv_bkg_total_val[iBin] * 100;
        Perc_bkg_0n0n_err[iBin] = Perc_bkg_0n0n_val[iBin] * TMath::Sqrt(TMath::Power(nEv_bkg_0n0n_err[iBin]/nEv_bkg_0n0n_val[iBin], 2) + TMath::Power(nEv_bkg_total_err[iBin]/nEv_bkg_total_val[iBin], 2));
        Perc_bkg_Xn0n_val[iBin] = nEv_bkg_Xn0n_val[iBin] / nEv_bkg_total_val[iBin] * 100;
        Perc_bkg_Xn0n_err[iBin] = Perc_bkg_Xn0n_val[iBin] * TMath::Sqrt(TMath::Power(nEv_bkg_Xn0n_err[iBin]/nEv_bkg_Xn0n_val[iBin], 2) + TMath::Power(nEv_bkg_total_err[iBin]/nEv_bkg_total_val[iBin], 2));
        Perc_bkg_0nXn_val[iBin] = nEv_bkg_0nXn_val[iBin] / nEv_bkg_total_val[iBin] * 100;
        Perc_bkg_0nXn_err[iBin] = Perc_bkg_0nXn_val[iBin] * TMath::Sqrt(TMath::Power(nEv_bkg_0nXn_err[iBin]/nEv_bkg_0nXn_val[iBin], 2) + TMath::Power(nEv_bkg_total_err[iBin]/nEv_bkg_total_val[iBin], 2));
        Perc_bkg_XnXn_val[iBin] = nEv_bkg_XnXn_val[iBin] / nEv_bkg_total_val[iBin] * 100;
        Perc_bkg_XnXn_err[iBin] = Perc_bkg_XnXn_val[iBin] * TMath::Sqrt(TMath::Power(nEv_bkg_XnXn_err[iBin]/nEv_bkg_XnXn_val[iBin], 2) + TMath::Power(nEv_bkg_total_err[iBin]/nEv_bkg_total_val[iBin], 2));
    } 

    // Print the results to the console
    // Numbers
    Printf("*******************************************************************");
    Printf("Numbers of events:");
    Printf("Bin \tTotal \tFired \tZNA \tZNC \tOnlyZNA\tOnlyZNC\tBoth \tNone\n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        nEv_all_sums[0] += nEv_all_total[iBin];
        nEv_all_sums[1] += nEv_all_fired[iBin];
        nEv_all_sums[2] += nEv_all_ZNA[iBin];
        nEv_all_sums[3] += nEv_all_ZNC[iBin];
        nEv_all_sums[4] += nEv_all_Xn0n[iBin];
        nEv_all_sums[5] += nEv_all_0nXn[iBin];
        nEv_all_sums[6] += nEv_all_XnXn[iBin];
        nEv_all_sums[7] += nEv_all_0n0n[iBin];
        Printf("%i \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f", 
            iBin+1, 
            nEv_all_total[iBin], 
            nEv_all_fired[iBin], 
            nEv_all_ZNA[iBin], 
            nEv_all_ZNC[iBin], 
            nEv_all_Xn0n[iBin], 
            nEv_all_0nXn[iBin], 
            nEv_all_XnXn[iBin],
            nEv_all_0n0n[iBin]);
    }
    Printf("sum \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f", nEv_all_sums[0], nEv_all_sums[1], nEv_all_sums[2], nEv_all_sums[3], nEv_all_sums[4], nEv_all_sums[5], nEv_all_sums[6], nEv_all_sums[7]);
    // Percentages
    Printf("*******************************************************************");
    Printf("Percentages:");
    Printf("Bin \t0n0n[%%]\terr \tXn0n[%%]\terr \t0nXn[%%]\terr \t0n0n[%%]\terr");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Printf("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f", 
            iBin+1,
            Perc_all_0n0n_val[iBin],
            Perc_all_0n0n_err[iBin],
            Perc_all_Xn0n_val[iBin],
            Perc_all_Xn0n_err[iBin],
            Perc_all_0nXn_val[iBin],
            Perc_all_0nXn_err[iBin],
            Perc_all_XnXn_val[iBin],
            Perc_all_XnXn_err[iBin]);
    }
    Printf("*******************************************************************");

    // Print numbers to txt files
    TString str_out_mass = Form("Results/ZNClasses/%ibins_mass%i.txt", nPtBins, iMassCut);
    ofstream fout_mass(str_out_mass.Data());
    fout_mass << Form("No. of events:\n");
    fout_mass << Form("Bin \ttotal \tfired \tfZNA \tfZNC \t0n0n\tXn0n\t0nXn\tXnXn\n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_mass << Form("%i \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \n", 
            iBin+1,
            nEv_all_total[iBin], 
            nEv_all_fired[iBin], 
            nEv_all_ZNA[iBin], 
            nEv_all_ZNC[iBin],
            nEv_all_0n0n[iBin],
            nEv_all_Xn0n[iBin], 
            nEv_all_0nXn[iBin], 
            nEv_all_XnXn[iBin]);
    }
    fout_mass << Form("sum \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \n", nEv_all_sums[0], nEv_all_sums[1], nEv_all_sums[2], nEv_all_sums[3], nEv_all_sums[4], nEv_all_sums[5], nEv_all_sums[6], nEv_all_sums[7]);
    fout_mass << Form("\n\n");
    fout_mass << Form("Percentages (all events):\n");
    fout_mass << Form("Bin \t0n0n[%%]\terr \tXn0n[%%]\terr \t0nXn[%%]\terr \tXnXn[%%]\terr \n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_mass << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            Perc_all_0n0n_val[iBin],
            Perc_all_0n0n_err[iBin],
            Perc_all_Xn0n_val[iBin],
            Perc_all_Xn0n_err[iBin],
            Perc_all_0nXn_val[iBin],
            Perc_all_0nXn_err[iBin],
            Perc_all_XnXn_val[iBin],
            Perc_all_XnXn_err[iBin]);
    }
    fout_mass.close();
    Printf("Results printed to %s.", str_out_mass.Data()); 

    TString str_out_J_bkg = Form("Results/ZNClasses/%ibins_perc_J_bkg.txt", nPtBins);
    ofstream fout_J_bkg(str_out_J_bkg.Data());
    fout_J_bkg << Form("Percentages (J/psi events):\n");
    fout_J_bkg << Form("Bin \t0n0n[%%]\terr \tXn0n[%%]\terr \t0nXn[%%]\terr \tXnXn[%%]\terr \n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_J_bkg << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            Perc_J_0n0n_val[iBin],
            Perc_J_0n0n_err[iBin],
            Perc_J_Xn0n_val[iBin],
            Perc_J_Xn0n_err[iBin],
            Perc_J_0nXn_val[iBin],
            Perc_J_0nXn_err[iBin],
            Perc_J_XnXn_val[iBin],
            Perc_J_XnXn_err[iBin]);
    }
    fout_J_bkg << Form("\n\n");
    fout_J_bkg << Form("Percentages (bkgr events):\n");
    fout_J_bkg << Form("Bin \t0n0n[%%]\terr \tXn0n[%%]\terr \t0nXn[%%]\terr \tXnXn[%%]\terr \n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_J_bkg << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            Perc_bkg_0n0n_val[iBin],
            Perc_bkg_0n0n_err[iBin],
            Perc_bkg_Xn0n_val[iBin],
            Perc_bkg_Xn0n_err[iBin],
            Perc_bkg_0nXn_val[iBin],
            Perc_bkg_0nXn_err[iBin],
            Perc_bkg_XnXn_val[iBin],
            Perc_bkg_XnXn_err[iBin]);
    }
    fout_J_bkg.close();
    Printf("Results printed to %s.", str_out_J_bkg.Data()); 

    // Print to TeX table (all events)
    TString str_out_mass_TeX = Form("Results/ZNClasses/%ibins_TeX_mass%i.txt", nPtBins, iMassCut);
    ofstream fout_mass_TeX(str_out_mass_TeX.Data());
    fout_mass_TeX << Form("Bin \tTotal \t0n0n[-]\tXn0n[-]\t0nXn[-]\tXnXn[-]\t0n0n[%%]\tXn0n[%%]\t0nXn[%%]\tXnXn[%%]\n");    
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_mass_TeX << std::fixed << std::setprecision(3)
                    << "$(" << ptBoundaries[iBin] << "," << ptBoundaries[iBin+1] << ")$ &\t"
                    << std::fixed << std::setprecision(0)
                    << nEv_all_total[iBin] << " &\t"
                    << nEv_all_0n0n[iBin] << " &\t"
                    << nEv_all_Xn0n[iBin] << " &\t"
                    << nEv_all_0nXn[iBin] << " &\t"
                    << nEv_all_XnXn[iBin] << " &\t$"
                    << std::fixed << std::setprecision(1)
                    << Perc_all_0n0n_val[iBin] << R"( \pm )" << Perc_all_0n0n_err[iBin] << "$ &\t$"
                    << Perc_all_Xn0n_val[iBin] << R"( \pm )" << Perc_all_Xn0n_err[iBin] << "$ &\t$"
                    << Perc_all_0nXn_val[iBin] << R"( \pm )" << Perc_all_0nXn_err[iBin] << "$ &\t$"
                    << Perc_all_XnXn_val[iBin] << R"( \pm )" << Perc_all_XnXn_err[iBin] << R"($ \\)" << "\n";
    }
    fout_mass_TeX.close();
    Printf("Results printed to %s.", str_out_mass_TeX.Data()); 

    // Print to TeX table (J/psi)
    TString str_out_J_TeX = Form("Results/ZNClasses/%ibins_TeX_J.txt", nPtBins);
    ofstream fout_J_TeX(str_out_J_TeX.Data()); 
    fout_J_TeX << Form("Bin\tTotal\t0n0n[-]\tXn0n[-]\t0nXn[-]\tXnXn[-]\n");    
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_J_TeX << std::fixed << std::setprecision(3)
                    << "$(" << ptBoundaries[iBin] << "," << ptBoundaries[iBin+1] << ")$ &\t$"
                    << std::fixed << std::setprecision(1)
                    << nEv_J_total_val[iBin] << R"( \pm )" << nEv_J_total_err[iBin] << "$ &\t$"
                    << nEv_J_0n0n_val[iBin] << R"( \pm )" << nEv_J_0n0n_err[iBin] << "$ &\t$"
                    << nEv_J_Xn0n_val[iBin] << R"( \pm )" << nEv_J_Xn0n_err[iBin] << "$ &\t$"
                    << nEv_J_0nXn_val[iBin] << R"( \pm )" << nEv_J_0nXn_err[iBin] << "$ &\t$"
                    << nEv_J_XnXn_val[iBin] << R"( \pm )" << nEv_J_XnXn_err[iBin] << R"($ \\)" << "\n";
    } 
    fout_J_TeX << Form("\n\nBin\t0n0n[%%]\tXn0n[%%]\t0nXn[%%]\tXnXn[%%]\n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_J_TeX << std::fixed << std::setprecision(3)
                    << "$(" << ptBoundaries[iBin] << "," << ptBoundaries[iBin+1] << ")$ &\t$"
                    << std::fixed << std::setprecision(1)
                    << Perc_J_0n0n_val[iBin] << R"( \pm )" << Perc_J_0n0n_err[iBin] << "$ &\t$"
                    << Perc_J_Xn0n_val[iBin] << R"( \pm )" << Perc_J_Xn0n_err[iBin] << "$ &\t$"
                    << Perc_J_0nXn_val[iBin] << R"( \pm )" << Perc_J_0nXn_err[iBin] << "$ &\t$"
                    << Perc_J_XnXn_val[iBin] << R"( \pm )" << Perc_J_XnXn_err[iBin] << R"($ \\)" << "\n";
    }
    fout_J_TeX.close();
    Printf("Results printed to %s.", str_out_J_TeX.Data()); 

    // Print to TeX table (background)
    TString str_out_bkg_TeX = Form("Results/ZNClasses/%ibins_TeX_bkg.txt", nPtBins);
    ofstream fout_bkg_TeX(str_out_bkg_TeX.Data()); 
    fout_bkg_TeX << Form("Bin\tTotal\t0n0n[-]\tXn0n[-]\t0nXn[-]\tXnXn[-]\n");    
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_bkg_TeX << std::fixed << std::setprecision(3)
                    << "$(" << ptBoundaries[iBin] << "," << ptBoundaries[iBin+1] << ")$ &\t$"
                    << std::fixed << std::setprecision(1)
                    << nEv_bkg_total_val[iBin] << R"( \pm )" << nEv_bkg_total_err[iBin] << "$ &\t$"
                    << nEv_bkg_0n0n_val[iBin] << R"( \pm )" << nEv_bkg_0n0n_err[iBin] << "$ &\t$"
                    << nEv_bkg_Xn0n_val[iBin] << R"( \pm )" << nEv_bkg_Xn0n_err[iBin] << "$ &\t$"
                    << nEv_bkg_0nXn_val[iBin] << R"( \pm )" << nEv_bkg_0nXn_err[iBin] << "$ &\t$"
                    << nEv_bkg_XnXn_val[iBin] << R"( \pm )" << nEv_bkg_XnXn_err[iBin] << R"($ \\)" << "\n";
    } 
    fout_bkg_TeX << Form("\n\nBin\t0n0n[%%]\tXn0n[%%]\t0nXn[%%]\tXnXn[%%]\n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        fout_bkg_TeX << std::fixed << std::setprecision(3)
                    << "$(" << ptBoundaries[iBin] << "," << ptBoundaries[iBin+1] << ")$ &\t$"
                    << std::fixed << std::setprecision(1)
                    << Perc_bkg_0n0n_val[iBin] << R"( \pm )" << Perc_bkg_0n0n_err[iBin] << "$ &\t$"
                    << Perc_bkg_Xn0n_val[iBin] << R"( \pm )" << Perc_bkg_Xn0n_err[iBin] << "$ &\t$"
                    << Perc_bkg_0nXn_val[iBin] << R"( \pm )" << Perc_bkg_0nXn_err[iBin] << "$ &\t$"
                    << Perc_bkg_XnXn_val[iBin] << R"( \pm )" << Perc_bkg_XnXn_err[iBin] << R"($ \\)" << "\n";
    }
    fout_bkg_TeX.close();
    Printf("Results printed to %s.", str_out_bkg_TeX.Data()); 

    return;
}

void InvMassFitsInClasses(){

    SetPtBinning();

    for(Int_t iBin = 0; iBin < UpToBin; iBin++){
        for(Int_t iZNClass = 0; iZNClass < UpToClass; iZNClass++){
            DoInvMassFitMain(iZNClass, ptBoundaries[iBin], ptBoundaries[iBin+1], iBin+1);
            if(iZNClass == 0){
                nEv_J_0n0n_val[iBin] = N_Jpsi_out[0];
                nEv_J_0n0n_err[iBin] = N_Jpsi_out[1];
                nEv_bkg_0n0n_val[iBin] = N_bkgr_out[0];
                nEv_bkg_0n0n_err[iBin] = N_bkgr_out[1];
            } else if(iZNClass == 1){
                nEv_J_Xn0n_val[iBin] = N_Jpsi_out[0];
                nEv_J_Xn0n_err[iBin] = N_Jpsi_out[1];
                nEv_bkg_Xn0n_val[iBin] = N_bkgr_out[0];
                nEv_bkg_Xn0n_err[iBin] = N_bkgr_out[1];
            } else if(iZNClass == 2){
                nEv_J_0nXn_val[iBin] = N_Jpsi_out[0];
                nEv_J_0nXn_err[iBin] = N_Jpsi_out[1];
                nEv_bkg_0nXn_val[iBin] = N_bkgr_out[0];
                nEv_bkg_0nXn_err[iBin] = N_bkgr_out[1];
            } else if(iZNClass == 3){
                nEv_J_XnXn_val[iBin] = N_Jpsi_out[0];
                nEv_J_XnXn_err[iBin] = N_Jpsi_out[1];
                nEv_bkg_XnXn_val[iBin] = N_bkgr_out[0];
                nEv_bkg_XnXn_err[iBin] = N_bkgr_out[1];
            } 
        }
        nEv_J_total_val[iBin] = nEv_J_0n0n_val[iBin] + nEv_J_Xn0n_val[iBin] + nEv_J_0nXn_val[iBin] + nEv_J_XnXn_val[iBin];
        nEv_J_total_err[iBin] = TMath::Sqrt(TMath::Power(nEv_J_0n0n_err[iBin],2) 
                                            + TMath::Power(nEv_J_Xn0n_err[iBin],2)
                                            + TMath::Power(nEv_J_0nXn_err[iBin],2)
                                            + TMath::Power(nEv_J_XnXn_err[iBin],2));
        nEv_bkg_total_val[iBin] = nEv_bkg_0n0n_val[iBin] + nEv_bkg_Xn0n_val[iBin] + nEv_bkg_0nXn_val[iBin] + nEv_bkg_XnXn_val[iBin];
        nEv_bkg_total_err[iBin] = TMath::Sqrt(TMath::Power(nEv_bkg_0n0n_err[iBin],2) 
                                            + TMath::Power(nEv_bkg_Xn0n_err[iBin],2)
                                            + TMath::Power(nEv_bkg_0nXn_err[iBin],2)
                                            + TMath::Power(nEv_bkg_XnXn_err[iBin],2));
    }

    // Print the results to text file
    // J/psi
    TString str_out_J = Form("Results/ZNClasses/%ibins_J.txt", nPtBins);
    ofstream f_out_J(str_out_J.Data());
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        f_out_J << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            nEv_J_total_val[iBin], 
            nEv_J_total_err[iBin],
            nEv_J_0n0n_val[iBin],
            nEv_J_0n0n_err[iBin],
            nEv_J_Xn0n_val[iBin], 
            nEv_J_Xn0n_err[iBin],
            nEv_J_0nXn_val[iBin],
            nEv_J_0nXn_err[iBin], 
            nEv_J_XnXn_val[iBin],
            nEv_J_XnXn_err[iBin]);
    }
    f_out_J.close();
    Printf("Results printed to %s.", str_out_J.Data()); 
    // background
    TString str_out_bkg = Form("Results/ZNClasses/%ibins_bkg.txt", nPtBins);
    ofstream f_out_bkg(str_out_bkg.Data());
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        f_out_bkg << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            nEv_bkg_total_val[iBin], 
            nEv_bkg_total_err[iBin],
            nEv_bkg_0n0n_val[iBin],
            nEv_bkg_0n0n_err[iBin],
            nEv_bkg_Xn0n_val[iBin], 
            nEv_bkg_Xn0n_err[iBin],
            nEv_bkg_0nXn_val[iBin],
            nEv_bkg_0nXn_err[iBin], 
            nEv_bkg_XnXn_val[iBin],
            nEv_bkg_XnXn_err[iBin]);
    }
    f_out_bkg.close();
    Printf("Results printed to %s.", str_out_bkg.Data()); 

    return;
}

void PlotEnergyDistribution(){

    TFile *fFileIn = TFile::Open("Trees/ZNClasses/ZNClasses.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tZNClasses"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("fPt", &fPt);
    fTreeIn->SetBranchAddress("fM", &fM);
    fTreeIn->SetBranchAddress("fZNA_energy", &fZNA_energy);
    fTreeIn->SetBranchAddress("fZNC_energy", &fZNC_energy);
    fTreeIn->SetBranchAddress("fZNA_time", &fZNA_time);
    fTreeIn->SetBranchAddress("fZNC_time", &fZNC_time);

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;     

    SetPtBinning();

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird); // https://root.cern/doc/master/classTColor.html#C05
    gStyle->SetPaintTextFormat("4.2f");

    // Variables for histograms
    Int_t nBins1D = 64;
    Int_t nBins2D = 48;
    Double_t EnLow = -4;
    Double_t EnUpp2D = 8;
    Double_t EnUpp1D = 12;

    for(Int_t iBin = 0; iBin < nPtBins; iBin++){

        TH2D *hZN_energies = new TH2D("hZN_energies", "correlation between ZNA and ZNC energies", nBins2D, EnLow, EnUpp2D, nBins2D, EnLow, EnUpp2D);
        TH1D *hZNA_energy = new TH1D("hZNA_energy", "hist of ZNA energies", nBins1D, EnLow, EnUpp1D);
        TH1D *hZNC_energy = new TH1D("hZNC_energy", "hist of ZNC energies", nBins1D, EnLow, EnUpp1D);    

        for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
            fTreeIn->GetEntry(iEntry);

            Bool_t SignalZN = kFALSE;
            for(Int_t i = 0; i < 4; i++){
                if(TMath::Abs(fZNA_time[i]) < 2 || TMath::Abs(fZNC_time[i]) < 2) SignalZN = kTRUE;
            }
            
            if(EventPassedZN(0,1,iBin+1) && SignalZN){
                hZN_energies->Fill(fZNA_energy/1000, fZNC_energy/1000);
                hZNA_energy->Fill(fZNA_energy/1000);
                hZNC_energy->Fill(fZNC_energy/1000);
            }
        }   
        // Print 2D histogram
        TCanvas *c1 = new TCanvas("c1","c1",700,700);
        c1->SetTopMargin(0.08);
        c1->SetBottomMargin(0.11);
        c1->SetLeftMargin(0.09);
        c1->SetRightMargin(0.11);
        c1->SetTitle("ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"); 
        // X-axis
        hZN_energies->GetXaxis()->SetTitle("ZNA energy (TeV)");
        hZN_energies->GetXaxis()->SetTitleSize(0.05);
        hZN_energies->GetXaxis()->SetLabelSize(0.05);
        // Y-axis
        hZN_energies->GetYaxis()->SetTitle("ZNC energy (TeV)");
        hZN_energies->GetYaxis()->SetTitleSize(0.05);
        hZN_energies->GetYaxis()->SetLabelSize(0.05);
        hZN_energies->GetYaxis()->SetTitleOffset(0.7);
        // Z-axis
        hZN_energies->GetZaxis()->SetLabelSize(0.042);
        hZN_energies->GetZaxis()->SetDecimals(1);
        // Draw    
        hZN_energies->Draw("COLZ");
        // Legend 1
        TLegend *l1 = new TLegend(0.12,0.93,0.6,1.0);
        l1->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
        l1->SetTextSize(0.05);
        l1->SetBorderSize(0);
        l1->SetFillStyle(0);
        l1->Draw();
        // Legend 2
        TLegend *l2 = new TLegend(0.10,0.12,0.22,0.22);
        l2->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", ptBoundaries[iBin],ptBoundaries[iBin+1]),"");
        l2->SetTextSize(0.048);
        l2->SetBorderSize(0);
        l2->SetFillStyle(0);
        l2->Draw();
        // Save canvas
        c1->Print(Form("Results/ZNClasses/EnergyDistribution/%ibins/bin%i_EnergyZN_2D.pdf", nPtBins, iBin+1));
        c1->Print(Form("Results/ZNClasses/EnergyDistribution/%ibins/bin%i_EnergyZN_2D.png", nPtBins, iBin+1));

        // Print 1D histograms
        TCanvas *c2 = new TCanvas("c2","c2",700,600);
        c2->SetTopMargin(0.02);
        c2->SetBottomMargin(0.11);
        c2->SetLeftMargin(0.13);
        c2->SetRightMargin(0.04);
        c2->SetLogy();
        // X-axis
        hZNA_energy->GetXaxis()->SetTitle("ZN energy (TeV)");
        hZNA_energy->GetXaxis()->SetTitleSize(0.05);
        hZNA_energy->GetXaxis()->SetLabelSize(0.05);
        // Y-axis
        hZNA_energy->GetYaxis()->SetTitle("Counts per 250 GeV");
        hZNA_energy->GetYaxis()->SetTitleSize(0.05);
        hZNA_energy->GetYaxis()->SetLabelSize(0.05);
        hZNA_energy->GetYaxis()->SetTitleOffset(1.2);
        // Style hist ZNA
        hZNA_energy->SetLineColor(kRed);
        hZNA_energy->SetLineWidth(1);
        hZNA_energy->SetMarkerStyle(21);
        hZNA_energy->SetMarkerColor(kRed);
        hZNA_energy->SetMarkerSize(0.5);
        // Style hist ZNC
        hZNC_energy->SetLineColor(kBlue);
        hZNC_energy->SetLineWidth(1);
        hZNC_energy->SetMarkerStyle(21);
        hZNC_energy->SetMarkerColor(kBlue);
        hZNC_energy->SetMarkerSize(0.5);
        // Draw
        hZNA_energy->Draw("E0");
        hZNC_energy->Draw("SAME E0");
        // Legend 3
        TLegend *l3 = new TLegend(0.215,0.91,0.99,0.96);
        l3->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
        l3->SetTextSize(0.05);
        l3->SetBorderSize(0);
        l3->SetFillStyle(0);
        l3->Draw();
        // Legend 4
        TLegend *l4 = new TLegend(0.78,0.71,1.08,0.83);
        l4->AddEntry(hZNA_energy,"ZNA","PL");
        l4->AddEntry(hZNC_energy,"ZNC","PL");
        l4->SetTextSize(0.05);
        l4->SetBorderSize(0);
        l4->SetFillStyle(0);
        l4->Draw();
        // Legend 5
        TLegend *l5 = new TLegend(0.41,0.84,0.99,0.90);
        l5->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", ptBoundaries[iBin],ptBoundaries[iBin+1]),"");
        l5->SetTextSize(0.05);
        l5->SetBorderSize(0);
        l5->SetFillStyle(0);
        l5->Draw();
        // Save canvas
        c2->Print(Form("Results/ZNClasses/EnergyDistribution/%ibins/bin%i_EnergyZN_1D.pdf", nPtBins, iBin+1));
        c2->Print(Form("Results/ZNClasses/EnergyDistribution/%ibins/bin%i_EnergyZN_1D.png", nPtBins, iBin+1));
    }

    return;
}

void PrepareData(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");  

    ConnectTreeVariables(fTreeIn);  

    TFile fFileOut("Trees/ZNClasses/ZNClasses.root","RECREATE");

    TTree *tZNClasses = new TTree("tZNClasses", "tZNClasses");
    tZNClasses->Branch("fPt", &fPt, "fPt/D");
    tZNClasses->Branch("fM", &fM, "fM/D");
    tZNClasses->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    tZNClasses->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    tZNClasses->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    tZNClasses->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");

    TTree *t0n0n = new TTree("t0n0n", "t0n0n");
    t0n0n->Branch("fPt", &fPt, "fPt/D");
    t0n0n->Branch("fM", &fM, "fM/D");
    t0n0n->Branch("fY", &fY, "fY/D");

    TTree *tXn0n = new TTree("tXn0n", "tXn0n");
    tXn0n->Branch("fPt", &fPt, "fPt/D");
    tXn0n->Branch("fM", &fM, "fM/D");
    tXn0n->Branch("fY", &fY, "fY/D");

    TTree *t0nXn = new TTree("t0nXn", "t0nXn");
    t0nXn->Branch("fPt", &fPt, "fPt/D");
    t0nXn->Branch("fM", &fM, "fM/D");
    t0nXn->Branch("fY", &fY, "fY/D");

    TTree *tXnXn = new TTree("tXnXn", "tXnXn");
    tXnXn->Branch("fPt", &fPt, "fPt/D");
    tXnXn->Branch("fM", &fM, "fM/D");
    tXnXn->Branch("fY", &fY, "fY/D");

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        // iMassCut == 0    => 2.2 < m < 4.5 GeV
        // iPtCut == 3      => 0.2 < pT < 1.0 GeV
        if(EventPassed(0,3)) tZNClasses->Fill();
        // fill the trees corresponding to ZN classes
        if(EventPassed(0,3)){
            Bool_t ZNA_fired = kFALSE;
            Bool_t ZNC_fired = kFALSE;
            for(Int_t i = 0; i < 4; i++){
                // Check ZN times: at least one time has to fall within pm 2 ns
                if(TMath::Abs(fZNA_time[i]) < 2) ZNA_fired = kTRUE;
                if(TMath::Abs(fZNC_time[i]) < 2) ZNC_fired = kTRUE;
            }
            if(!ZNA_fired && !ZNC_fired) t0n0n->Fill();
            if(ZNA_fired && !ZNC_fired) tXn0n->Fill();
            if(ZNC_fired && !ZNA_fired) t0nXn->Fill();
            if(ZNA_fired && ZNC_fired) tXnXn->Fill();
        }

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}

Bool_t EventPassedZN(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin){

    // Inv. mass cut
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case 0:
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
        case 1:
            if((fM > 2.5 && fM < 2.9) || (fM > 3.3 && fM < 4.0)) bMassCut = kTRUE;
            break;       
    }
    if(!bMassCut) return kFALSE;

    // Transverse momentu cut
    Bool_t bPtCut = kFALSE;
    switch(iPtCut){
        case 0:
            if(fPt > 0.2 && fPt < 1.0) bPtCut = kTRUE;
            break;
        case 1:
            if(fPt > ptBoundaries[iPtBin-1] && fPt <= ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = nominator, n = denominator

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}

void DoInvMassFitMain(Int_t ZNClass, Double_t fPtCutLow, Double_t fPtCutUpp, Int_t iPtBin){
    // class 0 = 0n0n
    // class 1 = Xn0n (ZNA hit, ZNC no hit)
    // class 2 = 0nXn (ZNC hit, ZNA no hit)
    // class 3 = XnXn (hit in both ZNs)

    // Cuts:
    char fStrReduce[120];
    Double_t fYCut      = 0.8;
    Double_t fMCutLow   = 2.5;
    Double_t fMCutUpp   = 4.0;

    sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);

    // Binning:
    Int_t nBins = 50; // so that each bin between 2.5 and 4.0 GeV is 20 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);
  
    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    //fM.setBinning(binM);

    // Get the data trees
    TFile *fFileIn = new TFile("Trees/ZNClasses/ZNClasses.root"); 
    TTree *fTreeIn = NULL;
    if(ZNClass == 0){
        fFileIn->GetObject("t0n0n",fTreeIn);
    } else if(ZNClass == 1){
        fFileIn->GetObject("tXn0n",fTreeIn);
    } else if(ZNClass == 2){
        fFileIn->GetObject("t0nXn",fTreeIn);
    } else if(ZNClass == 3){
        fFileIn->GetObject("tXnXn",fTreeIn);
    }

    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*fTreeIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    TString pathMC = "Results/InvMassFitMC/allbins/allbins.txt";

    ifstream ifs;
    ifs.open(pathMC.Data());
    char name[8];
    Double_t err;
    for(Int_t i = 0; i < 4; i++){
        if(i == 0) ifs >> name >> fAlpha_L >> err;
        if(i == 1) ifs >> name >> fAlpha_R >> err;
        if(i == 2) ifs >> name >> fN_L >> err;
        if(i == 3) ifs >> name >> fN_R >> err;
    }
    ifs.close();
    Printf("MC values of tail parameters loaded.");

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,30.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,30.);
    alpha_L.setConstant(kTRUE);
    alpha_R.setConstant(kTRUE);
    n_L.setConstant(kTRUE);
    n_R.setConstant(kTRUE);

    // Crystal Ball for J/Psi
    RooRealVar mass_Jpsi("mass_Jpsi","J/psi mass",3.097,3.00,3.20); 
    //mass_Jpsi.setConstant(kTRUE);
    RooRealVar sigma_Jpsi("sigma_Jpsi","J/psi resolution",0.021,0.01,0.1);
    if(SigmaFixed) sigma_Jpsi.setConstant(kTRUE);
    RooGenericPdf mean_R("mean_R","J/psi mass","mass_Jpsi",RooArgSet(mass_Jpsi));
    RooGenericPdf sigma_R("sigma_R","J/psi resolution","sigma_Jpsi",RooArgSet(sigma_Jpsi));
    RooRealVar N_Jpsi("N_Jpsi","number of J/psi events",0.4*nEvents,0,nEvents);

    // Background
    RooRealVar lambda("lambda","background exp",-1.2,-10.,0.);
    RooRealVar N_bkg("N_bkg","number of background events",0.6*nEvents,0,nEvents);

    // Functions for fitting
    // J/psi:
    RooCBShape CB_left("CB_left","CB_left",fM,mass_Jpsi,sigma_Jpsi,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));
    // Background:
    RooGenericPdf BkgPdf("BkgPdf","exp(fM*lambda)",RooArgSet(fM,lambda));

    // Create model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of J/psi events
    fM.setRange("WholeMassRange",fMCutLow,fMCutUpp);
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *iDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("WholeMassRange"));
    RooAbsReal *iBkgr = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_out[0] = iDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = iDSCB->getVal()*N_Jpsi.getError();
    N_bkgr_out[0] = iBkgr->getVal()*N_bkg.getVal();
    N_bkgr_out[1] = iBkgr->getVal()*N_bkg.getError();

    // ##########################################################
    // Plot the results
    // Draw Correlation Matrix
    TCanvas *cCorrMat = new TCanvas("cCorrMat","cCorrMat",700,600);
    DrawCorrelationMatrix(cCorrMat,fResFit);

    // Draw histogram with fit results
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    SetCanvas(cHist,kFALSE);

    RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(215),LineWidth(3));
    // Vertical axis
    fFrameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    fFrameM->GetYaxis()->SetTitleSize(0.05);
    fFrameM->GetYaxis()->SetTitleOffset(1.1);
    fFrameM->GetYaxis()->SetLabelSize(0.05);
    fFrameM->GetYaxis()->SetLabelOffset(0.01);
    fFrameM->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    fFrameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    fFrameM->GetXaxis()->SetTitleSize(0.05);
    fFrameM->GetXaxis()->SetLabelSize(0.05);
    fFrameM->GetXaxis()->SetDecimals(1);
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","data",fResFit->floatParsFinal().getSize());

    // -------------------------------------------------------------------------------- 
    // Legend1
    TLegend *l1 = new TLegend(0.09,0.70,0.3,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    l1->AddEntry((TObject*)0,Form("ZN class: %s", classes[ZNClass].Data()),"");
    l1->SetTextSize(0.042);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    TLegend *lTitle = new TLegend(0.325,0.88,0.95,0.935);
    lTitle->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    lTitle->SetTextSize(0.05);
    lTitle->SetBorderSize(0);
    lTitle->SetFillStyle(0);
    lTitle->Draw();

    // Legend2
    TLegend *l2 = new TLegend(0.50,0.42,0.95,0.87);
    //l2->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l2->AddEntry("DSCBAndBkgPdf","sum","L");
    //l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l2->AddEntry("DoubleSidedCB","J/#psi signal","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.1f #pm %.1f",N_Jpsi_out[0],N_Jpsi_out[1]),"");
    l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
    if(SigmaFixed) l2->AddEntry((TObject*)0,Form("#sigma = %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal()),"");
    else l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    //l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.3f, #alpha_{R} = %.3f", alpha_L.getVal(), (-1)*(alpha_R.getVal())),"");
    //l2->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f, #it{n}_{R} = %.2f", n_L.getVal(), n_R.getVal()),"");
    l2->AddEntry("BkgPdf","background","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{bkg} = %.1f #pm %.1f",N_bkgr_out[0],N_bkgr_out[1]),"");
    l2->AddEntry((TObject*)0,Form("#lambda = %.3f #pm %.3f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
    l2->SetTextSize(0.042);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // Prepare path
    TString pathOut = Form("Results/ZNClasses/InvMassFits/%ibins/", nPtBins);

    // Print the plots
    cHist->Print((pathOut + Form("bin%i_%s.pdf", iPtBin, classes[ZNClass].Data())).Data());
    cHist->Print((pathOut + Form("bin%i_%s.png", iPtBin, classes[ZNClass].Data())).Data());
    cCorrMat->Print((pathOut + Form("CorrM/bin%i_%s_cm.pdf", iPtBin, classes[ZNClass].Data())).Data());
    cCorrMat->Print((pathOut + Form("CorrM/bin%i_%s_cm.png", iPtBin, classes[ZNClass].Data())).Data());

    return;
}

void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit){

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    cCorrMat->SetTopMargin(0.03);
    cCorrMat->SetBottomMargin(0.12);
    cCorrMat->SetRightMargin(0.19);
    cCorrMat->SetLeftMargin(0.13);

    TH2* hCorr = fResFit->correlationHist();
    hCorr->SetMarkerSize(3.);

    hCorr->GetXaxis()->SetBinLabel(1,"#it{M}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{bkg}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(4,"#lambda");
    hCorr->GetXaxis()->SetBinLabel(5,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(1,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(2,"#lambda");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{bkg}");
    hCorr->GetYaxis()->SetBinLabel(5,"#it{M}_{J/#psi}");

    hCorr->GetXaxis()->SetLabelSize(0.1);
    hCorr->GetXaxis()->SetLabelOffset(0.012);
    hCorr->GetYaxis()->SetLabelSize(0.1);
    hCorr->GetZaxis()->SetLabelSize(0.07);
    // https://root-forum.cern.ch/t/colz-color-palette-font-and-size/15263
    hCorr->Draw("colz,text");

    return;
}

void SetCanvas(TCanvas *c, Bool_t bLogScale){

    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}