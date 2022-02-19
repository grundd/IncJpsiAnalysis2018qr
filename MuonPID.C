// MuonPID.C
// David Grund, Feb 18, 2022

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"

void PlotAndFitHistograms(Bool_t MC);
void ShiftPIDSignal(Bool_t pass3);

void MuonPID(){

    // pass1
    //PlotAndFitHistograms(kFALSE);
    
    // pass3
    //PlotAndFitHistograms(kTRUE);

    // pass1
    //ShiftPIDSignal(kFALSE);

    // pass3
    ShiftPIDSignal(kTRUE);

    return;
}

void ShiftPIDSignal(Bool_t pass3){

    // pass1
    TString str_file = "";
    TString str_tree = "";
    if(!pass3){
        str_file = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree = "AnalysisOutput/fTreeJPsiMCRec";
    } else {
        str_file = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree = "AnalysisOutput/fTreeJpsi";
    }  
        
    TFile *f_in = TFile::Open(str_file.Data(), "read");
    if(f_in) Printf("Input data loaded.");

    TTree *t_in = dynamic_cast<TTree*> (f_in->Get(str_tree.Data()));
    if(t_in) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(t_in, pass3);

    // create new file
    TString str_f_out = "";
    if(!pass3)  str_f_out = "Trees/AnalysisDataMC_pass1/kIncohJpsiToMu_calibSigmasTPC.root";
    else        str_f_out = "Trees/AnalysisDataMC_pass3/kIncohJpsiToMu_calibSigmasTPC.root";
    TFile *f_out = new TFile(str_f_out.Data(),"RECREATE");

    // create new tree
    TTree *fTreeJpsi = new TTree("fTreeJpsi", "fTreeJpsi");
    // Basic things:
    fTreeJpsi->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeJpsi->Branch("fTriggerName", &fTriggerName);
    // PID, sigmas:
    if(pass3){
        fTreeJpsi->Branch("fTrk1dEdx", &fTrk1dEdx, "fTrk1dEdx/D");
        fTreeJpsi->Branch("fTrk2dEdx", &fTrk2dEdx, "fTrk2dEdx/D");
    }
    fTreeJpsi->Branch("fTrk1SigIfMu", &fTrk1SigIfMu, "fTrk1SigIfMu/D");
    fTreeJpsi->Branch("fTrk1SigIfEl", &fTrk1SigIfEl, "fTrk1SigIfEl/D");
    fTreeJpsi->Branch("fTrk2SigIfMu", &fTrk2SigIfMu, "fTrk2SigIfMu/D");
    fTreeJpsi->Branch("fTrk2SigIfEl", &fTrk2SigIfEl, "fTrk2SigIfEl/D");
    // Kinematics:
    fTreeJpsi->Branch("fPt", &fPt, "fPt/D");
    fTreeJpsi->Branch("fPhi", &fPhi, "fPhi/D");
    fTreeJpsi->Branch("fY", &fY, "fY/D");
    fTreeJpsi->Branch("fM", &fM, "fM/D");
    // Two tracks:
    fTreeJpsi->Branch("fPt1", &fPt1, "fPt1/D");
    fTreeJpsi->Branch("fPt2", &fPt2, "fPt2/D");
    fTreeJpsi->Branch("fEta1", &fEta1, "fEta1/D");
    fTreeJpsi->Branch("fEta2", &fEta2, "fEta2/D");
    fTreeJpsi->Branch("fPhi1", &fPhi1, "fPhi1/D");
    fTreeJpsi->Branch("fPhi2", &fPhi2, "fPhi2/D");
    fTreeJpsi->Branch("fQ1", &fQ1, "fQ1/D");
    fTreeJpsi->Branch("fQ2", &fQ2, "fQ2/D");
    // Vertex info:
    if(pass3){
        fTreeJpsi->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
        fTreeJpsi->Branch("fVertexContrib", &fVertexContrib, "fVertexContrib/I"); 
    }
    // Info from the detectors:
    // ZDC:
    fTreeJpsi->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    fTreeJpsi->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    fTreeJpsi->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    fTreeJpsi->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    // V0, AD:
    fTreeJpsi->Branch("fV0A_dec", &fV0A_dec, "fV0A_dec/I");
    fTreeJpsi->Branch("fV0C_dec", &fV0C_dec, "fV0C_dec/I");
    fTreeJpsi->Branch("fADA_dec", &fADA_dec, "fADA_dec/I");
    fTreeJpsi->Branch("fADC_dec", &fADC_dec, "fADC_dec/I");
    if(!pass3){
        fTreeJpsi->Branch("fV0A_time", &fV0A_time, "fV0A_time/D");
        fTreeJpsi->Branch("fV0C_time", &fV0C_time, "fV0C_time/D");
        fTreeJpsi->Branch("fADA_time", &fADA_time, "fADA_time/D");
        fTreeJpsi->Branch("fADC_time", &fADC_time, "fADC_time/D");        
    }
    // Matching SPD clusters with FOhits:
    fTreeJpsi->Branch("fMatchingSPD", &fMatchingSPD, "fMatchingSPD/O");
    // Replayed trigger inputs:
    fTreeJpsi->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[11]/O"); 
    // Kinematics, MC gen:
    fTreeJpsi->Branch("fPtGen", &fPtGen, "fPtGen/D");
    fTreeJpsi->Branch("fMGen", &fMGen, "fMGen/D");
    fTreeJpsi->Branch("fYGen", &fYGen, "fYGen/D");
    fTreeJpsi->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");

    Double_t fShiftMu = 0;
    Double_t fShiftEl = 0;
    if(!pass3){
        fShiftMu = 0.350;
        fShiftEl = -2.628;
    } else {
        fShiftMu = 1.404;
        fShiftEl = -1.797;
    }

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        fTrk1SigIfMu = fTrk1SigIfMu - fShiftMu;
        fTrk1SigIfEl = fTrk1SigIfEl - fShiftEl;
        fTrk2SigIfMu = fTrk2SigIfMu - fShiftMu;
        fTrk2SigIfEl = fTrk2SigIfEl - fShiftEl;

        fTreeJpsi->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

    }

    f_out->Write("",TObject::kWriteDelete);

    return;
}

void PlotAndFitHistograms(Bool_t MC){
    // if MC, fit histograms with gaussian

    // pass1
    TString str_file_pass1 = "";
    TString str_tree_pass1 = "";
    if(MC){
        str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree_pass1 = "AnalysisOutput/fTreeJPsiMCRec";
    } else {
        str_file_pass1 = "Trees/AnalysisData_pass1/AnalysisResultsLHC18qrMerged.root";
        str_tree_pass1 = "AnalysisOutput/fTreeJPsi";
    }  
        
    TFile *f_pass1 = TFile::Open(str_file_pass1.Data(), "read");
    if(f_pass1) Printf("Input data loaded.");

    TTree *t_pass1 = dynamic_cast<TTree*> (f_pass1->Get(str_tree_pass1.Data()));
    if(t_pass1) Printf("Input tree loaded.");

    if(MC)  ConnectTreeVariablesMCRec(t_pass1, kFALSE);
    else    ConnectTreeVariables(t_pass1, kFALSE);

    // pass3
    TString str_file_pass3 = "";
    TString str_tree_pass3 = "";
    if(MC){
        str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree_pass3 = "AnalysisOutput/fTreeJpsi";
    } else {
        str_file_pass3 = "Trees/AnalysisData_pass3/AnalysisResults.root";
        str_tree_pass3 = "AnalysisOutput/fTreeJpsi";
    }  

    TFile *f_pass3 = TFile::Open(str_file_pass3.Data(), "read");
    if(f_pass3) Printf("Input data loaded.");

    TTree *t_pass3 = dynamic_cast<TTree*> (f_pass3->Get(str_tree_pass3.Data()));
    if(t_pass3) Printf("Input tree loaded.");

    if(MC)  ConnectTreeVariablesMCRec(t_pass3, kTRUE);
    else    ConnectTreeVariables(t_pass3, kTRUE);

    // histograms
    Int_t nBins = 100;
    Double_t low(-5.), high(5.);
    TH1D *hSigmaTPC[4] = { NULL };
    hSigmaTPC[0] = new TH1D("pass1_SigIfMu", "pass1_SigIfMu", nBins, low, high);
    hSigmaTPC[1] = new TH1D("pass1_SigIfEl", "pass1_SigIfEl", nBins, low, high);
    hSigmaTPC[2] = new TH1D("pass3_SigIfMu", "pass3_SigIfMu", nBins, low, high);
    hSigmaTPC[3] = new TH1D("pass3_SigIfEl", "pass3_SigIfEl", nBins, low, high);

    Printf("%lli entries found in the tree.", t_pass1->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
        t_pass1->GetEntry(iEntry);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        // Cuts 0-4 offline
        // if MC, trigger has to be replayed:
        if(MC){
            Bool_t CCUP31 = kFALSE;
            if(
                !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
                !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
                !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
                !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
                fTriggerInputsMC[10] &&  //  0STG (SPD topological)
                fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
            ) CCUP31 = kTRUE;
            if(!CCUP31) continue;    
        }
        // Cut 5a: ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;
        // Cut 5b: ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;
        // Cut 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;
        // Cut 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;
        // Cut 7: SPD cluster matches FOhits
        if(fMatchingSPD == kFALSE) continue;
        // Cut 8: Dilepton rapidity |y| < 0.8
        if(!(abs(fY) < 0.8)) continue;
        // Cut 9: Pseudorapidity of both tracks |eta| < 0.8
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // Cut 10: Tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;

        hSigmaTPC[0]->Fill(fTrk1SigIfMu);
        hSigmaTPC[0]->Fill(fTrk2SigIfMu);
        hSigmaTPC[1]->Fill(fTrk1SigIfEl);
        hSigmaTPC[1]->Fill(fTrk2SigIfEl);  
    }

    Printf("%lli entries found in the tree.", t_pass3->GetEntries());
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_pass3->GetEntries(); iEntry++){
        t_pass3->GetEntry(iEntry);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

        // Cuts 0-2 offline
        // if MC, trigger has to be replayed:
        if(MC){
            Bool_t CCUP31 = kFALSE;
            if(
                !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
                !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
                !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
                !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
                fTriggerInputsMC[10] &&  //  0STG (SPD topological)
                fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
            ) CCUP31 = kTRUE;
            if(!CCUP31) continue;    
        }
        // Cut 3: more than 2 tracks associated with the primary vertex
        if(fVertexContrib < 2) continue;
        // Cut 4: z-distance from the IP lower than 15 cm
        if(fVertexZ > 15) continue;
        // Cut 5a: ADA offline veto (no effect on MC)
        if(!(fADA_dec == 0)) continue;
        // Cut 5b: ADC offline veto (no effect on MC)
        if(!(fADC_dec == 0)) continue;
        // Cut 6a) V0A offline veto (no effect on MC)
        if(!(fV0A_dec == 0)) continue;
        // Cut 6b) V0C offline veto (no effect on MC)
        if(!(fV0C_dec == 0)) continue;
        // Cut 7: SPD cluster matches FOhits
        if(fMatchingSPD == kFALSE) continue;
        // Cut 8: Dilepton rapidity |y| < 0.8
        if(!(abs(fY) < 0.8)) continue;
        // Cut 9: Pseudorapidity of both tracks |eta| < 0.8
        if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) continue;
        // Cut 10: Tracks have opposite charges
        if(!(fQ1 * fQ2 < 0)) continue;

        hSigmaTPC[2]->Fill(fTrk1SigIfMu);
        hSigmaTPC[2]->Fill(fTrk2SigIfMu);
        hSigmaTPC[3]->Fill(fTrk1SigIfEl);
        hSigmaTPC[3]->Fill(fTrk2SigIfEl);  
    }

    // fit histograms
    TF1 *fGauss[4] = { 0 };
    if(MC){
        for(Int_t i = 0; i < 4; i++){
            fGauss[i] = new TF1(Form("fGauss%i", i+1), "gaus", low, high);
            hSigmaTPC[i]->Fit(fGauss[i]);
        }
    }

    // plot histograms
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(2,2,0,0);

    for(Int_t i = 0; i < 4; i++){
        c->cd(i+1);
        hSigmaTPC[i]->Draw();
    }

    // print the means
    if(MC){
        TString str_f_out = "Results/MuonPID/fitted_means.txt";
        ofstream f_out(str_f_out.Data());
        f_out << std::fixed << std::setprecision(3);
        for(Int_t i = 0; i < 4; i++){
            f_out << fGauss[i]->GetParameter(1) << endl;
        }
        f_out.close();
        Printf("*** Results printed to %s.***", str_f_out.Data());
    }


    TString str_file_out = "";
    if(MC)  str_file_out = "Results/MuonPID/MC_muons_pass1_vs_pass3";
    else    str_file_out = "Results/MuonPID/data_muons_pass1_vs_pass3";
    c->Print((str_file_out + ".png").Data());
    c->Print((str_file_out + ".pdf").Data());

    return;
}