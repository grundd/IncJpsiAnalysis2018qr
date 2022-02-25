// MuonPID.C
// David Grund, Feb 18, 2022

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"

void PlotAndFitHistograms(Bool_t MC, Int_t iMC = -1);
void ShiftPIDSignal(Bool_t pass3);

void MuonPID(){

    // data
    PlotAndFitHistograms(kFALSE);
    
    // MC
    for(Int_t i = 0; i < 5; i++) PlotAndFitHistograms(kTRUE, i);

    // pass1
    //ShiftPIDSignal(kFALSE);

    // pass3
    //ShiftPIDSignal(kTRUE);

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

    // histograms where the differences will be stored
    Double_t range = 0.02;
    TH1D *hPt = new TH1D("hPt", "hPt", 100, -range, range);
    TH1D *hPhi = new TH1D("hPhi", "hPhi", 100, -range, range);
    TH1D *hY = new TH1D("hY", "hY", 100, -range, range);
    TH1D *hM = new TH1D("hM", "hM", 100, -range, range);

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
    } else {
        fShiftMu = 1.404;
    }

    Printf("%lli entries found in the tree.", t_in->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in->GetEntries(); iEntry++){
        t_in->GetEntry(iEntry);

        fTrk1SigIfMu = fTrk1SigIfMu - fShiftMu;
        //fTrk1SigIfEl = fTrk1SigIfEl;
        fTrk2SigIfMu = fTrk2SigIfMu - fShiftMu;
        //fTrk2SigIfEl = fTrk2SigIfEl;

        // store the old values of J/psi kinematic variables
        Double_t fPt_old = fPt;
        Double_t fPhi_old = fPhi;
        Double_t fY_old = fY;
        Double_t fM_old = fM;

        // recalculate J/psi kinematics after assigning a proper mass to tracks
        Double_t isMuonPair = fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu;
        Double_t isElectronPair = fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl;
        Double_t massTracks = -1;
        if(isMuonPair < isElectronPair) massTracks = 0.105658; // GeV/c^2
        else                            massTracks = 0.000511; // GeV/c^2

        TLorentzVector vTrk1, vTrk2;
        vTrk1.SetPtEtaPhiM(fPt1, fEta1, fPhi1, massTracks);
        vTrk2.SetPtEtaPhiM(fPt2, fEta2, fPhi2, massTracks);
        TLorentzVector vTrkTrk = vTrk1 + vTrk2;

        // set tree variables
        fPt = vTrkTrk.Pt(); 
        fPhi = vTrkTrk.Phi();
        fY = vTrkTrk.Rapidity(); 
        fM = vTrkTrk.M();

        // calculate the differences
        Double_t Delta_pt = fPt_old - fPt;
        Double_t Delta_phi = fPhi_old - fPhi;
        Double_t Delta_y = fY_old - fY;
        Double_t Delta_m = fM_old - fM;
        // fill the histograms
        hPt->Fill(Delta_pt);
        hPhi->Fill(Delta_phi);
        hY->Fill(Delta_y);
        hM->Fill(Delta_m);

        fTreeJpsi->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }

    }

    f_out->Write("",TObject::kWriteDelete);

    // plot histograms
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(2,2,0,0);
    c->cd(1); hPt->Draw();
    c->cd(2); hPhi->Draw();
    c->cd(3); hY->Draw();
    c->cd(4); hM->Draw();

    TString str_diff = "";
    if(!pass3)  str_diff = "Results/MuonPID/differences_pass1";
    else        str_diff = "Results/MuonPID/differences_pass3";
    c->Print((str_diff + ".png").Data());
    c->Print((str_diff + ".pdf").Data());

    return;
}

void PlotAndFitHistograms(Bool_t MC, Int_t iMC){
    // if MC, fit histograms with gaussian

    // pass1
    TString str_file_pass1 = "";
    TString str_tree_pass1 = "";
    if(MC){
        if(iMC == 0) str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kCohJpsiToMu.root";
        if(iMC == 1) str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kCohPsi2sToMuPi.root";
        if(iMC == 2) str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohJpsiToMu.root";
        if(iMC == 3) str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohPsi2sToMuPi.root";
        if(iMC == 4) str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kTwoGammaToMuMedium.root";
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
        if(iMC == 0) str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohJpsiToMu.root";
        if(iMC == 1) str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohPsi2sToMuPi.root";
        if(iMC == 2) str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root";
        if(iMC == 3) str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohPsi2sToMuPi.root";
        if(iMC == 4) str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kTwoGammaToMuMedium.root";
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
    Double_t range(15.);
    TH1D *hSigmaTPC[4] = { NULL };
    hSigmaTPC[0] = new TH1D("pass1_SigIfMu", "pass1_SigIfMu", nBins, -range, range);
    hSigmaTPC[1] = new TH1D("pass1_SigIfEl", "pass1_SigIfEl", nBins, -range, range);
    hSigmaTPC[2] = new TH1D("pass3_SigIfMu", "pass3_SigIfMu", nBins, -range, range);
    hSigmaTPC[3] = new TH1D("pass3_SigIfEl", "pass3_SigIfEl", nBins, -range, range);

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

    // plot histograms
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(2,2,0,0);
    for(Int_t i = 0; i < 4; i++){
        c->cd(i+1);
        hSigmaTPC[i]->Draw();
    }

    // fit histograms
    TF1 *fGauss[4] = { 0 };
    if(MC){
        for(Int_t i = 0; i < 4; i++){
            fGauss[i] = new TF1(Form("fGauss%i", i+1), "gaus", -range, range);
            hSigmaTPC[i]->Fit(fGauss[i]);
        }
    }

    // print the plots and (if MC) print the values of the fitted means
    TString str_f_out = "";
    if(MC){
        if(iMC == 0) str_f_out = "MC_kCohJpsiToMu";
        if(iMC == 1) str_f_out = "MC_kCohPsi2sToMuPi";
        if(iMC == 2) str_f_out = "MC_kIncohJpsiToMu";
        if(iMC == 3) str_f_out = "MC_kIncohPsi2sToMuPi";
        if(iMC == 4) str_f_out = "MC_kTwoGammaToMuMedium";

        ofstream f_out(("Results/MuonPID/means_" + str_f_out + ".txt").Data());
        f_out << std::fixed << std::setprecision(3);
        for(Int_t i = 0; i < 4; i++){
            f_out << fGauss[i]->GetParameter(1) << endl;
        }
        f_out.close();
        Printf("*** Results printed to %s.***", ("means_" + str_f_out + ".txt").Data());

    } else str_f_out = "data_pass1_vs_pass3";

    c->Print(("Results/MuonPID/" + str_f_out + ".png").Data());
    c->Print(("Results/MuonPID/" + str_f_out + ".pdf").Data());

    return;
}