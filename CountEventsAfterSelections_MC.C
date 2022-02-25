// CountEventsAfterSelections_MC.c
// David Grund, Sep 19, 2021
// To compare the remaining number of MC events after applying the selections

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
// my headers
#include "AnalysisManager.h"

Int_t counter[16] = { 0 };
Int_t counter_1[16] = { 0 };
Int_t counter_2[16] = { 0 };

Bool_t cuts_AOD_vs_ESD = kFALSE;
Bool_t cuts_pass1_vs_pass3 = kTRUE;
Bool_t PID_calibrated = kFALSE;

Bool_t EventPassed_MC_noSPDmatch(Bool_t isESD);
Bool_t EventPassed_MC_pass1();
Bool_t EventPassed_MC_pass3();

void CountEventsAfterSelections_MC(){

    // pass1
    TString str_file_pass1, str_file_pass3, str_tree_pass1, str_tree_pass3;
    if(!PID_calibrated){
        str_file_pass1 = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree_pass1 = "AnalysisOutput/fTreeJPsiMCRec";
        str_file_pass3 = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree_pass3 = "AnalysisOutput/fTreeJpsi";
    } else {
        str_file_pass1 = "Trees/AnalysisDataMC_pass1/kIncohJpsiToMu_calibSigmasTPC.root";
        str_tree_pass1 = "fTreeJpsi";
        str_file_pass3 = "Trees/AnalysisDataMC_pass3/kIncohJpsiToMu_calibSigmasTPC.root";
        str_tree_pass3 = "fTreeJpsi";
    }

    TFile *f_pass1 = TFile::Open(str_file_pass1.Data(), "read");
    if(f_pass1) Printf("Input data loaded.");

    TTree *t_pass1 = dynamic_cast<TTree*> (f_pass1->Get(str_tree_pass1.Data()));
    if(t_pass1) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(t_pass1, kFALSE);

    TList *l_pass1 = NULL;
    TH1F *hCounterCuts_pass1 = NULL;
    if(!PID_calibrated){
        l_pass1 = dynamic_cast<TList*> (f_pass1->Get("AnalysisOutput/fOutputList"));
        if(l_pass1) Printf("Input list loaded.");

        hCounterCuts_pass1 = (TH1F*)l_pass1->FindObject("hCounterCuts");
        if(hCounterCuts_pass1) Printf("Input histogram loaded.");
    }

    // pass3
    TFile *f_pass3 = TFile::Open(str_file_pass3.Data(), "read");
    if(f_pass3) Printf("Input data loaded.");

    TTree *t_pass3 = dynamic_cast<TTree*> (f_pass3->Get(str_tree_pass3.Data()));
    if(t_pass3) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(t_pass3, kTRUE);

    TList *l_pass3 = NULL;
    TH1F *hCounterCuts_pass3 = NULL;
    if(!PID_calibrated){
        l_pass3 = dynamic_cast<TList*> (f_pass3->Get("AnalysisOutput/fOutputList"));
        if(l_pass3) Printf("Input list loaded.");

        hCounterCuts_pass3 = (TH1F*)l_pass3->FindObject("hCounterCuts");
        if(hCounterCuts_pass3) Printf("Input histogram loaded.");
    }

    // AOD
    TFile *f_AOD = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
    if(f_AOD) Printf("Input data loaded.");

    TTree *t_AOD = dynamic_cast<TTree*> (f_AOD->Get("analysisTree"));
    if(t_AOD) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec_AOD(t_AOD);

    if(cuts_AOD_vs_ESD){

        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);
            EventPassed_MC_noSPDmatch(kTRUE);

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Save the values to counter_1 and set counter to zeros
        for(Int_t i = 0; i < 13; i++){
            counter_1[i] = counter[i];
            counter[i] = 0;
        }        

        Printf("%lli entries found in the tree.", t_AOD->GetEntries());
        nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_AOD->GetEntries(); iEntry++){
            t_AOD->GetEntry(iEntry);
            EventPassed_MC_noSPDmatch(kFALSE);

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Save the values to counter_2
        for(Int_t i = 0; i < 12; i++){
            counter_2[i] = counter[i];
        }
        counter_2[12] = counter_2[11];

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections_MC/kIncohJpsiToMu_AOD_vs_ESD.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "kIncohJpsiToMu:      \tESD \tAOD \n";
        outfile << "After previous cuts: \t" << counter_1[0] << "\t" << counter_2[0] << "\n";
        outfile << "4) CCUP31 trigger:   \t" << counter_1[1] << "\t" << counter_2[1] << "\n";
        outfile << "5a) ADA offline veto:\t" << counter_1[2] << "\t" << counter_2[2] << "\n";
        outfile << "5b) ADC offline veto:\t" << counter_1[3] << "\t" << counter_2[3] << "\n";
        outfile << "6a) V0A offline veto:\t" << counter_1[4] << "\t" << counter_2[4] << "\n";
        outfile << "6b) V0C offline veto:\t" << counter_1[5] << "\t" << counter_2[5] << "\n";
        outfile << "7) muon pairs only:  \t" << counter_1[6] << "\t" << counter_2[6] << "\n";
        outfile << "8) dilept |y| < 0.8: \t" << counter_1[7] << "\t" << counter_2[7] << "\n";
        outfile << "9) trks |eta| < 0.8: \t" << counter_1[8] << "\t" << counter_2[8] << "\n";
        outfile << "10) opposite charges:\t" << counter_1[9] << "\t" << counter_2[9] << "\n";
        outfile << "11) mass 2.2 to 4.5: \t" << counter_1[10] << "\t" << counter_2[10] << "\n";
        outfile << "12) dimuon pt > 0.2: \t" << counter_1[11] << "\t" << counter_2[11] << "\n";

        outfile.close();
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 16; i++){
        counter[i] = 0;
        counter_1[i] = 0;
        counter_2[i] = 0;
    }    

    if(cuts_pass1_vs_pass3){

        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);
            EventPassed_MC_pass1();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Save the values to counter_1 and set counter to zeros
        for(Int_t i = 0; i < 14; i++){
            counter_1[i] = counter[i];
            counter[i] = 0;
        }        

        Printf("%lli entries found in the tree.", t_pass3->GetEntries());
        nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass3->GetEntries(); iEntry++){
            t_pass3->GetEntry(iEntry);
            EventPassed_MC_pass3();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Save the values to counter_2
        for(Int_t i = 0; i < 16; i++){
            counter_2[i] = counter[i];
        }

        // Print the numbers:
        TString name;
        if(!PID_calibrated) name = "Results/CountEventsAfterSelections_MC/kIncohJpsiToMu_pass1_vs_pass3.txt";
        else                name = "Results/CountEventsAfterSelections_MC/kIncohJpsiToMu_pass1_vs_pass3_PIDcalibrated.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "kIncohJpsiToMu:      \tpass1 \tpass3 \tratio \n";
        outfile << "After previous cuts: \t" << counter_1[0] << "\t" << counter_2[2] << "\t" << std::fixed << std::setprecision(3) << counter_1[0]/(Double_t)counter_2[2] << "\n";
        outfile << "4) CCUP31 trigger:   \t" << counter_1[1] << "\t" << counter_2[3] << "\t" << std::fixed << std::setprecision(3) << counter_1[1]/(Double_t)counter_2[3] << "\n";
        outfile << "5a) ADA offline veto:\t" << counter_1[2] << "\t" << counter_2[4] << "\t" << std::fixed << std::setprecision(3) << counter_1[2]/(Double_t)counter_2[4] << "\n";
        outfile << "5b) ADC offline veto:\t" << counter_1[3] << "\t" << counter_2[5] << "\t" << std::fixed << std::setprecision(3) << counter_1[3]/(Double_t)counter_2[5] << "\n";
        outfile << "6a) V0A offline veto:\t" << counter_1[4] << "\t" << counter_2[6] << "\t" << std::fixed << std::setprecision(3) << counter_1[4]/(Double_t)counter_2[6] << "\n";
        outfile << "6b) V0C offline veto:\t" << counter_1[5] << "\t" << counter_2[7] << "\t" << std::fixed << std::setprecision(3) << counter_1[5]/(Double_t)counter_2[7] << "\n";
        outfile << "7) SPD match FOhits: \t" << counter_1[6] << "\t" << counter_2[8] << "\t" << std::fixed << std::setprecision(3) << counter_1[6]/(Double_t)counter_2[8] << "\n";
        outfile << "8) muon pairs only:  \t" << counter_1[7] << "\t" << counter_2[9] << "\t" << std::fixed << std::setprecision(3) << counter_1[7]/(Double_t)counter_2[9] << "\n";
        outfile << "9) dilept |y| < 0.8: \t" << counter_1[8] << "\t" << counter_2[10] << "\t" << std::fixed << std::setprecision(3) << counter_1[8]/(Double_t)counter_2[10] << "\n";
        outfile << "10) trks |eta| < 0.8: \t" << counter_1[9] << "\t" << counter_2[11] << "\t" << std::fixed << std::setprecision(3) << counter_1[9]/(Double_t)counter_2[11] << "\n";
        outfile << "11) opposite charges:\t" << counter_1[10] << "\t" << counter_2[12] << "\t" << std::fixed << std::setprecision(3) << counter_1[10]/(Double_t)counter_2[12] << "\n";
        outfile << "12) mass 2.2 to 4.5: \t" << counter_1[11] << "\t" << counter_2[13] << "\t" << std::fixed << std::setprecision(3) << counter_1[11]/(Double_t)counter_2[13] << "\n";
        outfile << "13) dimuon pt > 0.2: \t" << counter_1[12] << "\t" << counter_2[14] << "\t" << std::fixed << std::setprecision(3) << counter_1[12]/(Double_t)counter_2[14] << "\n";
        outfile << "14) mass 3.0 to 3.2: \t" << counter_1[13] << "\t" << counter_2[15] << "\t" << std::fixed << std::setprecision(3) << counter_1[13]/(Double_t)counter_2[15] << "\n";

        outfile.close();
    }

    return;
}

Bool_t EventPassed_MC_noSPDmatch(Bool_t isESD){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    counter[0]++;

    // 4) Central UPC trigger CCUP31
    Bool_t CCUP31 = kFALSE;
    if(
        !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
        !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
        !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
        !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
        fTriggerInputsMC[10] &&  //  0STG (SPD topological)
        fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
    ) CCUP31 = kTRUE;
    if(!CCUP31) return kFALSE;
    counter[1]++;

    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[2]++;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[3]++;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[4]++;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[5]++;

    // 7) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[6]++;

    // 8) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[7]++;

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[8]++;

    // 10) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[9]++;

    // 11) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[10]++;

    // 12) Transverse momentum cut
    if(!(fPt > 0.20)) return kFALSE;
    counter[11]++;

    // 13) SPD cluster matches FOhits
    if(isESD){
        if(!(fMatchingSPD == kTRUE)) return kFALSE;
        counter[12]++;
    }

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassed_MC_pass1(){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    counter[0]++;

    // 4) Central UPC trigger CCUP31
    Bool_t CCUP31 = kFALSE;
    if(
        !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
        !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
        !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
        !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
        fTriggerInputsMC[10] &&  //  0STG (SPD topological)
        fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
    ) CCUP31 = kTRUE;
    if(!CCUP31) return kFALSE;
    counter[1]++;

    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[2]++;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[3]++;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[4]++;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[5]++;

    // 7) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[6]++;

    // 8) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[7]++;

    // 9) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[8]++;

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[9]++;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[10]++;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[11]++;

    // 13) Transverse momentum cut
    if(!(fPt > 0.20)) return kFALSE;
    counter[12]++;

    // 14) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counter[13]++;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassed_MC_pass3(){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    counter[0]++;

    // 2) At least two tracks associated with the vertex
    if(fVertexContrib < 2) return kFALSE;
    counter[1]++;

    // 3) Distance from the IP lower than 15 cm
    if(fVertexZ > 15) return kFALSE;
    counter[2]++;

    // 4) Central UPC trigger CCUP31
    Bool_t CCUP31 = kFALSE;
    if(
        !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
        !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
        !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
        !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
        fTriggerInputsMC[10] &&  //  0STG (SPD topological)
        fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
    ) CCUP31 = kTRUE;
    if(!CCUP31) return kFALSE;
    counter[3]++;

    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[4]++;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[5]++;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[6]++;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[7]++;

    // 7) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[8]++;

    // 8) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[9]++;

    // 9) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[10]++;

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[11]++;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[12]++;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[13]++;

    // 13) Transverse momentum cut
    if(!(fPt > 0.20)) return kFALSE;
    counter[14]++;

    // 14) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counter[15]++;

    // Event passed all the selections =>
    return kTRUE;
}