// CountEventsAfterSelectionsMC.c
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

Int_t counter[13] = {0};
Int_t counterESD[13] = {0};
Int_t counterAOD[13] = {0};

Bool_t EventPassedNoMatchingMC(Bool_t isESD);

void CountEventsAfterSelectionsMC(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(fTreeIn);

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        EventPassedNoMatchingMC(kTRUE);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Save the values to counterESD and set counter to zeros
    for(Int_t i = 0; i < 13; i++){
        counterESD[i] = counter[i];
        counter[i] = 0;
    }

    TFile *fFileIn2 = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
    if(fFileIn2) Printf("Input data loaded.");

    TTree *fTreeIn2 = dynamic_cast<TTree*> (fFileIn2->Get("analysisTree"));
    if(fTreeIn2) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec_AOD(fTreeIn2);

    Printf("%lli entries found in the tree.", fTreeIn2->GetEntries());
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn2->GetEntries(); iEntry++){
        fTreeIn2->GetEntry(iEntry);
        EventPassedNoMatchingMC(kFALSE);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Save the values to counterAOD
    for(Int_t i = 0; i < 12; i++){
        counterAOD[i] = counter[i];
    }
    counterAOD[12] = counterAOD[11];

    // Print the numbers:
    TString name = "Results/CountEventsAfterSelectionsMC/kIncohJpsiToMu.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
    outfile << "kIncohJpsiToMu:      \tESD \tAOD \n";
    outfile << "After previous cuts: \t" << counterESD[0] << "\t" << counterAOD[0] << "\n";
    outfile << "4) CCUP31 trigger:   \t" << counterESD[1] << "\t" << counterAOD[1] << "\n";
    outfile << "5a) ADA offln veto:  \t" << counterESD[2] << "\t" << counterAOD[2] << "\n";
    outfile << "5b) ADC offln veto:  \t" << counterESD[3] << "\t" << counterAOD[3] << "\n";
    outfile << "6a) V0A offln veto:  \t" << counterESD[4] << "\t" << counterAOD[4] << "\n";
    outfile << "6b) V0C offln veto:  \t" << counterESD[5] << "\t" << counterAOD[5] << "\n";
    outfile << "7) dilept rapidity:  \t" << counterESD[6] << "\t" << counterAOD[6] << "\n";
    outfile << "8) eta of both trks: \t" << counterESD[7] << "\t" << counterAOD[7] << "\n";
    outfile << "9) oposite charges:  \t" << counterESD[8] << "\t" << counterAOD[8] << "\n";
    outfile << "10) muon pairs only: \t" << counterESD[9] << "\t" << counterAOD[9] << "\n";
    outfile << "11) mass 2.2 to 4.5: \t" << counterESD[10] << "\t" << counterAOD[10] << "\n";
    outfile << "12) dimuon pt > 0.2: \t" << counterESD[11] << "\t" << counterAOD[11] << "\n";
    outfile << "13) SPD match FOhits:\t" << counterESD[12] << "\t" << counterAOD[12] << "\n";

    outfile.close();

    return;
}

Bool_t EventPassedNoMatchingMC(Bool_t isESD){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    counter[0]++;

    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
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

    // 7) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[6]++;

    // 8) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[7]++;

    // 9) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[8]++;

    // 10) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
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