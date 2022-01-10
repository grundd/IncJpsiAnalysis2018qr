// CountEventsAfterSelections.c
// David Grund, Sep 7, 2021
// To calculate the remaining number of events after applying the selections

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

Bool_t CutsMine = kTRUE;
Bool_t CutsRoman = kTRUE;
Bool_t CutsNoMatching = kTRUE;

Bool_t EventPassedMyCuts();
Bool_t EventPassedRoman();
Bool_t EventPassedNoMatching();

void CountEventsAfterSelections(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    TList *fListIn = dynamic_cast<TList*> (fFileIn->Get("AnalysisOutput/fOutputList"));
    if(fListIn) Printf("Input list loaded.");

    TH1F *hCounterCuts = (TH1F*)fListIn->FindObject("hCounterCuts");
    if(hCounterCuts) Printf("Input histogram loaded.");

    hCounterCuts->Draw();

    if(CutsMine){

        Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
            fTreeIn->GetEntry(iEntry);
            EventPassedMyCuts();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/CutsMine.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0) event non-empty: \t" << hCounterCuts->GetBinContent(1) << "\n";
        outfile << "1) vert # contribs: \t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "2) vert Z distance: \t" << hCounterCuts->GetBinContent(3) << "\n";
        outfile << "3) two good tracks: \t" << hCounterCuts->GetBinContent(4) << "\n";
        outfile << "4) CCUP31 trigger:  \t" << hCounterCuts->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "5) SPD match FOhits:\t" << counter[1] << "\n";
        outfile << "6a) ADA offln veto: \t" << counter[2] << "\n";
        outfile << "6b) ADC offln veto: \t" << counter[3] << "\n";
        outfile << "7a) V0A offln veto: \t" << counter[4] << "\n";
        outfile << "7b) V0C offln veto: \t" << counter[5] << "\n";
        outfile << "8) dilept rapidity: \t" << counter[6] << "\n";
        outfile << "9) eta of both trks:\t" << counter[7] << "\n";
        outfile << "10) oposite charges:\t" << counter[8] << "\n";
        outfile << "11) muon pairs only:\t" << counter[9] << "\n";
        outfile << "12) mass 2.2 to 4.5:\t" << counter[10] << "\n";
        outfile << "13) p_T 0.2 to 1.0: \t" << counter[11] << "\n";
        outfile << "14) mass 3.0 to 3.2:\t" << counter[12] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 13; i++){
        counter[i] = 0;
    }

    if(CutsRoman){

        Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
            fTreeIn->GetEntry(iEntry);
            EventPassedRoman();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/CutsRoman.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0a) event non-empty:\t" << hCounterCuts->GetBinContent(1) << "\n";
        outfile << "0b) vert # contribs:\t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "0c) vert Z distance:\t" << hCounterCuts->GetBinContent(3) << "\n";
        outfile << "0d) two good tracks:\t" << hCounterCuts->GetBinContent(4) << "\n";
        outfile << "0e) CCUP31 trigger: \t" << hCounterCuts->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "1) SPD match FOhits:\t" << counter[1] << "\n";
        outfile << "2) muon pairs only: \t" << counter[2] << "\n";
        outfile << "3) mass 2.2 to 4.5: \t" << counter[3] << "\n";
        outfile << "4a) ADA offln veto: \t" << counter[4] << "\n";
        outfile << "4b) ADC offln veto: \t" << counter[5] << "\n";
        outfile << "5a) V0A offln veto: \t" << counter[6] << "\n";
        outfile << "5b) V0C offln veto: \t" << counter[7] << "\n";
        outfile << "6) dilept rapidity: \t" << counter[8] << "\n";
        outfile << "7) eta of both trks:\t" << counter[9] << "\n";
        outfile << "8) oposite charges: \t" << counter[10] << "\n";
        outfile << "9) dimuon pt < 0.11:\t" << counter[11] << "\n";
        outfile << "10) mass 3.0 to 3.2:\t" << counter[12] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 13; i++){
        counter[i] = 0;
    }

    if(CutsNoMatching){

        Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
            fTreeIn->GetEntry(iEntry);
            EventPassedNoMatching();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/CutsNoMatching.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0) event non-empty:  \t" << hCounterCuts->GetBinContent(1) << "\n";
        outfile << "1) vert # contribs:  \t" << hCounterCuts->GetBinContent(2) << "\n";
        outfile << "2) vert Z distance:  \t" << hCounterCuts->GetBinContent(3) << "\n";
        outfile << "3) two good tracks:  \t" << hCounterCuts->GetBinContent(4) << "\n";
        outfile << "4) CCUP31 trigger:   \t" << hCounterCuts->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "5a) ADA offln veto:  \t" << counter[1] << "\n";
        outfile << "5b) ADC offln veto:  \t" << counter[2] << "\n";
        outfile << "6a) V0A offln veto:  \t" << counter[3] << "\n";
        outfile << "6b) V0C offln veto:  \t" << counter[4] << "\n";
        outfile << "7) dilept rapidity:  \t" << counter[5] << "\n";
        outfile << "8) eta of both trks: \t" << counter[6] << "\n";
        outfile << "9) oposite charges:  \t" << counter[7] << "\n";
        outfile << "10) muon pairs only: \t" << counter[8] << "\n";
        outfile << "11) mass 2.2 to 4.5: \t" << counter[9] << "\n";
        outfile << "12) dimuon pt > 0.2: \t" << counter[10] << "\n";
        outfile << "13) SPD match FOhits:\t" << counter[11] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    return;
}

Bool_t EventPassedMyCuts(){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    counter[0]++;

    // 5) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[1]++;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[2]++;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[3]++;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[4]++;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[5]++;

    // 8) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[6]++;

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[7]++;

    // 10) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[8]++;

    // 11) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[9]++;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[10]++;

    // 13) Transverse momentum cut
    if(!(fPt > 0.20 && fPt < 1.00)) return kFALSE;
    counter[11]++;

    // 14) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counter[12]++;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedRoman(){
    // To compare with Roman's numbers

    // 0) Selections applied on the GRID:
    // 0a) fEvent non-empty
    // 0b) At least two tracks associated with the vertex
    // 0c) Distance from the IP lower than 15 cm
    // 0d) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    // 0e) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    counter[0]++;

    // 1) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[1]++;

    // 2) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[2]++;

    // 3) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[3]++;

    // 4a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[4]++;

    // 4b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[5]++;

    // 5a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[6]++;

    // 5b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[7]++;

    // 6) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[8]++;

    // 7) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[9]++;

    // 8) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[10]++;

    // 9) Transverse momentum cut
    if(!(fPt < 0.11)) return kFALSE;
    counter[11]++;

    // 10) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counter[12]++;
    
    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedNoMatching(){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    counter[0]++;

    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[1]++;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[2]++;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[3]++;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[4]++;

    // 7) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[5]++;

    // 8) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[6]++;

    // 9) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[7]++;

    // 10) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[8]++;

    // 11) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[9]++;

    // 12) Transverse momentum cut
    if(!(fPt > 0.20)) return kFALSE;
    counter[10]++;

    // 13) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[11]++;

    // Event passed all the selections =>
    return kTRUE;
}