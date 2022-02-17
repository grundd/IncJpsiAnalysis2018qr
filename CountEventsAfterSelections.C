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

Int_t counter[15] = { 0 };

Bool_t cuts_pass1 = kTRUE;
Bool_t cuts_pass3 = kTRUE;
Bool_t cuts_Roman = kTRUE;
Bool_t cuts_noMatching = kTRUE;

Bool_t EventPassed_pass1();
Bool_t EventPassed_pass3();
Bool_t EventPassed_Roman();
Bool_t EventPassed_noMatching();

void CountEventsAfterSelections(){

    // pass1 data
    TFile *f_pass1 = TFile::Open("Trees/AnalysisData_pass1/AnalysisResultsLHC18qrMerged.root", "read");
    if(f_pass1) Printf("Input data loaded.");

    TTree *t_pass1 = dynamic_cast<TTree*> (f_pass1->Get("AnalysisOutput/fTreeJPsi"));
    if(t_pass1) Printf("Input tree loaded.");

    ConnectTreeVariables(t_pass1, kFALSE);

    TList *l_pass1 = dynamic_cast<TList*> (f_pass1->Get("AnalysisOutput/fOutputList"));
    if(l_pass1) Printf("Input list loaded.");

    TH1F *hCounterCuts_pass1 = (TH1F*)l_pass1->FindObject("hCounterCuts");
    if(hCounterCuts_pass1) Printf("Input histogram loaded.");

    //hCounterCuts_pass1->Draw();

    // pass3 data
    TFile *f_pass3 = TFile::Open("Trees/AnalysisData_pass3/AnalysisResults.root", "read");
    if(f_pass3) Printf("Input data loaded.");

    TTree *t_pass3 = dynamic_cast<TTree*> (f_pass3->Get("AnalysisOutput/fTreeJpsi"));
    if(t_pass3) Printf("Input tree loaded.");

    ConnectTreeVariables(t_pass3, kTRUE);

    TList *l_pass3 = dynamic_cast<TList*> (f_pass3->Get("AnalysisOutput/fOutputList"));
    if(l_pass3) Printf("Input list loaded.");

    TH1F *hCounterCuts_pass3 = (TH1F*)l_pass3->FindObject("hCounterCuts");
    if(hCounterCuts_pass3) Printf("Input histogram loaded.");

    //hCounterCuts_pass3->Draw();

    if(cuts_pass1){

        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);
            EventPassed_pass1();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/cuts_pass1.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0) event non-empty:  \t" << hCounterCuts_pass1->GetBinContent(1) << "\n";
        outfile << "1) vertex # contribs:\t" << hCounterCuts_pass1->GetBinContent(2) << "\n";
        outfile << "2) vertex Z distance:\t" << hCounterCuts_pass1->GetBinContent(3) << "\n";
        outfile << "3) two good tracks:  \t" << hCounterCuts_pass1->GetBinContent(4) << "\n";
        outfile << "4) CCUP31 trigger:   \t" << hCounterCuts_pass1->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "5a) ADA offline veto:\t" << counter[1] << "\n";
        outfile << "5b) ADC offline veto:\t" << counter[2] << "\n";
        outfile << "6a) V0A offline veto:\t" << counter[3] << "\n";
        outfile << "6b) V0C offline veto:\t" << counter[4] << "\n";
        outfile << "7) SPD match FOhits: \t" << counter[5] << "\n";
        outfile << "8) muon pairs only:  \t" << counter[6] << "\n";
        outfile << "9) dilept |y| < 0.8: \t" << counter[7] << "\n";
        outfile << "10) trks |eta| < 0.8:\t" << counter[8] << "\n";
        outfile << "11) opposite charges:\t" << counter[9] << "\n";
        outfile << "12) mass 2.2 to 4.5: \t" << counter[10] << "\n";
        outfile << "13) p_T 0.2 to 1.0:  \t" << counter[11] << "\n";
        outfile << "14) mass 3.0 to 3.2: \t" << counter[12] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 15; i++){
        counter[i] = 0;
    }

    if(cuts_pass3){

        Printf("%lli entries found in the tree.", t_pass3->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass3->GetEntries(); iEntry++){
            t_pass3->GetEntry(iEntry);
            EventPassed_pass3();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/cuts_pass3.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0) event non-empty:  \t" << hCounterCuts_pass3->GetBinContent(1) << "\n";
        outfile << "1) two good tracks:  \t" << hCounterCuts_pass3->GetBinContent(2) << "\n";
        outfile << "2) CCUP31 trigger:   \t" << hCounterCuts_pass3->GetBinContent(3) << " (check: " << counter[0] << ")\n";
        outfile << "3) vertex # contribs:\t" << counter[1] << "\n";
        outfile << "4) vertex Z distance:\t" << counter[2] << "\n";
        outfile << "5a) ADA offline veto:\t" << counter[3] << "\n";
        outfile << "5b) ADC offline veto:\t" << counter[4] << "\n";
        outfile << "6a) V0A offline veto:\t" << counter[5] << "\n";
        outfile << "6b) V0C offline veto:\t" << counter[6] << "\n";
        outfile << "7) SPD match FOhits: \t" << counter[7] << "\n";
        outfile << "8) muon pairs only:  \t" << counter[8] << "\n";
        outfile << "9) dilept |y| < 0.8: \t" << counter[9] << "\n";
        outfile << "10) trks |eta| < 0.8:\t" << counter[10] << "\n";
        outfile << "11) opposite charges:\t" << counter[11] << "\n";
        outfile << "12) mass 2.2 to 4.5: \t" << counter[12] << "\n";
        outfile << "13) p_T 0.2 to 1.0:  \t" << counter[13] << "\n";
        outfile << "14) mass 3.0 to 3.2: \t" << counter[14] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 15; i++){
        counter[i] = 0;
    }

    if(cuts_Roman){

        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);
            EventPassed_Roman();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/cuts_Roman.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0a) event non-empty:  \t" << hCounterCuts_pass1->GetBinContent(1) << "\n";
        outfile << "0b) vertex # contribs:\t" << hCounterCuts_pass1->GetBinContent(2) << "\n";
        outfile << "0c) vertex Z distance:\t" << hCounterCuts_pass1->GetBinContent(3) << "\n";
        outfile << "0d) two good tracks:  \t" << hCounterCuts_pass1->GetBinContent(4) << "\n";
        outfile << "0e) CCUP31 trigger:   \t" << hCounterCuts_pass1->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "1) SPD match FOhits:  \t" << counter[1] << "\n";
        outfile << "2) muon pairs only:   \t" << counter[2] << "\n";
        outfile << "3) mass 2.2 to 4.5:   \t" << counter[3] << "\n";
        outfile << "4a) ADA offline veto: \t" << counter[4] << "\n";
        outfile << "4b) ADC offline veto: \t" << counter[5] << "\n";
        outfile << "5a) V0A offline veto: \t" << counter[6] << "\n";
        outfile << "5b) V0C offline veto: \t" << counter[7] << "\n";
        outfile << "6) dilepton |y| < 0.8:\t" << counter[8] << "\n";
        outfile << "7) tracks |eta| < 0.8:\t" << counter[9] << "\n";
        outfile << "8) opposite charges:  \t" << counter[10] << "\n";
        outfile << "9) dimuon pt < 0.11:  \t" << counter[11] << "\n";
        outfile << "10) mass 3.0 to 3.2:  \t" << counter[12] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    // Set all counters to zeros
    for(Int_t i = 0; i < 15; i++){
        counter[i] = 0;
    }

    if(cuts_noMatching){

        Printf("%lli entries found in the tree.", t_pass1->GetEntries());
        Int_t nEntriesAnalysed = 0;

        for(Int_t iEntry = 0; iEntry < t_pass1->GetEntries(); iEntry++){
            t_pass1->GetEntry(iEntry);
            EventPassed_noMatching();

            if((iEntry+1) % 100000 == 0){
                nEntriesAnalysed += 100000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }

        // Print the numbers:
        TString name = "Results/CountEventsAfterSelections/cuts_noMatching.txt";
        ofstream outfile (name.Data());
        outfile << std::fixed << std::setprecision(0); // Set the precision to 0 dec places
        outfile << "0) event non-empty:  \t" << hCounterCuts_pass1->GetBinContent(1) << "\n";
        outfile << "1) vert # contribs:  \t" << hCounterCuts_pass1->GetBinContent(2) << "\n";
        outfile << "2) vert Z distance:  \t" << hCounterCuts_pass1->GetBinContent(3) << "\n";
        outfile << "3) two good tracks:  \t" << hCounterCuts_pass1->GetBinContent(4) << "\n";
        outfile << "4) CCUP31 trigger:   \t" << hCounterCuts_pass1->GetBinContent(5) << " (check: " << counter[0] << ")\n";
        outfile << "5a) ADA offline veto:\t" << counter[1] << "\n";
        outfile << "5b) ADC offline veto:\t" << counter[2] << "\n";
        outfile << "6a) V0A offline veto:\t" << counter[3] << "\n";
        outfile << "6b) V0C offline veto:\t" << counter[4] << "\n";
        outfile << "7) muon pairs only:  \t" << counter[5] << "\n";
        outfile << "8) dilept |y| < 0.8: \t" << counter[6] << "\n";
        outfile << "9) trks |eta| < 0.8: \t" << counter[7] << "\n";
        outfile << "10) opposite charges:\t" << counter[8] << "\n";
        outfile << "11) mass 2.2 to 4.5: \t" << counter[9] << "\n";
        outfile << "12) dimuon pt > 0.2: \t" << counter[10] << "\n";
        outfile << "13) SPD match FOhits:\t" << counter[11] << "\n";

        outfile.close();
        Printf("*** Results printed to %s.***", name.Data());
    }

    return;
}

Bool_t EventPassed_pass1(){

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

    // 7) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[5]++;

    // 8) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[6]++;

    // 9) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[7]++;

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[8]++;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
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

Bool_t EventPassed_pass3(){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) nGoodTracksTPC == 2 && nGoodTracksSPD == 2
    // 2) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    counter[0]++;

    // 3) At least two tracks associated with the vertex
    if(fVertexContrib < 2) return kFALSE;
    counter[1]++;

    // 4) Distance from the IP lower than 15 cm
    if(fVertexZ > 15) return kFALSE;
    counter[2]++;

    // 5a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;
    counter[3]++;

    // 5b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;
    counter[4]++;

    // 6a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;
    counter[5]++;

    // 6b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;
    counter[6]++;

    // 7) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;
    counter[7]++;

    // 8) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
    counter[8]++;

    // 9) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;
    counter[9]++;

    // 10) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;
    counter[10]++;

    // 11) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;
    counter[11]++;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;
    counter[12]++;

    // 13) Transverse momentum cut
    if(!(fPt > 0.20 && fPt < 1.00)) return kFALSE;
    counter[13]++;

    // 14) Invariant mass between 3.0 and 3.2 GeV/c^2
    if(!(fM > 3.0 && fM < 3.2)) return kFALSE;
    counter[14]++;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassed_Roman(){
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

Bool_t EventPassed_noMatching(){

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

    // 7) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;
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