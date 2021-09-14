// RunListCheck.c
// David Grund, Sep 7, 2021
// To check if the run numbers of analysed data match the official list

#include "TFile.h"
#include <fstream> // print output to txt file

#include "TreesManager.h"

Int_t RunList18q[123] = {
    295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 295668, 295673, 295675, 
    295676, 295712, 295714, 295717, 295718, 295719, 295721, 295723, 295725, 295754, 295755, 295758, 295759, 
    295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 
    295853, 295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 295937, 295941, 
    295942, 296016, 296060, 296062, 296063, 296065, 296066, 296123, 296132, 296133, 296134, 296135, 296142, 
    296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 
    296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312, 296377, 296378, 
    296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 
    296510, 296511, 296512, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296594, 296615, 
    296616, 296618, 296619, 296621, 296622, 296623
};

Int_t RunList18r[96] = {
    296690, 296691, 296693, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787, 296790, 296793, 
    296794, 296799, 296835, 296836, 296838, 296839, 296848, 296849, 296850, 296851, 296852, 296890, 296894, 
    296899, 296900, 296903, 296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 297029, 297031, 
    297035, 297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129, 297132, 297133, 297193, 297194, 
    297195, 297196, 297218, 297219, 297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 
    297336, 297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 297415, 297441, 
    297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540, 297541, 297542, 
    297544, 297558, 297588, 297590, 297595
};

Int_t RunList18qr[219] = {
    295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 295668, 295673, 295675, 
    295676, 295712, 295714, 295717, 295718, 295719, 295721, 295723, 295725, 295754, 295755, 295758, 295759, 
    295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 
    295853, 295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 295937, 295941, 
    295942, 296016, 296060, 296062, 296063, 296065, 296066, 296123, 296132, 296133, 296134, 296135, 296142, 
    296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 
    296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312, 296377, 296378, 
    296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 
    296510, 296511, 296512, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296594, 296615, 
    296616, 296618, 296619, 296621, 296622, 296623, 296690, 296691, 296693, 296694, 296749, 296750, 296781, 
    296784, 296785, 296786, 296787, 296790, 296793, 296794, 296799, 296835, 296836, 296838, 296839, 296848, 
    296849, 296850, 296851, 296852, 296890, 296894, 296899, 296900, 296903, 296930, 296931, 296932, 296934, 
    296935, 296938, 296941, 296966, 297029, 297031, 297035, 297085, 297117, 297118, 297119, 297123, 297124, 
    297128, 297129, 297132, 297133, 297193, 297194, 297195, 297196, 297218, 297219, 297221, 297222, 297278, 
    297310, 297311, 297317, 297332, 297333, 297335, 297336, 297363, 297366, 297367, 297372, 297379, 297380, 
    297405, 297406, 297413, 297414, 297415, 297441, 297442, 297446, 297450, 297451, 297452, 297479, 297481, 
    297483, 297512, 297537, 297540, 297541, 297542, 297544, 297558, 297588, 297590, 297595
};

Bool_t FoundRunNumbers[219] = { kFALSE };

void SortRunNumbers();
void PrepareTotalRunList();
void DoRunListCheck();

void RunListCheck(){

    //SortRunNumbers();

    //PrepareTotalRunList();

    DoRunListCheck();

    return;
}

void DoRunListCheck(){
    
    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    Int_t RunNumbersNotFound = 0;
    Int_t n = sizeof(RunList18qr) / sizeof(RunList18qr[0]);

    // Check if the current run number is on the list
    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        Int_t iRun = 0;
        Bool_t RunNumberFound = kFALSE;
        while(iRun < n){
            if(fRunNumber == RunList18qr[iRun]){
                FoundRunNumbers[iRun] = kTRUE;
                RunNumberFound = kTRUE;
            }
            iRun++;
        }
        if(RunNumberFound == kFALSE){
            RunNumbersNotFound++;
            Printf("Problematic run number: %i", fRunNumber);
        } 
        RunNumberFound = kFALSE;

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    if(RunNumbersNotFound == 0){
        Printf("No problematic run numbers found.");
    } else {
        Printf("%i problematic run numbers found.", RunNumbersNotFound);
    }

    // Check if all the run numbers were found
    Bool_t AllNumbersFound = kTRUE;
    for(Int_t j = 0; j < n; j++){
        if(FoundRunNumbers[j] == kFALSE) AllNumbersFound = kFALSE;
    }
    if(AllNumbersFound == kTRUE){
        Printf("All run numbers appear in the tree.");
    } else {
        Printf("Some run number is missing!");
    }

    return;
}

void PrepareTotalRunList(){

    // Fill RunList18qr:
    for(Int_t i = 0; i < 219; i++){
        if(i < 123){
            RunList18qr[i] = RunList18q[i];
        } else {
            RunList18qr[i] = RunList18r[i - 123];
        }
    } 

    // Sort the run numbers in the ascending order:
    Int_t n = sizeof(RunList18qr) / sizeof(RunList18qr[0]);
    std::sort(RunList18qr, RunList18qr + n);

    // Print the results:
    TString name = "Results/RunListCheck/TotalRunList.txt";
    ofstream outfile (name.Data());
    Int_t RunsPerLine = 13;
    for(Int_t i = 0; i < n; i++){
        if((i+1) % RunsPerLine == 0){
            outfile << RunList18qr[i] << ", \n";
        } else if(i == n-1){
            outfile << RunList18qr[i];
        } else {
            outfile << RunList18qr[i] << ", ";
        }
    }
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}

void SortRunNumbers(){

    Int_t nq = sizeof(RunList18q) / sizeof(RunList18q[0]);
    std::sort(RunList18q, RunList18q + nq);

    Int_t nr = sizeof(RunList18r) / sizeof(RunList18r[0]);
    std::sort(RunList18r, RunList18r + nr);

    // Print the results:
    TString name = "Results/RunListCheck/18qrRunLists.txt";
    ofstream outfile (name.Data());
    Int_t RunsPerLine = 13;

    // LHC18q
    outfile << "LHC18q (" << nq << " runs):\n";
    for(Int_t i = 0; i < nq; i++){
        if((i+1) % RunsPerLine == 0){
            outfile << RunList18q[i] << ", \n";
        } else if(i == nq-1){
            outfile << RunList18q[i] << "\n\n";
        } else {
            outfile << RunList18q[i] << ", ";
        }
    }

    // LHC18r
    outfile << "LHC18r (" << nr << " runs):\n"; 
    for(Int_t i = 0; i < nr; i++){
        if((i+1) % RunsPerLine == 0){
            outfile << RunList18r[i] << ", \n";
        } else if(i == nr-1){
            outfile << RunList18r[i] << "\n\n";
        } else {
            outfile << RunList18r[i] << ", ";
        }
    }    

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}