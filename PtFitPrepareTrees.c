// PtFitPrepareTrees.c
// David Grund, Sep 23, 2021
// Applied selections to MC events (defined in AnalysisManager.h):
// 3.0 < m < 3.2 GeV
// 0.0 < pt < 2.0 GeV

// root headers
#include "TFile.h"
// my headers
#include "AnalysisManager.h"

TString DatasetsMC[5] = {"kCohJpsiToMu","kIncohJpsiToMu","kCohPsi2sToMuPi","kIncohPsi2sToMuPi","kTwoGammaToMuMedium"};

void PrepareTreeMC(Int_t iMC, TTree *t);
// iMC == 0 => kCohJpsiToMu
// iMC == 1 => kIncohJpsiToMu
// iMC == 2 => kCohPsi2sToMuPi
// iMC == 3 => kIncohPsi2sToMuPi
// iMC == 4 => kTwoGammaToMuMedium

void PtFitPrepareTrees(){

    // ***************************************************************
    // 1) Go over real data
    TFile *file = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *fTreeIn = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    Printf("Tree %s has %lli entries.", fTreeIn->GetName(), fTreeIn->GetEntries());

    // Create the output file
    TFile fTreesPt("Trees/PtFit/PtFitTrees.root","RECREATE");

    // Define output tree to store pt of events that pass criteria
    TTree *TreeData = new TTree("fDataTree", "fDataTree");
    TreeData->Branch("fPt", &fPt, "fPt/D");

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        // iMassCut == 1 => 3.0 < m < 3.2 GeV
        // iPtCut == 2 => pt < 2.0 GeV
        if(EventPassed(1, 2)){
            nEvPassed++;
            TreeData->Fill();
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);

    // ***************************************************************
    // 2) Go over MC data
    // Define output tree to store pt of events that pass criteria
    TTree *TreeMC[5] = { NULL };
    TreeMC[0] = new TTree(DatasetsMC[0].Data(), DatasetsMC[0].Data());
    TreeMC[1] = new TTree(DatasetsMC[1].Data(), DatasetsMC[1].Data());
    TreeMC[2] = new TTree(DatasetsMC[2].Data(), DatasetsMC[2].Data());
    TreeMC[3] = new TTree(DatasetsMC[3].Data(), DatasetsMC[3].Data());
    TreeMC[4] = new TTree(DatasetsMC[4].Data(), DatasetsMC[4].Data());

    for(Int_t i = 0; i < 5; i++){
        // Create pt branch in each tree
        TreeMC[i]->Branch("fPt", &fPt, "fPt/D");

        PrepareTreeMC(i, TreeMC[i]);
    }    
    // ***************************************************************
    // 3) Create tree with the data from sidebands
    // (...)
    // ***************************************************************

    // Save results to the output file
    fTreesPt.Write("",TObject::kWriteDelete);

    return;
}

void PrepareTreeMC(Int_t iMC, TTree *t){

    // Load the data
    TFile *file = NULL;
    switch(iMC){
        case 0:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            break;
        case 1:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;
        case 2:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
        case 4:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kTwoGammaToMuMedium.root", "read");
            break;
    }
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);

        // iMassCut == 1 => 3.0 < m < 3.2 GeV
        // iPtCut == 2 => pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)){
            nEvPassed++;
            t->Fill();
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);

    return;
}