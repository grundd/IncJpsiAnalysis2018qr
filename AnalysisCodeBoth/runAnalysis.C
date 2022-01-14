// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskCentralJpsi_DG.h"

Bool_t MC = kFALSE;
Int_t iMCDataset = 3;
// 0 => kIncohJpsiToMu
// 1 => kIncohPsi2sToMuPi
// 2 => kCohJpsiToMu
// 3 => kCohPsi2sToMuPi
// 4 => kTwoGammaToMuMedium

void runAnalysis()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kTRUE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kTRUE;

    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskCentralJpsi_DG");
    // ESD input handler
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    // MC handler 
    AliMCEventHandler *MCHandler = NULL;
    if(MC){
        MCHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(MCHandler);
        if (!MCHandler) Printf("Could not retrieve MC event handler.");
    }

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    gInterpreter->LoadMacro("AliAnalysisTaskCentralJpsi_DG.cxx++g");
    char txt_cmd[120];
    sprintf(txt_cmd,"AddTaskCentralJpsi_DG.C()");
    AliAnalysisTaskCentralJpsi_DG *task = reinterpret_cast<AliAnalysisTaskCentralJpsi_DG*>(gInterpreter->ExecuteMacro(txt_cmd));
#else
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse();
    gROOT->LoadMacro("AliAnalysisTaskCentralJpsi_DG.cxx++g");
    gROOT->LoadMacro("AddTaskCentralJpsi_DG.C");
    AliAnalysisTaskCentralJpsi_DG *task = AddTaskCentralJpsi_DG();
#endif

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree"); // !!!
        // add a few files to the chain (change this so that your local files are added)
        if(!MC) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/Data_000296623/101_AliESDs.root");
        else{
            if(iMCDataset == 0) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohJpsiToMu_295585_001/AliESDs.root"); 
            if(iMCDataset == 1) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohPsi2sToMuPi_295585_001/AliESDs.root"); 
            if(iMCDataset == 2) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kCohJpsiToMu_295585_001/AliESDs.root"); 
            if(iMCDataset == 3) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kCohPsi2sToMuPi_295585_001/AliESDs.root"); 
            if(iMCDataset == 4) chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kTwoGammaToMuMedium_295585_001/AliESDs.root"); 
        }
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);

    // ================================================================
    // ======================== Running on GRID =======================
    // ================================================================

    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskCentralJpsi_DG.cxx AliAnalysisTaskCentralJpsi_DG.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskCentralJpsi_DG.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20200415_ROOT6-1"); // Guillermo: new version might solve a problem with duplicity of some runs
        // Guillermo: vAN-20200221_ROOT6-1
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        
        // select the input files
        // Data
        if(!MC){
            alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");   
            alienHandler->SetDataPattern("/*pass1/*AliESDs.root");
            alienHandler->SetRunPrefix("000");
            //alienHandler->AddRunNumber(295585); 
            alienHandler->AddRunNumber(295937); 
            //alienHandler->AddRunNumber(296623); 
            // working directory
            alienHandler->SetGridWorkingDir("AnalysisTrial_data");
            alienHandler->SetExecutable("AnalysisTrial_data.sh");
            alienHandler->SetJDLName("AnalysisTrial_data.jdl");
        }
        // MC
        if(MC){
            if(iMCDataset == 0) alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohJpsiToMu");   
            if(iMCDataset == 1) alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohPsi2sToMuPi");   
            if(iMCDataset == 2) alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohJpsiToMu");   
            if(iMCDataset == 3) alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohPsi2sToMuPi");  
            if(iMCDataset == 4) alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kTwoGammaToMuMedium");  
            alienHandler->SetDataPattern("/*AliESDs.root");
            // no run prefix for MC!
            alienHandler->AddRunNumber(295585); 
            //alienHandler->AddRunNumber(296062);
            // working directory
            alienHandler->SetGridWorkingDir(Form("trial_MC0%i", iMCDataset));
            alienHandler->SetExecutable(Form("trial_MC0%i.sh", iMCDataset));
            alienHandler->SetJDLName(Form("trial_MC0%i.jdl", iMCDataset));
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE); // set to kFALSE to download merged results
        // define the output folders
        alienHandler->SetGridOutputDir("myOutputDir");
        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test"); // "test" or "terminate"
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full"); // "full" or "terminate"
            mgr->StartAnalysis("grid");
        }
    }
    /*
    if(local && !MC){
        gSystem->Exec("mkdir -p local_data");
        gSystem->Exec("mv AnalysisResults.root track_cuts.root local_data");
    }
    if(local && MC){
        TString folder = TString::Format("local_MC0%i", iMCDataset);
        gSystem->Exec("mkdir -p " + folder);
        gSystem->Exec("mv AnalysisResults.root track_cuts.root " + folder);
    }
    */
    return;
}
