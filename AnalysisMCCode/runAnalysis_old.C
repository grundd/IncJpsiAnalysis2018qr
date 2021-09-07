// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskJPsiMC_DG.h"

void runAnalysis_old()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;

    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskJPsiMC_DG");
    // ESD handler
    AliESDInputHandler *ESDHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(ESDHandler);
    // MC handler 
    AliMCEventHandler *MCHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(MCHandler);
    if (!MCHandler) Printf("Could not retrieve MC event handler.");

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    gInterpreter->LoadMacro("AliAnalysisTaskJPsiMC_DG.cxx++g");
    char txt_cmd[120];
    sprintf(txt_cmd,"AddTaskJPsiMC_DG.C()");
    AliAnalysisTaskJPsiMC_DG *task = reinterpret_cast<AliAnalysisTaskJPsiMC_DG*>(gInterpreter->ExecuteMacro(txt_cmd));
#else
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse();
    gROOT->LoadMacro("AliAnalysisTaskJPsiMC_DG.cxx++g");
    gROOT->LoadMacro("AddTaskJPsiMC_DG.C");
    AliAnalysisTaskJPsiMC_DG *task = AddTaskJPsiMC_DG();
#endif

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree"); // !!!
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohJpsiToMu_295585_001/AliESDs.root"); 

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
        alienHandler->SetAdditionalLibs("AliAnalysisTaskJPsiMC_DG.cxx AliAnalysisTaskJPsiMC_DG.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskJPsiMC_DG.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20200415_ROOT6-1"); // Guillermo: new version might solve a problem with duplicity of some runs
        // Guillermo: vAN-20200221_ROOT6-1
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data

        alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohJpsiToMu");
        alienHandler->SetDataPattern("/*AliESDs.root");
        // no run prefix for MC!
        //alienHandler->AddRunNumber(295937); 
        alienHandler->AddRunNumber(296623); 

        // working directory
        alienHandler->SetGridWorkingDir("ESDs_MC_07-09-2021_test");
        alienHandler->SetExecutable(    "ESDs_MC_07-09-2021_test.sh");
        alienHandler->SetJDLName(       "ESDs_MC_07-09-2021_test.jdl");

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
        alienHandler->SetMergeViaJDL(kTRUE);
        // define the output folders
        alienHandler->SetGridOutputDir("myOutputDir");
        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test"); // "test" or "terminate" ::
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full"); // "full" or "terminate"
            mgr->StartAnalysis("grid");
        }
    }
}
