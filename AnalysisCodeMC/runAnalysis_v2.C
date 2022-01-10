// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskJPsiMC_DG.h"

void runAnalysis_v2(Bool_t Neutral = kFALSE)
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;

    Int_t DatasetMC = 5;
    // DatasetMC == 1   => kCohJpsiToMu
    // DatasetMC == 2   => kIncohJpsiToMu
    // DatasetMC == 3   => kCohPsi2sToMuPi
    // DatasetMC == 4   => kIncohPsi2sToMuPi
    // DatasetMC == 5   => kTwoGammaToMuMedium // inv mass from 1.8 to 15 GeV

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
    sprintf(txt_cmd,"AddTaskJPsiMC_DG.C(%d)", Neutral);
    AliAnalysisTaskJPsiMC_DG *task = reinterpret_cast<AliAnalysisTaskJPsiMC_DG*>(gInterpreter->ExecuteMacro(txt_cmd));
#else
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse();
    gROOT->LoadMacro("AliAnalysisTaskJPsiMC_DG.cxx++g");
    gROOT->LoadMacro("AddTaskJPsiMC_DG.C");
    AliAnalysisTaskJPsiMC_DG *task = AddTaskJPsiMC_DG(Neutral);
#endif

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree"); // !!!
        // add a few files to the chain (change this so that your local files are added)
        //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohJpsiToMu_295585_001/AliESDs.root"); 
        chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohPsi2sToMuPi_295585_001/AliESDs.root"); 
        //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kCohJpsiToMu_295585_001/AliESDs.root"); 
        //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kCohPsi2sToMuPi_295585_001/AliESDs.root"); 
        //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kTwoGammaToMuMedium_295585_001/AliESDs.root"); 

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

        //alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohJpsiToMu");
        //alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohJpsiToMu");
        //alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohPsi2sToMuPi");
        //alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohPsi2sToMuPi");
        //alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kTwoGammaToMuMedium");

        if(DatasetMC == 1 && Neutral == kFALSE){
            // kCohJpsiToMu
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohJpsiToMu");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kCohJpsiToMu");
            alienHandler->SetExecutable(    "ESDs_MC_kCohJpsiToMu.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kCohJpsiToMu.jdl");                
        } else if(DatasetMC == 2 && Neutral == kFALSE){
            // kIncohJpsiToMu
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohJpsiToMu");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kIncohJpsiToMu");
            alienHandler->SetExecutable(    "ESDs_MC_kIncohJpsiToMu.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kIncohJpsiToMu.jdl");                
        } else if(DatasetMC == 3 && Neutral == kFALSE){
            // kCohPsi2sToMuPi, no neutral
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohPsi2sToMuPi");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kCohPsi2sToMuPi");
            alienHandler->SetExecutable(    "ESDs_MC_kCohPsi2sToMuPi.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kCohPsi2sToMuPi.jdl");
        } else if(DatasetMC == 4 && Neutral == kFALSE){
            // kIncohPsi2sToMuPi, no neutral
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohPsi2sToMuPi");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kIncohPsi2sToMuPi");
            alienHandler->SetExecutable(    "ESDs_MC_kIncohPsi2sToMuPi.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kIncohPsi2sToMuPi.jdl");
        } else if(DatasetMC == 5 && Neutral == kFALSE){
            // kTwoGammaToMuMedium
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kTwoGammaToMuMedium");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kTwoGammaToMuMedium");
            alienHandler->SetExecutable(    "ESDs_MC_kTwoGammaToMuMedium.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kTwoGammaToMuMedium.jdl");
        } else if(DatasetMC == 3 && Neutral == kTRUE){
            // kCohPsi2sToMuPi, neutral
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kCohPsi2sToMuPi");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kCohPsi2sToMuPi_Neutral");
            alienHandler->SetExecutable(    "ESDs_MC_kCohPsi2sToMuPi_Neutral.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kCohPsi2sToMuPi_Neutral.jdl");
        } else if(DatasetMC == 4 && Neutral == kTRUE){
            // kIncohPsi2sToMuPi, neutral
            alienHandler->SetGridDataDir("/alice/sim/2019/LHC19k1/kIncohPsi2sToMuPi");
            // working directory
            alienHandler->SetGridWorkingDir("ESDs_MC_kIncohPsi2sToMuPi_Neutral");
            alienHandler->SetExecutable(    "ESDs_MC_kIncohPsi2sToMuPi_Neutral.sh");
            alienHandler->SetJDLName(       "ESDs_MC_kIncohPsi2sToMuPi_Neutral.jdl");
        } else {
            cout << "Not a valid option. Terminating..." << endl;
            return;
        }

        alienHandler->SetDataPattern("/*AliESDs.root");
        // no run prefix for MC!
        // run numbers
        // LHC18q
        alienHandler->AddRunNumber(	295585	); // 1
        alienHandler->AddRunNumber(	295586	); // 2
        alienHandler->AddRunNumber(	295588	); // 3
        alienHandler->AddRunNumber(	295589	); // 4
        alienHandler->AddRunNumber(	295610	); // 5
        alienHandler->AddRunNumber(	295611	); // 6
        alienHandler->AddRunNumber(	295612	); // 7
        alienHandler->AddRunNumber(	295615	); // 8
        alienHandler->AddRunNumber(	295666	); // 9
        alienHandler->AddRunNumber(	295667	); // 10
        alienHandler->AddRunNumber(	295668	); // 11
        alienHandler->AddRunNumber(	295673	); // 12
        alienHandler->AddRunNumber(	295675	); // 13
        alienHandler->AddRunNumber(	295676	); // 14
        alienHandler->AddRunNumber(	295712	); // 15
        alienHandler->AddRunNumber(	295714	); // 16
        alienHandler->AddRunNumber(	295717	); // 17
        alienHandler->AddRunNumber(	295718	); // 18
        alienHandler->AddRunNumber(	295719	); // 19
        alienHandler->AddRunNumber(	295721	); // 20
        alienHandler->AddRunNumber(	295723	); // 21
        alienHandler->AddRunNumber(	295725	); // 22
        alienHandler->AddRunNumber(	295754	); // 23
        alienHandler->AddRunNumber(	295755	); // 24
        alienHandler->AddRunNumber(	295758	); // 25
        alienHandler->AddRunNumber(	295759	); // 26
        alienHandler->AddRunNumber(	295762	); // 27
        alienHandler->AddRunNumber(	295763	); // 28
        alienHandler->AddRunNumber(	295786	); // 29
        alienHandler->AddRunNumber(	295788	); // 30
        alienHandler->AddRunNumber(	295791	); // 31
        alienHandler->AddRunNumber(	295816	); // 32
        alienHandler->AddRunNumber(	295818	); // 33
        alienHandler->AddRunNumber(	295819	); // 34
        alienHandler->AddRunNumber(	295822	); // 35
        alienHandler->AddRunNumber(	295825	); // 36
        alienHandler->AddRunNumber(	295826	); // 37
        alienHandler->AddRunNumber(	295829	); // 38
        alienHandler->AddRunNumber(	295831	); // 39
        alienHandler->AddRunNumber(	295853	); // 40
        alienHandler->AddRunNumber(	295854	); // 41
        alienHandler->AddRunNumber(	295855	); // 42
        alienHandler->AddRunNumber(	295856	); // 43
        alienHandler->AddRunNumber(	295859	); // 44
        alienHandler->AddRunNumber(	295860	); // 45
        alienHandler->AddRunNumber(	295861	); // 46
        alienHandler->AddRunNumber(	295909	); // 47
        alienHandler->AddRunNumber(	295910	); // 48
        alienHandler->AddRunNumber(	295913	); // 49
        alienHandler->AddRunNumber(	295936	); // 50
        alienHandler->AddRunNumber(	295937	); // 51
        alienHandler->AddRunNumber(	295941	); // 52
        alienHandler->AddRunNumber(	295942	); // 53
        alienHandler->AddRunNumber(	296016	); // 54
        alienHandler->AddRunNumber(	296060	); // 55
        alienHandler->AddRunNumber(	296062	); // 56
        alienHandler->AddRunNumber(	296063	); // 57
        alienHandler->AddRunNumber(	296065	); // 58
        alienHandler->AddRunNumber(	296066	); // 59
        alienHandler->AddRunNumber(	296123	); // 60
        alienHandler->AddRunNumber(	296132	); // 61
        alienHandler->AddRunNumber(	296133	); // 62
        alienHandler->AddRunNumber(	296134	); // 63
        alienHandler->AddRunNumber(	296135	); // 64
        alienHandler->AddRunNumber(	296142	); // 65
        alienHandler->AddRunNumber(	296143	); // 66
        alienHandler->AddRunNumber(	296191	); // 67
        alienHandler->AddRunNumber(	296192	); // 68
        alienHandler->AddRunNumber(	296194	); // 69
        alienHandler->AddRunNumber(	296195	); // 70
        alienHandler->AddRunNumber(	296196	); // 71
        alienHandler->AddRunNumber(	296197	); // 72
        alienHandler->AddRunNumber(	296198	); // 73
        alienHandler->AddRunNumber(	296240	); // 74
        alienHandler->AddRunNumber(	296241	); // 75
        alienHandler->AddRunNumber(	296242	); // 76
        alienHandler->AddRunNumber(	296243	); // 77
        alienHandler->AddRunNumber(	296244	); // 78
        alienHandler->AddRunNumber(	296246	); // 79
        alienHandler->AddRunNumber(	296247	); // 80
        alienHandler->AddRunNumber(	296269	); // 81
        alienHandler->AddRunNumber(	296270	); // 82
        alienHandler->AddRunNumber(	296273	); // 83
        alienHandler->AddRunNumber(	296279	); // 84
        alienHandler->AddRunNumber(	296280	); // 85
        alienHandler->AddRunNumber(	296303	); // 86
        alienHandler->AddRunNumber(	296304	); // 87
        alienHandler->AddRunNumber(	296309	); // 88
        alienHandler->AddRunNumber(	296312	); // 89
        alienHandler->AddRunNumber(	296377	); // 90
        alienHandler->AddRunNumber(	296378	); // 91
        alienHandler->AddRunNumber(	296379	); // 92
        alienHandler->AddRunNumber(	296380	); // 93
        alienHandler->AddRunNumber(	296381	); // 94
        alienHandler->AddRunNumber(	296383	); // 95
        alienHandler->AddRunNumber(	296414	); // 96
        alienHandler->AddRunNumber(	296415	); // 97
        alienHandler->AddRunNumber(	296419	); // 98
        alienHandler->AddRunNumber(	296420	); // 99
        alienHandler->AddRunNumber(	296423	); // 100
        alienHandler->AddRunNumber(	296424	); // 101
        alienHandler->AddRunNumber(	296433	); // 102
        alienHandler->AddRunNumber(	296472	); // 103
        alienHandler->AddRunNumber(	296509	); // 104
        alienHandler->AddRunNumber(	296510	); // 105
        alienHandler->AddRunNumber(	296511	); // 106
        alienHandler->AddRunNumber(	296512	); // 107
        alienHandler->AddRunNumber(	296516	); // 108
        alienHandler->AddRunNumber(	296547	); // 109
        alienHandler->AddRunNumber(	296548	); // 110
        alienHandler->AddRunNumber(	296549	); // 111
        alienHandler->AddRunNumber(	296550	); // 112
        alienHandler->AddRunNumber(	296551	); // 113
        alienHandler->AddRunNumber(	296552	); // 114
        alienHandler->AddRunNumber(	296553	); // 115
        alienHandler->AddRunNumber(	296594	); // 116
        alienHandler->AddRunNumber(	296615	); // 117
        alienHandler->AddRunNumber(	296616	); // 118
        alienHandler->AddRunNumber(	296618	); // 119
        alienHandler->AddRunNumber(	296619	); // 120
        alienHandler->AddRunNumber(	296621	); // 121
        alienHandler->AddRunNumber(	296622	); // 122
        alienHandler->AddRunNumber(	296623	); // 123

        // LHC18r
        alienHandler->AddRunNumber(	296690	); // 1
        alienHandler->AddRunNumber(	296691	); // 2
        alienHandler->AddRunNumber(	296693	); // 3
        alienHandler->AddRunNumber(	296694	); // 4
        alienHandler->AddRunNumber(	296749	); // 5
        alienHandler->AddRunNumber(	296750	); // 6
        alienHandler->AddRunNumber(	296781	); // 7
        alienHandler->AddRunNumber(	296784	); // 8
        alienHandler->AddRunNumber(	296785	); // 9
        alienHandler->AddRunNumber(	296786	); // 10
        alienHandler->AddRunNumber(	296787	); // 11
        alienHandler->AddRunNumber(	296790	); // 12
        alienHandler->AddRunNumber(	296793	); // 13
        alienHandler->AddRunNumber(	296794	); // 14
        alienHandler->AddRunNumber(	296799	); // 15
        alienHandler->AddRunNumber(	296835	); // 16
        alienHandler->AddRunNumber(	296836	); // 17
        alienHandler->AddRunNumber(	296838	); // 18
        alienHandler->AddRunNumber(	296839	); // 19
        alienHandler->AddRunNumber(	296848	); // 20
        alienHandler->AddRunNumber(	296849	); // 21
        alienHandler->AddRunNumber(	296850	); // 22
        alienHandler->AddRunNumber(	296851	); // 23
        alienHandler->AddRunNumber(	296852	); // 24
        alienHandler->AddRunNumber(	296890	); // 25
        alienHandler->AddRunNumber(	296894	); // 26
        alienHandler->AddRunNumber(	296899	); // 27
        alienHandler->AddRunNumber(	296900	); // 28
        alienHandler->AddRunNumber(	296903	); // 29
        alienHandler->AddRunNumber(	296930	); // 30
        alienHandler->AddRunNumber(	296931	); // 31
        alienHandler->AddRunNumber(	296932	); // 32
        alienHandler->AddRunNumber(	296934	); // 33
        alienHandler->AddRunNumber(	296935	); // 34
        alienHandler->AddRunNumber(	296938	); // 35
        alienHandler->AddRunNumber(	296941	); // 36
        alienHandler->AddRunNumber(	296966	); // 37
        alienHandler->AddRunNumber(	297029	); // 38
        alienHandler->AddRunNumber(	297031	); // 39
        alienHandler->AddRunNumber(	297035	); // 40
        alienHandler->AddRunNumber(	297085	); // 41
        alienHandler->AddRunNumber(	297117	); // 42
        alienHandler->AddRunNumber(	297118	); // 43
        alienHandler->AddRunNumber(	297119	); // 44
        alienHandler->AddRunNumber(	297123	); // 45
        alienHandler->AddRunNumber(	297124	); // 46
        alienHandler->AddRunNumber(	297128	); // 47
        alienHandler->AddRunNumber(	297129	); // 48
        alienHandler->AddRunNumber(	297132	); // 49
        alienHandler->AddRunNumber(	297133	); // 50
        alienHandler->AddRunNumber(	297193	); // 51
        alienHandler->AddRunNumber(	297194	); // 52
        alienHandler->AddRunNumber(	297195	); // 53
        alienHandler->AddRunNumber(	297196	); // 54
        alienHandler->AddRunNumber(	297218	); // 55
        alienHandler->AddRunNumber(	297219	); // 56
        alienHandler->AddRunNumber(	297221	); // 57
        alienHandler->AddRunNumber(	297222	); // 58
        alienHandler->AddRunNumber(	297278	); // 59
        alienHandler->AddRunNumber(	297310	); // 60
        alienHandler->AddRunNumber(	297311	); // 61
        alienHandler->AddRunNumber(	297317	); // 62
        alienHandler->AddRunNumber(	297332	); // 63
        alienHandler->AddRunNumber(	297333	); // 64
        alienHandler->AddRunNumber(	297335	); // 65
        alienHandler->AddRunNumber(	297336	); // 66
        alienHandler->AddRunNumber(	297363	); // 67
        alienHandler->AddRunNumber(	297366	); // 68
        alienHandler->AddRunNumber(	297367	); // 69
        alienHandler->AddRunNumber(	297372	); // 70
        alienHandler->AddRunNumber(	297379	); // 71
        alienHandler->AddRunNumber(	297380	); // 72
        alienHandler->AddRunNumber(	297405	); // 73
        alienHandler->AddRunNumber(	297406	); // 74
        alienHandler->AddRunNumber(	297413	); // 75
        alienHandler->AddRunNumber(	297414	); // 76
        alienHandler->AddRunNumber(	297415	); // 77
        alienHandler->AddRunNumber(	297441	); // 78
        alienHandler->AddRunNumber(	297442	); // 79
        alienHandler->AddRunNumber(	297446	); // 80
        alienHandler->AddRunNumber(	297450	); // 81
        alienHandler->AddRunNumber(	297451	); // 82
        alienHandler->AddRunNumber(	297452	); // 83
        alienHandler->AddRunNumber(	297479	); // 84
        alienHandler->AddRunNumber(	297481	); // 85
        alienHandler->AddRunNumber(	297483	); // 86
        alienHandler->AddRunNumber(	297512	); // 87
        alienHandler->AddRunNumber(	297537	); // 88
        alienHandler->AddRunNumber(	297540	); // 89
        alienHandler->AddRunNumber(	297541	); // 90
        alienHandler->AddRunNumber(	297542	); // 91
        alienHandler->AddRunNumber(	297544	); // 92
        alienHandler->AddRunNumber(	297558	); // 93
        alienHandler->AddRunNumber(	297588	); // 94
        alienHandler->AddRunNumber(	297590	); // 95
        alienHandler->AddRunNumber(	297595	); // 96

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
