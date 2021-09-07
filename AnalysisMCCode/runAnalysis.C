// include the header of your analysis task here
#include "AliAnalysisTaskJPsiMC_DG.h"

void runAnalysis() {
    // header location
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisMyTask");
    // ESD input handler
    AliESDInputHandler *ESDHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(ESDHandler);
    // AOD handler
    //AliAODInputHandler *AODHandler = new AliAODInputHandler();
    //mgr->SetInputEventHandler(AODHandler);
    // MC handler 
    //AliMCEventHandler *MCHandler = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    AliMCEventHandler *MCHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(MCHandler);
    if (!MCHandler) Printf("Could not retrieve MC event handler.");

    // compile the class (locally) with debug symbols
    gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    gInterpreter->LoadMacro("AliAnalysisTaskJPsiMC_DG.cxx++g");
    // load the addtask macro and create the task
    AliAnalysisTaskJPsiMC_DG *task = reinterpret_cast<AliAnalysisTaskJPsiMC_DG*>(gInterpreter->ExecuteMacro("AddTaskJPsiMC_DG.C()"));
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("esdTree");
    //chain->Add("Data/AliESDs_kIncohJpsiToMu_295937_001.root"); 
    //chain->Add("Data/AliESDs_kIncohPsi2sToMuPi_295585_001.root");
    chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohJpsiToMu_295585_001/AliESDs.root"); 
    //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kIncohPsi2sToMuPi_295585_001/AliESDs.root"); 
    //chain->Add("/home/david/alice/IncJpsiAnalysis2018qr/Data/MC_kTwoGammaToMuMedium_295585_001/AliESDs.root"); 
    
    // Print the number of entries
    Printf("Number of events: %lld", chain->GetEntries());
    
    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(0);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    // start the analysis locally
    mgr->StartAnalysis("local", chain);
}