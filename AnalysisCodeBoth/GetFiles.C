// header files
// c++ headers
#include <iostream>

// root headers
#include "TGrid.h"
#include "TGridResult.h"
#include "TFileMerger.h"

// main program
void GetFiles(Int_t iMCDataset)
// get a specific set of files from GRID
// 0 => kIncohJpsiToMu
// 1 => kIncohPsi2sToMuPi
// 2 => kCohJpsiToMu
// 3 => kCohPsi2sToMuPi
// 4 => kTwoGammaToMuMedium
{
    // connect to the GRID
    TGrid::Connect("alien://");

    TGridResult* result = NULL;
    TFileMerger m;

    if(iMCDataset == 0){
        result = gGrid->Query("/alice/cern.ch/user/d/dgrund/trial_MC00/myOutputDir","AnalysisResults.root");
        m.OutputFile("AnalysisResults_MC00.root");
    } else if(iMCDataset == 4){
        result = gGrid->Query("/alice/cern.ch/user/d/dgrund/trial_MC04/myOutputDir","AnalysisResults.root");
        m.OutputFile("AnalysisResults_MC04.root");
    }

    Int_t i = 0;
    // Loop over the TGridResult entries and add them to the TFileMerger
    while(result->GetKey(i,"turl")) {
    cout << endl << endl;
    cout << " adding " << result->GetKey(i,"turl") << endl;
    m.AddFile(result->GetKey(i,"turl"));
    i++;
    }

    cout << endl << endl;
    cout << endl << endl;
    cout << endl << endl;

    cout << " ********** MERGING ************ " << endl << endl << endl;
    // Merge
    m.Merge();

    cout << " ********** DONE ************ " << endl << endl << endl;
}

	
	
	
