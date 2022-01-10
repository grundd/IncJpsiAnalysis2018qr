// header files
// c++ headers
#include <iostream>

// root headers
#include "TGrid.h"
#include "TGridResult.h"
#include "TFileMerger.h"

// main program
void GetFiles(Int_t option)
// get a specific set of files from GRID
// 1 => MC kCohJpsiToMu
// 2 => MC kIncohJpsiToMu
// 3 => MC kCohPsi2sToMuPi charged
// 4 => MC kIncohPsi2sToMuPi charged
// 5 => MC kTwoGammaToMuMedium
// 6 => MC kCohPsi2sToMuPi neutral
// 7 => MC kIncohPsi2sToMuPi neutral
{
    // connect to the GRID
    TGrid::Connect("alien://");

    TGridResult* result = NULL;
    TFileMerger m;

    if(option == 2){
        result = gGrid->Query("/alice/cern.ch/user/d/dgrund/ESDs_MC_kIncohJpsiToMu/myOutputDir","AnalysisResults.root");
        m.OutputFile("AnalysisResults_kIncohJpsiToMu.root");
    }
    if(option == 3){
        result = gGrid->Query("/alice/cern.ch/user/d/dgrund/ESDs_MC_kIncohJpsiToMu_CCUP31correct/myOutputDir","AnalysisResults.root");
        m.OutputFile("AnalysisResults_kIncohJpsiToMu_CCUP31correct.root");
    }
    if(option == 4){
        result = gGrid->Query("/alice/cern.ch/user/d/dgrund/ESDs_MC_kIncohJpsiToMu_CCUP31wrong/myOutputDir","AnalysisResults.root");
        m.OutputFile("AnalysisResults_kIncohJpsiToMu_CCUP31wrong.root");
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

	
	
	
