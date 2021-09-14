// MakePtBins.c
// David Grund, Sep 14, 2021
// To calculate bin edges to investigate pt dependence of AxE

// cpp headers
#include <vector>
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TMath.h"

void MakePtBins(){

    const Int_t nBinTypes = 3;
    // in pt, GeV/c
    Double_t bin_widths[nBinTypes] = {0.04, 0.08, 0.15}; 
    Double_t bin_max[nBinTypes] = {0.60, 1.00, 1.60}; 

    vector<Double_t> edges;

    Double_t edge = 0;
    Int_t curr_bin_type = 0;

    // numeric problems here without "- 0.0001"... 
    while(edge < bin_max[nBinTypes-1] - 0.0001){
        edge += bin_widths[curr_bin_type];
        edges.push_back(edge);
        Printf("bin type %i, edge %.2f", curr_bin_type, edge);
        if(edge > bin_max[curr_bin_type] - 0.0001) curr_bin_type++;
    }

    TString name = "PtBinEdges.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(0); 
    outfile << "Double_t edges[" << edges.size() << "] = {";
    outfile << std::fixed << std::setprecision(2); 
    for(Int_t i = 0; i < (edges.size()-1); i++){
        outfile << edges[i] << ", ";
    }
    outfile << edges[edges.size()-1] << "};";

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}