// STARlight_PrepareOutputTrees.c
// David Grund, Nov 21, 2021

#include "STARlight_Utilities.h"

void PrepareOutputTrees_ChooseDataset(Int_t iDS);

void STARlight_PrepareOutputTrees(){

    //PrepareOutputTrees_ChooseDataset(1);

    // To find the optimal value of R_A
    TString str = "";
    /*
    str = "Trees/STARlight/OptimalRA/coh_modRA_0_6.624/";
    ConvertStarlightAsciiToTree(6e6, str);
    */
    ///*
    Double_t R_A = 0.;
    for(Int_t i = 11; i < 13; i++){
        R_A = 6.60 + i * 0.10;
        str = Form("Trees/STARlight/OptimalRA/coh_modRA_%i_%.3f/",i+1,R_A);
        ConvertStarlightAsciiToTree(6e6, str);
        //CompareTrees(str);        
    }
    //*/

}

void PrepareOutputTrees_ChooseDataset(Int_t iDS){

    Int_t nGenEv = 0;
    TString str = "";

    if(iDS == 1){
        // |t|-Dependence of incoherent photonuclear cross section.
        nGenEv = 5000000;
        str = "Trees/STARlight/inc_tDep/";
    } else if(iDS == 2){
        // R_A set to 7.53 in src/nucleus.cpp => to study the effect
        // of the change on the shape of coherent template in pT fit
        nGenEv = 4880000;
        str = "Trees/STARlight/coh_4880000_modRA/";
    } else if(iDS == 3){
        // R_A set to 6.62 in src/nucleus.cpp 
        nGenEv = 6000000;
        str = "Trees/STARlight/coh_6000000_stdRA/";
    } else if(iDS == 4){
        // R_A set to 7.53 in src/nucleus.cpp 
        nGenEv = 6000000;
        str = "Trees/STARlight/coh_6000000_modRA/";
    } else if(iDS == 5){
        // Sudakov corrections included in gammagammaleptonpair.cpp
        // (R_A set to 6.62 in src/nucleus.cpp)
        nGenEv = 20000;
        str = "Trees/STARlight/bkg_sudakov/";
    } 

    PrepareTreesPtGammaVMPom(nGenEv, str);
    ConvertStarlightAsciiToTree(nGenEv, str);
    //CompareTrees(str);

    return;
}