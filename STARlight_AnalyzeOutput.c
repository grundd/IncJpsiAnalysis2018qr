// STARlight_AnalyzeOutput.c
// David Grund, Nov 21, 2021

#include "STARlight_Utilities.h"

Int_t nGenEv = 0;
TString str = "";

Int_t iProd = 4;

void STARlight_AnalyzeOutput(){

    if(iProd == 1){
        // |t|-Dependence of incoherent photonuclear cross section.
        nGenEv = 5000000;
        str = "Trees/STARlight/inc_tDep/";
    } else if(iProd == 2){
        // R_A set to 7.53 in src/nucleus.cpp => to study the effect
        // of the change on the shape of coherent template in pT fit
        nGenEv = 4880000;
        str = "Trees/STARlight/coh_modRA/";
    } else if(iProd == 3){
        // R_A set to 6.62 in src/nucleus.cpp 
        nGenEv = 6000000;
        str = "Trees/STARlight/coh_stdRA/";
    } else if(iProd == 4){
        // R_A set to 7.53 in src/nucleus.cpp 
        nGenEv = 6000000;
        str = "Trees/STARlight/coh_modRA_2/";
    }

    PrepareTreesPtGammaVMPom(nGenEv, str);
    ConvertStarlightAsciiToTree(str);
    //CompareTrees(str);

    return;
}