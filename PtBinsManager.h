// PtBinsManager.h
// David Grund, Sep 11, 2021

const Int_t nPtBins = 4;
Int_t method = 3;

Double_t *ptBoundaries = NULL;
Double_t ptBoundaries1[5] = {0.200, 0.331, 0.455, 0.608, 1.000};
Double_t ptBoundaries2[5] = {0.200, 0.288, 0.412, 0.621, 1.000};
Double_t ptBoundaries3_4[5] = {0.200, 0.280, 0.380, 0.580, 1.000};
Double_t ptBoundaries3_5[6] = { 0 };
Double_t ptBoundaries4[5] = {0.200, 0.286, 0.402, 0.582, 1.000};

void SetPtBinning(){
    switch(method){
        case 1:
            // PtBinning "Method 1": CalcBinning() in tSlope.c
            if(nPtBins == 4) ptBoundaries = ptBoundaries1;
            break;
        case 2:
            // PtBinning "Method 2": CalcBinning2() in tSlope.c
            if(nPtBins == 4) ptBoundaries = ptBoundaries2;
            break;
        case 3:
            // PtBinning "Method 3": BinsThroughMassFit.c
            if(nPtBins == 4) ptBoundaries = ptBoundaries3_4;
            if(nPtBins == 5) ptBoundaries = ptBoundaries3_5;
            break;
        case 4:
            // PtBinning "Method 4": tSlopeWithoutBkg.c
            if(nPtBins == 4) ptBoundaries = ptBoundaries4;
            break;
    }
}


