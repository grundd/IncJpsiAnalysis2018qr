// PtFit.c
// David Grund, Sep 23, 2021

// cpp headers
#include <stdio.h> // printf
#include <fstream> // print output to txt file
// root headers
# include "TFile.h"
# include "TTree.h"
// roofit headers
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"

using namespace RooFit;

Int_t nBins = 200;

Double_t fPtLow = 0.0;
Double_t fPtUpp = 2.0;
Double_t fPt;

Double_t BinSizeDouble = (fPtUpp - fPtLow) * 1000 / nBins + 0.5;
Int_t BinSize = (Int_t)BinSizeDouble;

void DoPtFit(Bool_t isBkgData = kFALSE); // isBkgData == kTRUE if sideband data are used for bkg
void DoPtFitNoBkg();
void ConnectTreeVariablesPt(TTree *t);

void PtFit(){

    DoPtFit();

    DoPtFitNoBkg(); // after subtracting the background first

    return;
}

void DoPtFit(Bool_t isBkgData){

    // Load the data file
    TFile *file = TFile::Open("Trees/PtFit/PtFitTrees.root", "read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    // Load trees   
    TTree *tData = dynamic_cast<TTree*> (file->Get("fDataTree")); 
    if(tData) Printf("Tree %s loaded.", tData->GetName());

    TTree *tCohJpsi = dynamic_cast<TTree*> (file->Get("kCohJpsiToMu")); 
    if(tCohJpsi) Printf("Tree %s loaded.", tCohJpsi->GetName());

    TTree *tIncJpsi = dynamic_cast<TTree*> (file->Get("kIncohJpsiToMu")); 
    if(tIncJpsi) Printf("Tree %s loaded.", tIncJpsi->GetName());

    TTree *tCohPsi2s = dynamic_cast<TTree*> (file->Get("kCohPsi2sToMuPi")); 
    if(tCohPsi2s) Printf("Tree %s loaded.", tCohPsi2s->GetName());

    TTree *tIncPsi2s = dynamic_cast<TTree*> (file->Get("kIncohPsi2sToMuPi")); 
    if(tIncPsi2s) Printf("Tree %s loaded.", tIncPsi2s->GetName());

    TTree *tBkg = NULL;
    if(isBkgData == kFALSE) tBkg = dynamic_cast<TTree*> (file->Get("kTwoGammaToMuMedium")); 
    if(isBkgData == kTRUE) return; // (!)
    if(tBkg) Printf("Tree %s loaded.", tBkg->GetName());

    // Cuts
    char fStrReduce[120];
    sprintf(fStrReduce,"fPt > %f && fPt < %f", fPtLow, fPtUpp);

    // Binning
    // (!) is binning necessary or not?
    RooBinning bin(fPtLow, fPtUpp);
    bin.addUniform(nBins, fPtLow, fPtUpp);

    // Definition of roofit variables
    RooRealVar fPt("fPt", "fPt", fPtLow, fPtUpp);
    fPt.setBinning(bin);
    RooArgSet fSetOfVariables(fPt);

    // 1) Create PDF for kCohJpsiToMu
    ConnectTreeVariablesPt(tCohJpsi);
    RooDataSet *DSetCohJ = new RooDataSet("DSetCohJ","DSetCohJ",fSetOfVariables,Import(*tCohJpsi));
    Printf("Tree %s contains %lli entries.", tCohJpsi->GetName(), tCohJpsi->GetEntries());
    Printf("DataSet %s contains %i entries.", tCohJpsi->GetName(), DSetCohJ->numEntries());
    RooAbsData *AbsDCohJ = DSetCohJ->reduce(Cut(fStrReduce),EventRange(0,DSetCohJ->numEntries()));
    RooDataHist DHisCohJ("DHisCohJ","DHisCohJ",fPt,*AbsDCohJ,1);
    RooHistPdf  hPDFCohJ("hPDFCohJ","hPDFCohJ",fSetOfVariables,DHisCohJ,0);
    Printf("1) Template for %s created.", tCohJpsi->GetName());

    // 2) Create PDF for kIncohJpsiToMu
    ConnectTreeVariablesPt(tIncJpsi);
    RooDataSet *DSetIncJ = new RooDataSet("DSetIncJ","DSetIncJ",fSetOfVariables,Import(*tIncJpsi));
    Printf("Tree %s contains %lli entries.", tIncJpsi->GetName(), tIncJpsi->GetEntries());
    Printf("DataSet %s contains %i entries.", tIncJpsi->GetName(), DSetIncJ->numEntries());
    RooAbsData *AbsDIncJ = DSetIncJ->reduce(Cut(fStrReduce),EventRange(0,DSetIncJ->numEntries()));
    RooDataHist DHisIncJ("DHisIncJ","DHisIncJ",fPt,*AbsDIncJ,1);
    RooHistPdf  hPDFIncJ("hPDFIncJ","hPDFIncJ",fSetOfVariables,DHisIncJ,0);
    Printf("2) Template for %s created.", tIncJpsi->GetName());

    // 3) Create PDF for kCohPsi2sToMuPi
    ConnectTreeVariablesPt(tCohPsi2s);
    RooDataSet *DSetCohP = new RooDataSet("DSetCohP","DSetCohP",fSetOfVariables,Import(*tCohPsi2s));
    Printf("Tree %s contains %lli entries.", tCohPsi2s->GetName(), tCohPsi2s->GetEntries());
    Printf("DataSet %s contains %i entries.", tCohPsi2s->GetName(), DSetCohP->numEntries());
    RooAbsData *AbsDCohP = DSetCohP->reduce(Cut(fStrReduce),EventRange(0,DSetCohP->numEntries()));
    RooDataHist DHisCohP("DHisCohP","DHisCohP",fPt,*AbsDCohP,1);
    RooHistPdf  hPDFCohP("hPDFCohP","hPDFCohP",fSetOfVariables,DHisCohP,0);
    Printf("3) Template for %s created.", tCohPsi2s->GetName());

    // 4) Create PDF for kIncohPsi2sToMuPi
    ConnectTreeVariablesPt(tIncPsi2s);
    RooDataSet *DSetIncP = new RooDataSet("DSetIncP","DSetIncP",fSetOfVariables,Import(*tIncPsi2s));
    Printf("Tree %s contains %lli entries.", tIncPsi2s->GetName(), tIncPsi2s->GetEntries());
    Printf("DataSet %s contains %i entries.", tIncPsi2s->GetName(), DSetIncP->numEntries());
    RooAbsData *AbsDIncP = DSetIncP->reduce(Cut(fStrReduce),EventRange(0,DSetIncP->numEntries()));
    RooDataHist DHisIncP("DHisIncP","DHisIncP",fPt,*AbsDIncP,1);
    RooHistPdf  hPDFIncP("hPDFIncP","hPDFIncP",fSetOfVariables,DHisIncP,0);
    Printf("4) Template for %s created.", tIncPsi2s->GetName());

    // 5) Create PDF for kTwoGammaToMuMedium
    ConnectTreeVariablesPt(tBkg);
    RooDataSet *DSetBkg= new RooDataSet("DSetBkg","DSetBkg",fSetOfVariables,Import(*tBkg));
    Printf("Tree %s contains %lli entries.", tBkg->GetName(), tBkg->GetEntries());
    Printf("DataSet %s contains %i entries.", tBkg->GetName(), DSetBkg->numEntries());
    RooAbsData *AbsDBkg = DSetBkg->reduce(Cut(fStrReduce),EventRange(0,DSetBkg->numEntries()));
    RooDataHist DHisBkg("DHisBkg","DHisBkg",fPt,*AbsDBkg,1);
    RooHistPdf  hPDFBkg("hPDFBkg","hPDFBkg",fSetOfVariables,DHisBkg,0);
    Printf("5) Template for %s created.", tBkg->GetName());

    // 6) H1 parametrization of inc with nucleon dissociation
    // (...)

    // 7) Get the dataset
    ConnectTreeVariablesPt(tData);
    RooDataSet *DSetData = new RooDataSet("DSetData","DSetData",fSetOfVariables,Import(*tData));
    Printf("Tree %s contains %lli entries.", tData->GetName(), tData->GetEntries());
    Printf("DataSet %s contains %i entries.", tData->GetName(), DSetData->numEntries());
    RooAbsData *AbsDData = DSetData->reduce(Cut(fStrReduce),EventRange(0,DSetData->numEntries()));
    RooDataHist DHisData("DHisData","DHisData",fPt,*AbsDData,1);
    Printf("7) LHC18qr dataset created.");

    // 8) Create the model for fitting
    // 8.1) Normalizations:
    RooRealVar NCohJ("NCohJ","Number of coh J/psi events", 0.90*AbsDData->numEntries(),0.50*AbsDData->numEntries(),1.0*AbsDData->numEntries());
    RooRealVar NIncJ("NIncJ","Number of inc J/psi events", 0.05*AbsDData->numEntries(),0.01*AbsDData->numEntries(),0.3*AbsDData->numEntries());
    RooRealVar NDiss("NDiss","Number of dis J/psi events", 0.05*AbsDData->numEntries(),0.01*AbsDData->numEntries(),0.3*AbsDData->numEntries());
    // Load the number of bkg events
    Double_t N_bkg_val, N_bkg_err;
    ifstream file_in;
    file_in.open("Results/InvMassFit/all/all_bkg.txt");
    if(!(file_in.fail())){
        while(!file_in.eof()) file_in >> N_bkg_val >> N_bkg_err;
    } else {
        Printf("Results/InvMassFit/all/all_bkg.txt not found. Terminating...");
        return;
    }
    file_in.close();
    RooRealVar NBkg("NBkg","Number of background events",N_bkg_val,N_bkg_val-10,N_bkg_val+10);
    NBkg.setVal(N_bkg_val);
    NBkg.setConstant(kTRUE);
    Printf("Number of background events with 3.0 < m < 3.2 GeV: %.0f pm %.0f", N_bkg_val, N_bkg_err);

    RooGenericPdf NCohP("NCohP","Number of coh FD events","NCohJ*0.0714",RooArgSet(NCohJ));
    RooGenericPdf NIncP("NIncP","Number of inc FD events","NIncJ*0.0692",RooArgSet(NIncJ));
    // calculate the coefficients (!)

    // 8.2) The model:
    RooAddPdf Mod("Mod","Sum of all PDFs",
        RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohJ, hPDFIncP, hPDFBkg), // dissociative part missing (!)
        RooArgList(NCohJ, NIncJ, NCohP, NIncP, NBkg) // dissociative part missing (!)
    );

    // 9) Perform fitting
    RooFitResult* ResFit = Mod.fitTo(*AbsDData,Extended(kTRUE),Save(),Range(fPtLow,fPtUpp));

    // 10) Output to text file
    // (...)

    // 11) Plot the results
    // (...)

    return;

}

void DoPtFitNoBkg(){

    return;
}

void ConnectTreeVariablesPt(TTree *t){

    t->SetBranchAddress("fPt", &fPt);

    Printf("Variables from %s connected.", t->GetName());

    return;
}