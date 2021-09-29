// PtFit.c
// David Grund, Sep 23, 2021

// cpp headers
#include <stdio.h> // printf
#include <fstream> // print output to txt file
// root headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"
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
#include "RooPlot.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

Double_t fPtLow = 0.0;
Double_t fPtUpp = 2.0;

void DoPtFit(Bool_t isBkgData); // isBkgData == kTRUE if sideband data are used for bkg
void ConnectTreeVariablesPt(TTree *t);
void SetStyle();
void DrawCorrelationMatrix(TCanvas *cCM, RooFitResult* ResFit);

void PtFit(){

    DoPtFit(kFALSE);

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

    // Definition of roofit variables
    RooRealVar fPt("fPt", "fPt", fPtLow, fPtUpp);
    //fPt.setBinning(bin);
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
    // https://arxiv.org/pdf/1304.5162.pdf, p.16, HE
    RooRealVar b_pd("b_pd","b_pd",1.79, 1.70, 1.90);
    RooRealVar n_pd("n_pd","n_pd",3.58, 3.50, 3.70);
    b_pd.setConstant(kTRUE);
    n_pd.setConstant(kTRUE); 
    RooGenericPdf PDFDiss("PDFDiss","fPt*pow((1+pow(fPt,2)*b_pd/n_pd),-n_pd)",RooArgSet(fPt, b_pd, n_pd));
    Printf("5) PDF for diss part created.");

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

    // Load the values of fD coefficients
    Double_t fDCohCh, fDCohChErr, fDIncCh, fDIncChErr, fDCohNe, fDCohNeErr, fDIncNe, fDIncNeErr;
    char ch[8];
    file_in.open("Results/FeedDown/FeedDown_PtFit.txt");
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> fDCohCh >> fDCohChErr >> fDIncCh >> fDIncChErr >> fDCohNe >> fDCohNeErr >> fDIncNe >> fDIncNeErr;
            i++;   
        } 
        file_in.close();
        Printf("fD coefficients loaded.");
    } else {
        Printf("fD coefficients missing. Terminating...");
        return;
    }

    Double_t fDCoh = (fDCohCh + fDCohNe) / 100;
    Double_t fDInc = (fDIncCh + fDIncNe) / 100;
    Printf("********************");
    Printf("fD_coh = %.4f", fDCoh);
    Printf("fD_inc = %.4f", fDInc);
    Printf("********************");
    
    RooGenericPdf NCohP("NCohP","Number of coh FD events",Form("NCohJ*%.4f", fDCoh),RooArgSet(NCohJ));
    RooGenericPdf NIncP("NIncP","Number of inc FD events",Form("NIncJ*%.4f", fDInc),RooArgSet(NIncJ));

    // 8.2) The model:
    RooAddPdf Mod("Mod","Sum of all PDFs",
        RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, PDFDiss, hPDFBkg),
        RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss, NBkg)
    );

    // 9) Perform fitting
    RooFitResult* ResFit = Mod.fitTo(*AbsDData,Extended(kTRUE),Save(),Range(fPtLow,fPtUpp));

    // 10) Output to text file
    // (...)

    // 11) Plot the results
    SetStyle();

    // 11.1) Draw the Correlation Matrix
    TCanvas *cCM = new TCanvas("cCM","cCM",600,500);
    DrawCorrelationMatrix(cCM,ResFit);

    // 11.2) Draw the pt fit
    TCanvas *cPt = new TCanvas("cPt","cPt",900,600);
    cPt->SetLogy();
    cPt->SetTopMargin(0.05);
    cPt->SetBottomMargin(0.12);
    cPt->SetRightMargin(0.02);  

    // 11.3) Set the binning (for plotting only, fit is unbinned)
    RooBinning bin(fPtLow, fPtUpp);
    Int_t nBins = 200;
    bin.addUniform(nBins,fPtLow,fPtUpp);

    RooPlot* PtFrame = fPt.frame(Title("Pt fit"));
    AbsDData->plotOn(PtFrame,Name("AbsDData"),Binning(bin), MarkerStyle(20), MarkerSize(1.));
    Mod.plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222), LineStyle(1),LineWidth(3));
    Mod.plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),LineStyle(1),LineWidth(3));
    Mod.plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222), LineStyle(7),LineWidth(3));
    Mod.plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),LineStyle(7),LineWidth(3));
    Mod.plotOn(PtFrame,Name("PDFDiss"), Components(PDFDiss),  LineColor(15),  LineStyle(1),LineWidth(3));
    Mod.plotOn(PtFrame,Name("Mod"),LineColor(215),LineWidth(3));
    Mod.plotOn(PtFrame,Name("hPDFBkg"), Components(hPDFBkg),LineColor(kBlack),LineStyle(1),LineWidth(3));

    PtFrame->SetAxisRange(0,2,"X");
    // Set X axis
    PtFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    PtFrame->GetXaxis()->SetTitleSize(0.05);
    PtFrame->GetXaxis()->SetLabelSize(0.05);
    // Set Y axis
    PtFrame->GetYaxis()->SetTitle(Form("Counts per %.1f MeV/#it{c}", bin.binWidth(1)*1000));
    PtFrame->GetYaxis()->SetTitleSize(0.05);
    PtFrame->GetYaxis()->SetTitleOffset(0.95);
    PtFrame->GetYaxis()->SetLabelSize(0.05);
    PtFrame->GetYaxis()->SetLabelOffset(0.01);
    PtFrame->Draw("][");

    return;

}

void ConnectTreeVariablesPt(TTree *t){

    t->SetBranchAddress("fPt", &fPt);

    Printf("Variables from %s connected.", t->GetName());

    return;
}

void SetStyle(){

    gStyle->SetOptTitle(0); // suppress title
    gStyle->SetOptStat(0);  // the type of information printed in the histogram statistics box
                            // 0 = no information
    gStyle->SetPalette(1);  // set color map
    gStyle->SetPaintTextFormat("4.2f"); // precision if plotted with "TEXT"

    return;
}

void DrawCorrelationMatrix(TCanvas *cCM, RooFitResult* ResFit){

    // Set margins
    cCM->SetTopMargin(0.03);
    cCM->SetBottomMargin(0.11);
    cCM->SetRightMargin(0.17);
    cCM->SetLeftMargin(0.15);
    // Get 2D corr hist
    TH2* hCorr = ResFit->correlationHist();
    // Set X axis
    hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{coh}");
    // Set Y axis
    hCorr->GetYaxis()->SetBinLabel(1,"#it{N}_{inc}");
    hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{coh}");
    // Set corr hist and draw it
    hCorr->SetMarkerSize(3.6);
    hCorr->GetXaxis()->SetLabelSize(0.13);
    hCorr->GetYaxis()->SetLabelSize(0.13);
    hCorr->GetZaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}