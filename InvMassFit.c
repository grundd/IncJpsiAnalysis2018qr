// InvMassFit.c
// David Grund, Sep 8, 2021
// To perform fit of the invariant mass distribution of measured data

using namespace RooFit;

#include "fTreeJPsiManager.h"

// Main functions
void DoInvMassFitMain(Int_t opt = 0);
// Support functions
void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit);
void SetCanvas(TCanvas *c, Bool_t bLogScale);
void PrepareDataTree();

void InvMassFit(){

    //PrepareDataTree();

    DoInvMassFitMain(2);

    return;
}

void DoInvMassFitMain(Int_t opt = 0){
    // Fit the invariant mass distribution using Double-sided CB function
    // Fix the values of the tail parameters to MC values
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fPtCut;
    Double_t fPtCutLow;
    Double_t fPtCutUpp;
    Double_t fYCut = 0.80;
    Double_t fMCutLow = 2.2;
    Double_t fMCutUpp = 4.5;

    switch(opt){
        case 0: // Incoherent-enriched sample
            fPtCut = 0.20;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 1: // Coherent-enriched sample
            fPtCut = 0.11;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 2: // Total sample (pt < 2.0 GeV/c)
            fPtCut = 2.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 3: // pt bin 1
            break;
        case 4: // pt bin 2
            break;
        case 5: // pt bin 3
            break;
        case 6: // pt bin 4
            break;
    }

    // Binning:
    Int_t nBins = 115; // so that each bin between 2.2 and 4.5 GeV is 20 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);
  
    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    //fM.setBinning(binM);

    // Get the data trees
    TFile *fFileIn = new TFile("Trees/InvMassFit/InvMassFit.root"); 
    TTree *fTreeIn = NULL;
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6){
        fFileIn->GetObject("tIncEnrSample",fTreeIn);
    } else if(opt == 1){
        fFileIn->GetObject("tCohEnrSample",fTreeIn);
    } else if(opt == 2){
        fFileIn->GetObject("tMixedSample",fTreeIn);
    }

    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*fTreeIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("\n");
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    char name[20];
    Double_t values[4];
    Double_t errors[4];

    TString* path = new TString("Results/InvMassFitMC/");
    switch(opt){
        case 0: // Incoherent-enriched sample
            path->Append("inc_doubleCB.txt");
            break;
        case 1: // Coherent-enriched sample
            path->Append("coh_doubleCB.txt");
            break;
        case 2: // Total sample (pt < 2.0 GeV/c)
            path->Append("all_doubleCB.txt");
            break;
        case 3: // pt bin 1
            break;
        case 4: // pt bin 2
            break;
        case 5: // pt bin 3
            break;
        case 6: // pt bin 4
            break;
    }

    ifstream fTxtFileIn;
    fTxtFileIn.open(path->Data());
    if(fTxtFileIn.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!fTxtFileIn.eof()){
            fTxtFileIn >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
        fTxtFileIn.close();
    }
    fAlpha_L = values[0];
    fAlpha_R = values[1];
    fN_L = values[2];
    fN_R = values[3];

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,20.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,20.);
    alpha_L.setConstant(kTRUE);
    alpha_R.setConstant(kTRUE);
    n_L.setConstant(kTRUE);
    n_R.setConstant(kTRUE);

    // Crystal Ball for J/Psi
    RooRealVar mass_Jpsi("mass_Jpsi","J/psi mass",3.097,3.00,3.20); 
    //mass_Jpsi.setConstant(kTRUE);
    RooRealVar sigma_Jpsi("sigma_Jpsi","J/psi resolution",0.08,0.01,0.1);
    RooGenericPdf mean_R("mean_R","J/psi mass","mass_Jpsi",RooArgSet(mass_Jpsi));
    RooGenericPdf sigma_R("sigma_R","J/psi resolution","sigma_Jpsi",RooArgSet(sigma_Jpsi));
    RooRealVar N_Jpsi("N_Jpsi","number of J/psi events",0.4*nEvents,0,nEvents);

    // Background
    RooRealVar lambda("lambda","background exp",-1.2,-10.,0.);
    RooRealVar N_bkg("N_bkg","number of background events",0.6*nEvents,0,nEvents);

    // Functions for fitting
    // J/psi:
    RooCBShape CB_left("CB_left","CB_left",fM,mass_Jpsi,sigma_Jpsi,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));
    // Background:
    RooGenericPdf BkgPdf("BkgPdf","exp(fM*lambda)",RooArgSet(fM,lambda));

    // Create Model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of J/psi events
    Double_t N_Jpsi_out[2];
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *intDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_out[0] = intDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = intDSCB->getVal()*N_Jpsi.getError();

    // ##########################################################
    // Plot the results
    // Draw Correlation Matrix
    TCanvas *cCorrMat = new TCanvas("cCorrMat","cCorrMat",700,600);
    DrawCorrelationMatrix(cCorrMat,fResFit);

    // Draw histogram with fit results
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    SetCanvas(cHist,kFALSE);

    RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(215),LineWidth(3));
    // Vertical axis
    fFrameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    fFrameM->GetYaxis()->SetTitleSize(0.05);
    fFrameM->GetYaxis()->SetTitleOffset(1.1);
    fFrameM->GetYaxis()->SetLabelSize(0.05);
    fFrameM->GetYaxis()->SetLabelOffset(0.01);
    fFrameM->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    fFrameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    fFrameM->GetXaxis()->SetTitleSize(0.05);
    fFrameM->GetXaxis()->SetLabelSize(0.05);
    fFrameM->GetXaxis()->SetDecimals(1);
    if(opt == 1 || opt == 2){
        //cHist->SetRightMargin(0.02);
        fFrameM->GetYaxis()->SetTitleOffset(1.12);
    }
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","data",fResFit->floatParsFinal().getSize());

    // -------------------------------------------------------------------------------- 
    // Legend1
    Double_t x_pos;
    if(opt == 0 || opt == 1 || opt == 2){
        x_pos = 0.12;
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6){
        x_pos = 0.09;
    }
    TLegend *l1 = new TLegend(x_pos,0.76,x_pos+0.2,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    if(opt == 0){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} > %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 1 || opt == 2){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} < %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    }
    // Text size    
    if(opt == 0 || opt == 1 || opt == 2){
        l1->SetTextSize(0.05);
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6){
        l1->SetTextSize(0.042);
    }
    
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    TLegend *lTitle = new TLegend(0.325,0.88,0.95,0.935);
    lTitle->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    lTitle->SetTextSize(0.05);
    lTitle->SetBorderSize(0);
    lTitle->SetFillStyle(0);
    lTitle->Draw();

    // Legend2
    TLegend *l2 = new TLegend(0.465,0.29,0.95,0.87);
    //l2->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l2->AddEntry("DSCBAndBkgPdf","sum","L");
    //l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l2->AddEntry("DoubleSidedCB","J/#psi signal","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",N_Jpsi_out[0],N_Jpsi_out[1]),"");
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6){
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    } else if(opt == 1 || opt == 2){
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.4f #pm %.4f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.4f #pm %.4f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    }
    l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.3f", alpha_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.3f", (-1)*(alpha_R.getVal())),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f", n_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f", n_R.getVal()),"");

    l2->AddEntry("BkgPdf","background","L");
    l2->AddEntry((TObject*)0,Form("#lambda = %.3f #pm %.3f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
    l2->SetTextSize(0.041);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // Legend3
    Double_t y_pos;
    if(opt == 0 || opt == 1 || opt == 2){
        y_pos = 0.66;
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6){
        y_pos = 0.65;
    }
    TLegend *l3 = new TLegend(0.12,y_pos,0.3,y_pos + 0.1);
    l3->AddEntry((TObject*)0,"#bf{This thesis}","");
    l3->SetTextSize(0.05);
    l3->SetBorderSize(0);
    l3->SetFillStyle(0); //legend is transparent
    l3->Draw();

    // Print results
    switch(opt){
        case 0:
            cHist->Print("Results/InvMassFit/inc/inc_DoubleCB_noPsi.pdf"); 
            cHist->Print("Results/InvMassFit/inc/inc_DoubleCB_noPsi.png"); 
            cCorrMat->Print("Results/InvMassFit/inc/inc_DoubleCB_noPsi_CorrM.pdf");
            cCorrMat->Print("Results/InvMassFit/inc/inc_DoubleCB_noPsi_CorrM.png");
            break;
        case 1:
            cHist->Print("Results/InvMassFit/coh/coh_DoubleCB_noPsi.pdf"); 
            cHist->Print("Results/InvMassFit/coh/coh_DoubleCB_noPsi.png"); 
            cCorrMat->Print("Results/InvMassFit/coh/coh_DoubleCB_noPsi_CorrM.pdf");
            cCorrMat->Print("Results/InvMassFit/coh/coh_DoubleCB_noPsi_CorrM.png");
            break;
        case 2:
            cHist->Print("Results/InvMassFit/all/all_DoubleCB_noPsi.pdf"); 
            cHist->Print("Results/InvMassFit/all/all_DoubleCB_noPsi.png"); 
            cCorrMat->Print("Results/InvMassFit/all/all_DoubleCB_noPsi_CorrM.pdf");
            cCorrMat->Print("Results/InvMassFit/all/all_DoubleCB_noPsi_CorrM.png");
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            break;
        case 6:
            break;
    }

    // Number of background events in the mass range 3.0 to 3.2 GeV/c^2
    // For the transverse momentum fit
    if(opt == 0 || opt == 2){
        Double_t N_bkg_out[2];
        RooAbsReal *iBkg = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));

        N_bkg_out[0] = iBkg->getVal()*N_bkg.getVal();
        N_bkg_out[1] = iBkg->getVal()*N_bkg.getError();

        Printf("\n*** Number of background events for pt fit (if opt == 2): ***");
        Printf("*** N_bkg = %f pm %f \n",N_bkg_out[0],N_bkg_out[1]);
        Printf("*** Finished.");
    }

    return;
}

void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit){

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    cCorrMat->SetTopMargin(0.03);
    cCorrMat->SetBottomMargin(0.12);
    cCorrMat->SetRightMargin(0.19);
    cCorrMat->SetLeftMargin(0.13);

    TH2* hCorr = fResFit->correlationHist();
    hCorr->SetMarkerSize(3.);

    hCorr->GetXaxis()->SetBinLabel(1,"#it{M}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{bkg}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetXaxis()->SetBinLabel(4,"#lambda");
    hCorr->GetXaxis()->SetBinLabel(5,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(1,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(2,"#lambda");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{J/#psi}");
    hCorr->GetYaxis()->SetBinLabel(4,"#it{N}_{bkg}");
    hCorr->GetYaxis()->SetBinLabel(5,"#it{M}_{J/#psi}");

    hCorr->GetXaxis()->SetLabelSize(0.1);
    hCorr->GetXaxis()->SetLabelOffset(0.012);
    hCorr->GetYaxis()->SetLabelSize(0.1);
    hCorr->GetZaxis()->SetLabelSize(0.07);
    // https://root-forum.cern.ch/t/colz-color-palette-font-and-size/15263
    hCorr->Draw("colz,text");

    return;
}

void SetCanvas(TCanvas *c, Bool_t bLogScale){

    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}

void PrepareDataTree(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    // Create new data tree with applied cuts
    TFile fFileOut("Trees/InvMassFit/InvMassFit.root","RECREATE");

    TTree *tIncEnrSample = new TTree("tIncEnrSample", "tIncEnrSample");
    tIncEnrSample->Branch("fPt", &fPt, "fPt/D");
    tIncEnrSample->Branch("fM", &fM, "fM/D");
    tIncEnrSample->Branch("fY", &fY, "fY/D");

    TTree *tCohEnrSample = new TTree("tCohEnrSample", "tCohEnrSample");
    tCohEnrSample->Branch("fPt", &fPt, "fPt/D");
    tCohEnrSample->Branch("fM", &fM, "fM/D");
    tCohEnrSample->Branch("fY", &fY, "fY/D");

    TTree *tMixedSample = new TTree("tMixedSample", "tMixedSample");
    tMixedSample->Branch("fPt", &fPt, "fPt/D");
    tMixedSample->Branch("fM", &fM, "fM/D");
    tMixedSample->Branch("fY", &fY, "fY/D");

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(0,0)) tIncEnrSample->Fill();
        if(EventPassed(0,1)) tCohEnrSample->Fill();
        if(EventPassed(0,2)) tMixedSample->Fill();
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}