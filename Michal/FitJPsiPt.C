static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

using namespace RooFit;

void FitJPsiPt(Int_t iPoint = -1,Int_t channel = -1,Int_t ZDCcut = -1,Float_t minY = -0.8, Float_t maxY = 0.8,Float_t yieldGG = 1070,Float_t minM = 2.9,Float_t maxM = 3.2) {
//void FitJPsiPt(Int_t iPoint = -1,Int_t channel = 1,Int_t ZDCcut = -1,Float_t minY = -0.8, Float_t maxY = 0.8,Float_t yieldGG = 1061,Float_t minM = 3.0,Float_t maxM = 3.2) {

if(iPoint != -1){
	TFile *yieldGG_file;
	if(channel ==  1)yieldGG_file = TFile::Open("Yield_JPsiMu_Pt.root", "read");
	if(channel == -1)yieldGG_file = TFile::Open("Yield_JPsiEl_Pt.root", "read");
	TH1F *nGGhist;
	if(ZDCcut == -1) nGGhist = (TH1F*)yieldGG_file->Get("nGGhistRapStore");
	else nGGhist = (TH1F*)yieldGG_file->Get("nGGhistZDCStore");
	yieldGG = nGGhist->GetBinContent(iPoint+1);
	} 

TFile *fData_file = TFile::Open("AnalysisResults.4.root", "read");
TFile *fTemplates_file;
if(channel == -1)fTemplates_file = TFile::Open("./MC/templatesEl.29.02.root", "read");
if(channel ==  1)fTemplates_file = TFile::Open("./MC/templatesMu.root", "read");

//Int_t nBins = 200;
Float_t maxPt = 2.0;

Int_t nBins = 43;
Float_t bins[44];
for(Int_t i = 0; i<20; i++)bins[i] = 0.01*i;
for(Int_t i = 20; i<30; i++)bins[i] = 0.2+(0.05*(i-20));
for(Int_t i = 30; i<44; i++)bins[i] = 0.7+(0.1*(i-30));


TH1D *hData_Pt = new TH1D("hData_Pt","p_{T} of Data J/#psi",nBins,bins);
hData_Pt->SetLineWidth(2);
hData_Pt->SetMarkerStyle(kFullCircle);
hData_Pt->SetMarkerColor(kBlack);
hData_Pt->SetLineColor(kBlack);
hData_Pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
hData_Pt->GetXaxis()->SetTitleSize(0.045);
hData_Pt->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}",(maxPt*1000/nBins)));
hData_Pt->GetYaxis()->SetTitleOffset(1.4);
hData_Pt->GetYaxis()->SetTitleSize(0.045);
hData_Pt->SetTitle("TFractionFitter");
hData_Pt->SetStats(kFALSE);

TH1D *hGammaGammaDT_Pt = new TH1D("hGammaGammaDT_Pt"," ",nBins,0,maxPt);
hGammaGammaDT_Pt->SetLineWidth(3);
hGammaGammaDT_Pt->SetLineStyle(kDashed);
hGammaGammaDT_Pt->SetLineColor(2);
hGammaGammaDT_Pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
hGammaGammaDT_Pt->GetXaxis()->SetTitleSize(0.045);
hGammaGammaDT_Pt->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}",(maxPt*1000/nBins)));
hGammaGammaDT_Pt->GetYaxis()->SetTitleOffset(1.4);
hGammaGammaDT_Pt->GetYaxis()->SetTitleSize(0.045);
hGammaGammaDT_Pt->SetStats(kFALSE);

Float_t fPt, fY, fM, fDiLeptonM, fDiLeptonPt, fZNAenergy, fZNCenergy, fPIDsigma, fZNAtime[4], fZNCtime[4],fPIDsigma;
Int_t fChannel, fSign, fRunNumber, fADAdecision, fADCdecision,fV0Adecision, fV0Cdecision, fNGoodTracksITS, fNGoodTracksLoose, fNGoodTracksDCA;
Bool_t fInEtaRec, fTriggers[10];

Int_t nEntries = 0;

TTree *fTreeDataPt = new TTree("fTreeDataPt", "fTreeDataPt");
fTreeDataPt->Branch("rPt", &fPt, "rPt/F");

TList *fData_list = (TList*) fData_file->Get("Upc/ListHist");
fData_Tree  = dynamic_cast<TTree*> (fData_list->FindObject("fTreeJPsi"));
fData_Tree->SetBranchAddress("fY", &fY);
fData_Tree->SetBranchAddress("fPt", &fPt);
fData_Tree->SetBranchAddress("fM", &fM);
fData_Tree->SetBranchAddress("fChannel", &fChannel);
fData_Tree->SetBranchAddress("fSign", &fSign);
fData_Tree->SetBranchAddress("fZNCenergy", &fZNCenergy);
fData_Tree->SetBranchAddress("fZNAenergy", &fZNAenergy);
fData_Tree->SetBranchAddress("fZNCtime", &fZNCtime);
fData_Tree->SetBranchAddress("fZNAtime", &fZNAtime);
fData_Tree->SetBranchAddress("fRunNumber", &fRunNumber);
fData_Tree->SetBranchAddress("fInEtaRec", &fInEtaRec);
fData_Tree->SetBranchAddress("fTriggers", &fTriggers);
fData_Tree->SetBranchAddress("fPIDsigma", &fPIDsigma);

fData_Tree->SetBranchAddress("fADAdecision", &fADAdecision);
fData_Tree->SetBranchAddress("fADCdecision", &fADCdecision);
fData_Tree->SetBranchAddress("fV0Adecision", &fV0Adecision);
fData_Tree->SetBranchAddress("fV0Cdecision", &fV0Cdecision);
fData_Tree->SetBranchAddress("fNGoodTracksITS", &fNGoodTracksITS);
fData_Tree->SetBranchAddress("fNGoodTracksLoose", &fNGoodTracksLoose);
fData_Tree->SetBranchAddress("fNGoodTracksDCA", &fNGoodTracksDCA);

Bool_t ZNA0n, ZNC0n;

nEntries = fData_Tree->GetEntries();
for(Int_t i=0; i<nEntries; i++){
	fData_Tree->GetEntry(i);
	if(!fTriggers[5] && !fTriggers[6])continue;
	if(!fInEtaRec)continue;
	if(!fTriggers[9])continue;
	
	if(fChannel != channel) continue;
	
	if(fSign != -1) continue;
	if(TMath::Abs(fY)>maxY || TMath::Abs(fY)<minY) continue;
	
	//if(fNGoodTracksLoose>2)continue;
	//if(fNGoodTracksDCA<2)continue;
	//if(fNGoodTracksITS>0)continue;
	if(fADAdecision != 0)continue;
	if(fADCdecision != 0)continue;
	if(fV0Adecision != 0)continue;
	if(fV0Cdecision != 0)continue;
	
	ZNA0n = kFALSE;
	ZNC0n = kFALSE;
	
	Bool_t isTDCAfired = 0;
        Bool_t isTDCCfired = 0;
        for (Int_t k=0;k<4;k++) isTDCAfired|=TMath::Abs(fZNAtime[k])<5;
        for (Int_t k=0;k<4;k++) isTDCCfired|=TMath::Abs(fZNCtime[k])<5;
	
	if(!isTDCAfired)ZNA0n = kTRUE;
	if(!isTDCCfired)ZNC0n = kTRUE;
	
	if(ZDCcut == 0 && !(ZNA0n && ZNC0n))continue;
	if(ZDCcut == 1 && !(!ZNA0n && ZNC0n) && !(ZNA0n && !ZNC0n))continue;
	if(ZDCcut == 2 && !(!ZNA0n && !ZNC0n))continue;
	
	//if((fM > 2.7 && fM < 2.9)) hGammaGammaDT_Pt->Fill(fPt);
	if((fM > 3.2 && fM < 3.5)) hGammaGammaDT_Pt->Fill(fPt);
	
	if(fM > maxM || fM < minM) continue;
	hData_Pt->Fill(fPt);
	fTreeDataPt->Fill();	
	}

TH1D *hGammaGamma_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hGammaGamma_Pt_%d",iPoint));
TH1D *hCoherent_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hCoherent_Pt_%d",iPoint));
TH1D *hIncoherent_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hIncoherent_Pt_%d",iPoint));
TH1D *hCoherentFD_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hCoherentFD_Pt_%d",iPoint));
TH1D *hIncoherentFD_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hIncoherentFD_Pt_%d",iPoint));
TH1D *hDissociative_Pt = (TH1D*)fTemplates_file->Get(TString::Format("hDissociative_Pt_%d",iPoint));

TH1D *hCoherent_M = (TH1D*)fTemplates_file->Get(TString::Format("hCoherent_M_%d",iPoint));
hCoherent_M->Scale(1.0/hCoherent_M->Integral());
TH1D *hIncoherent_M = (TH1D*)fTemplates_file->Get(TString::Format("hIncoherent_M_%d",iPoint));
hIncoherent_M->Scale(1.0/hIncoherent_M->Integral());

Float_t correctionfI = (hIncoherent_M->Integral()/hCoherent_M->Integral())		      /(hIncoherent_M->Integral(hIncoherent_M->FindBin(minM),hIncoherent_M->FindBin(maxM))/hCoherent_M->Integral(hCoherent_M->FindBin(minM),hCoherent_M->FindBin(maxM)));

Float_t correctionfC = (hCoherent_M->Integral()/hIncoherent_M->Integral())		      /(hCoherent_M->Integral(hCoherent_M->FindBin(minM),hCoherent_M->FindBin(maxM))/hIncoherent_M->Integral(hIncoherent_M->FindBin(minM),hIncoherent_M->FindBin(maxM)));

//=============================RooFit=================================
Int_t order = 0; // order of the interpolation between bins. Zero is no interpolation

RooRealVar rPt("rPt","#it{p}_{T} (GeV/#it{c})",0,maxPt);
rPt.setRange("coherent_cut",0.0,0.2);
rPt.setRange("incoherent_cut",0.2,2.0);
RooDataSet rData_UnBin("rData_UnBin","rData_UnBin",RooArgSet(rPt),Import(*fTreeDataPt));

RooHistPdf rGammaGamma_Pdf("rGammaGamma_Pdf","rGammaGamma_Pdf",rPt,RooDataHist("rGammaGamma_DH","rGammaGamma_DH",RooArgList(rPt),hGammaGamma_Pt),order);
RooHistPdf rCoherent_Pdf("rCoherent_Pdf","rCoherent_Pdf",rPt,RooDataHist("rCoherent_DH","rCoherent_DH",RooArgList(rPt),hCoherent_Pt),order);
RooHistPdf rIncoherent_Pdf("rIncoherent_Pdf","rIncoherent_Pdf",rPt,RooDataHist("rIncoherent_DH","rIncoherent_DH",RooArgList(rPt),hIncoherent_Pt),order);
RooHistPdf rCoherentFD_Pdf("rCoherentFD_Pdf","rCoherentFD_Pdf",rPt,RooDataHist("rCoherentFD_DH","rCoherentFD_DH",RooArgList(rPt),hCoherentFD_Pt),order);
RooHistPdf rIncoherentFD_Pdf("rIncoherentFD_Pdf","rIncoherentFD_Pdf",rPt,RooDataHist("rIncoherentFD_DH","rIncoherentFD_DH",RooArgList(rPt),hIncoherentFD_Pt),order);
RooHistPdf rDissociative_Pdf("rDissociative_Pdf","rDissociative_Pdf",rPt,RooDataHist("rDissociative_DH","rDissociative_DH",RooArgList(rPt),hDissociative_Pt),order);

//pT*exp(-b*pT^2)
RooRealVar rCohFnc_B("rCohFnc_B","B",200,100,500);
//RooGenericPdf rCoherent_Pdf("rCoherent_Pdf","rCoherent_Pdf","rPt*exp(-rCohFnc_B*pow(rPt,2))",RooArgSet(rPt,rCohFnc_B));
/*/
RooRealVar rDissParB("rDissParB","rDissParB",1.79, 1.79, 1.79);
RooRealVar rDissParN("rDissParN","rDissParN",3.58, 3.58,3.58);
rDissParB.setConstant(kTRUE);
rDissParN.setConstant(kTRUE);

RooGenericPdf rDissociative_Pdf("rDissociative_Pdf","rPt*pow((1+pow(rPt,2)*rDissParB/rDissParN),-rDissParN)",RooArgSet(rPt, rDissParB, rDissParN));
/*/

RooRealVar rGammaGamma_Frc("rGammaGamma_Frc","Number of gg",100,0,10000); 
RooRealVar rCoherent_Frc("rCoherent_Frc","Number of coherent",100,0,10000);
RooRealVar rIncoherent_Frc("rIncoherent_Frc","Number of incoherent",100,0,10000);
RooRealVar rDissociative_Frc("rDissociative_Frc","Number of dissociative",100,0,10000);

if(channel ==  1){
	RooGenericPdf rCoherentFD_Frc("rCoherentFD_Frc","Number of coherent FD","rCoherent_Frc*0.070",RooArgSet(rCoherent_Frc));
	RooGenericPdf rIncoherentFD_Frc("rIncoherentFD_Frc","Number of incoherent FD","rIncoherent_Frc*0.074",RooArgSet(rIncoherent_Frc));
	}
if(channel == -1){
	RooGenericPdf rCoherentFD_Frc("rCoherentFD_Frc","Number of coherent FD","rCoherent_Frc*0.073",RooArgSet(rCoherent_Frc));
	RooGenericPdf rIncoherentFD_Frc("rIncoherentFD_Frc","Number of incoherent FD","rIncoherent_Frc*0.078",RooArgSet(rIncoherent_Frc));
	}

// Create the model as the sum of templates

RooAddPdf Model("Model","Extended sum templates",RooArgList(rGammaGamma_Pdf,rCoherent_Pdf,rIncoherent_Pdf,rCoherentFD_Pdf,rIncoherentFD_Pdf,rDissociative_Pdf),
						 RooArgList(RooConst(yieldGG),rCoherent_Frc,rIncoherent_Frc,rCoherentFD_Frc,rIncoherentFD_Frc,rDissociative_Frc));
						 //RooArgList(rGammaGamma_Frc,rCoherent_Frc,rIncoherent_Frc,rCoherentFD_Frc,rIncoherentFD_Frc,rDissociative_Frc));
// perform fit
RooFitResult* r = Model.fitTo(rData_UnBin,Extended(kTRUE),Save());

    RooBinning ptBins(0, 2);
    ptBins.addUniform(20, 0, 0.2);
    ptBins.addUniform(10, 0.2, 0.7);
    ptBins.addUniform(13, 0.7, 2.0);


   const char* ZDCcutnames[4] = {"0n0n","0nXn","XnXn"};
// plot
  TCanvas *cRooFit = new TCanvas("cRooFit","cRooFit",1600,1100);
  cRooFit->Divide(2,1);
  cRooFit->cd(1);

  RooPlot* frame = rPt.frame(Binning(ptBins),Title(" "));
  frame->GetYaxis()->SetTitleOffset(1.4) ;
  frame->GetYaxis()->SetTitle("J/#psi candidates (counts/1 MeV/#it{c})");
  rData_UnBin.plotOn(frame,Name("rData_UnBin"),MarkerStyle(20),MarkerSize(1.),Binning(ptBins));
  Model.plotOn(frame,Name("Model"),LineColor(kBlack),LineWidth(iPoint == -1 ? 2:1),Binning(ptBins));
  
  // Model.paramOn(frame,data);
  Model.plotOn(frame,Components(rGammaGamma_Pdf), LineColor(kGreen+2),LineWidth(2),Binning(ptBins));
  Model.plotOn(frame,Components(rCoherent_Pdf), LineColor(kBlue),LineWidth(2),Binning(ptBins));
  Model.plotOn(frame,Components(rIncoherent_Pdf), LineColor(kRed),LineWidth(2),Binning(ptBins));
  Model.plotOn(frame,Components(rCoherentFD_Pdf), LineColor(kOrange-6),LineWidth(2),Binning(ptBins));
  Model.plotOn(frame,Components(rIncoherentFD_Pdf), LineColor(kYellow+1),LineWidth(2),Binning(ptBins));
  Model.plotOn(frame,Components(rDissociative_Pdf), LineColor(kMagenta),LineWidth(2),Binning(ptBins));
  frame->GetXaxis()->SetRangeUser(0,0.7);
  
  TPad *myPad1 = new TPad("myPad1", "The pad",0,0.3,1,1);
  myPad1->Draw();
  myPad1->cd();
  gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.0); gPad->SetRightMargin(0.02);  
  frame->Draw();
  
  TLatex *system = new TLatex(0.17,0.844,"Pb-Pb  #sqrt{s_{NN}} = 5.02 TeV");
  system->SetNDC();
  system->SetTextFont(42);
  system->SetTextSize(0.04);
  system->Draw();

  TLegend *myLegend1 = new TLegend(0.37,0.29,0.84,0.72);
  myLegendSetUp(myLegend1,0.04,1);
  
  myLegend1->AddEntry((TObject*)0,"UPC, #it{L}_{int} #approx 229.4 #mub^{-1}","");
  if(channel == 1)myLegend1->AddEntry((TObject*)0,TString::Format("%.1f < M_{ #mu^{+} #mu^{-}} (GeV/#it{c}^{2}) < %.1f",minM,maxM),"");
  if(channel == -1)myLegend1->AddEntry((TObject*)0,TString::Format("%.1f < M_{ e^{+} e^{-}} (GeV/#it{c}^{2}) < %.1f",minM,maxM),"");
  if(ZDCcut != -1)myLegend1->AddEntry((TObject*)0,TString::Format("Neutron emission: %s",ZDCcutnames[ZDCcut]),"");
  myLegend1->AddEntry((TObject*)0,TString::Format("%1.2f < y < %1.2f",minY,maxY),"");
  myLegend1->AddEntry(hData_Pt,"ALICE data","p");
  myLegend1->AddEntry(hCoherent_Pt," Coherent J/#psi","l");
  myLegend1->AddEntry(hIncoherent_Pt," Incoherent J/#psi","l");
  myLegend1->AddEntry(hCoherentFD_Pt,"J/#psi from Coherent #psi(2S) decay","l");
  myLegend1->AddEntry(hIncoherentFD_Pt,"J/#psi from Incoherent #psi(2S) decay","l");
  myLegend1->AddEntry(hDissociative_Pt,"Incoherent J/#psi with dissociation","l");
  if(channel == -1)myLegend1->AddEntry(hGammaGamma_Pt,"Continuum #gamma #gamma to e^{+} e^{-}","l");
  if(channel ==  1)myLegend1->AddEntry(hGammaGamma_Pt,"Continuum #gamma #gamma to #mu^{+} #mu^{-}","l");
  myLegend1->Draw();
  
  TLegend *myLegend2 = new TLegend(0.26,0.51,0.74,0.65);
  myLegendSetUp(myLegend2,0.04,1);
  if(channel == 1)myLegend2->AddEntry((TObject*)0,TString::Format("%.1f < M_{ #mu^{+} #mu^{-}} (GeV/#it{c}^{2}) < %.1f",minM,maxM),"");
  if(channel == -1)myLegend2->AddEntry((TObject*)0,TString::Format("%.1f < M_{ e^{+} e^{-}} (GeV/#it{c}^{2}) < %.1f",minM,maxM),"");
  //myLegend2->Draw();
  
  
  cRooFit->cd(2);
  //frame->GetYaxis()->SetRangeUser(0.9,1.2*hData_Pt->GetMaximum());
  gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.066);gPad->SetRightMargin(0.02);gPad->SetTopMargin(0.066);
  gPad->SetLogy();
  
  TH1D *histData = (TH1D*)rData_UnBin.createHistogram("histData",rPt,Binning(ptBins));
  histData->Sumw2();
  histData->Scale(1.0/histData->Integral());
  
  RooPlot* frameLong = rPt.frame(Title(" "));
  frameLong->GetYaxis()->SetTitleOffset(1.4) ;
  frameLong->GetYaxis()->SetTitle("dN/d#it{p}_{T} (GeV/#it{c})^{-1}");
  frameLong->GetYaxis()->SetTitleOffset(1.4);
  frameLong->GetXaxis()->SetTitleOffset(1.4);
  frameLong->GetYaxis()->SetTitleSize(0.04);
  rData_UnBin.plotOn(frameLong,MarkerStyle(20),MarkerSize(0.8),DrawOption("PZ"),Binning(ptBins),Rescale(1.0/rData_UnBin.numEntries()));
  Model.plotOn(frameLong,Name("Model"),LineColor(kBlack),LineWidth(2), Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  
  // Model.paramOn(frameLong,data);
  Model.plotOn(frameLong,Components(rGammaGamma_Pdf), LineColor(kGreen+2),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  Model.plotOn(frameLong,Components(rCoherent_Pdf), LineColor(kBlue),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  Model.plotOn(frameLong,Components(rIncoherent_Pdf), LineColor(kRed),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  Model.plotOn(frameLong,Components(rCoherentFD_Pdf), LineColor(kOrange-6),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  Model.plotOn(frameLong,Components(rIncoherentFD_Pdf), LineColor(kYellow+1),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  Model.plotOn(frameLong,Components(rDissociative_Pdf), LineColor(kMagenta),LineWidth(2),Binning(ptBins),Normalization(1,RooAbsReal::NumEvent));
  frameLong->GetXaxis()->SetRangeUser(0,2.0);
  frameLong->GetYaxis()->SetRangeUser(0.0002,0.7);
  frameLong->Draw();
  //myLegend1->Draw();
  //myLegend2->Draw();
  
  
  RooCurve *curveModel = (RooCurve*)frameLong->findObject("Model");
  for(Int_t i = 0;i<histData->GetNbinsX();i++)if(curveModel->Eval(histData->GetBinCenter(i+1))!=0){
    histData->SetBinContent(i+1,histData->GetBinContent(i+1)/curveModel->Eval(histData->GetBinCenter(i+1))/histData->GetBinWidth(i+1)*2.0/43.0);
    histData->SetBinError(i+1,histData->GetBinError(i+1)/curveModel->Eval(histData->GetBinCenter(i+1)));
    }
  
  histData->SetStats(kFALSE);
  histData->SetLineColor(kBlack);
  histData->SetLineWidth(2);
  histData->SetTitle(" ");
  histData->GetYaxis()->SetRangeUser(0,2);
  histData->GetXaxis()->SetRangeUser(0,0.7);
  histData->GetYaxis()->SetTitle("Data/Model");
  
  cRooFit->cd(1);
  TPad *myPad2 = new TPad("myPad2", "Pad for Data/MC",0,0,1,0.3);
  myPad2->Draw();
  myPad2->cd();
  gPad->SetTopMargin(0.0); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.2); gPad->SetRightMargin(0.02);
  histData->GetYaxis()->SetTitleSize(0.08);
  histData->GetYaxis()->SetTitleOffset(0.6); 
  histData->GetYaxis()->SetLabelSize(0.08);
  histData->GetXaxis()->SetTitleSize(0.08); 
  histData->GetXaxis()->SetLabelSize(0.08);
  
  RooPlot* framePull = rPt.frame(Title(" "));
  RooHist* hpull = frame->pullHist("rData_UnBin","Model");
  //hpull->plotOn(framePull,MarkerStyle(20),MarkerSize(0.8),DrawOption("PZ"));
  framePull->GetYaxis()->SetTitle("(Data - Model) / #sigma_{Data}");
  framePull->addPlotable(hpull,"P") ;
  framePull->GetXaxis()->SetRangeUser(0,0.7);
  framePull->GetYaxis()->SetTitleSize(0.08);
  framePull->GetYaxis()->SetTitleOffset(0.6); 
  framePull->GetYaxis()->SetLabelSize(0.08);
  framePull->GetXaxis()->SetTitleSize(0.08); 
  framePull->GetXaxis()->SetLabelSize(0.08);
  framePull->Draw();
  //histData->Draw("E0");
  cRooFit->cd(2);
  
  RooAbsReal* rDissociative_Prob;
  RooAbsReal* rIncoherent_Prob;
  RooAbsReal* rCoherent_Prob;

  rDissociative_Prob = rDissociative_Pdf.createIntegral(rPt,NormSet(rPt),Range("coherent_cut"));
  rIncoherent_Prob = rIncoherent_Pdf.createIntegral(rPt,NormSet(rPt),Range("coherent_cut"));
  rCoherent_Prob = rCoherent_Pdf.createIntegral(rPt,NormSet(rPt),Range("coherent_cut"));

  Float_t result_fI = (rIncoherent_Prob->getVal()*rIncoherent_Frc.getVal() +  rDissociative_Prob->getVal()*rDissociative_Frc.getVal())/(rCoherent_Prob->getVal()*rCoherent_Frc.getVal()); 
  Float_t error_fI = TMath::Sqrt(TMath::Power((rIncoherent_Prob->getVal()*rIncoherent_Frc.getError())/(rCoherent_Prob->getVal()*rCoherent_Frc.getVal()),2) + 
		      TMath::Power((rIncoherent_Prob->getVal()*rIncoherent_Frc.getVal()*rCoherent_Frc.getError())/(rCoherent_Prob->getVal()*rCoherent_Frc.getVal()*rCoherent_Frc.getVal()),2));  


  cout<<"RooFit f_I ="<<result_fI<<" +- "<<error_fI<<endl;
  cout<<"Correction = "<<correctionfI<<endl;
  cout<<"Corrected f_I ="<<result_fI*correctionfI<<" +- "<<error_fI*correctionfI<<endl;

  rDissociative_Prob = rDissociative_Pdf.createIntegral(rPt,NormSet(rPt),Range("incoherent_cut"));
  rIncoherent_Prob = rIncoherent_Pdf.createIntegral(rPt,NormSet(rPt),Range("incoherent_cut"));
  rCoherent_Prob = rCoherent_Pdf.createIntegral(rPt,NormSet(rPt),Range("incoherent_cut")); 

  Float_t result_fC = (rCoherent_Prob->getVal()*rCoherent_Frc.getVal())/(rIncoherent_Prob->getVal()*rIncoherent_Frc.getVal() + rDissociative_Prob->getVal()*rDissociative_Frc.getVal());
  Float_t error_fC = TMath::Sqrt(TMath::Power((rCoherent_Prob->getVal()*rCoherent_Frc.getError())/(rIncoherent_Prob->getVal()*rIncoherent_Frc.getVal()),2) + 
		      TMath::Power((rCoherent_Prob->getVal()*rCoherent_Frc.getVal()*rIncoherent_Frc.getError())/(rIncoherent_Prob->getVal()*rIncoherent_Frc.getVal()*rIncoherent_Frc.getVal()),2)); 

  cout<<"RooFit f_C ="<<result_fC<<endl;
  
  TLatex *resultFI = new TLatex(0.65,0.844,TString::Format("f_{I} = %1.3f #pm %1.3f",result_fI*correctionfI,error_fI*correctionfI));
  resultFI->SetNDC();
  resultFI->SetTextFont(42);
  resultFI->SetTextSize(0.04);
  resultFI->Draw();
  
  
  TLatex *resultChi = new TLatex(0.65,0.79,TString::Format("#chi^{2} = %5.3f",frame->chiSquare("Model","rData_UnBin", 3)));
  resultChi->SetNDC();
  resultChi->SetTextFont(42);
  resultChi->SetTextSize(0.04);
  resultChi->Draw();

  
  TLatex *resultFC = new TLatex(0.57,0.744,TString::Format("f_{C} = %1.3f #pm %1.3f",result_fC*correctionfC,error_fC*correctionfC));
  resultFC->SetNDC();
  resultFC->SetTextFont(42);
  resultFC->SetTextSize(0.04);
  //if(iPoint != -1)resultFC->Draw();
  
  
  TCanvas *cRooFitSingle = new TCanvas("cRooFitSingle","cRooFitSingle",1000,800);
  cRooFitSingle->cd();
  TPad *myPad3 = new TPad("myPad3", "The pad",0,0.3,1,1);
  myPad3->Draw();
  myPad3->cd();
  gPad->SetBottomMargin(0.0);gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02);
  frame->Draw();
  myLegend1->Draw();
  resultFI->Draw();
  resultChi->Draw();
  system->Draw();
  cRooFitSingle->cd();
  TPad *myPad4 = new TPad("myPad4", "Pad for Data/MC",0,0,1,0.3);
  myPad4->Draw();
  myPad4->cd();
  gPad->SetTopMargin(0.0);gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.2); gPad->SetRightMargin(0.02);
  //histData->Draw();
  framePull->Draw();
  
  TCanvas *cRooFitPaper = new TCanvas("cRooFitPaper","cRooFitPaper",285*3,750);
  cRooFitPaper->cd();
  gPad->SetBottomMargin(0.11);gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.01);gPad->SetTopMargin(0.01);
  frameLong->Draw();
    
  Float_t lumi = 233.0;
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42);
  latex->SetTextAlign(22);
  latex->DrawLatex(0.56,0.95,Form("ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
  latex->SetTextSize(0.033);
  latex->SetTextAlign(12);
  latex->DrawLatex(0.48,0.89,Form("UPC, L_{#lower[-0.3]{int}} = %.0f #pm %.0f #mub^{-1}",lumi,lumi*0.027));
  if(channel == 1)latex->DrawLatex(0.48,0.84,Form("%.2f < #it{m}_{#mu#mu} < %.2f GeV/#it{c}^{2}",minM,maxM));
  if(channel == -1)latex->DrawLatex(0.48,0.84,Form("%.2f < #it{m}_{ee} < %.2f GeV/#it{c}^{2}",minM,maxM));
  latex->DrawLatex(0.30,0.84,Form("|#it{y}| < 0.8"));

  TLegend* myLegend8 = new TLegend(0.47,0.50,0.63,0.80);
  myLegendSetUp(myLegend8,0.033,1);
  myLegend8->AddEntry(hData_Pt,"ALICE data","p");
  myLegend8->AddEntry(hCoherent_Pt,"Coherent J/#psi");
  myLegend8->AddEntry(hIncoherent_Pt,"Incoherent J/#psi");
  myLegend8->AddEntry(hDissociative_Pt,"Incoherent J/#psi with nucleon dissociation");
  myLegend8->AddEntry(hCoherentFD_Pt,"Coherent J/#psi from #psi' decay");
  myLegend8->AddEntry(hIncoherentFD_Pt,"Incoherent J/#psi from #psi' decay");  
  myLegend8->AddEntry(hGammaGamma_Pt,"Continuum #gamma#gamma #rightarrow ll","l");
  myLegend8->AddEntry(histData,Form("Fit: #chi^{2}/#it{dof}=%.2f",frame->chiSquare("Model","rData_UnBin", 3)));
  myLegend8->Draw();
  gPad->SetLogy();
  
  
  TLatex *resultFunction = new TLatex(0.47,0.7,TString::Format("#splitline{#frac{dN^{coh}}{d#it{p}_{T}} = #it{p}_{T}.e^{-B.#it{p}_{T}^{2}}}{B = %3.1f #pm %3.1f}",rCohFnc_B.getVal(),rCohFnc_B.getError()));
  resultFunction->SetNDC();
  resultFunction->SetTextFont(42);
  resultFunction->SetTextSize(0.04);
  //resultFunction->Draw();
  
  if(iPoint == -1){
  	if(channel == -1){
		cRooFit->SaveAs("./NoteFigures/PtFit_El.pdf");
		cRooFit->SaveAs("./NoteFigures/PtFit_El.png");
		cRooFitPaper->SaveAs("./NoteFigures/PtFit_ElPP.pdf");
		cRooFitSingle->SaveAs("./NoteFigures/PtFitShort_El.png");
		}
  	if(channel ==  1){
		cRooFit->SaveAs("./NoteFigures/PtFit_Mu.pdf");
		cRooFit->SaveAs("./NoteFigures/PtFit_Mu.png");
		cRooFitPaper->SaveAs("./NoteFigures/PtFit_MuPP.pdf");
		cRooFitSingle->SaveAs("./NoteFigures/PtFitShort_Mu.png");
		//cRooFitLog->SaveAs("./NoteFigures/PtFitLog_Mu.png");
		//cRooFit->SaveAs("./NoteFigures/PtFitAlt_Mu.pdf");
		//cRooFit->SaveAs("./NoteFigures/PtFitAlt_Mu.png");
		}
	}



if(iPoint != -1){
	TFile *fractionFile =0x0;
	if(channel ==  1) fractionFile = new TFile("Fraction_JPsiMu.root","UPDATE");
	if(channel == -1) fractionFile = new TFile("Fraction_JPsiEl.root","UPDATE");

	if(ZDCcut == -1){
		TH1F *fractionRap = 0x0;
		fractionRap = (TH1F*)fractionFile->Get("fractionRapStore");

		if(!fractionRap){ 
			cout<<"Creating new histo"<<endl;
			Float_t binsRap[6] = {-0.8, -0.35, -0.15, 0.15, 0.35, 0.8};
			TH1F *fractionRap = new TH1F("fractionRap","fractionRap",5,binsRap);
			if(channel ==  1)fractionRap->SetTitle("J/#psi #rightarrow #mu^{+} #mu^{-}");
			if(channel == -1)fractionRap->SetTitle("J/#psi #rightarrow e^{+} e^{-}");
			fractionRap->GetYaxis()->SetTitle("f_{I}");
			fractionRap->GetXaxis()->SetTitle("y");
			}
		fractionRap->SetName("fractionRap");
		
		fractionRap->SetBinContent(iPoint+1,result_fI*correctionfI);
		fractionRap->SetBinError(iPoint+1,error_fI*correctionfI); 
		fractionRap->SetBinContent(5-iPoint,result_fI*correctionfI);
		fractionRap->SetBinError(5-iPoint,error_fI*correctionfI);
		fractionFile->Delete("fractionRapStore;*");
		fractionRap->SetName("fractionRapStore");
		fractionRap->Write();

		if(channel ==  1)cRooFitSingle->SetName(TString::Format("FitDiMuon_%1.2f_%1.2f",minY,maxY));
		if(channel == -1)cRooFitSingle->SetName(TString::Format("FitDiElectron_%1.2f_%1.2f",minY,maxY));
		cRooFitSingle->Write();
		fractionFile->Close();
		}
	if(ZDCcut != -1){
		TH1F *fractionZDC = 0x0;
		fractionZDC = (TH1F*)fractionFile->Get("fractionZDCStore");
		TH1F *nGGhistZDC = 0x0;
		nGGhistZDC = (TH1F*)fractionFile->Get("nGGhistZDCStore");

		if(!fractionZDC){ 
			cout<<"Creating new histo"<<endl;

			TH1F *fractionZDC = new TH1F("fractionZDC","fractionZDC",3,0.5,3.5);
			if(channel ==  1)fractionZDC->SetTitle("J/#psi #rightarrow #mu^{+} #mu^{-}");
			if(channel == -1)fractionZDC->SetTitle("J/#psi #rightarrow e^{+} e^{-}");
			fractionZDC->GetYaxis()->SetTitle("f_{I}");
			fractionZDC->GetXaxis()->SetTitle("Neutron emission");
			}
		fractionZDC->SetName("fractionZDC");

		fractionZDC->SetBinContent(iPoint+1,result_fI*correctionfI);
		fractionZDC->SetBinError(iPoint+1,error_fI*correctionfI); 
		fractionFile->Delete("fractionZDCStore;*");
		fractionZDC->SetName("fractionZDCStore");
		fractionZDC->Write();

		if(channel ==  1)cRooFitSingle->SetName(TString::Format("FitDiMuon_ZDC%d",ZDCcut));
		if(channel == -1)cRooFitSingle->SetName(TString::Format("FitDiElectron_ZDC%d",ZDCcut));
		cRooFitSingle->Write();
		fractionFile->Close();
	
	}

  }

}

void myLegendSetUp(TLegend *currentLegend,float currentTextSize,int columns){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}

/*/==================Fraction Fitter===================================

TObjArray *fMCTemplates = new TObjArray(5); // MC histograms are put in this array
TH1F *hTemplateFit = new TH1F();
TH1F *hGammaGamma_Fit = new TH1F();
TH1F *hCoherent_Fit = new TH1F();
TH1F *hIncoherent_Fit = new TH1F();
TH1F *hCoherentFD_Fit = new TH1F();
TH1F *hIncoherentFD_Fit = new TH1F();


fMCTemplates->Add(hGammaGamma_Pt);
fMCTemplates->Add(hCoherent_Pt);
fMCTemplates->Add(hIncoherent_Pt);
//fMCTemplates->Add(hCoherentFD_Pt);
//fMCTemplates->Add(hIncoherentFD_Pt);

TFractionFitter* fit = new TFractionFitter(hData_Pt, fMCTemplates); // initialise
fit->Constrain(1,0.0,0.8);	       
fit->Constrain(2,0.0,0.8);	       
fit->Constrain(3,0.0,0.5); 
//fit->Constrain(4,0.0,0.01);
//fit->Constrain(5,0.0,0.01);	    

Int_t status = fit->Fit();	       
cout << "fit status: " << status << endl;
if (status == 0) {			 // check on fit status
        TH1F *hTemplateFit = (TH1F*) fit->GetPlot();
        hTemplateFit->SetLineColor(kBlack);
	hTemplateFit->SetLineWidth(3);
	
	TCanvas *cTemplates = new TCanvas("cTemplates"," ",800,600);
	cTemplates->cd(); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15) ;
	
	hData_Pt->Draw("E1");
	hTemplateFit->Draw("same");
	
	Double_t p0, p1, p2, p3, p4, errP0, errP1, errP2, errP3, errP4;
	fit->GetResult( 0, p0, errP0);
	fit->GetResult( 1, p1, errP1);
	fit->GetResult( 2, p2, errP2);
	//fit->GetResult( 3, p3, errP3);
	//fit->GetResult( 4, p4, errP4);
	Float_t Ndata = hData_Pt->Integral();
	
	TH1F *mcp0 = (TH1F*)fit->GetMCPrediction(0);
	mcp0->SetLineColor(2);
	mcp0->Scale(Ndata*p0/mcp0->Integral());
	mcp0->Draw("same");
	TH1F *mcp1 = (TH1F*)fit->GetMCPrediction(1);
	mcp1->SetLineColor(3);
	mcp1->Scale(Ndata*p1/mcp1->Integral());
	mcp1->Draw("same");
	TH1F *mcp2 = (TH1F*)fit->GetMCPrediction(2);
	mcp2->SetLineColor(4);
	mcp2->Scale(Ndata*p2/mcp2->Integral());
	mcp2->Draw("same");
	
	mcp3 = (TH1F*)fit->GetMCPrediction(3);
	mcp3->SetLineColor(6);
	mcp3->Scale(Ndata*p3/mcp3->Integral());
	mcp3->Draw("same");
	mcp4 = (TH1F*)fit->GetMCPrediction(4);
	mcp4->SetLineColor(7);
	mcp4->Scale(Ndata*p4/mcp4->Integral());
	mcp4->Draw("same");
	
	
	myLegend1->Draw();
	myLegend2->Draw();
	
	TCanvas *cTemplatesLog = new TCanvas("cTemplatesLog"," ",800,600);
	cTemplatesLog->cd(); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15) ;
	gPad->SetLogy();
	
	hData_Pt->Draw("E1");
	hTemplateFit->Draw("same");
	mcp0->Draw("same");
	mcp1->Draw("same");
	mcp2->Draw("same");
	myLegend1->Draw();
	myLegend2->Draw();
	
	if(channel == -1){
		cTemplates->SaveAs("PtFit.El.FractionFitter.png");
		cTemplatesLog->SaveAs("PtFitLog.El.FractionFitter.png");
		}
	if(channel ==  1){
		cTemplates->SaveAs("PtFit.Mu.FractionFitter.png");
		cTemplatesLog->SaveAs("PtFitLog.Mu.FractionFitter.png");
		}
	
	Float_t fI = mcp2->Integral(mcp2->FindBin(0.0),mcp2->FindBin(0.2))/mcp1->Integral(mcp1->FindBin(0.0),mcp1->FindBin(0.2));
	Float_t fC = mcp1->Integral(mcp1->FindBin(0.2),mcp1->FindBin(1.0))/mcp2->Integral(mcp2->FindBin(0.2),mcp2->FindBin(1.0));
	cout<<"FractionFitter f_I = "<<fI<<endl;
	cout<<"FractionFitter f_C = "<<fC<<endl;
	
}
fMCTemplates->Clear();
/*/
