// LuminosityCalculationStandAlone.C
// Updated: David Grund, Sep 14, 2021
// To calculate the integrated luminosity of the analyzed samples using the trending files

// C/C++
#include "iostream"
#include "iomanip"
#include "fstream"
#include "sstream"
#include "stdio.h"
#include "string.h"

// ROOT
#include "TCanvas.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TEfficiency.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitter.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"
#include "TParticle.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector.h"

// AliROOT
#include "AliCDBManager.h"
#include "AliTriggerScalers.h"
#include "AliTriggerRunScalers.h"
#include "AliTimeStamp.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerConfiguration.h"
#include "AliLHCData.h"
#include "AliTriggerClass.h"
#include "AliTriggerBCMask.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

Int_t opt_global = -1;

using namespace std;

Int_t RunList18q[123] = {
    //2018q, 123 runs
    295585, 295586, 295588, 295589, 295610, 295611, 295612, 295615, 295666, 295667, 295668, 295673, 295675, 
    295676, 295712, 295714, 295717, 295718, 295719, 295721, 295723, 295725, 295754, 295755, 295758, 295759, 
    295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829, 295831, 
    295853, 295854, 295855, 295856, 295859, 295860, 295861, 295909, 295910, 295913, 295936, 295937, 295941, 
    295942, 296016, 296060, 296062, 296063, 296065, 296066, 296123, 296132, 296133, 296134, 296135, 296142, 
    296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242, 296243, 296244, 
    296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303, 296304, 296309, 296312, 296377, 296378, 
    296379, 296380, 296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472, 296509, 
    296510, 296511, 296512, 296516, 296547, 296548, 296549, 296550, 296551, 296552, 296553, 296594, 296615, 
    296616, 296618, 296619, 296621, 296622, 296623
};

Int_t RunList18r[96] = {
    //2018r, 96 runs
    296690, 296691, 296693, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787, 296790, 296793, 
    296794, 296799, 296835, 296836, 296838, 296839, 296848, 296849, 296850, 296851, 296852, 296890, 296894, 
    296899, 296900, 296903, 296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 297029, 297031, 
    297035, 297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129, 297132, 297133, 297193, 297194, 
    297195, 297196, 297218, 297219, 297221, 297222, 297278, 297310, 297311, 297317, 297332, 297333, 297335, 
    297336, 297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414, 297415, 297441, 
    297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540, 297541, 297542, 
    297544, 297558, 297588, 297590, 297595
};

Int_t nRunsInList = -1;
Int_t *RunList = NULL;
Int_t *CountsList = NULL;

void SetLumiHisto(TH1D* h, Color_t color = kBlue){
  // A function to set the properties of the final histograms
  gStyle->SetOptStat(0);
  h->SetTitleFont(43);
  h->SetTitleSize(25);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(28);
  h->GetXaxis()->SetLabelSize(13);
  h->GetYaxis()->SetLabelSize(28);
  h->GetYaxis()->SetTickLength(0.01);
  if(opt_global == 0){
    h->GetYaxis()->SetTitleOffset(0.45);
  } else if(opt_global == 1) h->GetYaxis()->SetTitleOffset(0.3);
  h->GetYaxis()->SetDecimals(1);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetFillColor(color);
  h->SetMarkerSize(0.7);
  h->SetMarkerStyle(kFullCross);
  h->LabelsOption("v");
  h->SetMinimum(0);
  h->SetLineWidth(2);
  h->Sumw2(kFALSE);
}

// Plot lumi per run
void SumLumi(){
  // In here the total integrated luminosity is computed:
  // open input file and read lumi per run
  TFile *fInputLumiHisto = TFile::Open(Form("LumiHisto.root"), "read");
  TH1D *hLumi = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumi"));
  TH1D *hLumiS = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumiS"));
  Double_t integrated_lumi_ana(0.), integrated_lumi_rec(0.);
  for(Int_t i(0); i < hLumiS->GetNbinsX(); i++){
    integrated_lumi_ana += hLumiS->GetBinContent(i+1);
    integrated_lumi_rec += hLumi->GetBinContent(i+1);
  }

  SetLumiHisto(hLumi);
  SetLumiHisto(hLumiS,kRed);

  TCanvas *cLumi = new TCanvas("cLumi","cLumi",1300,400);
  cLumi->SetTopMargin(0.02);
  cLumi->SetRightMargin(0.01);
  //cLumi->SetLeftMargin(0.07);
  cLumi->SetBottomMargin(0.15);
  hLumi->GetYaxis()->SetTitle("L_{int} [#mub^{-1}]");
  hLumi->Draw();
  hLumiS->Draw("sameP0");
  TLegend* legLumi = new TLegend(0.28,0.70,0.52,0.96);
  legLumi->SetFillColor(kWhite);
  if(opt_global == 0){
    legLumi->SetHeader("CCUP31 trigger class, LHC18q","l");
    cLumi->SetLeftMargin(0.07);
  } else if(opt_global == 1){
    legLumi->SetHeader("CCUP31 trigger class, LHC18r","l");
    cLumi->SetLeftMargin(0.05);
  }
  legLumi->AddEntry(hLumi,Form("Total lumi rec.: %.3f #mub^{-1}",integrated_lumi_rec),"l");
  legLumi->AddEntry(hLumiS,Form("Total lumi ana.: %.3f #mub^{-1}",integrated_lumi_ana),"l");
    legLumi->AddEntry((TObject*)0,Form("#bf{This thesis}"),"");
  legLumi->SetTextSize(0.055);
  legLumi->SetBorderSize(0);
  legLumi->SetFillStyle(0);
  legLumi->Draw();
  if(opt_global == 0){
    // 2018q
    cLumi->SaveAs(Form("lumi_ccup31_18q.pdf"));
  } else if(opt_global == 1){
    cLumi->SaveAs(Form("lumi_ccup31_18r.pdf"));    
  }


}//end sumlumi

void LuminosityCalculationStandAlone(Int_t opt){
  opt_global = opt;
  
  // Arrays with fired triggers per run in my dataset
  // 2018q:
  // Increasing order (from 295585 to 296623)
  Int_t Counts18q[123] = {
    1437, 9412, 8871, 20469, 1194, 1552, 27268, 3754, 3587, 1627, 6494, 54341, 22194, 
    83847, 19649, 33473, 1424, 7941, 6451, 4050, 7976, 18441, 50909, 53413, 102486, 24467, 
    24079, 78993, 26100, 170833, 67443, 60581, 7861, 151223, 158432, 18147, 203468, 86928, 74886, 
    12423, 92099, 109025, 105882, 80091, 56133, 72294, 20309, 283366, 154593, 87766, 19324, 97752, 
    104581, 16855, 15839, 93031, 147048, 130988, 31833, 22272, 117295, 88909, 202591, 109475, 73312, 
    19463, 248482, 21027, 162120, 40184, 103878, 74517, 26328, 10041, 38644, 44455, 60796, 438079, 
    80228, 46970, 197663, 77608, 374508, 31408, 65267, 115394, 311491, 87954, 80068, 317138, 299847, 
    80530, 107767, 57698, 65936, 263008, 158927, 142479, 57342, 55723, 17954, 178499, 32183, 144623, 
    447963, 127687, 20080, 24370, 54637, 77211, 284437, 222096, 90738, 23139, 31230, 16106, 109986, 
    38995, 167669, 113039, 14201, 49971, 108904
  };
  // 2018r:
  // Increasing order (from 296690 to 297595)
  Int_t Counts18r[96] = {
    445822, 40468, 19060, 295954, 591829, 568073, 44130, 166915, 90397, 39207, 331402, 14923, 71113, 
    155519, 141597, 15307, 152917, 24253, 171664, 117416, 565813, 151481, 29071, 66309, 521260, 281767, 
    143134, 294328, 59652, 278584, 37110, 100410, 174076, 274603, 195054, 202052, 197906, 379260, 300011, 
    68256, 62412, 125344, 120861, 233226, 164489, 54190, 202393, 242382, 123389, 43871, 404369, 530537, 
    20778, 101409, 346548, 514265, 129032, 75880, 33768, 31449, 11574, 189079, 23496, 8359, 32374, 
    9717, 105290, 102229, 297634, 327571, 377723, 71205, 32472, 20938, 173432, 136227, 387223, 283883, 
    106716, 423901, 104321, 72747, 56314, 424440, 554439, 111360, 88212, 98061, 23936, 208110, 76480, 
    377664, 22174, 294363, 199435, 82380
  };
  
  if(opt == 0){
    // 2018q
    RunList = RunList18q;
    CountsList = Counts18q;
    nRunsInList = 123;

  } else if(opt == 1){
    // 2018r
    RunList = RunList18r;
    CountsList = Counts18r;
    nRunsInList = 96;
  }

  const Int_t nrunsmax = 130;

  TString className1 = "CCUP31-B-NOPF-CENTNOTRD";
  TString className2 = "CCUP31-B-SPD2-CENTNOTRD";
  TChain* t = new TChain("trending");
  t->AddFile("trending_merged_PbPb_2018.root");
  t->LoadTree(0);
  TObjArray* classes = new TObjArray();
  Double_t  lumi_seen[nrunsmax] = {0};
  Double_t  class_lumi[nrunsmax] = {0};
  Double_t  class_ds[nrunsmax] = {0};
  ULong64_t  class_l2a[nrunsmax] = {0};
  Int_t run;
  Double_t mu = 0;
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("class_ds",&class_ds);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->BuildIndex("run");

  TH1D* hLumi = new TH1D("hLumi","",nRunsInList,0,nRunsInList); // Lumi per event from the trending file for a corresponding trigger class
  TH1D* hLumiS = new TH1D("hLumiS","",nRunsInList,0,nRunsInList); // My recorded (seen) lumi
  TH1D* hScale = new TH1D("hScale","",nRunsInList,0,nRunsInList); // Scale between my lumi < official lumi
  TH1D* hCCUP31ds = new TH1D("hCCUP31ds","CCUP31 downscaling",nRunsInList,0,nRunsInList); // Downscale of the trigger class in each run?

  // Calculate seen luminosity for specific trigger class
  //for (Int_t i=0;i<nRunsInList;i++){  // ordering of runs in hist
  Int_t counter = 0;

  for (Int_t i = 0; i < nRunsInList; i++){
    Int_t r = RunList[i];
    char* srun = Form("%i",r); // Int run number to char/string run number
    //Printf("%s %i %i",srun, i, r);
    t->GetEntryWithIndex(r);

    // Check trigger class name
    AliTriggerClass* cl;
    if(r < 295881){
      cl = (AliTriggerClass*) classes->FindObject(className1.Data());
    }
    if(r >= 295881){
      cl = (AliTriggerClass*) classes->FindObject(className2.Data());
    }    
    if (!cl) continue;

    Int_t iclass = classes->IndexOf(cl);
    //Printf("%i %i %s",run, iclass, cl->GetName());
    Double_t l2a = (Double_t) class_l2a[iclass];
    //Printf("%.llu",class_l2a[iclass]);
    //Printf("%.10f",class_lumi[iclass]);
    Double_t scale = CountsList[i]/l2a;

    if(scale > 1.0){
      cout << "In run " << r << " the scale is " << scale << endl;
      counter++;
    }

    hScale->Fill(srun,scale);
    hLumiS->Fill(srun,scale*class_lumi[iclass]);
    hLumi->Fill(srun,class_lumi[iclass]);
    hCCUP31ds->Fill(srun,class_ds[iclass]);
  }

  // Write the results to the output root file:
  TFile* fOutputLumiHisto = new TFile(Form("LumiHisto.root"),"recreate");
  hLumi->Write();
  hLumiS->Write();
  hScale->Write();
  hCCUP31ds->Write();
  fOutputLumiHisto->Close();

  if(counter == 0){
    cout << "No scale factors above 1.00." << endl;
  } else if(counter > 0){
    cout << counter << " suspicious runs." << endl;
  }

  // Calculate the integrated luminosity (seen and official)
  SumLumi();
}
