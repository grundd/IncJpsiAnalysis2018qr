// AccAndEffMC.c
// David Grund, 14-09-2021
// To calculate the acceptance x efficiency from MC data

#include <iostream>
#include <fstream>  

using namespace std;

Int_t runNumber;
Bool_t triggerInputsMC[11];
Double_t xPt, xM, xY, xEta1, xEta2, xQ1, xQ2;
Double_t xPtGen, xYGen, xMGen;
Int_t xV0A_dec, xV0C_dec, xADA_dec, xADC_dec;
Double_t xSigIfEl1, xSigIfEl2, xSigIfMu1, xSigIfMu2;

void ConnectTreeVariables(TTree *t);
void ConnectTreeVariablesMCgen(TTree *t);
Bool_t EventPassedMCrec(Int_t pt_bin, Bool_t *cuts);
Bool_t EventPassedMCgen(Int_t pt_bin);
void FillHistNrec(Bool_t *cuts);
void FillHistNrecPtGen(Bool_t *cuts);
void FillHistNgen();
void SaveToFile(TH1D* hist, TString name);
TString ConvertCutsToString(Bool_t *cuts);
void PlotHistogram(TH1D* hist, TString option, TString titleY, Int_t color, Double_t width);
void SetPad(TPad* pad);
void PlotLegendMain();
void PlotLegendThisWork();
void PlotNrecRatios(TCanvas* canvas);
void PlotLegendRatios(Int_t i_hist);

//const Int_t n_bins = 30;
const Int_t n_bins = 20;
//Double_t edges[n_bins + 1] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
Double_t edges[n_bins + 1] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.0,1.2,1.4,1.6,1.8,2.0};
TH1D* hist_Nrec = new TH1D("hist_Nrec","N rec per bin",n_bins,edges);
TH1D* hist_Ngen = new TH1D("hist_Ngen","N gen per bin",n_bins,edges);
//TH1D* hist_AxE = new TH1D("hist_AxE","AxE per bin",n_bins,edges);
// For the AxE_both
TH1D* hist_Nrec_PtGen = new TH1D("hist_Nrec","N rec per bin",n_bins,edges);

// When reading the columns from txt files
Int_t InBin;
Double_t InValue;

const Int_t n_cuts = 10;
Bool_t cuts[n_cuts] = {
    kFALSE, //0 pt cut rec (kFALSE) or gen (kTRUE)
    kTRUE, //1 CCUP31
    kTRUE, //2 AD veto => negligible effect
    kTRUE, //3 V0 veto => negligible effect
    kTRUE, //4 Rap
    kTRUE, //5 Eta
    kTRUE, //6 Charge 
    kTRUE, //7 Muons
    kTRUE, //8 Mass
    kTRUE  //9 Pt
};

Bool_t AxE = kTRUE;
Bool_t ratios = kFALSE;
Bool_t AxE_both = kFALSE;

void AccAndEff_new(){
    Printf("*** This is AccAndEff_new.c ***");

    if(AxE == kTRUE){
        // Acceptance and efficiency:
        TCanvas *canvas1 = new TCanvas("canvas1","canvas1",800,600);
        TPad *pad = new TPad("pad", "pad",0.0,0.0,1.0,1.0);
        pad->Draw();
        pad->cd();
        // Calculate: 
        FillHistNrec(cuts);
        FillHistNgen();
        hist_Nrec->Sumw2();
        hist_Nrec->Divide(hist_Ngen);
        //hist_Nrec->Scale(1.,"width");
        // Calculate the mean value as weighted average over bins
        Int_t N_rec_tot = 0;
        Int_t N_gen_tot = 0;
        // Bin 5 starts at pt = 0.2 GeV/c
        for(Int_t i = 1; i <= n_bins; i++){
            N_gen_tot += hist_Ngen->GetBinContent(i);
        }
        for(Int_t i = 5; i <= n_bins; i++){
            N_rec_tot += hist_Nrec->GetBinContent(i)*hist_Ngen->GetBinContent(i);
        }
        Double_t AxE_tot = N_rec_tot/(Double_t)N_gen_tot;
        Printf("Total N rec: %i", N_rec_tot);
        Printf("Total N gen: %i", N_gen_tot);
        Printf("Total AxE calculated: %f", AxE_tot);
        // Make the plot:
        SetPad(pad);
        PlotHistogram(hist_Nrec,"P E1","#it{N}_{rec}/#it{N}_{gen}", 4, 1);
        PlotLegendMain();
        PlotLegendThisWork();
        // Save the plot as pdf and png
        TString* path1 = new TString("AxE_files_new/output/");
        path1->Append(ConvertCutsToString(cuts));
        path1->Append("_AxE.");
        TString* path2 = new TString(*path1);
        canvas1->Print(path1->Append("pdf")); 
        canvas1->Print(path2->Append("png"));
    }

    if(ratios == kTRUE){
        // N_rec ratios:
        TCanvas *canvas2 = new TCanvas("canvas2","canvas2",1000,600);
        //SetCanvas(canvas2);
        PlotNrecRatios(canvas2);
        // Save the plot as pdf and png
        TString* path3 = new TString("AxE_files_new/output/");
        TString* path4 = new TString("AxE_files_new/output/");
        canvas2->Print(path3->Append("N_rec_ratios.pdf")); 
        canvas2->Print(path4->Append("N_rec_ratios.png"));
    }

    if(AxE_both == kTRUE){
        // Acceptance and efficiency:
        TCanvas *canvas1 = new TCanvas("canvas1","canvas1",800,600);
        TPad *pad = new TPad("pad", "pad",0.0,0.0,1.0,1.0);
        pad->Draw();
        pad->cd();
        // Calculate with bins of pt rec: 
        cuts[0] = kFALSE;
        FillHistNrec(cuts);
        FillHistNgen();
        hist_Nrec->Sumw2();
        hist_Nrec->Divide(hist_Ngen);
        // Calculate with bins of pt gen: 
        cuts[0] = kTRUE;
        FillHistNrecPtGen(cuts);
        hist_Nrec_PtGen->Sumw2();
        hist_Nrec_PtGen->Divide(hist_Ngen);

        // Make the plot:
        SetPad(pad);
        PlotHistogram(hist_Nrec,"P E1","#it{N}_{rec}/#it{N}_{gen}", 4, 1);
        hist_Nrec_PtGen->SetLineColor(kRed);
        hist_Nrec_PtGen->SetLineWidth(1);
        hist_Nrec_PtGen->SetMarkerStyle(21);
        hist_Nrec_PtGen->SetMarkerColor(kRed);
        hist_Nrec_PtGen->SetMarkerSize(1.0);
        hist_Nrec_PtGen->Draw("SAME E1 P");
        // Draw the legend
        PlotLegendMain();
        PlotLegendThisWork();
        TLegend *leg2 = new TLegend(0.19,0.44,0.4,0.60); 
        leg2->AddEntry(hist_Nrec, Form("vs #it{p}_{T} rec"));
        leg2->AddEntry(hist_Nrec_PtGen, Form("vs #it{p}_{T} gen"));
        leg2->SetTextSize(0.058);
        leg2->SetBorderSize(0); // no border
        leg2->SetFillStyle(0);  // legend is transparent
        leg2->Draw();
        // Save the plot as pdf and png
        TString* path1 = new TString("AxE_files_new/output/");
        path1->Append(ConvertCutsToString(cuts));
        path1->Append("_AxE_both.");
        TString* path2 = new TString(*path1);
        canvas1->Print(path1->Append("pdf")); 
        canvas1->Print(path2->Append("png"));
    }

    return;
}

void FillHistNrec(Bool_t *cuts){
    // Check if the corresponding file already exists
    TString path("AxE_files_new/");
    TString cut = ConvertCutsToString(cuts);
    cut.Append(".txt");
    path.Append(cut);

    ifstream InFile;
    InFile.open(path);
    if(!(InFile.fail())){
        // This sequence of cuts has been already analysed
        Printf("*** The file %s already exists. ***", cut.Data());
        // Fill the hist_Ngen from the txt file
        Int_t i_line = 0;
        while(!InFile.eof()){
            InFile >> InBin >> InValue; // fist and second column
            hist_Nrec->SetBinContent(InBin, InValue);
            i_line++;
        }
        InFile.close(); 

        return;
    } else {
        // Reconstructed level:
        Printf("*** Calculating N rec per bin for %s... ***", cut.Data());

        TFile *f_rec = TFile::Open("MC_rec_gen_results/merged_LHC18qr/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
        //if(f_rec) Printf("MC rec file loaded.");

        TTree *t_rec = dynamic_cast<TTree*> (f_rec->Get("analysisTree"));
        //if(t_rec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariables(t_rec);

        // Loop over all pt bins
        for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
            Int_t N_rec = 0;
            for(Int_t iEntry = 0; iEntry < t_rec->GetEntries(); iEntry++){
                t_rec->GetEntry(iEntry);
                if(EventPassedMCrec(i_bin, cuts)) N_rec++;
            }
            hist_Nrec->SetBinContent(i_bin, N_rec);
            Printf("*** Bin %i done. ***", i_bin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hist_Nrec, path);
        return;
    }
}

void FillHistNrecPtGen(Bool_t *cuts){
    // Check if the corresponding file already exists
    TString path("AxE_files_new/");
    TString cut = ConvertCutsToString(cuts);
    cut.Append(".txt");
    path.Append(cut);

    ifstream InFile;
    InFile.open(path);
    if(!(InFile.fail())){
        // This sequence of cuts has been already analysed
        Printf("*** The file %s already exists. ***", cut.Data());
        // Fill the hist_Ngen from the txt file
        Int_t i_line = 0;
        while(!InFile.eof()){
            InFile >> InBin >> InValue; // fist and second column
            hist_Nrec_PtGen->SetBinContent(InBin, InValue);
            i_line++;
        }
        InFile.close(); 

        return;
    } else {
        // Reconstructed level:
        Printf("*** Calculating N rec per bin for %s... ***", cut.Data());

        TFile *f_rec = TFile::Open("MC_rec_gen_results/merged_LHC18qr/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
        //if(f_rec) Printf("MC rec file loaded.");

        TTree *t_rec = dynamic_cast<TTree*> (f_rec->Get("analysisTree"));
        //if(t_rec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariables(t_rec);

        // Loop over all pt bins
        for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
            Int_t N_rec = 0;
            for(Int_t iEntry = 0; iEntry < t_rec->GetEntries(); iEntry++){
                t_rec->GetEntry(iEntry);
                if(EventPassedMCrec(i_bin, cuts)) N_rec++;
            }
            hist_Nrec_PtGen->SetBinContent(i_bin, N_rec);
            Printf("*** Bin %i done. ***", i_bin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hist_Nrec_PtGen, path);
        return;
    }
}

void FillHistNgen(){
    TString path("AxE_files_new/N_gen.txt");

    ifstream InFile;
    InFile.open(path);
    if(!(InFile.fail())){
        // This sequence of cuts has been already analysed
        Printf("*** The file N_gen.txt already exists. ***");
        // Fill the hist_Ngen from the txt file
        Int_t i_line = 0;
        while(!InFile.eof()){
            InFile >> InBin >> InValue; // fist and second column
            hist_Ngen->SetBinContent(InBin, InValue);
            i_line++;
        }
        InFile.close(); 

        return;
    } else {
        // Generated level:
        Printf("*** Calculating N gen per bin... ***");

        TFile *f_gen = TFile::Open("MC_rec_gen_results/merged_LHC18qr/MC_gen_18qr_kIncohJpsiToMu_migr.root", "read");
        //if(f_gen) Printf("MC gen file loaded.");

        TTree *t_gen = dynamic_cast<TTree*> (f_gen->Get("MCgenTree"));
        //if(t_gen) Printf("MC gen tree loaded.");

        ConnectTreeVariablesMCgen(t_gen);
        // Loop over all pt bins
        for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
            Int_t N_gen = 0;
            for(Int_t iEntry = 0; iEntry < t_gen->GetEntries(); iEntry++){
                t_gen->GetEntry(iEntry);
                if(EventPassedMCgen(i_bin)) N_gen++;
            }
            hist_Ngen->SetBinContent(i_bin, N_gen);
            Printf("*** Bin %i done. ***", i_bin);
        }
        Printf("*** Finished! ***");

        TString s("AxE_files_new/N_gen.txt");
        SaveToFile(hist_Ngen, s);

        return;
    }
}

void SetPad(TPad* pad){
    pad->SetTopMargin(0.02);
    pad->SetBottomMargin(0.14);
    pad->SetRightMargin(0.02);
    pad->SetLeftMargin(0.145);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    return;
}

void PlotHistogram(TH1D* hist, TString option, TString titleY, Int_t color, Double_t width){
    // Vertical axis
    //hist->GetYaxis()->SetTitle("(#it{A} #times #it{#varepsilon})_{MC}");
    hist->GetYaxis()->SetTitle(titleY.Data()); // "#it{N}_{rec}/#it{N}_{gen}"
    hist->GetYaxis()->SetTitleSize(0.058);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetLabelSize(0.058);
    hist->GetYaxis()->SetLabelOffset(0.01);
    if(ratios == kTRUE){
        hist->GetYaxis()->SetRangeUser(0.0,1.05);
    }
    
    // Horizontal axis
    hist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hist->GetXaxis()->SetTitleSize(0.058);
    hist->GetXaxis()->SetLabelSize(0.058);
    //hist->GetXaxis()->SetRange(0,6);

    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(kBlue);
    hist->SetMarkerSize(1.0);
    hist->SetLineColor(color);
    hist->SetLineWidth(width);
    hist->Draw(option.Data()); // P E1

    /*
    Double_t x = 0.2;
    TLine *l = new TLine(x,0.01,x,0.062);
    l->SetLineColor(222);
    l->SetLineWidth(2);
    l->SetLineStyle(9);
    l->Draw("");
    */

    return;
}

void PlotLegendMain(){
    // Legend1
    TLegend *leg1 = new TLegend(0.32,0.77,0.85,0.97);
    //leg1->SetHeader("ALICE Simulation, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    leg1->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    leg1->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    leg1->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    //leg1->AddEntry((TObject*)0,Form("#bf{This thesis}"),"");
    leg1->SetTextSize(0.058);
    leg1->SetBorderSize(0); // no border
    leg1->SetFillStyle(0);  // legend is transparent
    leg1->Draw();

    return;
}

void PlotLegendThisWork(){
    // Legend1
    TLegend *leg1 = new TLegend(0.12,0.17,0.4,0.36);
    leg1->AddEntry((TObject*)0,Form("#bf{This thesis}"),"");
    leg1->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),""); 
    leg1->AddEntry((TObject*)0,Form("3.0 < #it{m}_{#mu#mu} < 3.2 GeV/#it{c}^{2}"),"");
    //leg1->AddEntry((TObject*)0,Form("#bf{This work}"),"");
    leg1->SetTextSize(0.058);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw();

    return;
}

void PlotLegendRatios(Int_t n_hist, TH1D* hist[n_hist]){
    // Legend1
    TLegend *leg1 = new TLegend(0.0,0.1,1.0,0.98);
    leg1->AddEntry(hist[0],"only S9");
    leg1->AddEntry((TObject*)0,Form("(= 1.0)"),"");
    leg1->AddEntry((TObject*)0,Form(""),"");
    leg1->AddEntry((TObject*)0,Form("S9 and:"),"");
    for(Int_t i_hist = 1; i_hist < n_hist; i_hist++){
        leg1->AddEntry(hist[i_hist],Form("S%i", i_hist));
    }
    leg1->SetTextSize(0.14);
    leg1->SetBorderSize(0); // no border
    leg1->SetFillStyle(0);  // legend is transparent
    leg1->Draw();

    return;
}

void PlotNrecRatios(TCanvas* canvas){
    Bool_t cuts_ratios[n_cuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kTRUE};
    TH1D* hist_Nrec_ratios[n_cuts-1];
    // Calculate N_rec per bin with just the pt cut
    FillHistNrec(cuts_ratios);
    // Store the values in the zeroth component
    hist_Nrec_ratios[0] = new TH1D("hist0","only pt cut",n_bins,edges);
    for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
        hist_Nrec_ratios[0]->SetBinContent(i_bin, hist_Nrec->GetBinContent(i_bin));
    }
    // Calculate the N_rec per bin for the given i-th cut
    for(Int_t i_cut = 1; i_cut < n_cuts-1; i_cut++){
        cuts_ratios[i_cut] = kTRUE;
        // Create corresponding i-th histogram
        TString name("hist%i", i_cut);
        TString title("hist after the cut %i", i_cut);
        hist_Nrec_ratios[i_cut] = new TH1D(name,title,n_bins,edges);
        // Calculate N_rec
        FillHistNrec(cuts_ratios);
        // Copy the values
        for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
            hist_Nrec_ratios[i_cut]->SetBinContent(i_bin, hist_Nrec->GetBinContent(i_bin));
        }
        cuts_ratios[i_cut] = kFALSE;  
    }
    // Plot the results
    Int_t DrawUpToCut = 8; // maximum = 8
    //hist_Nrec_ratios[0]->Scale(1.,"width");
    TH1D* hist_one = new TH1D("hist_one","ones",n_bins,edges);
    for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
        hist_one->SetBinContent(i_bin, 1.);
    }
    TPad *pad_l = new TPad("pad_l", "pad_l",0.00,0.0,0.76,1.0);
    TPad *pad_r = new TPad("pad_r", "pad_r",0.76,0.0,1.00,1.0);
    pad_l->Draw();
    pad_r->Draw();
    SetPad(pad_l);
    pad_l->cd();
    PlotHistogram(hist_one,"HIST","#it{N}_{rec}[S9 + Si]/#it{N}_{rec}[S9]",1,2);
    for(Int_t i_cut = 1; i_cut <= DrawUpToCut; i_cut++){
        hist_Nrec_ratios[i_cut]->Sumw2();
        hist_Nrec_ratios[i_cut]->Divide(hist_Nrec_ratios[0]);
        //hist_Nrec_ratios[i_cut]->Scale(1.,"width");
        // no yellow:
        if(i_cut != 4) hist_Nrec_ratios[i_cut]->SetLineColor(i_cut+1);
        else hist_Nrec_ratios[i_cut]->SetLineColor(218);
        hist_Nrec_ratios[i_cut]->SetLineWidth(1);
        hist_Nrec_ratios[i_cut]->SetMarkerStyle(21);
        if(i_cut != 4) hist_Nrec_ratios[i_cut]->SetMarkerColor(i_cut+1);
        else hist_Nrec_ratios[i_cut]->SetMarkerColor(218);
        hist_Nrec_ratios[i_cut]->SetMarkerSize(1.);
        hist_Nrec_ratios[i_cut]->Draw("SAME E1 P");
    }
    // Plot the extra legend with this thesis
    TLegend *leg2 = new TLegend(0.15,0.2,0.4,0.28); 
    leg2->AddEntry((TObject*)0,Form("#bf{This thesis}"),"");
    leg2->SetTextSize(0.058);
    leg2->SetBorderSize(0); // no border
    leg2->SetFillStyle(0);  // legend is transparent
    leg2->Draw();
    // Right pad
    pad_r->cd();
    PlotLegendRatios(DrawUpToCut+1, hist_Nrec_ratios);
    return;
}

//***********************************************************************************************
//************************************ Supporting functions *************************************
//***********************************************************************************************

void SaveToFile(TH1D* hist, TString name){
    ofstream outfile (name.Data());
    for(Int_t i_bin = 1; i_bin <= hist->GetNbinsX(); i_bin++){
        outfile << i_bin << "\t" << hist->GetBinContent(i_bin) << "\n";
    }
    outfile.close();
    Printf("*** File saved in %s.***", name.Data());
}

TString ConvertCutsToString(Bool_t *cuts){
    TString s("cut_");
    for(Int_t i_cut = 0; i_cut < n_cuts; i_cut++){
        if(cuts[i_cut] == kTRUE) s.Append("1");
        else s.Append("0");
    }
    //cout << s << endl;
    return s;
}

void ConnectTreeVariables(TTree *t){
    // Set Branch Addresses
    // Triggers
    t->SetBranchAddress("runNumber", &runNumber);
    t->SetBranchAddress("triggerInputsMC", &triggerInputsMC[0]);
    // Kinematics
    t->SetBranchAddress("fPt", &xPt);
    t->SetBranchAddress("fM", &xM);
    t->SetBranchAddress("fY", &xY);
    t->SetBranchAddress("Eta_1", &xEta1);
    t->SetBranchAddress("Eta_2", &xEta2);
    t->SetBranchAddress("Q_1", &xQ1);
    t->SetBranchAddress("Q_2", &xQ2);
    // Generated values
    t->SetBranchAddress("fPtGen", &xPtGen);
    // Forward detectors
    t->SetBranchAddress("ADA_decision", &xADA_dec);
    t->SetBranchAddress("ADC_decision", &xADC_dec);
    t->SetBranchAddress("V0A_decision", &xV0A_dec);
    t->SetBranchAddress("V0C_decision", &xV0C_dec);
    // PID
    t->SetBranchAddress("trk1SigIfEl", &xSigIfEl1);
    t->SetBranchAddress("trk2SigIfEl", &xSigIfEl2);
    t->SetBranchAddress("trk1SigIfMu", &xSigIfMu1);
    t->SetBranchAddress("trk2SigIfMu", &xSigIfMu2);

    //Printf("MC rec tree variables connected.");
    return;
}

void ConnectTreeVariablesMCgen(TTree *t){
    // Set Branch Addresses
    t->SetBranchAddress("runNumber", &runNumber);
    // Kinematics
    t->SetBranchAddress("fPtGen", &xPtGen);
    t->SetBranchAddress("fMGen", &xMGen);
    t->SetBranchAddress("fYGen", &xYGen);

    //Printf("MC gen tree variables connected.");
    return;
}

Bool_t EventPassedMCrec(Int_t pt_bin, Bool_t *cuts){
    // 0) Two good central tracks: applied in the primary analysis

    // 1) CCUP31 trigger
    if(cuts[1] == kTRUE){
        Bool_t CCUP31 = kFALSE;
        if(
            //!triggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
            //!triggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
            //!triggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
            //!triggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
            triggerInputsMC[10] //&&  // 0STG (SPD topological)
            //triggerInputsMC[4]      // 0OMU (TOF two hits topology)
        ) CCUP31 = kTRUE;

        if(!CCUP31) return kFALSE;
    }
    // 2) AD offline veto (negligible effect)
    if(cuts[2] == kTRUE){
        if(!(xADA_dec == 0 && xADC_dec == 0)) return kFALSE;
    }
    // 3) V0 offline veto (negligible effect)
    if(cuts[3] == kTRUE){
        if(!(xV0A_dec == 0 && xV0C_dec == 0)) return kFALSE;
    }
    // 4) Dilepton ("Jpsi") rapidity
    if(cuts[4] == kTRUE){
        if(!(abs(xY) < 0.8)) return kFALSE;
    }
    // 5) Pseudorapidity of both tracks |eta| < 0.8
    if(cuts[5] == kTRUE){
        if(!(abs(xEta1) < 0.8 && abs(xEta2) < 0.8)) return kFALSE;
    }
    // 6) Tracks have opposite charges
    if(cuts[6] == kTRUE){
        if(!(xQ1 * xQ2 < 0)) return kFALSE;
    }
    // 7) Muon pairs only
    if(cuts[7] == kTRUE){
        if(!(xSigIfMu1*xSigIfMu1 + xSigIfMu2*xSigIfMu2 < xSigIfEl1*xSigIfEl1 + xSigIfEl2*xSigIfEl2)) return kFALSE;
    }
    // 8) Invariant mass cut
    if(cuts[8] == kTRUE){
        if(!(xM > 3.0 && xM < 3.2)) return kFALSE;
        //if(!(xM > 2.2 && xM < 4.5)) return kFALSE;
    }
    // 9) Transverse momentum cut 
    if(cuts[9] == kTRUE){
        if(cuts[0] == kFALSE){
            if(!(xPt > edges[pt_bin-1] && xPt < edges[pt_bin])) return kFALSE;
        } else if(cuts[0] == kTRUE){
            if(!(xPtGen > edges[pt_bin-1] && xPtGen < edges[pt_bin])) return kFALSE;
        }
    }
    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCgen(Int_t pt_bin){
    // 1) J/psi rapidity |y| < 0.8
    if(!(abs(xYGen) < 0.8)) return kFALSE;

    // 2) Transverse momentum cut 
    if(!(xPtGen > edges[pt_bin-1] && xPtGen < edges[pt_bin])) return kFALSE; 

    // Event passed the selection =>
    return kTRUE;
}

//******************************************
// Backup:

/*
// Fill the AxE histogram - loop over all bins:
Printf("*** AxE histogram: ***");
for(Int_t i_bin = 1; i_bin <= n_bins; i_bin++){
    hist_AxE->SetBinContent(i_bin, hist_Nrec->GetBinContent(i_bin)/hist_Ngen->GetBinContent(i_bin));
    Printf("%.0f %.0f %.4f", hist_Nrec->GetBinContent(i_bin), hist_Ngen->GetBinContent(i_bin), hist_AxE->GetBinContent(i_bin));
} 
*/