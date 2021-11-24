// PhenoPredictions.c
// David Grund, Oct 30, 2021

// c++ headers
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>   // getline

// root headers
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>

// my headers
#include "AnalysisManager.h"

Int_t iFeedDown = 1;
// 0 = feed-down from PtFitWithoutBkg.c
// 1 = feed-down from FeedDown.c (FeedDown_debug.c)
Bool_t plot1 = kTRUE;
Bool_t plot2 = kTRUE;
Bool_t plot3 = kTRUE;
Bool_t plot4 = kTRUE;
// 1 = HS model
// 2 = Guzey's model
// 3 = Heikki's model
// 4 = STARlight

Int_t lineWidth = 2;

// To read the values from the files:
Double_t Sigma_val[nPtBins] = { 0 };
Double_t Sigma_err_stat[nPtBins] = { 0 };
Double_t Sigma_err_syst[nPtBins] = { 0 };
Double_t t_avg_val[nPtBins] = { 0 };
Double_t t_boundaries[nPtBins+1] = { 0 };

// For TGraphAsymmErrors:
Double_t Sigma_err_stat_upp[nPtBins] = { 0 };
Double_t Sigma_err_stat_low[nPtBins] = { 0 };
Double_t Sigma_err_syst_upp[nPtBins] = { 0 };
Double_t Sigma_err_syst_low[nPtBins] = { 0 };
Double_t t_value[nPtBins] = { 0 };
Double_t t_err_low[nPtBins] = { 0 };
Double_t t_err_upp[nPtBins] = { 0 };

// HS predictions
// reserve space for data to be read
const Int_t nData1 = 75;
Double_t abs_t_1[nData1];
Double_t coh_n[nData1];
Double_t coh_n_err[nData1];
Double_t inc_n[nData1];
Double_t inc_n_err[nData1];
Double_t coh_hs[nData1];
Double_t coh_hs_err[nData1];
Double_t inc_hs[nData1];
Double_t inc_hs_err[nData1];

// Guzey predictions
const Int_t nData2 = 51;
Double_t abs_t_2[nData2];
Double_t sig_ela_min[nData2];
Double_t sig_ela_max[nData2];
Double_t sig_dis_min[nData2];
Double_t sig_dis_max[nData2];
Double_t sig_tot_min[nData2];
Double_t sig_tot_max[nData2];

// Heikki predictions
const Int_t nData3 = 183;
Double_t abs_t_3[nData3];
Double_t sig_inc_fluct[nData3];
Double_t sig_inc_noflu[nData3];

// STARlight predictions
const Int_t nData4 = 50;
Double_t abs_t_4[nData4];
Double_t sig_SL[nData4];

void ReadInputMeasurement();
void ReadInputHSModel();
void ReadInputGuzey();
void ReadInputHeikki();
void ReadInputSTARlight();

void PhenoPredictions()
{
    ReadInputMeasurement();

    ReadInputHSModel();

    ReadInputGuzey();

    ReadInputHeikki();

    ReadInputSTARlight();

    // Fill the histogram
    for(Int_t i = 0; i < nPtBins; i++){
        // from microbarns to milibarns
        Sigma_val[i] = Sigma_val[i] / 1000.;
        Sigma_err_stat[i] = Sigma_err_stat[i] / 1000.;
        Sigma_err_syst[i] = Sigma_err_syst[i] / 1000.;
        Sigma_err_stat_low[i] = Sigma_err_stat[i];
        Sigma_err_stat_upp[i] = Sigma_err_stat[i];
        Sigma_err_syst_low[i] = Sigma_err_syst[i];
        Sigma_err_syst_upp[i] = Sigma_err_syst[i];
        Printf("Cross section in bin %i: %.5f pm %.5f(stat.) pm %.5f(syst.).", i+1, Sigma_val[i], Sigma_err_stat_low[i], Sigma_err_syst_low[i]);
        t_value[i] = t_avg_val[i];
        t_err_low[i] = t_avg_val[i] - t_boundaries[i];
        t_err_upp[i] = t_boundaries[i+1] - t_avg_val[i];
    }
    TGraphAsymmErrors *grData_stat = new TGraphAsymmErrors(nPtBins,t_value,Sigma_val,t_err_low,t_err_upp,Sigma_err_stat_low,Sigma_err_stat_upp);
    TGraphAsymmErrors *grData_syst = new TGraphAsymmErrors(nPtBins,t_value,Sigma_val,t_err_low,t_err_upp,Sigma_err_syst_low,Sigma_err_syst_upp);
    // with stat errors
    grData_stat->SetLineStyle(1);
    grData_stat->SetLineColor(kBlack);
    grData_stat->SetLineWidth(1);
    grData_stat->SetMarkerSize(1);
    grData_stat->SetMarkerStyle(8);
    grData_stat->SetMarkerColor(kBlack);
    // with syst errors 
    grData_syst->SetFillColor(17);
    grData_syst->SetMarkerSize(0);
    grData_syst->SetMarkerStyle(1);
    grData_syst->SetMarkerColor(kBlack);

    // Define graphs for incoherent predictions of HS model
    // Without subnucleonic degrees of freedom:
    TGraphErrors *gr1_inc_n = new TGraphErrors(nData1,abs_t_1,inc_n,NULL,inc_n_err);
    // https://root.cern.ch/doc/master/classTGraphErrors.html#a0f51786d0f0e210869a53ab58c0a3ffb 
    // number of points; x-values; y-values; x-errors; y-errors
    gr1_inc_n->SetMarkerStyle(20);
    gr1_inc_n->SetMarkerColor(kBlue);
    gr1_inc_n->SetLineStyle(2);
    gr1_inc_n->SetLineColor(217);
    gr1_inc_n->SetLineWidth(lineWidth);
    // With hotspots:
    TGraphErrors *gr1_inc_hs = new TGraphErrors(nData1,abs_t_1,inc_hs,NULL,inc_hs_err);
    gr1_inc_hs->SetMarkerStyle(20);
    gr1_inc_hs->SetMarkerColor(kRed);  
    gr1_inc_hs->SetLineStyle(7);
    gr1_inc_hs->SetLineColor(kRed);
    gr1_inc_hs->SetLineWidth(lineWidth);    

    // Define graphs for incoherent predictions of Guzey's model
    // First scale the values (Guzey uses nb instead of mb)
    for (Int_t i = 0; i < nData2; i++){
        sig_tot_min[i] = sig_tot_min[i] / 1e6;
        sig_tot_max[i] = sig_tot_max[i] / 1e6;
    }
    // Then fill the graph
    // https://root.cern/doc/master/graphShade_8C.html
    TGraph *gr2_min = new TGraph(nData2, abs_t_2, sig_tot_min);
    TGraph *gr2_max = new TGraph(nData2, abs_t_2, sig_tot_max);
    TGraph *gr2_area = new TGraph(2*nData2);
    for (Int_t i = 0; i < nData2; i++){
        gr2_area->SetPoint(i, abs_t_2[i], sig_tot_max[i]);
        gr2_area->SetPoint(nData2+i, abs_t_2[nData2-i-1], sig_tot_min[nData2-i-1]);
    }
    gr2_min->SetLineStyle(10);
    gr2_min->SetLineColor(8);
    gr2_min->SetLineWidth(lineWidth);
    gr2_max->SetLineStyle(10);
    gr2_max->SetLineColor(8);
    gr2_max->SetLineWidth(lineWidth);
    gr2_area->SetFillStyle(3013);
    gr2_area->SetFillColor(212);

    // Define graphs for incoherent predictions of Heikki's model
    TGraph *gr3_fluct = new TGraph(nData3, abs_t_3, sig_inc_fluct);
    TGraph *gr3_noflu = new TGraph(nData3, abs_t_3, sig_inc_noflu);
    gr3_fluct->SetLineStyle(9);
    gr3_fluct->SetLineColor(1);
    gr3_fluct->SetLineWidth(lineWidth);
    gr3_noflu->SetLineStyle(2);
    gr3_noflu->SetLineColor(222);
    gr3_noflu->SetLineWidth(lineWidth);

    // Define graphs for incoherent predictions of STARlight
    TGraph *gr4_STARlight = new TGraph(nData4, abs_t_4, sig_SL);
    gr4_STARlight->SetLineStyle(1);
    gr4_STARlight->SetLineColor(215);
    gr4_STARlight->SetLineWidth(lineWidth);

    // TStyle settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Canvas
    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetLogy();  
    //c->SetLogx();
    // Margins
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.12);
    //Plot the graphs
    TH1 *h = (TH1*) gr2_area->GetHistogram();
    h->SetTitle(";|#it{t}| (GeV^{2}); d#sigma_{#gammaPb}/d|#it{t}| (mb/GeV^{2})");
    h->SetMinimum(1e-6);
    h->SetMaximum(0.1);
    // Vertical axis
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetRangeUser(0.04,1.0);
    // Draw
    // https://root.cern.ch/doc/master/classTGraphPainter.html
    gr2_area->Draw("AF");
    grData_syst->Draw("P2 SAME");
    if(plot2){
        gr2_min->Draw("L SAME");
        gr2_max->Draw("L SAME");
    }    
    if(plot3){
        gr3_fluct->Draw("L SAME");
        gr3_noflu->Draw("L SAME");
    }
    if(plot4){
        gr4_STARlight->Draw("L SAME");
    }
    if(plot1){
        gr1_inc_hs->Draw("CX SAME");
        gr1_inc_n->Draw("CX SAME");
    }
    grData_stat->Draw("P SAME");
    // Legend
    TLegend *l = new TLegend(0.15,0.18,0.48,0.56);
    if(plot4){
        l->AddEntry(gr4_STARlight,"STARlight","L");
    }
    if(plot3){
        l->AddEntry(gr3_fluct,"MS: IPsat flu.","L");
        l->AddEntry(gr3_noflu,"MS: IPsat no flu.","L");
    }
    if(plot2){
        l->AddEntry(gr2_area,"GSZ: el. + diss.","F");
    }
    if(plot1){
        l->AddEntry(gr1_inc_hs,"CCK: GG-hs","L");
        l->AddEntry(gr1_inc_n,"CCK: GG-n","L");
    }
    l->AddEntry(grData_stat,"ALICE measurement","EP");
    l->SetTextSize(0.048);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    //gr1_inc_hs->Print();
    //gr2_area->Print();
    //gr3_fluct->Print();
    //grData_stat->Print();

    c->Print(Form("PhenoPredictions/fig/Plot_FeedDown%i_%ibins.pdf", iFeedDown, nPtBins));
    c->Print(Form("PhenoPredictions/fig/Plot_FeedDown%i_%ibins.png", iFeedDown, nPtBins));

    return;
}

void ReadInputMeasurement()
{
    // read the input file for measured cross section
    ifstream ifs;
    t_boundaries[0] = 0.04;
    TString str = Form("Results/CrossSection/%ibins_FeedDown%i_photo.txt", nPtBins, iFeedDown);
    ifs.open(str.Data());
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            // skip first line
            Int_t bin;
            Double_t tLow;
            if(i > 0) istr >> bin >> tLow >> t_boundaries[i] >> Sigma_val[i-1] >> Sigma_err_stat[i-1] >> Sigma_err_syst[i-1];
            i++;   
        }
        ifs.close();
    }
    Printf("Values of the photonuclear cross section loaded.");

    str = Form("DependenceOnT/output_%ibins.txt", nPtBins);
    ifs.open(str.Data()); 
    if(!ifs.fail()){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            Int_t bin;
            istr >> bin >> t_avg_val[i];
            i++;   
        }
        ifs.close();
    }  
    Printf("Values of an avg |t| value per bin loaded.");

    return;
}

void ReadInputHSModel()
{
    // read the input file for hot-spot model predictions
    ifstream ifs;
    ifs.open("PhenoPredictions/HSModel/data-dtdy-y_0.6-Run1.txt");
    for(Int_t i = 0; i < nData1; i++){
        Double_t x;
        ifs >> x;
        Double_t tmp;
        ifs >> tmp;
        ifs >> tmp;
        ifs >> abs_t_1[i];
        ifs >> coh_n[i];
        ifs >> coh_n_err[i];        
        ifs >> inc_n[i];
        ifs >> inc_n_err[i];        
        ifs >> coh_hs[i];
        ifs >> coh_hs_err[i];        
        ifs >> inc_hs[i];
        ifs >> inc_hs_err[i];
        //std::cout << i << " " << abs_t[i] << " " <<inc_n[i]<< " " <<inc_hs[i] << endl;
    }
    ifs.close();
    Printf("Predictions of HS model loaded.");

    return;
}

void ReadInputGuzey()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("PhenoPredictions/Guzey/incoh_tdep_nuc_run2.dat");
    for(Int_t i = 0; i < nData2; i++){
        ifs >> abs_t_2[i];
        ifs >> sig_ela_min[i];
        ifs >> sig_ela_max[i];
        ifs >> sig_dis_min[i];
        ifs >> sig_dis_max[i];
        ifs >> sig_tot_min[i];
        ifs >> sig_tot_max[i];
        //std::cout << i << " " << abs_t_2[i] << " " << sig_dis_min[i]<< " " << sig_dis_max[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Guzey's model loaded.");

    return;
}

void ReadInputHeikki()
{
    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("PhenoPredictions/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_fluct");
    for(Int_t i = 0; i < nData3; i++){
        ifs >> abs_t_3[i];
        ifs >> sig_inc_fluct[i];
        //std::cout << i << " " << abs_t_3[i] << " " << sig_inc_fluct[i] << endl;
    }
    ifs.close();
    ifs.open("PhenoPredictions/Heikki/ipsat_hight_alice_112021/no_photon_flux/incoherent_nofluct");
    for(Int_t i = 0; i < nData3; i++){
        ifs >> abs_t_3[i];
        ifs >> sig_inc_noflu[i];
        //std::cout << i << " " << abs_t_3[i] << " " << sig_inc_noflu[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Heikki's model loaded.");

    return;
}

void ReadInputSTARlight(){

    // read the input file for Guzey's model predictions
    ifstream ifs;
    ifs.open("PhenoPredictions/STARlight/inc_tDep.txt");
    for(Int_t i = 0; i < nData4; i++){
        ifs >> abs_t_4[i];
        ifs >> sig_SL[i];
        //std::cout << i << " " << abs_t_4[i] << " " << sig_SL[i] << endl;
    }
    ifs.close();
    Printf("STARlight predictions loaded.");

    return;
}
