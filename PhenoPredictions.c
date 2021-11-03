// PhenoPredictions.c
// David Grund, Oct 30, 2021

// c++ headers
#include <fstream>
#include <sstream> 
#include <string>   // getline

// root headers
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

// my headers
#include "AnalysisManager.h"

// 1 = HS model
// 2 = Guzey's model
// 3 = Heikki's model

Bool_t plot1 = kTRUE;
Bool_t plot2 = kTRUE;
Bool_t plot3 = kTRUE;

Double_t PhotonFlux = 84.9;
// To read the values from the files:
Double_t sig_val[nPtBins] = { 0 };
Double_t sig_err[nPtBins] = { 0 };
Double_t t_avg_val[nPtBins] = { 0 };
Double_t t_boundaries[nPtBins+1] = { 0 };

// For TGraphAsymmErrors:
Double_t sig_value[nPtBins] = { 0 };
Double_t sig_errLo[nPtBins] = { 0 };
Double_t sig_errUp[nPtBins] = { 0 };
Double_t t_value[nPtBins] = { 0 };
Double_t t_errLo[nPtBins] = { 0 };
Double_t t_errUp[nPtBins] = { 0 };

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

void ReadInputMeasurement();
void ReadInputHSModel();
void ReadInputGuzey();
void ReadInputHeikki();

void PhenoPredictions()
{
    ReadInputMeasurement();

    ReadInputHSModel();

    ReadInputGuzey();

    // Scale the measured results with photon flux and fill the histogram
    for(Int_t i = 0; i < nPtBins; i++){
        sig_value[i] = sig_val[i] / 2. / PhotonFlux;
        sig_errLo[i] = sig_err[i] / 2. / PhotonFlux;
        sig_errUp[i] = sig_err[i] / 2. / PhotonFlux;
        Printf("Cross section in bin %i: %.5f pm %.5f.", i+1, sig_value[i], sig_errLo[i]);
        t_value[i] = t_avg_val[i];
        t_errLo[i] = t_avg_val[i] - t_boundaries[i];
        t_errUp[i] = t_boundaries[i+1] - t_avg_val[i];
    }
    TGraphAsymmErrors *grData = new TGraphAsymmErrors(nPtBins,t_value,sig_value,t_errLo,t_errUp,sig_errLo,sig_errUp);
    grData->SetMarkerStyle(21);
    grData->SetMarkerColor(kBlack);

    // Define graphs for incoherent predictions of HS model
    // Without subnucleonic degrees of freedom:
    TGraphErrors *gr1_inc_n = new TGraphErrors(nData1,abs_t_1,inc_n,NULL,inc_n_err);
    // https://root.cern.ch/doc/master/classTGraphErrors.html#a0f51786d0f0e210869a53ab58c0a3ffb 
    // number of points; x-values; y-values; x-errors; y-errors
    gr1_inc_n->SetMarkerStyle(20);
    gr1_inc_n->SetMarkerColor(kBlue);
    gr1_inc_n->SetLineStyle(2);
    gr1_inc_n->SetLineColor(kBlue);
    gr1_inc_n->SetLineWidth(4);
    // With hotspots:
    TGraphErrors *gr1_inc_hs = new TGraphErrors(nData1,abs_t_1,inc_hs,NULL,inc_hs_err);
    gr1_inc_hs->SetMarkerStyle(20);
    gr1_inc_hs->SetMarkerColor(kRed);  
    gr1_inc_hs->SetLineStyle(7);
    gr1_inc_hs->SetLineColor(kRed);
    gr1_inc_hs->SetLineWidth(4);    

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
    gr2_min->SetLineWidth(3);
    gr2_max->SetLineStyle(10);
    gr2_max->SetLineColor(8);
    gr2_max->SetLineWidth(3);
    gr2_area->SetFillStyle(3013);
    gr2_area->SetFillColor(212);

    // Define graphs for incoherent predictions of Heikki's model
    // (...)

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
    if(plot2){
        gr2_min->Draw("L SAME");
        gr2_max->Draw("L SAME");
    }    
    if(plot1){
        gr1_inc_hs->Draw("CX SAME");
        gr1_inc_n->Draw("CX SAME");
    }
    if(plot3){

    }
    grData->Draw("P SAME");
    // Legend
    TLegend *l = new TLegend(0.15,0.18,0.5,0.42);
    l->AddEntry(grData,"ALICE measurement","P");
    if(plot1){
        l->AddEntry(gr1_inc_hs,"GG-hs, incoherent","L");
        l->AddEntry(gr1_inc_n,"GG-n, incoherent","L");
    }
    if(plot2){
        l->AddEntry(gr2_min,"GSZ, el+diss","L");
    }
    l->SetTextSize(0.048);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    //grData->Print();
    //gr2_area->Print();

    c->Print("PhenoPredictions/fig/ComparisonWithPheno.pdf");
    c->Print("PhenoPredictions/fig/ComparisonWithPheno.png");

    return;
}

void ReadInputMeasurement()
{
    // read the input file for measured cross section
    ifstream ifs;
    t_boundaries[0] = 0.04;
    TString sPath = Form("Results/CrossSection/%ibins_values_plot.txt", nPtBins);
    ifs.open(sPath.Data());
    if(!(ifs.fail())){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            // skip first line
            Int_t bin;
            Double_t tLow;
            if(i > 0) istr >> bin >> tLow >> t_boundaries[i] >> sig_val[i-1] >> sig_err[i-1];
            i++;   
        }
        ifs.close();
    }
    Printf("Values of cross section loaded.");

    sPath = Form("DependenceOnT/output_%ibins.txt", nPtBins);
    ifs.open(sPath.Data()); 
    if(!(ifs.fail())){
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
    Printf("Values of avg |t| per bin loaded.");

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
        //cout << i << " " << abs_t[i] << " " <<inc_n[i]<< " " <<inc_hs[i] << endl;
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
        //cout << i << " " << abs_t_2[i] << " " << sig_dis_min[i]<< " " << sig_dis_max[i] << endl;
    }
    ifs.close();
    Printf("Predictions of Guzey's model loaded.");

    return;
}

void ReadInputHeikki(){

    return;
}
