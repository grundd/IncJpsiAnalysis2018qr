// PlotHotSpotPredictions.c
// Program to read and plot predictions of the hot-spot model for
// the |t|-dependence of J/psi production in Pb-Pb collisions at
// x = 0.000615802 (corresponding to midrapidity Run 2)

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

Double_t PhotonFlux = 84.9;
Double_t Sig_val[5] = { 0 };
Double_t Sig_err[5] = { 0 };
Double_t tAvgValues[5] = {0.054826, 0.090448, 0.152451, 0.296726, 0.616564};
Double_t tBoundaries[6] = { 0 };

// For TGraphAsymmErrors:
Double_t sig_value[5] = { 0 };
Double_t sig_errLo[5] = { 0 };
Double_t sig_errUp[5] = { 0 };
Double_t t_value[5] = { 0 };
Double_t t_errLo[5] = { 0 };
Double_t t_errUp[5] = { 0 };

void PlotHotSpotPredictions()
{
    // reserve space for data to be read
    const Int_t nData = 75;
    Double_t abs_t[nData];
    Double_t coh_n[nData];
    Double_t coh_n_err[nData];
    Double_t inc_n[nData];
    Double_t inc_n_err[nData];
    Double_t coh_hs[nData];
    Double_t coh_hs_err[nData];
    Double_t inc_hs[nData];
    Double_t inc_hs_err[nData];

    // read the input file for hot-spot model predictions
    ifstream ifs;
    ifs.open("data-dtdy-y_0.6-Run1.txt");
    for(Int_t i = 0; i < nData; i++){
        Double_t x;
        ifs >> x;
        Double_t tmp;
        ifs >> tmp;
        ifs >> tmp;
        ifs >> abs_t[i];
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
    // read the input file for measured cross section
    tBoundaries[0] = 0.04;
    ifs.open("5bins_values_plot.txt");
    if(!(ifs.fail())){
        Int_t i = 0;
        std::string str;
        while(std::getline(ifs,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream istr(str);
            // skip first line
            Int_t bin;
            Double_t tLow;
            if(i > 0) istr >> bin >> tLow >> tBoundaries[i] >> Sig_val[i-1] >> Sig_err[i-1];
            i++;   
        }
        ifs.close();
    }
    // Scale the measured results with photon flux and fill the histogram
    //TH1D *hData = new TH1D("hData","hData",5,tBoundaries);
    for(Int_t i = 0; i < 5; i++){
        sig_value[i] = Sig_val[i] / 2. / PhotonFlux;
        sig_errLo[i] = Sig_err[i] / 2. / PhotonFlux;
        sig_errUp[i] = Sig_err[i] / 2. / PhotonFlux;
        Printf("Cross section in bin %i: %.5f pm %.5f.", i+1, sig_value[i], sig_errLo[i]);
        t_value[i] = tAvgValues[i];
        t_errLo[i] = tAvgValues[i] - tBoundaries[i];
        t_errUp[i] = tBoundaries[i+1] - tAvgValues[i];
        //hData->SetBinContent(i+1, Sig_val[i]);
        //hData->SetBinError(i+1, Sig_err[i]);
    }
    TGraphAsymmErrors *grData = new TGraphAsymmErrors(5,t_value,sig_value,t_errLo,t_errUp,sig_errLo,sig_errUp);
    //hData->SetMarkerStyle(20);
    grData->SetMarkerStyle(21);
    grData->SetMarkerColor(kBlack);

    // Define graphs for incoherent
    // https://root.cern.ch/doc/master/classTGraphErrors.html#a0f51786d0f0e210869a53ab58c0a3ffb 
    // number of points; x-values; y-values; x-errors; y-errors
    // Without subnucleonic degrees of freedom:
    TGraphErrors *gr_inc_n = new TGraphErrors(nData,abs_t,inc_n,NULL,inc_n_err);
    gr_inc_n->SetMarkerStyle(20);
    gr_inc_n->SetMarkerColor(kBlue);
    gr_inc_n->SetLineStyle(2);
    gr_inc_n->SetLineColor(kBlue);
    gr_inc_n->SetLineWidth(4);
    // With hotspots:
    TGraphErrors *gr_inc_hs = new TGraphErrors(nData,abs_t,inc_hs,NULL,inc_hs_err);
    gr_inc_hs->SetMarkerStyle(20);
    gr_inc_hs->SetMarkerColor(kRed);  
    gr_inc_hs->SetLineStyle(7);
    gr_inc_hs->SetLineColor(kRed);
    gr_inc_hs->SetLineWidth(4);    
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
    //Plot the histograms
    TH1 *h = (TH1*) gr_inc_hs->GetHistogram();
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
    h->GetXaxis()->SetRangeUser(0.02,1.0);
    // Draw
    // https://root.cern.ch/doc/master/classTGraphPainter.html
    gr_inc_hs->Draw("ACX");
    gr_inc_n->Draw("CX SAME");
    //hData->Draw("E1 SAME");
    grData->Draw("P SAME");
    // Legend
    TLegend *l = new TLegend(0.2,0.2,0.5,0.4);
    l->AddEntry(gr_inc_hs,"GG-hs, incoherent","L");
    l->AddEntry(gr_inc_n,"GG-n, incoherent","L");
    //l->AddEntry(hData,"ALICE measurement","P");
    l->AddEntry(grData,"ALICE measurement","P");
    l->SetTextSize(0.05);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();

    c->Print("ComparisonWithHSModel.pdf");
    c->Print("ComparisonWithHSModel.png");

    return;
}
