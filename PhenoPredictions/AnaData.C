//
// program to read and plot predictions of the hot-spot model for
// the t-dependence of j/psi production in Pb-Pb collisions at
// x = 0.000615802 (corresponding to midrapidity Run 2
//

//-------------------------------------------------------
// Header files from C++ and ROOT
//-------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>
using namespace std;

// root headers
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>

//-------------------------------------------------------
// Main program
//-------------------------------------------------------

void AnaData()
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

    // read the input file
    ifstream ifs;
    ifs.open("data-dtdy-y_0.6-Run1.txt");
    for(Int_t i=0;i<nData;i++) {
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

    // define graph
    TGraphErrors *gr_inc_n = new TGraphErrors(nData,abs_t,inc_n,NULL,inc_n_err);
    gr_inc_n->SetMarkerStyle(20);
    gr_inc_n->SetMarkerColor(kRed);
    TGraphErrors *gr_inc_hs = new TGraphErrors(nData,abs_t,inc_hs,NULL,inc_hs_err);
    gr_inc_hs->SetMarkerStyle(20);
    gr_inc_hs->SetMarkerColor(kBlue);  
    
    // plot histo
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TCanvas *t_dep_C = new TCanvas("t_dep","t_dep",1200,800);
    TH1 *h = (TH1*) gr_inc_hs->GetHistogram();
    h->SetTitle(";|t| (GeV^{2}); d#sigma_{#gammaPb}/d|t| (mb/GeV^{2})");
    h->SetMinimum(1e-6);
    h->SetMaximum(0.1);
    h->GetXaxis()->SetTitleOffset(1.3);
    gPad->SetLogy();  
    gPad->SetLogx();
    gr_inc_hs->Draw("AP");
    gr_inc_n->Draw("p,same");

}
