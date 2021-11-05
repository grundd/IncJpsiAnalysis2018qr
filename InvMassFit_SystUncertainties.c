// InvMassFit_SystUncertainties.c
// David Grund, Nov 4, 2021
// To perform fit of the invariant mass distribution of measured data

// cpp headers
#include <iostream>
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <algorithm>
// root headers
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h" // gStyle
#include "TMath.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

// Main function
void PrepareDataTree();
void FillParametersAndYields();
void DoInvMassFitMain(Int_t opt, Double_t fMCutLow, Double_t fMCutUpp, Double_t fAlpha_L, Double_t fAlpha_R, Double_t fN_L, Double_t fN_R, TString path);
// Support functions
void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit);
void SetCanvas(TCanvas *c, Bool_t bLogScale);

Double_t N_Jpsi_out[2];
Double_t N_bkgr_out[2];
Double_t N_Jpsi_out2[2];

Double_t fAlpha_L_val[nPtBins] = { 0 };
Double_t fAlpha_L_err[nPtBins] = { 0 };
Double_t fAlpha_R_val[nPtBins] = { 0 };
Double_t fAlpha_R_err[nPtBins] = { 0 };
Double_t fN_L_val[nPtBins] = { 0 };
Double_t fN_L_err[nPtBins] = { 0 };
Double_t fN_R_val[nPtBins] = { 0 };
Double_t fN_R_err[nPtBins] = { 0 };

Double_t fN_yield_val[nPtBins] = { 0 };
Double_t fN_yield_err[nPtBins] = { 0 };

Bool_t low = kFALSE;
Bool_t upp = kFALSE;
Bool_t alphL = kFALSE;
Bool_t alphR = kFALSE;
Bool_t nL = kFALSE;
Bool_t nR = kFALSE;

void FillParametersAndYields(){

    TString pathMC = Form("Results/InvMassFitMC/%ibins/", nPtBins); 
    TString pathData = Form("Results/InvMassFit/%ibins/", nPtBins);
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        ifstream ifs;
        ifs.open(Form("%sbin%i.txt", pathMC.Data(), iBin+1));
        char name[8];
        for(Int_t i = 0; i < 4; i++){
            if(i == 0) ifs >> name >> fAlpha_L_val[iBin] >> fAlpha_L_err[iBin];
            if(i == 1) ifs >> name >> fAlpha_R_val[iBin] >> fAlpha_R_err[iBin];
            if(i == 2) ifs >> name >> fN_L_val[iBin] >> fN_L_err[iBin];
            if(i == 3) ifs >> name >> fN_R_val[iBin] >> fN_R_err[iBin];
        }
        ifs.close();
        ifs.open(Form("%sbin%i_signal.txt", pathData.Data(), iBin+1));
        ifs >> fN_yield_val[iBin] >> fN_yield_err[iBin];
        ifs.close();
        Printf("Values for bin %i loaded.", iBin+1);
    }

    // Print the results
    TString path_output = Form("Results/InvMassFit_SystUncertainties/%ibins/input_values_%ibins.txt", nPtBins, nPtBins);
    ofstream outfile(path_output.Data());
    outfile << "bin \talpha_L\terr\talpha_R\terr\tn_L\terr\tn_R\terr\tYield\terr\n";
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        outfile << iBin+1 << "\t" 
                << std::fixed << std::setprecision(4)
                << fAlpha_L_val[iBin] << "\t" 
                << fAlpha_L_err[iBin] << "\t" 
                << fAlpha_R_val[iBin] << "\t" 
                << fAlpha_R_err[iBin] << "\t" 
                << fN_L_val[iBin] << "\t" 
                << fN_L_err[iBin] << "\t" 
                << fN_R_val[iBin] << "\t" 
                << fN_R_err[iBin] << "\t"
                << std::fixed << std::setprecision(2)
                << fN_yield_val[iBin] << "\t"
                << fN_yield_err[iBin] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", path_output.Data()); 

    return;
}

void InvMassFit_SystUncertainties(){

    //PrepareDataTree();

    SetPtBinning();

    FillParametersAndYields();

    // 1) vary the lower boundary
    const Int_t n1 = 6;
    TString str1 = Form("Results/InvMassFit_SystUncertainties/%ibins/boundary_low/", nPtBins);
    Double_t signal_1_val[nPtBins][n1] = { 0 };
    Double_t signal_1_err[nPtBins][n1] = { 0 };
    Double_t signal_1_devAbs[nPtBins][n1] = { 0 };
    Double_t signal_1_devRel[nPtBins][n1] = { 0 };
    if(low){
        Double_t low_bound[n1] = {2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
        TString str1_all[n1];
        for(Int_t iVar = 0; iVar < n1; iVar++){
            str1_all[iVar] = str1 + Form("mass_%.1f/", low_bound[iVar]);
            // here add the inv. mass fit of allbins
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                DoInvMassFitMain(iBin+4,low_bound[iVar],4.5,fAlpha_L_val[iBin],fAlpha_R_val[iBin], fN_L_val[iBin], fN_R_val[iBin], str1_all[iVar]);
                signal_1_val[iBin][iVar] = N_Jpsi_out[0];
                signal_1_err[iBin][iVar] = N_Jpsi_out[1];
                signal_1_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_1_val[iBin][iVar]);
                signal_1_devRel[iBin][iVar] = signal_1_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str1 + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n1; iVar++) outfile << Form("low:%.1f\t", low_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n1; iVar++){
                outfile << signal_1_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n1; iVar++) outfile << Form("low:%.1f\t", low_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n1; iVar++){
                outfile << signal_1_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n1; iVar++) outfile << Form("low:%.1f\t", low_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n1; iVar++){
                outfile << signal_1_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n1; iVar++) outfile << Form("low:%.1f\t", low_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n1; iVar++){
                outfile << signal_1_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str1 + "output.txt").Data()); 
    }

    // 2) vary the upper boundary
    const Int_t n2 = 7;
    TString str2 = Form("Results/InvMassFit_SystUncertainties/%ibins/boundary_upp/", nPtBins);
    Double_t signal_2_val[nPtBins][n2] = { 0 };
    Double_t signal_2_err[nPtBins][n2] = { 0 };
    Double_t signal_2_devAbs[nPtBins][n2] = { 0 };
    Double_t signal_2_devRel[nPtBins][n2] = { 0 };
    if(upp){
        Double_t upp_bound[n2] = {4.0, 4.2, 4.4, 4.5, 4.6, 4.8, 5.0};
        TString str2_all[n2];
        for(Int_t iVar = 0; iVar < n2; iVar++){
            str2_all[iVar] = str2 + Form("mass_%.1f/", upp_bound[iVar]);
            // here add the inv. mass fit of allbins
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                DoInvMassFitMain(iBin+4,2.2,upp_bound[iVar],fAlpha_L_val[iBin],fAlpha_R_val[iBin], fN_L_val[iBin], fN_R_val[iBin], str2_all[iVar]);
                signal_2_val[iBin][iVar] = N_Jpsi_out[0];
                signal_2_err[iBin][iVar] = N_Jpsi_out[1];
                signal_2_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_2_val[iBin][iVar]);
                signal_2_devRel[iBin][iVar] = signal_2_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str2 + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n2; iVar++) outfile << Form("upp:%.1f\t", upp_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n2; iVar++){
                outfile << signal_2_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n2; iVar++) outfile << Form("upp:%.1f\t", upp_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n2; iVar++){
                outfile << signal_2_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n2; iVar++) outfile << Form("upp:%.1f\t", upp_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n2; iVar++){
                outfile << signal_2_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n2; iVar++) outfile << Form("upp:%.1f\t", upp_bound[iVar]);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n2; iVar++){
                outfile << signal_2_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str2 + "output.txt").Data()); 
    }

    // 3a) vary the alpha_L parameter
    const Int_t n3a = 6;
    TString str3a = Form("Results/InvMassFit_SystUncertainties/%ibins/alpha_L/", nPtBins);
    Double_t signal_3a_val[nPtBins][n3a] = { 0 };
    Double_t signal_3a_err[nPtBins][n3a] = { 0 };
    Double_t signal_3a_devAbs[nPtBins][n3a] = { 0 };
    Double_t signal_3a_devRel[nPtBins][n3a] = { 0 };
    if(alphL){
        Double_t alpha_L[nPtBins][n3a] = { 0 };
        TString str3a_all[nPtBins][n3a];
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                str3a_all[iBin][iVar] = str3a + Form("var_%i/", iVar+1);
                if(iVar == 0) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] - fAlpha_L_err[iBin];
                if(iVar == 1) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] - 2./3. * fAlpha_L_err[iBin];
                if(iVar == 2) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] - 1./3. * fAlpha_L_err[iBin];
                if(iVar == 3) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] + 1./3. * fAlpha_L_err[iBin];
                if(iVar == 4) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] + 2./3. * fAlpha_L_err[iBin];
                if(iVar == 5) alpha_L[iBin][iVar] = fAlpha_L_val[iBin] + fAlpha_L_err[iBin];
            }
        }
        // Print the values
        Bool_t debug = kTRUE;
        if(debug){
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                std::cout << iBin << "\t";
                for(Int_t iVar = 0; iVar < n3a; iVar++){
                   std::cout << alpha_L[iBin][iVar] << "\t";
                }
                std::cout << "\n";
            }
        }
        // Do invariant mass fits
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                DoInvMassFitMain(iBin+4,2.2,4.5,alpha_L[iBin][iVar],fAlpha_R_val[iBin], fN_L_val[iBin], fN_R_val[iBin], str3a_all[iBin][iVar]);    
                signal_3a_val[iBin][iVar] = N_Jpsi_out[0];
                signal_3a_err[iBin][iVar] = N_Jpsi_out[1];
                signal_3a_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_3a_val[iBin][iVar]);
                signal_3a_devRel[iBin][iVar] = signal_3a_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str3a + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "parameters val (alpha_L) \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3a; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                outfile << alpha_L[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3a; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                outfile << signal_3a_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3a; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                outfile << signal_3a_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3a; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                outfile << signal_3a_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3a; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3a; iVar++){
                outfile << signal_3a_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str3a + "output.txt").Data()); 
    }

    // 3b) vary the alpha_R parameter
    const Int_t n3b = 6;
    TString str3b = Form("Results/InvMassFit_SystUncertainties/%ibins/alpha_R/", nPtBins);
    Double_t signal_3b_val[nPtBins][n3b] = { 0 };
    Double_t signal_3b_err[nPtBins][n3b] = { 0 };
    Double_t signal_3b_devAbs[nPtBins][n3b] = { 0 };
    Double_t signal_3b_devRel[nPtBins][n3b] = { 0 };
    if(alphR){
        Double_t alpha_R[nPtBins][n3b] = { 0 };
        TString str3b_all[nPtBins][n3b];
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                str3b_all[iBin][iVar] = str3b + Form("var_%i/", iVar+1);
                if(iVar == 0) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] - fAlpha_R_err[iBin];
                if(iVar == 1) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] - 2./3. * fAlpha_R_err[iBin];
                if(iVar == 2) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] - 1./3. * fAlpha_R_err[iBin];
                if(iVar == 3) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] + 1./3. * fAlpha_R_err[iBin];
                if(iVar == 4) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] + 2./3. * fAlpha_R_err[iBin];
                if(iVar == 5) alpha_R[iBin][iVar] = fAlpha_R_val[iBin] + fAlpha_R_err[iBin];
            }
        }
        // Print the values
        Bool_t debug = kTRUE;
        if(debug){
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                std::cout << iBin << "\t";
                for(Int_t iVar = 0; iVar < n3b; iVar++){
                   std::cout << alpha_R[iBin][iVar] << "\t";
                }
                std::cout << "\n";
            }
        }
        // Do invariant mass fits
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                DoInvMassFitMain(iBin+4,2.2,4.5,fAlpha_L_val[iBin], alpha_R[iBin][iVar], fN_L_val[iBin], fN_R_val[iBin], str3b_all[iBin][iVar]);    
                signal_3b_val[iBin][iVar] = N_Jpsi_out[0];
                signal_3b_err[iBin][iVar] = N_Jpsi_out[1];
                signal_3b_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_3b_val[iBin][iVar]);
                signal_3b_devRel[iBin][iVar] = signal_3b_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str3b + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "parameters val (alpha_R) \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3b; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                outfile << alpha_R[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3b; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                outfile << signal_3b_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3b; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                outfile << signal_3b_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3b; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                outfile << signal_3b_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3b; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3b; iVar++){
                outfile << signal_3b_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str3b + "output.txt").Data()); 
    }

    // 3c) vary the n_L parameter
    const Int_t n3c = 6;
    TString str3c = Form("Results/InvMassFit_SystUncertainties/%ibins/n_L/", nPtBins);
    Double_t signal_3c_val[nPtBins][n3c] = { 0 };
    Double_t signal_3c_err[nPtBins][n3c] = { 0 };
    Double_t signal_3c_devAbs[nPtBins][n3c] = { 0 };
    Double_t signal_3c_devRel[nPtBins][n3c] = { 0 };
    if(nL){
        Double_t n_L[nPtBins][n3c] = { 0 };
        TString str3c_all[nPtBins][n3c];
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                str3c_all[iBin][iVar] = str3c + Form("var_%i/", iVar+1);
                if(iVar == 0) n_L[iBin][iVar] = fN_L_val[iBin] - fN_L_err[iBin];
                if(iVar == 1) n_L[iBin][iVar] = fN_L_val[iBin] - 2./3. * fN_L_err[iBin];
                if(iVar == 2) n_L[iBin][iVar] = fN_L_val[iBin] - 1./3. * fN_L_err[iBin];
                if(iVar == 3) n_L[iBin][iVar] = fN_L_val[iBin] + 1./3. * fN_L_err[iBin];
                if(iVar == 4) n_L[iBin][iVar] = fN_L_val[iBin] + 2./3. * fN_L_err[iBin];
                if(iVar == 5) n_L[iBin][iVar] = fN_L_val[iBin] + fN_L_err[iBin];
            }
        }
        // Print the values
        Bool_t debug = kTRUE;
        if(debug){
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                std::cout << iBin << "\t";
                for(Int_t iVar = 0; iVar < n3c; iVar++){
                   std::cout << n_L[iBin][iVar] << "\t";
                }
                std::cout << "\n";
            }
        }
        // Do invariant mass fits
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                DoInvMassFitMain(iBin+4,2.2,4.5,fAlpha_L_val[iBin], fAlpha_R_val[iBin], n_L[iBin][iVar], fN_R_val[iBin], str3c_all[iBin][iVar]);    
                signal_3c_val[iBin][iVar] = N_Jpsi_out[0];
                signal_3c_err[iBin][iVar] = N_Jpsi_out[1];
                signal_3c_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_3c_val[iBin][iVar]);
                signal_3c_devRel[iBin][iVar] = signal_3c_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str3c + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "parameters val (n_L) \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3c; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                outfile << n_L[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3c; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                outfile << signal_3c_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3c; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                outfile << signal_3c_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3c; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                outfile << signal_3c_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3c; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3c; iVar++){
                outfile << signal_3c_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str3c + "output.txt").Data()); 
    }

    // 3d) vary the n_R parameter
    const Int_t n3d = 6;
    TString str3d = Form("Results/InvMassFit_SystUncertainties/%ibins/n_R/", nPtBins);
    Double_t signal_3d_val[nPtBins][n3d] = { 0 };
    Double_t signal_3d_err[nPtBins][n3d] = { 0 };
    Double_t signal_3d_devAbs[nPtBins][n3d] = { 0 };
    Double_t signal_3d_devRel[nPtBins][n3d] = { 0 };
    if(nR){
        Double_t n_R[nPtBins][n3d] = { 0 };
        TString str3d_all[nPtBins][n3d];
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                str3d_all[iBin][iVar] = str3d + Form("var_%i/", iVar+1);
                if(iVar == 0) n_R[iBin][iVar] = fN_R_val[iBin] - fN_R_err[iBin];
                if(iVar == 1) n_R[iBin][iVar] = fN_R_val[iBin] - 2./3. * fN_R_err[iBin];
                if(iVar == 2) n_R[iBin][iVar] = fN_R_val[iBin] - 1./3. * fN_R_err[iBin];
                if(iVar == 3) n_R[iBin][iVar] = fN_R_val[iBin] + 1./3. * fN_R_err[iBin];
                if(iVar == 4) n_R[iBin][iVar] = fN_R_val[iBin] + 2./3. * fN_R_err[iBin];
                if(iVar == 5) n_R[iBin][iVar] = fN_R_val[iBin] + fN_R_err[iBin];
            }
        }
        // Print the values
        Bool_t debug = kTRUE;
        if(debug){
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            for(Int_t iBin = 0; iBin < nPtBins; iBin++){
                std::cout << iBin << "\t";
                for(Int_t iVar = 0; iVar < n3d; iVar++){
                   std::cout << n_R[iBin][iVar] << "\t";
                }
                std::cout << "\n";
            }
        }
        // Do invariant mass fits
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                DoInvMassFitMain(iBin+4,2.2,4.5,fAlpha_L_val[iBin], fAlpha_R_val[iBin], fN_L_val[iBin], n_R[iBin][iVar], str3d_all[iBin][iVar]);    
                signal_3d_val[iBin][iVar] = N_Jpsi_out[0];
                signal_3d_err[iBin][iVar] = N_Jpsi_out[1];
                signal_3d_devAbs[iBin][iVar] = TMath::Abs(fN_yield_val[iBin] - signal_3d_val[iBin][iVar]);
                signal_3d_devRel[iBin][iVar] = signal_3d_devAbs[iBin][iVar] / fN_yield_val[iBin] * 100.;
            }
        }
        // Print the output to txt file
        ofstream outfile((str3d + "output.txt").Data());
        outfile << std::fixed << std::setprecision(3);
        outfile << "parameters val (n_R) \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3d; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                outfile << n_R[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal val \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3d; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                outfile << signal_3d_val[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal err \n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3d; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                outfile << signal_3d_err[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << "signal deviation absolute\n";
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3d; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                outfile << signal_3d_devAbs[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile << Form("signal deviation relative [%%]\n");
        outfile << "bin \t";
        for(Int_t iVar = 0; iVar < n3d; iVar++) outfile << Form("var_%i\t", iVar+1);
        outfile << "\n";
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            outfile << iBin+1 << "\t";
            for(Int_t iVar = 0; iVar < n3d; iVar++){
                outfile << signal_3d_devRel[iBin][iVar] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
        Printf("Results printed to %s.", (str3d + "output.txt").Data()); 
    }

    // Calculate the systematic uncertainty for each bin and each source as a maximum relative deviation over bins
    Double_t SystUncr_low[nPtBins];
    Double_t SystUncr_upp[nPtBins];
    Double_t SystUncr_aL[nPtBins];
    Double_t SystUncr_aR[nPtBins];
    Double_t SystUncr_nL[nPtBins];
    Double_t SystUncr_nR[nPtBins];
    Double_t SystUncr_tot[nPtBins];
    // Find the maximum values
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Double_t arr1_row[n1];
        for(Int_t i = 0; i < n1; i++) arr1_row[i] = signal_1_devRel[iBin][i];
        SystUncr_low[iBin] = *max_element(arr1_row, arr1_row + n1);
        Double_t arr2_row[n2];
        for(Int_t i = 0; i < n2; i++) arr2_row[i] = signal_2_devRel[iBin][i];
        SystUncr_upp[iBin] = *max_element(arr2_row, arr2_row + n2);
        Double_t arr3a_row[n3a];
        for(Int_t i = 0; i < n3a; i++) arr3a_row[i] = signal_3a_devRel[iBin][i];
        SystUncr_aL[iBin] = *max_element(arr3a_row, arr3a_row + n3a);
        Double_t arr3b_row[n3b];
        for(Int_t i = 0; i < n3b; i++) arr3b_row[i] = signal_3b_devRel[iBin][i];
        SystUncr_aR[iBin] = *max_element(arr3b_row, arr3b_row + n3b);
        Double_t arr3c_row[n3c];
        for(Int_t i = 0; i < n3c; i++) arr3c_row[i] = signal_3c_devRel[iBin][i];
        SystUncr_nL[iBin] = *max_element(arr3c_row, arr3c_row + n3c);
        Double_t arr3d_row[n3d];
        for(Int_t i = 0; i < n3d; i++) arr3d_row[i] = signal_3d_devRel[iBin][i];
        SystUncr_nR[iBin] = *max_element(arr3d_row, arr3d_row + n3d);
    }
    ofstream outfile(Form("Results/InvMassFit_SystUncertainties/%ibins/syst_uncertainties_%ibins.txt", nPtBins, nPtBins));
    outfile << std::fixed << std::setprecision(3);
    outfile << "systematic uncertainties per bin [%]\n";
    outfile << "bin \tlow \tupp \taL \taR \tnL \tnR \ttotal \n";
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        SystUncr_tot[iBin] = TMath::Sqrt(TMath::Power(SystUncr_low[iBin], 2) 
                                        + TMath::Power(SystUncr_upp[iBin], 2) 
                                        + TMath::Power(SystUncr_aL[iBin], 2) 
                                        + TMath::Power(SystUncr_aR[iBin], 2) 
                                        + TMath::Power(SystUncr_nL[iBin], 2) 
                                        + TMath::Power(SystUncr_nR[iBin], 2));
        outfile << iBin+1 << "\t" 
                << SystUncr_low[iBin] << "\t" 
                << SystUncr_upp[iBin] << "\t" 
                << SystUncr_aL[iBin] << "\t" 
                << SystUncr_aR[iBin] << "\t" 
                << SystUncr_nL[iBin] << "\t" 
                << SystUncr_nR[iBin] << "\t"
                << SystUncr_tot[iBin] << "\n"; 
    }
    outfile.close();
    Printf("Results printed to the output file."); 

    Printf("Done.");

    return;
}

void DoInvMassFitMain(Int_t opt, Double_t fMCutLow, Double_t fMCutUpp, Double_t fAlpha_L, Double_t fAlpha_R, Double_t fN_L, Double_t fN_R, TString str){
    // Fit the invariant mass distribution using Double-sided CB function
    // Fix the values of the tail parameters to MC values
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fPtCut     = -999;
    Double_t fPtCutLow  = -999;
    Double_t fPtCutUpp  = -999;
    Double_t fYCut      = 0.80;

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
        case 3: // Sample with pt from 0.2 to 1 GeV/c 
            fPtCutLow = 0.20;
            fPtCutUpp = 1.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 4: // pt bin 1
            fPtCutLow = ptBoundaries[0];
            fPtCutUpp = ptBoundaries[1];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 5: // pt bin 2
            fPtCutLow = ptBoundaries[1];
            fPtCutUpp = ptBoundaries[2];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 6: // pt bin 3
            fPtCutLow = ptBoundaries[2];
            fPtCutUpp = ptBoundaries[3];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 7: // pt bin 4
            fPtCutLow = ptBoundaries[3];
            fPtCutUpp = ptBoundaries[4];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 8:
            fPtCutLow = ptBoundaries[4];
            fPtCutUpp = ptBoundaries[5];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
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
    TFile *fFileIn = new TFile("Trees/InvMassFit_SystUncertainties/InvMassFit_SystUncertainties.root"); 
    TTree *fTreeIn = NULL;
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
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
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,30.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,30.);
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

    // Create model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of J/psi events
    fM.setRange("WholeMassRange",fMCutLow,fMCutUpp);
    RooAbsReal *iDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("WholeMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_out[0] = iDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = iDSCB->getVal()*N_Jpsi.getError();

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
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","data",fResFit->floatParsFinal().getSize());

    // -------------------------------------------------------------------------------- 
    // Legend1
    TLegend *l1 = new TLegend(0.09,0.76,0.3,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    if(opt == 0){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} > %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 1 || opt == 2){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} < %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    }
    l1->SetTextSize(0.042);
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
    // Incoherent: lower precision:
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
        l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    // No incoherent: higher precision:
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
    l2->SetTextSize(0.042);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    TString *str_local = NULL;
    if(opt == 3) str_local = new TString((str + "allbins").Data());
    else if(opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8) str_local = new TString(Form("%sbin%i", str.Data(), opt-3));

    // Print the plots
    cHist->Print((*str_local + ".pdf").Data());
    cHist->Print((*str_local + ".png").Data());
    //cHist->Print((*str_local + "_cm.pdf").Data());
    //cHist->Print((*str_local + "_cm.png").Data());

    // Calculate the number of signal and bkg events with mass in 3.0 to 3.2 GeV/c^2
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *iBkg = BkgPdf.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_bkgr_out[0] = iBkg->getVal()*N_bkg.getVal();
    N_bkgr_out[1] = iBkg->getVal()*N_bkg.getError();
    RooAbsReal *iDSCB2 = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    N_Jpsi_out2[0] = iDSCB2->getVal()*N_Jpsi.getVal();
    N_Jpsi_out2[1] = iDSCB2->getVal()*N_Jpsi.getError();

    // Print the number of events to text file
    ofstream outfile((*str_local + ".txt").Data());
    outfile << "The whole mass region:" << endl;
    outfile << "N_J/psi:\t" << N_Jpsi_out[0] << " pm " << N_Jpsi_out[1] << endl;
    outfile << "J/psi peak (3.0 < m < 3.2 GeV):" << endl;
    outfile << "N_J/psi:\t" << N_Jpsi_out2[0] << " pm " << N_Jpsi_out2[1] << endl;    
    outfile << "N_bkg:  \t" << N_bkgr_out[0] << " pm " << N_bkgr_out[1] << endl;
    outfile.close();
    Printf("*** Results printed to %s.***", (*str_local + ".txt").Data());

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
    TFile fFileOut("Trees/InvMassFit_SystUncertainties/InvMassFit_SystUncertainties.root","RECREATE");

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

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(2,0)) tIncEnrSample->Fill();
        if(EventPassed(2,1)) tCohEnrSample->Fill();
        if(EventPassed(2,2)) tMixedSample->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}