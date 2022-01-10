// FeedDown_debug.c
// David Grund, Nov 10, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <string> // getline
// root headers
#include "TFile.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

TString path = Form("Results/FeedDown/debug/%ibins/", nPtBins);
TString DatasetsMC[6] = {"CohJ","IncJ","CohPCh","IncPCh","CohPNe","IncPNe"};
Double_t CrosSecSL[6] = { 0 };
Double_t BR_val[4] = {0.3468, 0.3468, 0.2538, 0.2538};
Double_t BR_err[4] = {0.0030, 0.0030, 0.0032, 0.0032};
Double_t fD_val[4] = { 0 };
Double_t fD_err[4] = { 0 };
Double_t nEvCoh = 0;
Double_t nEvInc = 0;
Double_t fD_tot_val = 0;
Double_t fD_tot_err = 0;

Double_t NRec;
Double_t NGen;
Double_t AxE;
Double_t AxE_err;

Double_t lumi_tot = 227.926; // 1/(mu barn)

// STARlight cross sections for |y| < 1.0
Double_t sig_SL_j_inc = 5.247; //mb
Double_t sig_SL_j_coh = 12.504;//mb
Double_t sig_SL_p_inc = 0.92;  //mb
Double_t sig_SL_p_coh = 2.52;  //mb
// Measured cross section for |y| < 1.0
// Michal:
Double_t sig_meas_j_coh = 8.20;//mb
Double_t sig_meas_p_coh = 1.52;//mb
// My "result":
Double_t sig_meas_j_inc = 2.36;//mb
Double_t sig_meas_p_inc = 0.41;//mb (sig_meas_j_inc / sig_SL_j_inc * sig_SL_p_inc)
// Branching ratios
Double_t BR_J = 0.05961;
Double_t BR_J_err = 0.00033;
Double_t BR_ch = 0.3468;
Double_t BR_ch_err = 0.0030;
Double_t BR_ne = 0.2538;
Double_t BR_ne_err = 0.0032;

void CalculateFeedDown(Int_t iBin, Int_t iMC);
void CalculateFeedDownAll();
void CalculateNGenTot(Int_t iMC);
void CalculateAxE(Int_t iMC, Int_t iRapCut, Int_t iPtBin);
Bool_t EventPassedMCRecLocal(Int_t iPtBin);
Bool_t EventPassedMCGenLocal(Int_t iRapCut, Int_t iPtBin);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

void FeedDown_debug(){

    Bool_t CSmeasured = kFALSE;

    if(CSmeasured){
        CrosSecSL[0] = sig_meas_j_coh;
        CrosSecSL[1] = sig_meas_j_inc;
        CrosSecSL[2] = sig_meas_p_coh;
        CrosSecSL[3] = sig_meas_p_inc;
        CrosSecSL[4] = sig_meas_p_coh;
        CrosSecSL[5] = sig_meas_p_inc;
    } else {
        CrosSecSL[0] = sig_SL_j_coh;
        CrosSecSL[1] = sig_SL_j_inc;
        CrosSecSL[2] = sig_SL_p_coh;
        CrosSecSL[3] = sig_SL_p_inc;
        CrosSecSL[4] = sig_SL_p_coh; 
        CrosSecSL[5] = sig_SL_p_inc;        
    }

    CalculateFeedDownAll();

    return;
}

void CalculateFeedDownAll(){

    Int_t iUpToBin = nPtBins;
    Int_t iUpToChannel = 4; // 1 for CohPCh, 2 for IncPCh, 3 for CohPNe, 4 for IncPNe

    SetPtBinning();

    // Calculate fD in bins
    ofstream outfile((path + "FD_coeff.txt").Data());
    outfile << std::fixed << std::setprecision(3);
    outfile << "Bin \tCohCh \terr \tIncCh \terr \tCohNe \terr \tIncNe \terr \tfD tot \terr\n";

    ofstream outTeX_fD((path + "TeX_fD.txt").Data());
    outTeX_fD << std::fixed << std::setprecision(2);

    ofstream outTeX_AxE((path + "TeX_AxE.txt").Data());
    outTeX_AxE << std::fixed << std::setprecision(2);

    ofstream out_nEv((path + "FD_nEv.txt").Data());
    out_nEv << std::fixed << std::setprecision(2);
    out_nEv << "Bin \tIncJ \tCohPCh \tIncPCh \tCohPNe \tIncPNe \tCohTot \tIncTot\n";

    ofstream outTeX_nEv((path + "TeX_nEv.txt").Data());
    outTeX_nEv << std::fixed << std::setprecision(2);

    for(Int_t iBin = 1; iBin <= iUpToBin; iBin++){
        outfile << iBin << "\t";
        out_nEv << iBin << "\t";
        // TeX files
        outTeX_fD << std::fixed << std::setprecision(3);
        outTeX_fD << "$(" << ptBoundaries[iBin-1] << "," << ptBoundaries[iBin] << ")$";
        outTeX_fD << std::fixed << std::setprecision(2);
        outTeX_AxE << std::fixed << std::setprecision(3);
        outTeX_AxE << "$(" << ptBoundaries[iBin-1] << "," << ptBoundaries[iBin] << ")$";
        outTeX_AxE << std::fixed << std::setprecision(2);
        outTeX_nEv << std::fixed << std::setprecision(3);
        outTeX_nEv << "$(" << ptBoundaries[iBin-1] << "," << ptBoundaries[iBin] << ")$";
        outTeX_nEv << std::fixed << std::setprecision(2);
        for(Int_t iChannel = 1; iChannel <= iUpToChannel; iChannel++){
            outTeX_fD << " &\t";
            outTeX_AxE << " &\t";
            outTeX_nEv << " &\t";
            Printf("******************************************************");
            Printf("Bin %i, channel %s:", iBin, DatasetsMC[iChannel+1].Data());
            Printf("******************************************************");
            Double_t NGen_tot_J, NGen_bin_J, NRec_bin_J, NGen_tot_P, NGen_bin_P, NRec_bin_P;
            Double_t AxE_J_val, AxE_J_err, AxE_P_val, AxE_P_err; 
            Double_t CrossSectionJ, CrossSectionP;
            Double_t fD_Roman_val, fD_Roman_err; // calculated using the formula from Roman's AN
            Double_t nEvJ, nEvP;
            // N_tot with |y| < 1.0 for IncJ
            CalculateNGenTot(1);
            NGen_tot_J = NGen;
            // N_bin with |y| < 0.8 for IncJ
            CalculateAxE(1, 1, iBin);
            NRec_bin_J = NRec;
            NGen_bin_J = NGen;
            // AxE for IncJ
            AxE_J_val = AxE;
            AxE_J_err = AxE_err;
            // N_tot with |y| < 1.0 for iChannel
            CalculateNGenTot(iChannel+1);
            NGen_tot_P = NGen;
            // N_bin with |y| < 0.8 for iChannel
            CalculateAxE(iChannel+1, 1, iBin);
            NRec_bin_P = NRec;
            NGen_bin_P = NGen;
            // AxE for iChannel
            AxE_P_val = AxE;
            AxE_P_err = AxE_err;
            // Calculate fD for iChannel (in percent)
            CrossSectionJ = CrosSecSL[1] * NGen_bin_J / NGen_tot_J;
            CrossSectionP = CrosSecSL[iChannel+1] * NGen_bin_P / NGen_tot_P;
            fD_val[iChannel-1] = CrossSectionP / CrossSectionJ * AxE_P_val / AxE_J_val * BR_val[iChannel-1] * 100;
            fD_err[iChannel-1] = fD_val[iChannel-1] * TMath::Sqrt(TMath::Power(AxE_J_err/AxE_J_val, 2) + TMath::Power(AxE_P_err/AxE_P_val, 2) + TMath::Power(BR_err[iChannel-1]/BR_val[iChannel-1], 2));
            fD_Roman_val = CrosSecSL[iChannel+1] / CrosSecSL[1] * NGen_tot_J / NGen_tot_P * NRec_bin_P / NRec_bin_J * BR_val[iChannel-1] * 100;
            // Calculate the predicted number of J and P events
            // multiply by 1000 to get lumi in 1/(mili barns)
            nEvJ = CrossSectionJ * AxE_J_val/100. * lumi_tot * 1000 * BR_J;
            nEvP = CrossSectionP * AxE_P_val/100. * lumi_tot * 1000 * BR_val[iChannel-1] * BR_J;
            // Print the results
            Printf("******************************************************");
            Printf("NGen J/psi bin: %.0f", NGen_bin_J);
            Printf("NGen J/psi tot: %.0f", NGen_tot_J);
            Printf("NGen Psi2S bin: %.0f", NGen_bin_P);
            Printf("NGen Psi2S tot: %.0f", NGen_tot_P);
            Printf("Cross section for %s: %.4f", DatasetsMC[1].Data(), CrossSectionJ);
            Printf("Cross section for %s: %.4f", DatasetsMC[iChannel+1].Data(), CrossSectionP);
            Printf("AxE for %s: (%.4f pm %.4f)%%", DatasetsMC[1].Data(), AxE_J_val, AxE_J_err);
            Printf("AxE for %s: (%.4f pm %.4f)%%", DatasetsMC[iChannel+1].Data(), AxE_P_val, AxE_P_err);
            Printf("fD = (%.4f pm %.4f)%%", fD_val[iChannel-1], fD_err[iChannel-1]);
            Printf("Roman's formula: fD = %.4f%%", fD_Roman_val);
            Printf("Number of %s events: %.3f", DatasetsMC[1].Data(), nEvJ);
            Printf("Number of %s events: %.3f", DatasetsMC[iChannel+1].Data(), nEvP);
            Printf("******************************************************");
            outfile << fD_val[iChannel-1] << "\t" << fD_err[iChannel-1] << "\t";
            outTeX_fD << "$" << fD_val[iChannel-1] << R"( \pm )" << fD_err[iChannel-1] << "$";
            outTeX_AxE << "$" << AxE_P_val << R"( \pm )" << AxE_P_err << "$";
            if(iChannel == 1){
                out_nEv << nEvJ << "\t" << nEvP << "\t";
                outTeX_nEv << nEvJ << " &\t" << nEvP;
            } 
            else{
                out_nEv << nEvP << "\t";
                outTeX_nEv << nEvP;
            } 
            // To calculate total fD correction per bin
            fD_tot_val += fD_val[iChannel-1];
            fD_tot_err += TMath::Power(fD_err[iChannel-1],2);
            // To calculate total fD for coh and inc channel
            if(iChannel % 2 == 1) nEvCoh += nEvP;
            if(iChannel % 2 == 0) nEvInc += nEvP;
        }
        outfile << fD_tot_val << "\t" << TMath::Sqrt(fD_tot_err) << "\n";
        outTeX_fD << " &\t$" << fD_tot_val << R"( \pm )" << TMath::Sqrt(fD_tot_err) << R"($ \\)" << "\n";
        outTeX_AxE << R"( \\)" << "\n";
        out_nEv << nEvCoh << "\t" << nEvInc << "\n";
        outTeX_nEv << " &\t" << nEvCoh << " &\t" << nEvInc << R"( \\)" << "\n";
        // Set all counters to zero
        fD_tot_val = 0;
        fD_tot_err = 0;
        nEvCoh = 0;
        nEvInc = 0;
    }
    outfile.close();
    Printf("*** Results printed to %s.***", (path + "FD_coeff.txt").Data());
    outTeX_fD.close();
    Printf("*** Results printed to %s.***", (path + "TeX_fD.txt").Data());
    outTeX_AxE.close();
    Printf("*** Results printed to %s.***", (path + "TeX_AxE.txt").Data());
    out_nEv.close();
    Printf("*** Results printed to %s.***", (path + "FD_nEv.txt").Data());
    outTeX_nEv.close();
    Printf("*** Results printed to %s.***", (path + "TeX_nEv.txt").Data());

    // Calculate fC in bins
    ofstream outfile_fC((path + "FC_coeff.txt").Data());
    outfile_fC << std::fixed << std::setprecision(3);
    outfile_fC << "Bin \tfC \terr \n";

    ofstream out_nEv_fC((path + "FC_nEv.txt").Data());
    out_nEv_fC << std::fixed << std::setprecision(2);
    out_nEv_fC << "Bin \tIncJ \tCohJ \n";

    for(Int_t iBin = 1; iBin <= iUpToBin; iBin++){
        outfile_fC << iBin << "\t";
        out_nEv_fC << iBin << "\t";
        Printf("******************************************************");
        Printf("Bin %i, channel %s:", iBin, DatasetsMC[0].Data());
        Printf("******************************************************");
        Double_t NGen_tot_coh, NGen_bin_coh, NRec_bin_coh, NGen_tot_inc, NGen_bin_inc, NRec_bin_inc;
        Double_t AxE_coh_val, AxE_coh_err, AxE_inc_val, AxE_inc_err; 
        Double_t CrossSectionCoh, CrossSectionInc;
        Double_t fC_val, fC_err;
        Double_t nEvCoh, nEvInc;
        // N_tot with |y| < 1.0 for inc J
        CalculateNGenTot(1);
        NGen_tot_inc = NGen;
        // N_bin with |y| < 0.8 for inc J
        CalculateAxE(1, 1, iBin);
        NRec_bin_inc = NRec;
        NGen_bin_inc = NGen;
        // AxE for IncJ
        AxE_inc_val = AxE;
        AxE_inc_err = AxE_err;
        // N_tot with |y| < 1.0 for coh J
        CalculateNGenTot(0);
        NGen_tot_coh = NGen;
        // N_bin with |y| < 0.8 for coh J
        CalculateAxE(0, 1, iBin);
        NRec_bin_coh = NRec;
        NGen_bin_coh = NGen;
        // AxE for iChannel
        AxE_coh_val = AxE;
        AxE_coh_err = AxE_err;
        // Calculate fD for iChannel (in percent)
        CrossSectionInc = CrosSecSL[1] * NGen_bin_inc / NGen_tot_inc;
        CrossSectionCoh = CrosSecSL[0] * NGen_bin_coh / NGen_tot_coh;
        fC_val = CrossSectionCoh / CrossSectionInc * AxE_coh_val / AxE_inc_val * 100;
        if(AxE_inc_val != 0 && AxE_coh_val != 0) fC_err = fC_val * TMath::Sqrt(TMath::Power(AxE_inc_err/AxE_inc_val, 2) + TMath::Power(AxE_coh_err/AxE_coh_val, 2));
        else fC_err = 0.;
        // Calculate the predicted number of J and P events
        // multiply by 1000 to get lumi in 1/(mili barns)
        nEvInc = CrossSectionInc * AxE_inc_val/100. * lumi_tot * 1000 * BR_J;
        nEvCoh = CrossSectionCoh * AxE_coh_val/100. * lumi_tot * 1000 * BR_J;
        // Print the results
        Printf("******************************************************");
        Printf("NGen J/psi bin: %.0f", NGen_bin_inc);
        Printf("NGen J/psi tot: %.0f", NGen_tot_inc);
        Printf("NGen Psi2S bin: %.0f", NGen_bin_coh);
        Printf("NGen Psi2S tot: %.0f", NGen_tot_coh);
        Printf("Cross section for %s: %.4f", DatasetsMC[1].Data(), CrossSectionInc);
        Printf("Cross section for %s: %.4f", DatasetsMC[0].Data(), CrossSectionCoh);
        Printf("AxE for %s: (%.4f pm %.4f)%%", DatasetsMC[1].Data(), AxE_inc_val, AxE_inc_err);
        Printf("AxE for %s: (%.4f pm %.4f)%%", DatasetsMC[0].Data(), AxE_coh_val, AxE_coh_err);
        Printf("fC = (%.4f pm %.4f)%%", fC_val, fC_err);
        Printf("Number of %s events: %.3f", DatasetsMC[1].Data(), nEvInc);
        Printf("Number of %s events: %.3f", DatasetsMC[0].Data(), nEvCoh);
        Printf("******************************************************");
        outfile_fC << fC_val << "\t" << fC_err << "\n";
        out_nEv_fC << nEvInc << "\t" << nEvCoh << "\n";
    }
    outfile_fC.close();
    Printf("*** Results printed to %s.***", (path + "FC_coeff.txt").Data());
    out_nEv_fC.close();
    Printf("*** Results printed to %s.***", (path + "FC_nEv.txt").Data());

    return;
}

void CalculateFeedDown(Int_t iBin, Int_t iMC){

    return;
}

void CalculateNGenTot(Int_t iMC){

    TString path_NGen_tot = Form("%s/NGen_tot/%s_NGen_tot.txt", path.Data(), DatasetsMC[iMC].Data());
    ifstream ifs(path_NGen_tot.Data());

    if(!ifs.fail()){

        ifs >> NGen;
        Printf("NGen_tot values loaded.");

    } else {

        TFile *f = NULL;
        switch(iMC){
            case 0:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
                break; 
            case 1:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
                break;        
            case 2:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
                break;
            case 3:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
                break;
            case 4:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi_neutral.root", "read");
                break;
            case 5:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi_neutral.root", "read");
                break;
        }
        if(f) Printf("MC file for %s loaded.", DatasetsMC[iMC].Data());

        TTree *tGen = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsiMCGen"));
        if(tGen) Printf("MC gen tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        NGen = 0;
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
            tGen->GetEntry(iEntry);
            if(EventPassedMCGenLocal(0, 0)) NGen++;
        }

        Printf("N_gen = %.0f", NGen);

        ofstream outfile(path_NGen_tot.Data());
        outfile << std::fixed << std::setprecision(0);
        outfile << NGen;
    }
    ifs.close();

    return;
}

void CalculateAxE(Int_t iMC, Int_t iRapCut, Int_t iPtBin){

    TString path_AxE = Form("%s/bin%i/%s_AxE", path.Data(), iPtBin, DatasetsMC[iMC].Data());
    ifstream ifs(path_AxE.Data());

    if(!ifs.fail()){

        ifs >> NRec;
        ifs >> NGen;
        ifs >> AxE >> AxE_err;
        Printf("AxE values loaded.");

    } else {

        TFile *f = NULL;
        switch(iMC){
            case 0:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
                break; 
            case 1:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
                break;        
            case 2:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
                break;
            case 3:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
                break;
            case 4:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi_neutral.root", "read");
                break;
            case 5:
                f = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi_neutral.root", "read");
                break;
        }
        if(f) Printf("MC file for %s loaded.", DatasetsMC[iMC].Data());

        TTree *tRec = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsiMCRec"));
        if(tRec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCRec(tRec);

        NRec = 0;
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
            tRec->GetEntry(iEntry);
            if(EventPassedMCRecLocal(iPtBin)) NRec++;
        }

        TTree *tGen = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJPsiMCGen"));
        if(tGen) Printf("MC gen tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        NGen = 0;
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
            tGen->GetEntry(iEntry);
            if(EventPassedMCGenLocal(iRapCut, iPtBin)) NGen++;
        }

        Printf("N_rec = %.0f", NRec);
        Printf("N_gen = %.0f", NGen);

        if(NGen != 0){
            AxE = NRec / NGen * 100; // in percent already
            AxE_err = CalculateErrorBayes(NRec, NGen) * 100;
        } else {
            AxE = 0;
            AxE_err = 0;
        }

        Printf("AxE = (%.4f pm %.4f)%%", AxE, AxE_err);

        ofstream outfile(path_AxE.Data());
        outfile << std::fixed << std::setprecision(0);
        outfile << NRec << "\n" << NGen << "\n";
        outfile << std::fixed << std::setprecision(2);
        outfile << AxE << "\t" << AxE_err << "\n";
    }
    ifs.close();

    return;
}

Bool_t EventPassedMCRecLocal(Int_t iPtBin){

    // Selections applied on the GRID:
    // 0) fEvent non-empty
    // 1) At least two tracks associated with the vertex
    // 2) Distance from the IP lower than 15 cm
    // 3) nGoodTracksTPC == 2 && nGoodTracksSPD == 2

    // 4) Central UPC trigger CCUP31:
    // for fRunNumber < 295881: CCUP31-B-NOPF-CENTNOTRD
    // for fRunNumber >= 295881: CCUP31-B-SPD2-CENTNOTRD
    Bool_t CCUP31 = kFALSE;
    if(
        !fTriggerInputsMC[0] &&  // !0VBA (no signal in the V0A)
        !fTriggerInputsMC[1] &&  // !0VBC (no signal in the V0C)
        !fTriggerInputsMC[2] &&  // !0UBA (no signal in the ADA)
        !fTriggerInputsMC[3] &&  // !0UBC (no signal in the ADC)
        fTriggerInputsMC[10] &&  //  0STG (SPD topological)
        fTriggerInputsMC[4]      //  0OMU (TOF two hits topology)
    ) CCUP31 = kTRUE;
    if(!CCUP31) return kFALSE;

    // 5) SPD cluster matches FOhits
    if(!(fMatchingSPD == kTRUE)) return kFALSE;

    // 6a) ADA offline veto (no effect on MC)
    if(!(fADA_dec == 0)) return kFALSE;

    // 6b) ADC offline veto (no effect on MC)
    if(!(fADC_dec == 0)) return kFALSE;

    // 7a) V0A offline veto (no effect on MC)
    if(!(fV0A_dec == 0)) return kFALSE;

    // 7b) V0C offline veto (no effect on MC)
    if(!(fV0C_dec == 0)) return kFALSE;

    // 8) Dilepton rapidity |y| < 0.8
    if(!(abs(fY) < 0.8)) return kFALSE;

    // 9) Pseudorapidity of both tracks |eta| < 0.8
    if(!(abs(fEta1) < 0.8 && abs(fEta2) < 0.8)) return kFALSE;

    // 10) Tracks have opposite charges
    if(!(fQ1 * fQ2 < 0)) return kFALSE;

    // 11) Muon pairs only
    if(!(fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu < fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl)) return kFALSE;

    // 12) Invariant mass between 2.2 and 4.5 GeV/c^2
    if(!(fM > 2.2 && fM < 4.5)) return kFALSE;

    // 13) Transverse momentum cut
    if(!(fPt > ptBoundaries[iPtBin-1] && fPt <= ptBoundaries[iPtBin])) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Bool_t EventPassedMCGenLocal(Int_t iRapCut, Int_t iPtBin){

    Double_t fYGenMax = -1.0;
    if(iRapCut == 0) fYGenMax = 1.0;
    else if(iRapCut == 1) fYGenMax = 0.8;
    else Printf("Not a valid option for iRapCut");

    // 1) Dilepton rapidity |y| < 0.8
    if(!(abs(fYGen) < fYGenMax)) return kFALSE;

    // 2) Transverse momentum cut
    if(iPtBin == 0) return kTRUE;
    else {
        if(!(fPtGen > ptBoundaries[iPtBin-1] && fPtGen <= ptBoundaries[iPtBin])) return kFALSE;
    }

    // Event passed all the selections =>
    return kTRUE;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}