/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJpsi_DG_H
#define AliAnalysisTaskJpsi_DG_H

class TH1;
class TTree;
class TList;
class TBits;
class AliESDtrackCuts;

#include "AliTimeRangeCut.h"
#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliTOFTriggerMask;

class AliAnalysisTaskJpsi_DG : public AliAnalysisTaskSE
{
    public:
                AliAnalysisTaskJpsi_DG(); // a constructor
                AliAnalysisTaskJpsi_DG(const char *name);
        virtual ~AliAnalysisTaskJpsi_DG();	// a destructor

        virtual void    UserCreateOutputObjects(); // here one can define output objects
        virtual void    UserExec(Option_t* option);	// called for each single event
        virtual void    Terminate(Option_t* option); // usually empty, called at the end

        void    SetCrossed(Int_t spd[4], TBits &crossed);
        Int_t   GetChipId(Int_t index, Int_t &chipId2, Bool_t debug = 0);
        Bool_t  IsSTGFired(TBits bits, Int_t dphiMin = 4, Int_t dphiMax = 10, Bool_t tolerance = 1);

    private:
        AliPIDResponse  *fPIDResponse;
        AliTimeRangeCut fTimeRangeCut;
        AliESDtrackCuts *fTrackCutsBit4;

        AliVEvent   *fEvent;
        TList       *fOutputList;   //! output list
        TTree       *fTreeJPsi;     //! analysis tree
        Int_t       fRunNumber;
        TString     fTriggerName;
        // Histograms:
        TH1F    *hCounterCuts;      //! to count the number of events passing each of the cuts
        TH1F    *hCounterTrigger;   //! to count the number of events per run passing trigger conditions
        TH1F    *hVertexContrib;
        TH1F    *hVertexZ;
        TH2I    *hADdecision;
        TH2I    *hV0decision;
        TH2D    *hTPCdEdx;
        TH2D    *hTPCdEdxMuon;
        TH2D    *hTPCdEdxElectron;
        // PID, sigmas:
        Double_t fTrk1dEdx;
        Double_t fTrk2dEdx;
        Double_t fTrk1SigIfMu;
        Double_t fTrk1SigIfEl;
        Double_t fTrk2SigIfMu;
        Double_t fTrk2SigIfEl;
        // Kinematics:
        Double_t fPt;   //! transverse momentum
        Double_t fPhi;  //! azimuthal angle
        Double_t fY;    //! rapidity
        Double_t fM;    //! invariant mass
        // Two tracks:
        Double_t fPt1;  //! transverse momenta
        Double_t fPt2;
        Double_t fEta1; //! pseudorapidities
        Double_t fEta2;
        Double_t fPhi1; //! azimuthal angles
        Double_t fPhi2;
        Double_t fQ1;   //! charges
        Double_t fQ2;
        // Vertex info:
        Double_t fVertexZ;
        Int_t    fVertexContrib;
        // Info from the detectors:
        // ZDC
        Double_t fZNA_energy;
        Double_t fZNC_energy;
        // TDC = Time-to-Digital Converter
        Double_t fZNA_time[4];
        Double_t fZNC_time[4];
        // V0:
        Int_t fV0A_dec;
        Int_t fV0C_dec;
        // AD:
        Int_t fADA_dec;
        Int_t fADC_dec;
        // Matching SPD clusters with FOhits
        Bool_t fMatchingSPD;
        TBits fFOCrossFiredChips;

        AliAnalysisTaskJpsi_DG(const AliAnalysisTaskJpsi_DG&); // not implemented
        AliAnalysisTaskJpsi_DG& operator=(const AliAnalysisTaskJpsi_DG&); // not implemented

        ClassDef(AliAnalysisTaskJpsi_DG, 2);
};

#endif
