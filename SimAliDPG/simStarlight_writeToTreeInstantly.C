// First in console:
// export ALIDPG_ROOT=/home/david/alice/AliDPG-master/
// echo $ALIDPG_ROOT

// -*- C++ -*-
TString comment;
TString processConfig;
TString systemConfig = "Pb-Pb";
Float_t energyConfig = 5022; // 2018 Pb-Pb, in GeV
Float_t yminConfig   = -1.0;
Float_t ymaxConfig   = +1.0;
Int_t   seedConfig   = 23659458;
// https://alimonitor.cern.ch/users/download.jsp?view=true&path=/alice/sim/2019/LHC19k1/kCohJpsiToMu/295585/001/sim.log

TString folder;
TString process_name;
Int_t kEvents; // number of events to be generated (in thousands!) 
TTree *tPtGen = NULL;
Double_t fPtGen = 0;
//Double_t fPtGen_443 = 0;

// ***************
// Options to set:
Int_t opt = 7;
// = 0 => kCohPsi2sToMuPi, 01-03-2022 (2e6 events, std RA = 6.624 fm)
// = 1 => kCohPsi2sToMuPi, 01-07-2022 (2e6 events, mod RA = 7.350 fm)
// ***************

void simStarlight_writeToTreeInstantly()
{
    // set parameters
    /*
    if(opt == 1){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohP_6.624";
        nEvTh = 2000;
    }
    if(opt == 2){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohP_7.350";
        nEvTh = 2000;
    }
    if(opt == 3){
        process_name = "kCohJpsiToMu";
        folder = "CohJ_6.624";
        nEvTh = 2000;
    }
    if(opt == 4){
        process_name = "kCohJpsiToMu";
        folder = "CohJ_7.350";
        nEvTh = 2000;
    }
    if(opt == 5){
        process_name = "kIncohJpsiToMu";
        folder = "IncJ_6.624";
        nEvTh = 2000;
    }
    if(opt == 6){
        process_name = "kIncohJpsiToMu";
        folder = "IncJ_7.350";
        nEvTh = 2000;
    }
    */
    if(opt == 7){
        process_name = "kIncohJpsiToMu";
        folder = "IncJ_1.000";
        nEvTh = 2000;
    }
    /*
    if(opt == 101){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohP_6.624_new";
        nEvTh = 2000;
    }
    if(opt == 102){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohP_7.350_new";
        nEvTh = 2000;
    }
    if(opt == 1000){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohFD_RA_6.624_10kEvents";
        nEvTh = 10;
    }
    if(opt == 1001){
        process_name = "kCohJpsiToMu";
        folder = "CohJpsi_RA_6.624_10kEvents";
        nEvTh = 10;
    }
    if(opt == 1002){
        process_name = "kCohPsi2sToMuPi";
        folder = "CohFD_RA_7.350_10kEvents";
        nEvTh = 10;
    }
    if(opt == 1003){
        process_name = "kCohJpsiToMu";
        folder = "CohJpsi_RA_7.350_10kEvents";
        nEvTh = 10;
    }
    */

    gROOT->LoadMacro("$ALIDPG_ROOT/MC/GeneratorConfig.C");

    processConfig = TString::Format("%s", process_name.Data());

    // create output text file
    ofstream outfile;
    outfile.open(folder + "/out.txt"); 
    // create output root file
    TFile fFileOut(folder + "/tPtGen.root","RECREATE");
    // create output tree
    tPtGen = new TTree("tPtGen", "tPtGen");
    tPtGen->Branch("fPtGen", &fPtGen, "fPtGen/D");
    //tPtGen->Branch("fPtGen_443", &fPtGen_443, "fPtGen_443/D");

    AliGenerator *g = GeneratorStarlight();
    g->Init();
    for(Int_t iEvTh = 1; iEvTh <= nEvTh; iEvTh++){
        Sim(g);
        outfile << iEvTh * 1e3 << " events analyzed" << endl;
    }
    delete g;

    fFileOut.Write("",TObject::kWriteDelete);
    outfile.close();
    
    return;
}

void Sim(AliGenerator *g, Int_t nEv = 1000)
{
    AliRunLoader* rl = AliRunLoader::Open("galice.root", "FASTRUN", "recreate");
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(nEv);
    rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");
    AliRun *mygAlice = new AliRun();
    mygAlice->SetRunLoader(rl);

    rl->MakeStack();
    AliStack*  stack  = rl->Stack();
    AliHeader* header = rl->GetHeader();

    //  Create and Initialize Generator
    g->SetStack(stack);

    TLorentzVector vGen;        // 4-vector of a current particle
    TLorentzVector vJpsiMuons;  // 4-vector of J/psi reconstructed from mu and antimu 

    for (Int_t iEv = 0; iEv < nEv; ++iEv){
        vJpsiMuons.SetXYZM(0.,0.,0.,0.);

        Printf("Event number %5d", iEv);

        header->Reset(0,iEv);
        rl->SetEventNumber(iEv);
        stack->Reset();
        rl->MakeTree("K");

        g->Generate();

        // Analysis
        Int_t nPart = stack->GetNprimary();
        printf("Analyse %d Particles\n", nPart);

        for (Int_t iPart = 0; iPart < nPart; iPart++){
            // get particle with the index iPart
            TParticle *part = stack->Particle(iPart);
            vGen.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
            // if muons
            if(TMath::Abs(part->GetPdgCode()) == 13) vJpsiMuons += vGen;
            // if J/psi
            //if(TMath::Abs(part->GetPdgCode()) == 443) fPtGen_443 = vGen.Pt();
        }
        // save p_T of reconstructed J/psi to the tree
        fPtGen = vJpsiMuons.Pt();
        tPtGen->Fill();

        header->SetNprimary(stack->GetNprimary());
        header->SetNtrack(stack->GetNtrack());

        stack->FinishEvent();
        header->SetStack(stack);
        rl->TreeE()->Fill();
        rl->WriteKinematics("OVERWRITE");
    } // event loop

    g->FinishRun();
    rl->WriteHeader("OVERWRITE");
    g->Write();
    rl->Write();

    delete rl;

    return;
}