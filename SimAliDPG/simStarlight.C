// -*- C++ -*-
TString comment;
TString processConfig;
TString systemConfig = "Pb-Pb";
Float_t energyConfig = 5022; // 2018 Pb-Pb, in GeV
Float_t yminConfig   = -1.0;
Float_t ymaxConfig   = +1.0;
Int_t   seedConfig   = 12345;

// ***************
// Options to set:
Int_t opt = 2;
// = 0 => kCohPsi2sToMuPi, 12-30-2021 (1000 events, std RA = 6.602 fm)
// = 1 => kCohPsi2sToMuPi, 01-03-2022 (1e4 events, std RA = 6.602 fm)
// ***************

void simStarlight()
{
    TString process_name = "";
    Int_t nEv = 0;
    // set parameters
    if(opt == 0){
        process_name = "kCohPsi2sToMuPi";
        nEv = 1000;
    } 
    if(opt == 1){
        process_name = "kCohPsi2sToMuPi";
        nEv = 10;        
    } 
    if(opt == 2){
        process_name = "kCohJpsiToMu";
        nEv = 10;        
    } 

    gSystem->Exec("export ALIDPG_ROOT=/home/david/alice/AliDPG-master/");
    //gSystem->Exec("echo $ALIDPG_ROOT");

    gROOT->LoadMacro("$ALIDPG_ROOT/MC/GeneratorConfig.C");

    processConfig = TString::Format("%s", process_name.Data());
    
    AliGenerator *g = GeneratorStarlight();
    Sim(g, nEv);
    delete g;
}

void Sim(AliGenerator *g, Int_t nev)
{
  TCanvas *c1    = new TCanvas;
  TH2 *hMassPt   = new TH2D("hMassPt", processConfig+";M_{inv} (GeV/#it{c}^{2});p_{T}(mother) (GeV/#it{c})",
			   1000,0,10,200,0,2);
  TH1 *hRapidity = new TH1D("hRapidity", processConfig+";Rapidity(mother)", 100, -10, 10);
  TH2 *hEta12    = new TH2D("hEta12",    processConfig+";#eta_{1};#eta_{2}", 100, -10, 10, 100,-10,10);
  TH2 *hRap12    = new TH2D("hRap12",    processConfig+";Y_{1};Y_{2}",       100, -10, 10, 100,-10,10);

  AliRunLoader* rl = AliRunLoader::Open("galice.root", "FASTRUN", "recreate");
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  AliRun *mygAlice = new AliRun();
  mygAlice->SetRunLoader(rl);

  rl->MakeStack();
  AliStack*  stack  = rl->Stack();
  AliHeader* header = rl->GetHeader();

  //  Create and Initialize Generator
  g->Init();
  g->SetStack(stack);

  for (Int_t iev=0; iev<nev; ++iev) {
    Printf("Event number %5d", iev);

    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    //      stack->ConnectTree();

    g->Generate();
    //      Analysis
    Int_t npart = stack->GetNprimary();
    printf("Analyse %d Particles\n", npart);
    TLorentzVector v, vi[10];
    Float_t eta[2], rap[2];
    for (Int_t part=0; part<npart; part++) {
      TParticle *MPart = stack->Particle(part);
      Int_t mpart  = MPart->GetPdgCode();
      MPart->Momentum(vi[part]);
      if (part < 2) {
	eta[part] = vi[part].Eta();
	rap[part] = vi[part].Rapidity();
      }
      v += vi[part];
    }
    hMassPt->Fill(v.M(), v.Perp());
    hRapidity->Fill(v.Rapidity());
    hEta12->Fill(eta[0], eta[1]);
    hRap12->Fill(rap[0], rap[1]);

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

  c1->Divide(2,2);
  c1->cd(1);
  hMassPt->Draw("COLZ");
  c1->cd(2);
  hRapidity->Draw();
  c1->cd(3)->SetGrid();
  hEta12->Draw("COLZ");
  c1->cd(4)->SetGrid();
  hRap12->Draw("COLZ");
  c1->SaveAs("QA.pdf");
  c1->SaveAs("QA.root");

  gSystem->Exec("mkdir -p " + processConfig);
  gSystem->Exec("mv slight.txt galice.root Kinematics.root QA*.{root,pdf} log.txt " + processConfig);

  delete hMassPt;
  delete hRapidity;
  delete c1;
  delete rl;
}