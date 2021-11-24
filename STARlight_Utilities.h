// STARlight_Utilities.h
// David Grund, Nov 21, 2021

// cpp headers
#include <iostream>
#include <fstream> 
#include <sstream>
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
// needed by STARlight macros
#include "TLorentzVector.h"
#include "TClonesArray.h"

Double_t fPtGm, fPtVM, fPtPm;
TLorentzVector *parent;
TClonesArray *daughters;

using namespace std;

//###############################################################################
void PrepareTreesPtGammaVMPom(Int_t nGenEv, TString str_folder){

    TTree *tPtGammaVMPom = new TTree("tPtGammaVMPom", "tPtGammaVMPom");
    tPtGammaVMPom->Branch("fPtGm", &fPtGm, "fPtGm/D");
    tPtGammaVMPom->Branch("fPtVM", &fPtVM, "fPtVM/D");
    tPtGammaVMPom->Branch("fPtPm", &fPtPm, "fPtPm/D");

    Int_t nEntriesAnalysed = 0;
    Int_t nEntriesProgress = (Double_t)nGenEv / 20.;
    Int_t nPercent = 0;

    ifstream ifs;
    ifs.open((str_folder + "PtGammaVMPom.txt").Data());
    if(!ifs.fail()){
        for(Int_t i = 0; i < nGenEv; i++){
            ifs >> fPtGm >> fPtVM >> fPtPm;
            tPtGammaVMPom->Fill();

            if((i+1) % nEntriesProgress == 0){
            nPercent += 5;
            nEntriesAnalysed += nEntriesProgress;
            Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        ifs.close();
    } else {
        Printf("File %s missing. Terminating.", (str_folder + "PtGammaVMPom.txt").Data());
        return;
    }

    TList *l = new TList();
    l->Add(tPtGammaVMPom);

    TFile *f = new TFile((str_folder + "trees_tPtGammaVMPom.root").Data(),"RECREATE");
    l->Write("TreeList", TObject::kSingleKey);
    f->ls();
    f->Close();

    Printf("\n\n");

    return;
}
//###############################################################################
void ConnectTreeVariables_tPtGammaVMPom(TTree *tPtGammaVMPom){

    tPtGammaVMPom->SetBranchAddress("fPtGm", &fPtGm);
    tPtGammaVMPom->SetBranchAddress("fPtVM", &fPtVM);
    tPtGammaVMPom->SetBranchAddress("fPtPm", &fPtPm);

    Printf("Variables from %s connected.", tPtGammaVMPom->GetName());

    return;
}
//###############################################################################
void ConnectTreeVariables_tSL(TTree *tSL){

    tSL->SetBranchAddress("parent", &parent);
    tSL->SetBranchAddress("daughters", &daughters);

    Printf("Variables from %s connected.", tSL->GetName());

    return;
}
//###############################################################################
void CompareTrees(TString str_folder){

    TFile *fSL = TFile::Open((str_folder + "trees_starlight.root").Data(), "read");
    if(fSL) Printf("Input file SL loaded.");

    TTree *tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
    if(tSL) Printf("Tree %s loaded.", tSL->GetName());
    ConnectTreeVariables_tSL(tSL);

    TFile *fPt = TFile::Open((str_folder + "trees_tPtGammaVMPom.root").Data(), "read");
    if(fPt) Printf("Input file Pt loaded.");

    TList *lPt = (TList*) fPt->Get("TreeList");
    if(lPt) Printf("List %s loaded.", lPt->GetName()); 

    TTree *tPtGammaVMPom = (TTree*)lPt->FindObject("tPtGammaVMPom");
    if(tPtGammaVMPom) Printf("Tree %s loaded.", tPtGammaVMPom->GetName());
    ConnectTreeVariables_tPtGammaVMPom(tPtGammaVMPom);


    Printf("STARlight tree contains %lli entries.", tSL->GetEntries());
    Printf("tPtGammaVMPom tree contains %lli entries.", tPtGammaVMPom->GetEntries());

    for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++){
        tSL->GetEntry(iEntry);
        tPtGammaVMPom->GetEntry(iEntry);

        Printf("%i) J/psi pT: %.5f, %.5f", iEntry+1, fPtVM, parent->Pt());
    }

    return;
}
//###############################################################################
// From STARlight:
double IDtoMass(int particleCode){
    double mass;
    if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
    else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
    else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
    else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
    else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
    else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
    else if (particleCode == 14)	{mass = 0.93827231;} // proton
    else {
        cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << endl;
        mass = -1.0;
        //exit(0); 
    } 

    return mass;
}

void ConvertStarlightAsciiToTree(TString str_folder){

	// create output tree
	TFile* outFile = new TFile((str_folder + "trees_starlight.root").Data(), "RECREATE");
	if (!outFile) {
		cerr << "    error: could not create output file '" << (str_folder + "trees_starlight.root").Data() << "'" << endl;
		return;
	}
	TTree*          outTree           = new TTree("starlightTree", "starlightTree");
	TLorentzVector* parentParticle    = new TLorentzVector();
  	TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
	outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
	outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);

	ifstream inFile;
	inFile.open((str_folder + "slight.out").Data());
	unsigned int countLines = 0;
	while (inFile.good()) {
		string       line;
		stringstream lineStream;
		
		// read EVENT
		string label;
		int    eventNmb, nmbTracks;
        // no more lines => end
		if (!getline(inFile, line))
			break;
		++countLines;
		lineStream.str(line);
		lineStream >> label >> eventNmb >> nmbTracks;
        //Printf("%s", line.data()); // DGRUND
        //cout << countLines << "\t" << eventNmb << "\t" << nmbTracks << endl; // DGRUND
		if (!(label == "EVENT:"))
			continue;
		
		// read VERTEX
        // no more lines => end
		if (!getline(inFile, line))
			break;
		++countLines;
		lineStream.str(line);
		lineStream >> label;
        //Printf("%s", line.data()); // DGRUND
		assert(label == "VERTEX:");
			
		*parentParticle = TLorentzVector(0, 0, 0, 0);
		for (int i = 0; i < nmbTracks; ++i) {
			// read tracks
			int    particleCode;
			double momentum[3];
            // no more lines => end
			if (!getline(inFile, line))
				break;
			++countLines;
			lineStream.str(line);
			lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2];
            //Printf("%s", line.data()); // DGRUND
			assert(label == "TRACK:");
			Double_t daughterMass = IDtoMass(particleCode);
			if (daughterMass < 0) {break;}
			const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
			                      + momentum[2] * momentum[2] + daughterMass * daughterMass);
			new ( (*daughterParticles)[i] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
			*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(i)));
		}
		daughterParticles->Compress();
		outTree->Fill();
	}

	outTree->Write("",TObject::kWriteDelete);
	if (outFile) {
		outFile->Close();
		delete outFile;
	}
}
//###############################################################################