#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"



void addJetActivityVars() {
	TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");

	// get variables needed from the input tree
	float pt, eta, phi, tag_pt, tag_eta, tag_phi, jet_pt, jet_eta, jet_phi, jet_mass, JEC, R_miniIsoCone;
    tIn->SetBranchAddress("pt", &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("phi", &phi);
    tIn->SetBranchAddress("tag_pt", &tag_pt);
    tIn->SetBranchAddress("tag_eta", &tag_eta);
    tIn->SetBranchAddress("tag_phi", &tag_phi);
    tIn->SetBranchAddress("RawJetPt", &jet_pt);
    tIn->SetBranchAddress("RawJetEta", &jet_eta);
    tIn->SetBranchAddress("RawJetPhi", &jet_phi);
    tIn->SetBranchAddress("RawJetMass", &jet_mass);
    tIn->SetBranchAddress("JEC_L1L2L3Res", &JEC);
    tIn->SetBranchAddress("R_miniIsoCone", &R_miniIsoCone);

    TFile *fOut = new TFile("tnpZ_withJetAct.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);

    // initialize variables to be added as branches to the output tree
    float cleanLSjetpt, cleanLSjetmass, dR_lepjml, dR_lepjet, dPhi_lepjet, dPhi_lepjml, ptratio, dR_norm;

    tOut->Branch("JetPt", &cleanLSjetpt, "JetPt/F");
    tOut->Branch("JetMass", &cleanLSjetmass, "JetMass/F");
    tOut->Branch("dR_lepjet", &dR_lepjet, "dR_lepjet/F");
    tOut->Branch("dR_lepjml", &dR_lepjml, "dR_lepjml/F");
    tOut->Branch("dPhi_lepjet", &dPhi_lepjet, "dPhi_lepjet/F");
    tOut->Branch("dPhi_lepjml", &dPhi_lepjml, "dPhi_lepjml/F");
    tOut->Branch("LSjetlep_ptratio", &ptratio, "LSjetlep_ptratio/F");
    tOut->Branch("dR_norm", &dR_norm, "dR_norm/F");

    float dRmax = 0.4;

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();

    // variable computations go in this loop over tree entries. Recall that jet variables are for the jet closest to the probe lepton.
    int debug1count = 0;
    int debug2count = 0;
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
    	tIn->GetEntry(i);

    	TLorentzVector mu, tag_mu, nearestRawJet, nearestJet;
    	mu.SetPtEtaPhiM(pt, eta, phi, 0.105);
    	tag_mu.SetPtEtaPhiM(tag_pt, tag_eta, tag_phi, 0.105);
    	nearestRawJet.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass);

    	// do the lepton subtraction from the jet if the lepton is within the dRmax cone of the jet, and apply the JEC
    	bool nearprobe = false;
    	bool neartag = false;
    	if (nearestRawJet.DeltaR(mu) < dRmax) nearprobe = true;
    	if (nearestRawJet.DeltaR(tag_mu) < dRmax) neartag = true;

    	nearestJet = nearestRawJet;
        bool nearmu = nearprobe || neartag;
        TLorentzVector lepmom;

        if (nearprobe) lepmom += mu;
        if (neartag) lepmom += tag_mu;
        if (nearmu) {
            if ((nearestJet-lepmom).Rho() < 0.0001) {
                cleanLSjetmass = 0;
                cleanLSjetpt = 0;
                dR_lepjet = 99;
                dR_lepjml = 99;
                dPhi_lepjet = 99;
                dPhi_lepjml = 99;
                debug1count += 1;
            }
            else if (lepmom.E() >= nearestRawJet.E()) {
                cleanLSjetmass = 0;
                cleanLSjetpt = 0;
                dR_lepjet = 99;
                dR_lepjml = 99;
                dPhi_lepjet = 99;
                dPhi_lepjml = 99;
             debug2count += 1;
            } else {
                nearestJet -= lepmom;
                nearestJet *= JEC;
                cleanLSjetmass = nearestJet.M();
                cleanLSjetpt = nearestJet.Pt();
                dR_lepjet = nearestRawJet.DeltaR(mu);
                dR_lepjml = nearestJet.DeltaR(mu);
                dPhi_lepjet = nearestRawJet.DeltaPhi(mu);
                dPhi_lepjml = nearestJet.DeltaPhi(mu);
            }
        } else {
            nearestJet *= JEC;
            cleanLSjetmass = nearestJet.M();
            cleanLSjetpt = nearestJet.Pt();
            dR_lepjet = nearestRawJet.DeltaR(mu);
            dR_lepjml = nearestJet.DeltaR(mu);
            dPhi_lepjet = nearestRawJet.DeltaPhi(mu);
            dPhi_lepjml = nearestJet.DeltaPhi(mu);
        }

		if (cleanLSjetmass < 0) {
		   printf("Probe {Px, Py, Pz, E} = {%4.2f, %4.2f, %4.2f, %4.2f} --> nearprobe = %d \n", mu.Px(), mu.Py(), mu.Pz(), mu.E(), nearprobe);
		   printf("Tag {Px, Py, Pz, E} = {%4.2f, %4.2f, %4.2f, %4.2f} --> neartag = %d \n", tag_mu.Px(), tag_mu.Py(), tag_mu.Pz(), tag_mu.E(), neartag);
		   printf("Jet {Px, Py, Pz, E} = {%4.2f, %4.2f, %4.2f, %4.2f} \n", nearestRawJet.Px(), nearestRawJet.Py(), nearestRawJet.Pz(), nearestRawJet.E());
		   printf("Clean Jet {Px, Py, Pz, E} = {%4.2f, %4.2f, %4.2f, %4.2f} \n", nearestJet.Px(), nearestJet.Py(), nearestJet.Pz(), nearestJet.E());
		   printf("Clean Jet mass: %4.2f \n", cleanLSjetmass);
           printf("\n");
		}

        ptratio = cleanLSjetpt / pt;
	    dR_norm = dR_lepjml / R_miniIsoCone;


    	tOut->Fill();
    /*	if ((i+1) % step == 0) { 
            double totalTime = timer.RealTime()/60.; timer.Continue();
            double fraction = double(i+1)/double(n+1), remaining = totalTime*(1-fraction)/fraction;
            printf("Done %9d/%9d   %5.1f%%   (elapsed %5.1f min, remaining %5.1f min)\n", i, n, i*evDenom, totalTime, remaining); 
            fflush(stdout); 
        } */
    }
    printf("debug1count = %d \n", debug1count);
    printf("debug2count = %d \n", debug2count);
    tOut->AutoSave();
    fOut->Close();
}
