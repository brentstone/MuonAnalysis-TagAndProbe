#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "MuonEffectiveArea.h"

void addNewJetActParams() {
    TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");
    float MDch_eta, MDnt_eta, MDpuch_eta, MDph_eta, MDch_pt, MDnt_pt, MDpuch_pt, MDph_pt, MDch_phi, MDnt_phi, MDpuch_phi, MDph_phi, MDch_mass, MDnt_mass, MDpuch_mass, MDph_mass;
    float SMch_eta, SMnt_eta, SMpuch_eta, SMph_eta, SMch_pt, SMnt_pt, SMpuch_pt, SMph_pt, SMch_phi, SMnt_phi, SMpuch_phi, SMph_phi, SMch_mass, SMnt_mass, SMpuch_mass, SMph_mass;
    float rho, pt, eta, phi, r_iso, tag_pt, tag_eta, tag_phi;

    tIn->SetBranchAddress("fixedGridRhoFastjetCentralNeutral",     &rho);
    tIn->SetBranchAddress("pt", &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("phi", &phi);
    tIn->SetBranchAddress("tag_pt", &tag_pt);
    tIn->SetBranchAddress("tag_eta", &tag_eta);
    tIn->SetBranchAddress("tag_phi", &tag_phi);    

    tIn->SetBranchAddress("R_miniIsoCone", &r_iso);

    tIn->SetBranchAddress("MiniToDeadConeEta_Charged", &MDch_eta);
    tIn->SetBranchAddress("MiniToDeadConeEta_Neutral", &MDnt_eta);
    tIn->SetBranchAddress("MiniToDeadConeEta_PUCharged", &MDpuch_eta);
    tIn->SetBranchAddress("MiniToDeadConeEta_Phot", &MDph_eta);

    tIn->SetBranchAddress("MiniToDeadConePt_Charged", &MDch_pt);
    tIn->SetBranchAddress("MiniToDeadConePt_Neutral", &MDnt_pt);
    tIn->SetBranchAddress("MiniToDeadConePt_PUCharged", &MDpuch_pt);
    tIn->SetBranchAddress("MiniToDeadConePt_Phot", &MDph_pt);

    tIn->SetBranchAddress("MiniToDeadConePhi_Charged", &MDch_phi);
    tIn->SetBranchAddress("MiniToDeadConePhi_Neutral", &MDnt_phi);
    tIn->SetBranchAddress("MiniToDeadConePhi_PUCharged", &MDpuch_phi);
    tIn->SetBranchAddress("MiniToDeadConePhi_Phot", &MDph_phi);

    tIn->SetBranchAddress("MiniToDeadConeMass_Charged", &MDch_mass);
    tIn->SetBranchAddress("MiniToDeadConeMass_Neutral", &MDnt_mass);
    tIn->SetBranchAddress("MiniToDeadConeMass_PUCharged", &MDpuch_mass);
    tIn->SetBranchAddress("MiniToDeadConeMass_Phot", &MDph_mass);

    tIn->SetBranchAddress("SAtoMiniConeEta_Charged", &SMch_eta);
    tIn->SetBranchAddress("SAtoMiniConeEta_Neutral", &SMnt_eta);
    tIn->SetBranchAddress("SAtoMiniConeEta_PUCharged", &SMpuch_eta);
    tIn->SetBranchAddress("SAtoMiniConeEta_Phot", &SMph_eta);

    tIn->SetBranchAddress("SAtoMiniConePt_Charged", &SMch_pt);
    tIn->SetBranchAddress("SAtoMiniConePt_Neutral", &SMnt_pt);
    tIn->SetBranchAddress("SAtoMiniConePt_PUCharged", &SMpuch_pt);
    tIn->SetBranchAddress("SAtoMiniConePt_Phot", &SMph_pt);

    tIn->SetBranchAddress("SAtoMiniConePhi_Charged", &SMch_phi);
    tIn->SetBranchAddress("SAtoMiniConePhi_Neutral", &SMnt_phi);
    tIn->SetBranchAddress("SAtoMiniConePhi_PUCharged", &SMpuch_phi);
    tIn->SetBranchAddress("SAtoMiniConePhi_Phot", &SMph_phi);

    tIn->SetBranchAddress("SAtoMiniConeMass_Charged", &SMch_mass);
    tIn->SetBranchAddress("SAtoMiniConeMass_Neutral", &SMnt_mass);
    tIn->SetBranchAddress("SAtoMiniConeMass_PUCharged", &SMpuch_mass);
    tIn->SetBranchAddress("SAtoMiniConeMass_Phot", &SMph_mass);


    TFile *fOut = new TFile("tnpZ_withNewJetVars.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Float_t dR_lepact, dR_lepactNORM, PtRatio_lepact, dR_tagprobe;
    tOut->Branch("dR_tagprobe", &dR_tagprobe, "dR_tagprobe/F");
    tOut->Branch("dR_lepact", &dR_lepact, "dR_lepact/F");
    tOut->Branch("dR_lepactNORM", &dR_lepactNORM, "dR_lepactNORM/F");
    tOut->Branch("PtRatio_lepact", &PtRatio_lepact, "PtRatio_lepact/F");

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEASpring15_25ns; // new 2015
    MuonEffectiveArea::MuonEffectiveAreaType   effAreaType   = MuonEffectiveArea::kMuMiniIso03;

    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        Float_t ea_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType, fabs(eta), effAreaTarget);

        TLorentzVector MDmom_ch, MDmom_puch, MDmom_nt, MDmom_ph, SMmom_ch, SMmom_puch, SMmom_nt, SMmom_ph, Mom_ch, Mom_puch, Mom_nt, Mom_ph;
        TLorentzVector MomForDR, probemom, tagmom;

        probemom.SetPtEtaPhiM(pt, eta, phi, 0.105);
	tagmom.SetPtEtaPhiM(tag_pt, tag_eta, tag_phi, 0.105);

        MDmom_ch.SetPtEtaPhiM(MDch_pt, MDch_eta, MDch_phi, MDch_mass);
        MDmom_puch.SetPtEtaPhiM(MDpuch_pt, MDpuch_eta, MDpuch_phi, MDpuch_mass);
        MDmom_nt.SetPtEtaPhiM(MDnt_pt, MDnt_eta, MDnt_phi, MDnt_mass);
        MDmom_ph.SetPtEtaPhiM(MDph_pt, MDph_eta, MDph_phi, MDph_mass);

        SMmom_ch.SetPtEtaPhiM(SMch_pt, SMch_eta, SMch_phi, SMch_mass);
        SMmom_puch.SetPtEtaPhiM(SMpuch_pt, SMpuch_eta, SMpuch_phi, SMpuch_mass);
        SMmom_nt.SetPtEtaPhiM(SMnt_pt, SMnt_eta, SMnt_phi, SMnt_mass);
        SMmom_ph.SetPtEtaPhiM(SMph_pt, SMph_eta, SMph_phi, SMph_mass);

        Mom_ch = MDmom_ch + SMmom_ch;
        Mom_puch = MDmom_puch + SMmom_puch;
        Mom_nt = MDmom_nt + SMmom_nt;
        Mom_ph = MDmom_ph + SMmom_ph;

        MomForDR = Mom_ch + Mom_ph + Mom_nt + Mom_puch;

	dR_tagprobe = probemom.DeltaR(tagmom);
        dR_lepact = probemom.DeltaR(MomForDR);
	dR_lepactNORM = dR_lepact / r_iso;
        PtRatio_lepact = (Mom_ch.Pt() + max(0.0, Mom_ph.Pt() + Mom_nt.Pt() - rho*ea_tot*pow(0.4/0.3, 2)))/pt;
	                
        tOut->Fill();
    }

    tOut->AutoSave(); // according to root tutorial this is the right thing to do
    fOut->Close();
}
