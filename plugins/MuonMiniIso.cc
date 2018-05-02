// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TLorentzVector.h"
//
// class declaration
//

class MuonMiniIso : public edm::EDProducer {
public:

  typedef std::vector< edm::FwdPtr<reco::PFCandidate> > PFCollection;

  explicit MuonMiniIso(const edm::ParameterSet&);
  ~MuonMiniIso();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const edm::EDGetTokenT<edm::View<reco::Muon>> probes_;    
  const edm::EDGetTokenT<PFCollection> pfCandidates_;
  double dRCandProbeVeto_;
  double dRCandSoftActivityCone_;
  double CandPtThreshold_;
  double ChargedPVdZ_;
  bool usePUcands_;

  /// Store extra information in a ValueMap
  template<typename Hand, typename T>
  void storeMap(edm::Event &iEvent, 
  const Hand & handle,
  const std::vector<T> & values,
  const std::string    & label) const ;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonMiniIso::MuonMiniIso(const edm::ParameterSet& iConfig):
probes_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("probes"))),
pfCandidates_(consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
dRCandProbeVeto_(iConfig.getParameter<double>("dRCandProbeVeto")),
dRCandSoftActivityCone_(iConfig.getParameter<double>("dRCandSoftActivityCone")),
CandPtThreshold_(iConfig.getParameter<double>("CandPtThreshold"))
{
  produces<edm::ValueMap<float> >("miniIso");
  produces<edm::ValueMap<float> >("activity");

  produces<edm::ValueMap<float> >("SAtoMiniConePt");
  produces<edm::ValueMap<float> >("SAtoMiniConeEta");
  produces<edm::ValueMap<float> >("SAtoMiniConePhi");
  produces<edm::ValueMap<float> >("SAtoMiniConeMass");

  produces<edm::ValueMap<float> >("DeadConePt");
  produces<edm::ValueMap<float> >("DeadConeEta");
  produces<edm::ValueMap<float> >("DeadConePhi");
  produces<edm::ValueMap<float> >("DeadConeMass");

  produces<edm::ValueMap<float> >("MiniToDeadConePt");
  produces<edm::ValueMap<float> >("MiniToDeadConeEta");
  produces<edm::ValueMap<float> >("MiniToDeadConePhi");
  produces<edm::ValueMap<float> >("MiniToDeadConeMass");
}


MuonMiniIso::~MuonMiniIso()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonMiniIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByToken(probes_, probes);

  //Handle<View<reco::PFCandidate> > pfCandidates;
  Handle<PFCollection> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  View<reco::Muon>::const_iterator probe, endprobes=probes->end();
  PFCollection::const_iterator iP, beginpf = pfCandidates->begin(), endpf=pfCandidates->end();
  unsigned int n = probes->size();
// printf("Num probes: %d", n); 
  std::vector<float> iso(n,0);
  std::vector<float> activity(n,0);
  // added by Brent
  // want to fill these variables for dR = deadcone radius --> dRCandSoftActivityCone

  std::vector<float> SAtoMiniConePt(n,0);
  std::vector<float> SAtoMiniConeEta(n,0);
  std::vector<float> SAtoMiniConePhi(n,0);
  std::vector<float> SAtoMiniConeMass(n,0);

  std::vector<float> DeadConePt(n,0);
  std::vector<float> DeadConeEta(n,0);
  std::vector<float> DeadConePhi(n,0);
  std::vector<float> DeadConeMass(n,0);

  std::vector<float> MiniToDeadConePt(n,0);
  std::vector<float> MiniToDeadConeEta(n,0);
  std::vector<float> MiniToDeadConePhi(n,0);
  std::vector<float> MiniToDeadConeMass(n,0);

  // loop on PROBES
  unsigned int imu = 0;
  for (probe = probes->begin(); probe != endprobes; ++probe, ++imu) {
    const reco::Muon &mu = *probe;
//std::cout<<"Probe Event "<< imu << std::endl;
//printf("\n"); 
    double dR_miniIso;
    if (mu.pt() <= 50) dR_miniIso = 0.2;
    else if (50 < mu.pt() && mu.pt() < 200) dR_miniIso = 10/mu.pt();
    else dR_miniIso = 0.05;

    TLorentzVector SAtoMiniConeMom;
    TLorentzVector DeadConeMom;
    TLorentzVector MiniToDeadConeMom;

    // loop on PF candidates
    int ipf=0;
    for (iP = beginpf; iP != endpf; ++iP, ++ipf) {
 
      // check pf candidate threshold
      if(iP->get()->pt() < CandPtThreshold_) continue;
      
      // get delta R of the PF Cand from the probe muon
      double dr = deltaR( *(iP->get() ) , mu );

      // if PF cand is outside soft activity cone, continue to next iteration
      if (dr > dRCandSoftActivityCone_) continue;
//      std::cout<<"PF Candidate #"<<ipf<<": "<<std::endl;
//        std::cout<<"    ID = " << iP->get()->pdgId() << std::endl;
//        std::cout<<"    Pt = " << iP->get()->pt() << std::endl;
//        std::cout<<"    Eta = " << iP->get()->eta() << std::endl;
//        std::cout<<"    Phi = " << iP->get()->phi() << std::endl;
//        std::cout<<"    Mass = " << iP->get()->mass() << std::endl;
//	std::cout<<"	dR = " << dr << std::endl;
//	printf("\n");
      TLorentzVector Candmom;
      Candmom.SetPtEtaPhiM(iP->get()->pt(), iP->get()->eta(), iP->get()->phi(), iP->get()->mass());

      // if the PFcand is within the Deadcone, add its p4 to the total Deadcone 4-momentum, then continue to next iteratio
      if (dr < dRCandProbeVeto_) {
        DeadConeMom += Candmom;
        continue;
      }
      
      if (dr <= dR_miniIso) {
        MiniToDeadConeMom += Candmom;
        iso[imu] += iP->get()->pt();
      } else {
        SAtoMiniConeMom += Candmom;
        activity[imu] += iP->get()->pt();
      }

    } // end loop on PF cands

    SAtoMiniConePt[imu] = SAtoMiniConeMom.Pt();
    SAtoMiniConeEta[imu] = SAtoMiniConeMom.Eta();
    SAtoMiniConePhi[imu] = SAtoMiniConeMom.Phi();
    SAtoMiniConeMass[imu] = SAtoMiniConeMom.M();

    DeadConePt[imu] = DeadConeMom.Pt();
    DeadConeEta[imu] = DeadConeMom.Eta();
    DeadConePhi[imu] = DeadConeMom.Phi();
    DeadConeMass[imu] = DeadConeMom.M();

    MiniToDeadConePt[imu] = MiniToDeadConeMom.Pt();
    MiniToDeadConeEta[imu] = MiniToDeadConeMom.Eta();
    MiniToDeadConePhi[imu] = MiniToDeadConeMom.Phi();
    MiniToDeadConeMass[imu] = MiniToDeadConeMom.M();

  }// end loop on probes

  storeMap(iEvent, probes, iso, "miniIso");
  storeMap(iEvent, probes, activity, "activity");

  storeMap(iEvent, probes, SAtoMiniConePt, "SAtoMiniConePt");
  storeMap(iEvent, probes, SAtoMiniConeEta, "SAtoMiniConeEta");  
  storeMap(iEvent, probes, SAtoMiniConePhi, "SAtoMiniConePhi");  
  storeMap(iEvent, probes, SAtoMiniConeMass, "SAtoMiniConeMass");

  storeMap(iEvent, probes, DeadConePt, "DeadConePt");  
  storeMap(iEvent, probes, DeadConeEta, "DeadConeEta");  
  storeMap(iEvent, probes, DeadConePhi, "DeadConePhi");  
  storeMap(iEvent, probes, DeadConeMass, "DeadConeMass");  

  storeMap(iEvent, probes, MiniToDeadConePt, "MiniToDeadConePt");  
  storeMap(iEvent, probes, MiniToDeadConeEta, "MiniToDeadConeEta");  
  storeMap(iEvent, probes, MiniToDeadConePhi, "MiniToDeadConePhi");  
  storeMap(iEvent, probes, MiniToDeadConeMass, "MiniToDeadConeMass"); 
}

template<typename Hand, typename T>
void
MuonMiniIso::storeMap(edm::Event &iEvent,
const Hand & handle,
const std::vector<T> & values,
const std::string    & label) const {
  using namespace edm; using namespace std;
  auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMiniIso);
