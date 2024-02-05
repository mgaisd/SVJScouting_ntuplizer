#include "PhysicsTools/SVJScouting/plugins/FixedGridRhoProducerFastjetScouting.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"

using namespace std;

FixedGridRhoProducerFastjetScouting::FixedGridRhoProducerFastjetScouting(const edm::ParameterSet& iConfig)
    : bge_(iConfig.getParameter<double>("maxRapidity"), iConfig.getParameter<double>("gridSpacing")) {
  pfCandidatesTag_ = iConfig.getParameter<edm::InputTag>("pfCandidatesTag");
  electronsTag_ = iConfig.getParameter<edm::InputTag>("electronsTag");
  produces<double>();

  input_pfcoll_token_ = consumes<std::vector<ScoutingParticle>>(pfCandidatesTag_);
  input_elecoll_token_ = consumes<std::vector<ScoutingElectron>>(electronsTag_);
}

FixedGridRhoProducerFastjetScouting::~FixedGridRhoProducerFastjetScouting() {}

void FixedGridRhoProducerFastjetScouting::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::vector<fastjet::PseudoJet> inputs;

  edm::Handle<std::vector<ScoutingParticle>> pfColl;
  iEvent.getByToken(input_pfcoll_token_, pfColl);
  for (std::vector<ScoutingParticle>::const_iterator ibegin = pfColl->begin(), iend = pfColl->end(), i = ibegin; i != iend;
       ++i) {
    math::PtEtaPhiMLorentzVector v{i->pt(), i->eta(), i->phi(), i->m()};
    inputs.push_back(fastjet::PseudoJet(v.px(), v.py(), v.pz(), v.energy()));
  }

  //electrons are not in PF candidate collection
  edm::Handle<std::vector<ScoutingElectron>> eleColl;
  iEvent.getByToken(input_elecoll_token_, eleColl);
  for (std::vector<ScoutingElectron>::const_iterator ibegin = eleColl->begin(), iend = eleColl->end(), i = ibegin; i != iend;
       ++i) {
    math::PtEtaPhiMLorentzVector v{i->pt(), i->eta(), i->phi(), i->m()};
    inputs.push_back(fastjet::PseudoJet(v.px(), v.py(), v.pz(), v.energy()));
  }
  bge_.set_particles(inputs);
  iEvent.put(std::make_unique<double>(bge_.rho()));
}

DEFINE_FWK_MODULE(FixedGridRhoProducerFastjetScouting);
