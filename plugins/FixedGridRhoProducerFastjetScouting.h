#ifndef PhysicsTools_SUEPScouting_plugins_FixedGridRhoProducerFastjetScouting_h
#define PhysicsTools_SUEPScouting_plugins_FixedGridRhoProducerFastjetScouting_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

class FixedGridRhoProducerFastjetScouting : public edm::stream::EDProducer<> {
public:
  explicit FixedGridRhoProducerFastjetScouting(const edm::ParameterSet& iConfig);
  ~FixedGridRhoProducerFastjetScouting() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::InputTag pfCandidatesTag_;
  edm::InputTag electronsTag_;
  fastjet::GridMedianBackgroundEstimator bge_;

  edm::EDGetTokenT<std::vector<ScoutingParticle>> input_pfcoll_token_;
  edm::EDGetTokenT<std::vector<ScoutingElectron>> input_elecoll_token_;
};

#endif
