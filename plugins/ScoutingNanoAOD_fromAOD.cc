// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>
#include <math.h>

#include "boost/algorithm/string.hpp"

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

//Added for offline jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//Added for MET
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"


//Added for offline electrons and muons
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Adding Gen jets
#include "DataFormats/JetReco/interface/GenJet.h"

//Adding Reco vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"

#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"

// for mini isolation
#include "ScoutingMiniIsolation.h"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//------------------------------------------------------------------------
 // the user information
 // 
 // To associate extra information to a PseudoJet, one first has to
 // create a class, derived from UserInfoBase, that contains
 // that information.
 //
 // In our simple example, we shall use 2 informations
 //  - the PDG id associated with the particle
 //  - the "vertex number" associated with the particle
 class PdgIdInfo : public PseudoJet::UserInfoBase{
 public:
   // default ctor
   //  - pdg_id        the PDG id of the particle
   PdgIdInfo(const int & pdg_id_in) :
     _pdg_id(pdg_id_in){}
  
   /// access to the PDG id
   int pdg_id() const { return _pdg_id;}
   
 protected:
   int _pdg_id;         // the associated pdg id
 };

//------------------------------------------------------------------------


class ScoutingNanoAOD_fromAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingNanoAOD_fromAOD(const edm::ParameterSet&);
  ~ScoutingNanoAOD_fromAOD();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  int getCharge(int pdgId);
  bool jetID(const ScoutingPFJet &pfjet);
  bool jetIDoff(const reco::PFJet &pfjet);

  //Scouting tokens
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;
  const edm::EDGetTokenT<std::vector<ScoutingMuon> >        muonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingElectron> >  	electronsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPhoton> >  	  photonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingParticle> >  	pfcandsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPFJet> >  		  pfjetsToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >  	  verticesToken;
  const edm::EDGetTokenT<double> metPtToken;
  const edm::EDGetTokenT<double> metPhiToken;

  double jetAK4ScoutPtMin = 0.;
  double jetAK8ScoutPtMin = 0.;
  
  //Offline tokens
  const edm::EDGetTokenT<reco::VertexCollection>              recoverticeToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet>> recoak4PFJetToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet>> recoak4PuppiJetToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet>> recoak8PuppiJetToken;
  const edm::EDGetTokenT<std::vector<reco::GsfElectron>> recoElectronToken;
  const edm::EDGetTokenT<edm::View<reco::Muon> > recoMuonToken;
  const edm::EDGetTokenT<std::vector<reco::PFCandidate>> recoPfCandidateToken;
  const edm::EDGetTokenT<std::vector<reco::PFMET>> recoMetToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

  double jetAK4PtMin = 0.;
  double jetAK8PtMin = 0.;

  //Gen info
  const edm::EDGetTokenT<std::vector<reco::GenJet> >            genak4jetsToken; 
  const edm::EDGetTokenT<std::vector<reco::GenJet> >            genak8jetsToken; 
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> >       pileupInfoToken;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> >       pileupInfoToken2;
  const edm::EDGetTokenT<GenEventInfoProduct>                  genEvtInfoToken;
  const edm::EDGetTokenT<GenLumiInfoHeader>  	genLumiInfoHeadTag_;
  const edm::EDGetTokenT<std::vector<reco::GenMET> >            genMetToken;


  const edm::EDGetTokenT<double>  	rhoToken2;
  const edm::EDGetTokenT<double>  	prefireToken;
  const edm::EDGetTokenT<double>  	prefireTokenup;
  const edm::EDGetTokenT<double>  	prefireTokendown;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
	
  // Generator-level information
  // Flags for the different types of triggers used in the analysis
  // For now we are interested in events passing either the single or double lepton triggers


  // Trigger information 
  bool doL1;       
  bool doData;       
  bool doSignal;       
  bool isMC;
  bool onlyScouting;
  //bool monitor;
  bool runScouting = false;
  bool runOffline =false;
  bool runGen = false;
  std::string label;

  HLTPrescaleProvider hltPSProv_;
  
  std::string hltProcess_; //name of HLT process, usually "HLT"

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::Handle<edm::TriggerResults> triggerBits;

  std::vector<std::string>     l1Seeds_;
  std::vector<std::string>     hltSeeds_;
  std::vector<bool>            l1Result_;
  std::vector<int>             l1Prescale_;
  std::vector<bool>            hltResult_;
  std::vector<std::string>     hltResultName_;
  vector<double>            PSweights;

  UInt_t scouting_trig; 
  UInt_t scouting_trig_prescaled;
  UInt_t scouting_trig_zero_bias;
  UInt_t offline_trig; 
  UInt_t veto_trig;


  //Photon
  UInt_t n_pho;
  vector<Float16_t> 	     Photon_pt;
  vector<Float16_t>        Photon_eta;
  vector<Float16_t>        Photon_phi;
  vector<Float16_t>	       Photon_m;
  vector<Float16_t>	       Photon_sigmaietaieta;
  vector<Float16_t>	       Photon_HoE;
  vector<Float16_t>        Photon_ecaliso;
  vector<Float16_t>	       Photon_hcaliso;

  //Scouting Electron
  UInt_t n_ele;
  vector<Float16_t> 	     Electron_pt;
  vector<Float16_t>        Electron_eta;
  vector<Float16_t>        Electron_phi;
  vector<Float16_t>	       Electron_m;
  vector<Float16_t>        Electron_d0;
  vector<Float16_t>	       Electron_dz;
  vector<Float16_t>	       Electron_detain;
  vector<Float16_t>	       Electron_dphiin;
  vector<Float16_t>	       Electron_sigmaietaieta;
  vector<Float16_t>	       Electron_HoE;
  vector<Float16_t>	       Electron_ooEMOop;
  vector<Float16_t>	       Electron_mHits;
  vector<Float16_t>        Electron_charge;
  vector<Float16_t>        Electron_ecaliso;
  vector<Float16_t>	       Electron_hcaliso;
  vector<Float16_t>        Electron_trkiso;
  vector<Float16_t>        Electron_combinediso;
  vector<bool>             Electron_ID;
  vector<Float16_t>        Electron_chargedMiniIso;
  vector<Float16_t>        Electron_combinedMiniIso;


  //Offline Electron
  UInt_t n_ele_off;
  vector<Float16_t> 	     OffElectron_pt;
  vector<Float16_t>        OffElectron_eta;
  vector<Float16_t>        OffElectron_phi;
  vector<Float16_t>	       OffElectron_m;
  vector<Float16_t>        OffElectron_d0;
  vector<Float16_t>	       OffElectron_dz;
  vector<Float16_t>	       OffElectron_detain;
  vector<Float16_t>	       OffElectron_dphiin;
  vector<Float16_t>	       OffElectron_sigmaietaieta;
  vector<Float16_t>	       OffElectron_HoE;
  vector<Float16_t>	       OffElectron_ooEMOop;
  vector<Float16_t>	       OffElectron_mHits;
  vector<Float16_t>        OffElectron_charge;
  vector<Float16_t>        OffElectron_ecaliso;
  vector<Float16_t>	       OffElectron_hcaliso;
  vector<Float16_t>        OffElectron_trkiso;
  vector<Float16_t>        OffElectron_combinediso;
  vector<bool>             OffElectron_ID;
  Float_t isoChargedHadrons_;
  Float_t isoNeutralHadrons_;
  Float_t isoPhotons_;
  Float_t isoChargedFromPU_;

  //Muon
  UInt_t n_mu;
  vector<Float16_t>            Muon_pt;
  vector<Float16_t>            Muon_eta;
  vector<Float16_t>            Muon_phi;
  vector<Float16_t>            Muon_m;
  vector<Float16_t>            Muon_ecaliso;
  vector<Float16_t>            Muon_hcaliso;
  vector<Float16_t>            Muon_trkiso;
  vector<Float16_t>            Muon_chi2;
  vector<bool>                 Muon_isGlobalMuon;
  vector<bool>                 Muon_isTrackerMuon;
  vector<Float16_t>            Muon_ndof;
  vector<Float16_t>            Muon_charge;
  vector<Float16_t>            Muon_dxy;
  vector<Float16_t>            Muon_dz;
  vector<Float16_t>            Muon_nvalidmuon_hits;
  vector<Float16_t>            Muon_nvalidpixelhits;
  vector<Float16_t>            Muon_nmatchedstations;
  vector<Float16_t>            Muon_type;
  vector<Float16_t>            Muon_nvalidstriphits;
  vector<Float16_t>            Muon_trkqoverp;
  vector<Float16_t>            Muon_trklambda;
  vector<Float16_t>            Muon_trkpt;
  vector<Float16_t>            Muon_trkphi;
  vector<Float16_t>            Muon_trketa;
  vector<Float16_t>            Muon_trkqoverperror;
  vector<Float16_t>            Muon_trklambdaerror;
  vector<Float16_t>            Muon_trkpterror;
  vector<Float16_t>            Muon_trkphierror;
  vector<Float16_t>            Muon_trketaerror;
  vector<Float16_t>            Muon_trkdszerror;
  vector<Float16_t>            Muon_trkdsz;
  vector<Float16_t>            Muon_chargedMiniIso;
  vector<Float16_t>            Muon_combinedMiniIso;

  //Offline Muon
  UInt_t n_mu_off;
  vector<Float16_t>            OffMuon_pt;
  vector<Float16_t>            OffMuon_eta;
  vector<Float16_t>            OffMuon_phi;
  vector<Float16_t>            OffMuon_m;
  vector<Float16_t>            OffMuon_charge;
  vector<bool>                 OffMuon_isGlobalMuon;
  vector<bool>                 OffMuon_isTrackerMuon;
  vector<bool>                 OffMuon_isStandAloneMuon;
  vector<bool>                 OffMuon_isCaloMuon;
  vector<bool>                 OffMuon_isPFMuon;

  UInt_t                       PU_num;

  //Scouting PFJets
  UInt_t                       n_jet;
  UInt_t                       n_jetId;
  UInt_t                       n_jetoff;
  UInt_t                       n_jetIdoff;
  UInt_t                       n_jetoffpuppi;
  UInt_t                       n_jetIdoffpuppi;
  float                        ht;
  float                        htoff;
  float                        htoffpuppi;
  bool                         passJetId;
  vector<Float16_t> 	         Jet_pt;
  vector<Float16_t>            Jet_eta;
  vector<Float16_t>            Jet_phi;
  vector<Float16_t>            Jet_m;
  vector<Float16_t>            Jet_area;
  vector<Float16_t>            Jet_chargedHadronEnergy;
  vector<Float16_t>            Jet_neutralHadronEnergy;
  vector<Float16_t>            Jet_photonEnergy;
  vector<Float16_t>            Jet_electronEnergy;
  vector<Float16_t>	       Jet_muonEnergy;
  vector<Float16_t>	       Jet_HFHadronEnergy;
  vector<Float16_t>	       Jet_HFEMEnergy;
  vector<Float16_t>	       Jet_HOEnergy;
  vector<Float16_t>	       Jet_chargedHadronMultiplicity;
  vector<Float16_t>            Jet_neutralHadronMultiplicity;
  vector<Float16_t>	       Jet_photonMultiplicity;
  vector<Float16_t>	       Jet_electronMultiplicity;
  vector<Float16_t>	       Jet_muonMultiplicity;
  vector<Float16_t>	       Jet_HFHadronMultiplicity;
  vector<Float16_t>	       Jet_HFEMMultiplicity;
  vector<Float16_t> 	       Jet_csv;
  vector<Float16_t> 	       Jet_mvaDiscriminator;
  vector<Float16_t>  	       Jet_nConstituents;
  // for jet id
  vector<Float16_t>        Jet_neutralHadronEnergyFraction;
  vector<Float16_t>	       Jet_photonEnergyFraction;
  vector<Float16_t>	       Jet_chargedHadronEnergyFraction;
  vector<Float16_t>	       Jet_muonEnergyFraction;
  vector<Float16_t>	       Jet_electronEnergyFraction;
  vector<Float16_t>	       Jet_neutralMultiplicity;
  vector<Float16_t>	       Jet_chargedMultiplicity;

  vector<bool>                 Jet_passId;

  //Offline AK4 PFJets (regular PF)
  vector<Float16_t> 	     OffJet_pt;
  vector<Float16_t>        OffJet_eta;
  vector<Float16_t>        OffJet_phi;
  vector<Float16_t>	       OffJet_m;
  vector<Float16_t>	       OffJet_area;
  vector<Float16_t>	       OffJet_nConstituents;
  vector<Float16_t>	       OffJet_chargedHadronEnergy;
  vector<Float16_t>        OffJet_neutralHadronEnergy;
  vector<Float16_t>	       OffJet_photonEnergy;
  vector<Float16_t>	       OffJet_electronEnergy;
  vector<Float16_t>	       OffJet_muonEnergy;
  vector<Float16_t>	       OffJet_HFHadronEnergy;
  vector<Float16_t>	       OffJet_HFEMEnergy;
  vector<Float16_t>	       OffJet_HOEnergy;
  // for jet id
  vector<Float16_t>        OffJet_neutralHadronEnergyFraction;
  vector<Float16_t>	       OffJet_neutralEmEnergyFraction;
  vector<Float16_t>	       OffJet_chargedHadronEnergyFraction;
  vector<Float16_t>	       OffJet_muonEnergyFraction;
  vector<Float16_t>	       OffJet_chargedEmEnergyFraction;
  vector<Float16_t>	       OffJet_electronEnergyFraction;
  vector<Float16_t>	       OffJet_photonEnergyFraction;
  vector<Float16_t>	       OffJet_neutralMultiplicity;
  vector<Float16_t>	       OffJet_chargedMultiplicity;

  vector<Float16_t>	       OffJet_chargedHadronMultiplicity;
  vector<Float16_t>        OffJet_neutralHadronMultiplicity;
  vector<Float16_t>	       OffJet_photonMultiplicity;
  vector<Float16_t>	       OffJet_electronMultiplicity;
  vector<Float16_t>	       OffJet_muonMultiplicity;
  vector<Float16_t>	       OffJet_HFHadronMultiplicity;
  vector<Float16_t>	       OffJet_HFEMMultiplicity;
  vector<bool>             OffJet_passId;

  //Offline AK4 Puppi Jets
  vector<Float16_t> 	     OffPuppiJet_pt;
  vector<Float16_t>        OffPuppiJet_eta;
  vector<Float16_t>        OffPuppiJet_phi;
  vector<Float16_t>	       OffPuppiJet_m;
  vector<Float16_t>	       OffPuppiJet_area;
  vector<Float16_t>	       OffPuppiJet_chargedHadronEnergy;
  vector<Float16_t>        OffPuppiJet_neutralHadronEnergy;
  vector<Float16_t>	       OffPuppiJet_photonEnergy;
  vector<Float16_t>	       OffPuppiJet_electronEnergy;
  vector<Float16_t>	       OffPuppiJet_muonEnergy;
  vector<Float16_t>	       OffPuppiJet_HFHadronEnergy;
  vector<Float16_t>	       OffPuppiJet_HFEMEnergy;
  vector<Float16_t>	       OffPuppiJet_HOEnergy;
  // for jet id
  vector<Float16_t>        OffPuppiJet_neutralHadronEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_neutralEmEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_chargedHadronEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_muonEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_chargedEmEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_electronEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_photonEnergyFraction;
  vector<Float16_t>	       OffPuppiJet_nConstituents;
  vector<Float16_t>	       OffPuppiJet_neutralMultiplicity;
  vector<Float16_t>	       OffPuppiJet_chargedMultiplicity;

  vector<Float16_t>	       OffPuppiJet_chargedHadronMultiplicity;
  vector<Float16_t>        OffPuppiJet_neutralHadronMultiplicity;
  vector<Float16_t>	       OffPuppiJet_photonMultiplicity;
  vector<Float16_t>	       OffPuppiJet_electronMultiplicity;
  vector<Float16_t>	       OffPuppiJet_muonMultiplicity;
  vector<Float16_t>	       OffPuppiJet_HFHadronMultiplicity;
  vector<Float16_t>	       OffPuppiJet_HFEMMultiplicity;
  vector<bool>             OffPuppiJet_passId;
  
  //GenJets Ak4
  UInt_t                       n_genjet;
  vector<Float16_t>            GenJet_pt;
  vector<Float16_t>            GenJet_eta;
  vector<Float16_t>            GenJet_phi;
  vector<Float16_t>            GenJet_mass;

  //GenJets Ak8
  UInt_t                       n_genfatjet;
  vector<Float16_t>            GenFatJet_pt;
  vector<Float16_t>            GenFatJet_eta;
  vector<Float16_t>            GenFatJet_phi;
  vector<Float16_t>            GenFatJet_mass;

  //GenJet matching indices
  vector<int>                  Jet_genJetIdx;        // index of matched gen jet for AK4
  vector<int>                  FatJet_genJetAK8Idx;  // index of matched gen jet for AK8

  
  //CZZ: to add OffPFCands
  UInt_t                       n_offpfcand;
  UInt_t                       n_offpfMu;
  UInt_t                       n_offpfEl;
  vector<Float16_t>            OffPFcand_pt;
  vector<Float16_t>            OffPFcand_eta;
  vector<Float16_t>            OffPFcand_phi;
  vector<Float16_t>            OffPFcand_m;
  vector<Float16_t>            OffPFcand_pdgid;
  vector<Float16_t>            OffPFcand_q;
  vector<Float16_t>            OfflineFatJetPFCands_jetIdx;
  vector<Float16_t>            OfflineFatJetPFCands_pFCandsIdx;
  


  // Scouting PFCand
  UInt_t                       n_pfcand;
  UInt_t                       n_pfMu;
  UInt_t                       n_pfEl;
  vector<Float16_t>            PFcand_pt;
  vector<Float16_t>            PFcand_eta;
  vector<Float16_t>            PFcand_phi;
  vector<Float16_t>            PFcand_m;
  vector<Float16_t>            PFcand_pdgid;
  vector<Float16_t>            PFcand_q;
  vector<Float16_t>            PFcand_vertex;
  vector<Float16_t>            FatJetPFCands_jetIdx;
  vector<Float16_t>            FatJetPFCands_pFCandsIdx;


  // Fatjets 
  UInt_t                       n_fatjet;
  //vector<Float16_t>            FatJet_jesc; 
  vector<Float16_t>            FatJet_area;
  vector<Float16_t>            FatJet_eta;
  //vector<Float16_t>            FatJet_n2b1;
  //vector<Float16_t>            FatJet_n3b1;
  vector<Float16_t>            FatJet_phi;
  vector<Float16_t>            FatJet_pt;
  //vector<Float16_t>            FatJet_tau1;
  //vector<Float16_t>            FatJet_tau2;
  //vector<Float16_t>            FatJet_tau3;
  //vector<Float16_t>            FatJet_tau4;
  //vector<Float16_t>            FatJet_tau21;
  //vector<Float16_t>            FatJet_tau32;
  vector<Float16_t>            FatJet_mass;
  //vector<Float16_t>            FatJet_msoftdrop;
  //vector<Float16_t>            FatJet_mtrim;
  vector<Float16_t>            FatJet_nConstituents;
  vector<bool>                 FatJet_passId;

  // for jet id
  vector<Float16_t>        FatJet_neutralHadronEnergyFraction;
  vector<Float16_t>	       FatJet_photonEnergyFraction;
  vector<Float16_t>	       FatJet_chargedHadronEnergyFraction;
  vector<Float16_t>	       FatJet_muonEnergyFraction;
  vector<Float16_t>	       FatJet_electronEnergyFraction;
  vector<Float16_t>	       FatJet_neutralMultiplicity;
  vector<Float16_t>	       FatJet_chargedMultiplicity;

  
  // Fatjets offline
  UInt_t                       n_fatjet_off;
  vector<Float16_t>            OfflineFatJet_area;
  vector<Float16_t>            OfflineFatJet_eta;
  // vector<Float16_t>            OfflineFatJet_n2b1;
  // vector<Float16_t>            OfflineFatJet_n3b1;
  vector<Float16_t>            OfflineFatJet_phi;
  vector<Float16_t>            OfflineFatJet_pt;
  // vector<Float16_t>            OfflineFatJet_tau1;
  // vector<Float16_t>            OfflineFatJet_tau2;
  // vector<Float16_t>            OfflineFatJet_tau3;
  // vector<Float16_t>            OfflineFatJet_tau4;
  // vector<Float16_t>            OfflineFatJet_tau21;
  // vector<Float16_t>            OfflineFatJet_tau32;
  vector<Float16_t>            OfflineFatJet_mass;
  // vector<Float16_t>            OfflineFatJet_msoftdrop;
  // vector<Float16_t>            OfflineFatJet_mtrim;
  vector<Float16_t>            OfflineFatJet_nConstituents;
  


  //Offline Puppi PFJets
  vector<Float16_t> 	     OffPuppiFatJet_pt;
  vector<Float16_t>        OffPuppiFatJet_eta;
  vector<Float16_t>        OffPuppiFatJet_phi;
  vector<Float16_t>	       OffPuppiFatJet_m;
  vector<Float16_t>	       OffPuppiFatJet_area;
  vector<Float16_t>	       OffPuppiFatJet_chargedHadronEnergy;
  vector<Float16_t>        OffPuppiFatJet_neutralHadronEnergy;
  vector<Float16_t>	       OffPuppiFatJet_photonEnergy;
  vector<Float16_t>	       OffPuppiFatJet_electronEnergy;
  vector<Float16_t>	       OffPuppiFatJet_muonEnergy;
  vector<Float16_t>	       OffPuppiFatJet_HFHadronEnergy;
  vector<Float16_t>	       OffPuppiFatJet_HFEMEnergy;
  vector<Float16_t>	       OffPuppiFatJet_HOEnergy;
  vector<Float16_t>	       OffPuppiFatJet_chargedHadronMultiplicity;
  vector<Float16_t>        OffPuppiFatJet_neutralHadronMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_photonMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_electronMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_muonMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_HFHadronMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_HFEMMultiplicity;

  // for jet id
  vector<Float16_t>        OffPuppiFatJet_neutralHadronEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_neutralEmEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_chargedHadronEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_muonEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_chargedEmEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_electronEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_photonEnergyFraction;
  vector<Float16_t>	       OffPuppiFatJet_nConstituents;
  vector<Float16_t>	       OffPuppiFatJet_neutralMultiplicity;
  vector<Float16_t>	       OffPuppiFatJet_chargedMultiplicity;


  vector<bool>             OffPuppiFatJet_passId;
  UInt_t                   n_fatjet_off_puppi;
  UInt_t                   n_fatjetIdoffpuppi;
  bool                     passOffPuppiFatJetId;

  //CZZ: to add OffPFCands
  UInt_t                       n_offpuppipfcand;
  vector<Float16_t>            OffPuppiPFcand_pt;
  vector<Float16_t>            OffPuppiPFcand_eta;
  vector<Float16_t>            OffPuppiPFcand_phi;
  vector<Float16_t>            OffPuppiPFcand_m;
  vector<Float16_t>            OffPuppiPFcand_pdgid;
  vector<Float16_t>            OffPuppiPFcand_q;
  vector<Float16_t>            OfflinePuppiFatJetPFCands_jetIdx;
  vector<Float16_t>            OfflinePuppiFatJetPFCands_pFCandsIdx;

  // Primary vertices
  UInt_t n_pvs;
  vector<Float16_t>            Vertex_x;
  vector<Float16_t>            Vertex_y;
  vector<Float16_t>            Vertex_z;
  vector<Float16_t>            Vertex_tracksSize;
  vector<Float16_t>            Vertex_chi2;
  vector<Float16_t>            Vertex_ndof;
  vector<Float16_t>            Vertex_isValidVtx;

  //prefire
  float                        rho2;
  float                        prefire;
  float                        prefireup;
  float                        prefiredown;



  // MET
  double met_pt_reco, met_phi_reco;

  //Scouting MET
  double met_pt, met_phi;

  //Gen MET
  double genMET_pt, genMET_phi;

  //reco vertices
  Int_t nPV_;        // number of reconsrtucted primary vertices
    
  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;
  int event_;

};

ScoutingNanoAOD_fromAOD::ScoutingNanoAOD_fromAOD(const edm::ParameterSet& iConfig): 
  
  //Scouting tokens
  muonsToken               (consumes<std::vector<ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken           (consumes<std::vector<ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken             (consumes<std::vector<ScoutingPhoton> >           (iConfig.getParameter<edm::InputTag>("photons"))), 
  pfcandsToken             (consumes<std::vector<ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
  pfjetsToken              (consumes<std::vector<ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))), 
  verticesToken            (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("vertices"))),
  metPtToken               (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("metPt"))),
  metPhiToken              (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("metPhi"))),

  jetAK4ScoutPtMin         (iConfig.getParameter<double>("jetAK4ScoutPtMin")),
  jetAK8ScoutPtMin         (iConfig.getParameter<double>("jetAK8ScoutPtMin")),

  //Offline tokens
  recoverticeToken     (consumes<reco::VertexCollection>          (iConfig.getParameter<edm::InputTag>("verticesReco"))),
  recoak4PFJetToken    (consumes<std::vector<reco::PFJet>>        (iConfig.getParameter<edm::InputTag>("ak4pfjetsReco"))),
  recoak4PuppiJetToken (consumes<std::vector<reco::PFJet>>        (iConfig.getParameter<edm::InputTag>("ak4puppijetsReco"))),
  recoak8PuppiJetToken (consumes<std::vector<reco::PFJet>>        (iConfig.getParameter<edm::InputTag>("ak8pfjetsReco"))),
  recoElectronToken    (consumes<std::vector<reco::GsfElectron> > (iConfig.getParameter<edm::InputTag>("electronsReco"))),
  recoMuonToken        (consumes<edm::View<reco::Muon>>           (iConfig.getParameter<edm::InputTag>("muonsReco"))),
  recoPfCandidateToken (consumes<std::vector<reco::PFCandidate>>  (iConfig.getParameter<edm::InputTag>("pfcandsReco"))),
  recoMetToken         (consumes<std::vector<reco::PFMET>>        (iConfig.getParameter<edm::InputTag>("metReco"))),


  jetAK4PtMin          (iConfig.getParameter<double>("jetAK4PtMin")),
  
  jetAK8PtMin          (iConfig.getParameter<double>("jetAK8PtMin")),

  //Gen info
  genak4jetsToken             (consumes<std::vector<reco::GenJet> >             (iConfig.getParameter<edm::InputTag>("genak4jets"))),
  genak8jetsToken             (consumes<std::vector<reco::GenJet> >             (iConfig.getParameter<edm::InputTag>("genak8jets"))),
  pileupInfoToken          (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo"))),
  pileupInfoToken2         (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo_sig"))),
  genEvtInfoToken          (consumes<GenEventInfoProduct>                    (iConfig.getParameter<edm::InputTag>("geneventinfo"))), 
  genLumiInfoHeadTag_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),   
  genMetToken              (consumes<std::vector<reco::GenMET> >             (iConfig.getParameter<edm::InputTag>("genMet"))),
  
  rhoToken2                (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho2"))),
  prefireToken             (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProb"))),
  prefireTokenup           (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
  prefireTokendown         (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProbDown"))),

  doL1                     (iConfig.existsAs<bool>("doL1")              ?    iConfig.getParameter<bool>  ("doL1")            : false),
  doData                   (iConfig.existsAs<bool>("doData")            ?    iConfig.getParameter<bool>  ("doData")            : false),
  doSignal                 (iConfig.existsAs<bool>("doSignal")          ?    iConfig.getParameter<bool>  ("doSignal")            : false),
  isMC                     (iConfig.existsAs<bool>("isMC")              ?    iConfig.getParameter<bool>  ("isMC")            : true),
  onlyScouting             (iConfig.existsAs<bool>("onlyScouting")      ?    iConfig.getParameter<bool>  ("onlyScouting")            : false),


  hltPSProv_(iConfig,consumesCollector(),*this), //it needs a referernce to the calling module for some reason, hence the *this   
  hltProcess_(iConfig.getParameter<std::string>("hltProcess")),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  l1Seeds_(iConfig.getParameter<std::vector<std::string> >("l1Seeds")),
  hltSeeds_(iConfig.getParameter<std::vector<std::string> >("hltSeeds"))


{
 // now do whatever initialization is needed
  usesResource("TFileService");


 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("Events"       , "tree");

  // Event weights 
  tree->Branch("lumSec"		                  ,&lumSec            ,"lumSec/i");
  tree->Branch("run"		                    ,&run                  ,"run/i");
  tree->Branch("event"		                    ,&event_                  ,"event/i");
  tree->Branch("PSweights"            	    ,&PSweights 	                 );
  tree->Branch("prefire"		                ,&prefire                      );
  tree->Branch("prefireup"		              ,&prefireup                    );
  tree->Branch("prefiredown"		            ,&prefiredown                  );
    
  // HLT and L1 Triggers results
  tree->Branch("hltResult"                      ,&hltResult_                    );              
  tree->Branch("hltResultName"                  ,&hltResultName_                );              
  tree->Branch("l1Result"		                    ,&l1Result_	                );		
  tree->Branch("l1Prescale"		                  ,&l1Prescale_                   );		

  
  //scouting, offline triggers
  tree->Branch("scouting_trig_prescaled"            	         ,&scouting_trig_prescaled 			,"scouting_trig_prescaled/i");
  tree->Branch("scouting_trig"            	                   ,&scouting_trig 			,"scounting_trig/i");
  tree->Branch("scouting_trig_zero_bias"            	        ,&scouting_trig_zero_bias 			,"scounting_trig_zero_bias/i");
  tree->Branch("offline_trig"            	        ,&offline_trig 			,"offline_trig/i");
  tree->Branch("veto_trig"            	        ,&veto_trig 			,"veto_trig/i");
  tree->Branch("genModel"            	        ,&label 			);

  //Scouting Electrons
  tree->Branch("nElectron"               	         ,&n_ele                        ,"nElectron/i");
  tree->Branch("Electron_pt"                    ,&Electron_pt                   );
  tree->Branch("Electron_eta"                   ,&Electron_eta 	                );
  tree->Branch("Electron_phi"                   ,&Electron_phi                  );
  tree->Branch("Electron_charge"                ,&Electron_charge               );
  tree->Branch("Electron_mass"            	        ,&Electron_m                    );
  tree->Branch("Electron_hoe"                   ,&Electron_HoE                  );
  tree->Branch("Electron_sieie"         ,&Electron_sigmaietaieta        );
  tree->Branch("Electron_dphiin"                ,&Electron_dphiin 	        );
  tree->Branch("Electron_detain"                ,&Electron_detain 	        );
  tree->Branch("Electron_mHits"                 ,&Electron_mHits 	        );
  tree->Branch("Electron_ooEMOop"               ,&Electron_ooEMOop              );            
  tree->Branch("Electron_trkiso"                 ,&Electron_trkiso 	        );
  tree->Branch("Electron_ecaliso"               ,&Electron_ecaliso              );
  tree->Branch("Electron_hcaliso"               ,&Electron_hcaliso              );
  tree->Branch("Electron_combinediso"               ,&Electron_combinediso   );
  tree->Branch("Electron_ID"               ,&Electron_ID   );
  tree->Branch("Electron_d0"               ,&Electron_d0              );
  tree->Branch("Electron_dz"               ,&Electron_dz              );
  tree->Branch("Electron_chargedMiniIso"               ,&Electron_chargedMiniIso              );
  tree->Branch("Electron_combinedMiniIso"               ,&Electron_combinedMiniIso              );


  //Scouting Photons
  tree->Branch("nPhotons"            	          ,&n_pho 			,"nPhotons/i");
  tree->Branch("Photon_pt"            	        ,&Photon_pt                     );
  tree->Branch("Photon_eta"            	        ,&Photon_eta                    );
  tree->Branch("Photon_phi"            	        ,&Photon_phi                    );	
  tree->Branch("Photon_mass"            	      ,&Photon_m 	                );
  tree->Branch("Photon_hcaliso"                 ,&Photon_hcaliso 		);
  tree->Branch("Photon_ecaliso"                 ,&Photon_ecaliso 		);
  tree->Branch("Photon_hoe"            	        ,&Photon_HoE                    );
  tree->Branch("Photon_sieie"                   ,&Photon_sigmaietaieta	        );

  //Scouting muons
  tree->Branch("nMuons"            	            ,&n_mu 	                        ,"nMuons/i");
  tree->Branch("Muon_pt"                        ,&Muon_pt                       );
  tree->Branch("Muon_eta"                       ,&Muon_eta                      );
  tree->Branch("Muon_phi"                       ,&Muon_phi                      );
  tree->Branch("Muon_mass"                      ,&Muon_m                        );
  tree->Branch("Muon_ecaliso"                   ,&Muon_ecaliso                  );
  tree->Branch("Muon_hcaliso"                   ,&Muon_hcaliso                  );
  tree->Branch("Muon_tkIso"                     ,&Muon_trkiso                   );
  tree->Branch("Muon_chi2"                      ,&Muon_chi2                     );
  tree->Branch("Muon_isGlobal"                  ,&Muon_isGlobalMuon             );
  tree->Branch("Muon_isTracker"                 ,&Muon_isTrackerMuon            );
  tree->Branch("Muon_ndof"                      ,&Muon_ndof                     );
  tree->Branch("Muon_charge"                    ,&Muon_charge	                  );
  tree->Branch("Muon_dxy"                       ,&Muon_dxy                      );
  tree->Branch("Muon_dz"                        ,&Muon_dz                       );
  tree->Branch("Muon_nvalidmuon_hits"           ,&Muon_nvalidmuon_hits          );
  tree->Branch("Muon_validpixelhits"            ,&Muon_nvalidpixelhits          );
  tree->Branch("Muon_nmatchedstations"          ,&Muon_nmatchedstations         );
  tree->Branch("Muon_type"                      ,&Muon_type                     );
  tree->Branch("Muon_nvalidstriphits"           ,&Muon_nvalidstriphits          );
  tree->Branch("Muon_tkqoverp"                 ,&Muon_trkqoverp                );
  tree->Branch("Muon_tklambda"                 ,&Muon_trklambda                );
  tree->Branch("Muon_tkpt"                     ,&Muon_trkpt                    );
  tree->Branch("Muon_tkphi"                    ,&Muon_trkphi                   );
  tree->Branch("Muon_tketa"                    ,&Muon_trketa                   );
  tree->Branch("Muon_tkqoverperror"            ,&Muon_trkqoverperror           );
  tree->Branch("Muon_tklambdaerror"            ,&Muon_trklambdaerror           );
  tree->Branch("Muon_tkpterror"                ,&Muon_trkpterror               );
  tree->Branch("Muon_tkphierror"               ,&Muon_trkphierror              );
  tree->Branch("Muon_tketaerror"               ,&Muon_trketaerror              );
  tree->Branch("Muon_tkdszerror"               ,&Muon_trkdszerror              );
  tree->Branch("Muon_tkdsz"                    ,&Muon_trkdsz                   );
  tree->Branch("Muon_chargedMiniIso"           ,&Muon_chargedMiniIso           );
  tree->Branch("Muon_combinedMiniIso"          ,&Muon_combinedMiniIso          );


  //Scouting AK4 PFJets
  tree->Branch("nJet"            	        ,&n_jet                         ,"nJet/i");
  tree->Branch("nJetId"            	        ,&n_jetId                       ,"nJetId/i");
  tree->Branch("Jet_pt"            	        ,&Jet_pt                        );
  tree->Branch("Jet_eta"            	        ,&Jet_eta                       );
  tree->Branch("Jet_phi"            	        ,&Jet_phi                       );
  tree->Branch("Jet_mass"            	        ,&Jet_m                         );
  tree->Branch("Jet_area"            	        ,&Jet_area                      );
  // Not needed for Type-1 MET
  // tree->Branch("Jet_chargedHadronEnergy"        ,&Jet_chargedHadronEnergy       );
  // tree->Branch("Jet_neutralHadronEnergy"        ,&Jet_neutralHadronEnergy       );
  // tree->Branch("Jet_photonEnergy"               ,&Jet_photonEnergy 	        );
  // tree->Branch("Jet_electronEnergy"             ,&Jet_electronEnergy            );
  // tree->Branch("Jet_muonEnergy"    	        ,&Jet_muonEnergy                );
  // tree->Branch("Jet_HFHadronEnergy"             ,&Jet_HFHadronEnergy            );
  // tree->Branch("Jet_HFEMEnergy"                 ,&Jet_HFEMEnergy                );
  // tree->Branch("Jet_HOEnergy"                   ,&Jet_HOEnergy                  );
  // tree->Branch("Jet_chargedHadronMultiplicity"  ,&Jet_chargedHadronMultiplicity );
  // tree->Branch("Jet_neutralHadronMultiplicity"  ,&Jet_neutralHadronMultiplicity );
  // tree->Branch("Jet_photonMultiplicity"         ,&Jet_photonMultiplicity        );
  // tree->Branch("Jet_electronMultiplicity"       ,&Jet_electronMultiplicity      );
  // tree->Branch("Jet_muonMultiplicity"           ,&Jet_muonMultiplicity          );
  // tree->Branch("Jet_HFHadronMultiplicity"       ,&Jet_HFHadronMultiplicity      );
  // tree->Branch("Jet_HFEMMultiplicity"           ,&Jet_HFEMMultiplicity          );
  // tree->Branch("Jet_csv"            	        ,&Jet_csv                       );
  // tree->Branch("Jet_mvaDiscriminator"           ,&Jet_mvaDiscriminator          );
  tree->Branch("Jet_nConstituents"              ,&Jet_nConstituents             );
  // tree->Branch("Jet_passId"                     ,&Jet_passId                    );

  //for jet id
  // tree->Branch("Jet_chargedHadronEnergyFraction"        ,&Jet_chargedHadronEnergyFraction       );
  // tree->Branch("Jet_neutralHadronEnergyFraction"        ,&Jet_neutralHadronEnergyFraction       );
  // tree->Branch("Jet_electronEnergyFraction"        ,&Jet_electronEnergyFraction       );
  // tree->Branch("Jet_muonEnergyFraction"        ,&Jet_muonEnergyFraction       );
  // tree->Branch("Jet_photonEnergyFraction"        ,&Jet_photonEnergyFraction       );
  // tree->Branch("Jet_neutralMultiplicity"        ,&Jet_neutralMultiplicity       );
  // tree->Branch("Jet_chargedMultiplicity"        ,&Jet_chargedMultiplicity       );

  //Scouting AK8 PFJets
  tree->Branch("nFatJet"                       ,&n_fatjet                      ,"nFatJet/i");
  tree->Branch("FatJet_area"                    ,&FatJet_area                   );
  tree->Branch("FatJet_eta"                     ,&FatJet_eta                    );
  //tree->Branch("FatJet_n2b1"                    ,&FatJet_n2b1                   );
  //tree->Branch("FatJet_n3b1"                    ,&FatJet_n3b1                   );
  tree->Branch("FatJet_phi"                     ,&FatJet_phi                    );
  tree->Branch("FatJet_pt"                      ,&FatJet_pt                     );
  //tree->Branch("FatJet_tau1"                    ,&FatJet_tau1                   );
  //tree->Branch("FatJet_tau2"                    ,&FatJet_tau2                   );
  //tree->Branch("FatJet_tau3"                    ,&FatJet_tau3                   );
  //tree->Branch("FatJet_tau4"                    ,&FatJet_tau4                   );
  //tree->Branch("FatJet_tau21"                   ,&FatJet_tau21                  );
  //tree->Branch("FatJet_tau32"                   ,&FatJet_tau32                  );
  tree->Branch("FatJet_mass"                    ,&FatJet_mass                   );
  //tree->Branch("FatJet_msoftdrop"               ,&FatJet_msoftdrop              );
  //tree->Branch("FatJet_mtrim"                   ,&FatJet_mtrim                  );
  tree->Branch("FatJet_nConstituents"                  ,&FatJet_nConstituents                 );
  tree->Branch("FatJet_passId"                  ,&FatJet_passId                 );

  //for jet id
  tree->Branch("FatJet_chargedHadronEnergyFraction"        ,&FatJet_chargedHadronEnergyFraction       );
  tree->Branch("FatJet_neutralHadronEnergyFraction"        ,&FatJet_neutralHadronEnergyFraction       );
  tree->Branch("FatJet_electronEnergyFraction"        ,&FatJet_electronEnergyFraction       );
  tree->Branch("FatJet_muonEnergyFraction"        ,&FatJet_muonEnergyFraction       );
  tree->Branch("FatJet_photonEnergyFraction"        ,&FatJet_photonEnergyFraction       );
  tree->Branch("FatJet_neutralMultiplicity"        ,&FatJet_neutralMultiplicity       );
  tree->Branch("FatJet_chargedMultiplicity"        ,&FatJet_chargedMultiplicity       );

  tree->Branch("rho"                            ,&rho2                           );


  //add gen info
  tree->Branch("n_genjet"                          ,&n_genjet                         ,"nGenJets/i");
  tree->Branch("GenJet_pt"                         ,&GenJet_pt                        );
  tree->Branch("GenJet_eta"                        ,&GenJet_eta                       );
  tree->Branch("GenJet_phi"                        ,&GenJet_phi                       );
  tree->Branch("GenJet_mass"                       ,&GenJet_mass                      );

  //add gen info
  tree->Branch("n_genfatjet"                          ,&n_genfatjet                         ,"nGenFatJets/i");
  tree->Branch("GenFatJet_pt"                         ,&GenFatJet_pt                        );
  tree->Branch("GenFatJet_eta"                        ,&GenFatJet_eta                       );
  tree->Branch("GenFatJet_phi"                        ,&GenFatJet_phi                       );
  tree->Branch("GenFatJet_mass"                       ,&GenFatJet_mass                      );

  //add genjet matching indices
  tree->Branch("Jet_genJetIdx"                        ,&Jet_genJetIdx                       );
  tree->Branch("FatJet_genJetAK8Idx"                  ,&FatJet_genJetAK8Idx                 );

  //Scouting PF Candidates
  tree->Branch("nPFCands"            	        ,&n_pfcand 		        ,"nPFCands/i");	
  tree->Branch("nPFMuons"            	        ,&n_pfMu 		        ,"nPFMuons/i");	
  tree->Branch("nPFElectrons"            	        ,&n_pfEl 		        ,"nPFElectrons/i");	
  tree->Branch("PFCands_pt"        	        ,&PFcand_pt 		        );
  tree->Branch("PFCands_eta"            	        ,&PFcand_eta 	                );
  tree->Branch("PFCands_phi"            	        ,&PFcand_phi		        );
  tree->Branch("PFCands_mass"            	        ,&PFcand_m 		        );
  tree->Branch("PFCands_pdgId"                   ,&PFcand_pdgid                  );
  tree->Branch("PFCands_charge"                       ,&PFcand_q                      );
  tree->Branch("PFCands_vertex"                  ,&PFcand_vertex                 );

  //CZZ: added jet indices for Scouting PF candidates
  tree->Branch("FatJetPFCands_jetIdx"                   ,&FatJetPFCands_jetIdx 	                );
  tree->Branch("FatJetPFCands_pFCandsIdx"                   ,&FatJetPFCands_pFCandsIdx 	                );

  //Primary vertices 
  tree->Branch("nPVs"            	        ,&n_pvs                         ,"nPVs/i");	
  tree->Branch("nOfflinePVs"            	        ,&nPV_                         ,"nOfflinePVs/i");	
  tree->Branch("PV_x"        	        ,&Vertex_x  		        );
  tree->Branch("PV_y"                       ,&Vertex_y   	                );
  tree->Branch("PV_z"                       ,&Vertex_z  		        );
  tree->Branch("PV_tracksSize"              ,&Vertex_tracksSize 	        );
  tree->Branch("PV_chi2"                    ,&Vertex_chi2	                );
  tree->Branch("PV_ndof"                    ,&Vertex_ndof	                );
  tree->Branch("PV_isValidVtx"              ,&Vertex_isValidVtx 	        );

  //Other variables
  tree->Branch("ht"                             ,&ht                            );
  tree->Branch("htoff"                             ,&htoff                            );
  tree->Branch("Pileup_nPU"            	        ,&PU_num                        ,"PU_num/i");

  if(!onlyScouting){
    //Offline Electrons
    tree->Branch("nOfflineElectron"               	         ,&n_ele_off                        ,"nOfflineElectron/i");
    tree->Branch("OfflineElectron_pt"                    ,&OffElectron_pt                   );
    tree->Branch("OfflineElectron_eta"                   ,&OffElectron_eta 	                );
    tree->Branch("OfflineElectron_phi"                   ,&OffElectron_phi                  );
    tree->Branch("OfflineElectron_charge"                ,&OffElectron_charge               );
    tree->Branch("OfflineElectron_mass"            	        ,&OffElectron_m                    );
    tree->Branch("OfflineElectron_hoe"                   ,&OffElectron_HoE                  );
    tree->Branch("OfflineElectron_sieie"           ,&OffElectron_sigmaietaieta        );
    tree->Branch("OfflineElectron_dphiin"                ,&OffElectron_dphiin 	        );
    tree->Branch("OfflineElectron_detain"                ,&OffElectron_detain 	        );
    tree->Branch("OfflineElectron_mHits"                 ,&OffElectron_mHits 	        );
    tree->Branch("OfflineElectron_ooEMOop"               ,&OffElectron_ooEMOop              );
    tree->Branch("OfflineElectron_trkiso"                 ,&OffElectron_trkiso 	        );
    tree->Branch("OfflineElectron_ecaliso"               ,&OffElectron_ecaliso              );
    tree->Branch("OfflineElectron_hcaliso"               ,&OffElectron_hcaliso              );
    tree->Branch("OfflineElectron_combinediso"               ,&OffElectron_combinediso   );
    tree->Branch("OfflineElectron_ID"               ,&OffElectron_ID   );
    tree->Branch("OfflineElectron_d0"               ,&OffElectron_d0              );
    tree->Branch("OfflineElectron_dz"               ,&OffElectron_dz              );

    //Offline muons
    tree->Branch("nOfflineMuons"            	        ,&n_mu_off 	                        ,"nOfflineMuons/i");
    tree->Branch("OfflineMuon_pt"                        ,&OffMuon_pt                       );
    tree->Branch("OfflineMuon_eta"                       ,&OffMuon_eta                      );
    tree->Branch("OfflineMuon_phi"                       ,&OffMuon_phi                      );
    tree->Branch("OfflineMuon_mass"                         ,&OffMuon_m                        );
    tree->Branch("OfflineMuon_charge"                    ,&OffMuon_charge	                  );
    tree->Branch("OfflineMuon_isGlobal"              ,&OffMuon_isGlobalMuon             );
    tree->Branch("OfflineMuon_isTracker"             ,&OffMuon_isTrackerMuon            );
    tree->Branch("OfflineMuon_isStandAlone"                   ,&OffMuon_isStandAloneMuon                  );
    tree->Branch("OfflineMuon_isPF"                   ,&OffMuon_isPFMuon                  );
    tree->Branch("OfflineMuon_isCalo"                   ,&OffMuon_isCaloMuon                 );
    

    //Offline AK4 PFJets
    tree->Branch("nJet"            	        ,&n_jetoff                         ,"nOfflineJet/i");
    tree->Branch("nJetId"            	        ,&n_jetIdoff                       ,"nOfflineJetId/i");
    tree->Branch("OfflineJet_pt"            	           ,&OffJet_pt                        );
    tree->Branch("OfflineJet_eta"            	           ,&OffJet_eta                       );
    tree->Branch("OfflineJet_phi"            	           ,&OffJet_phi                       );
    tree->Branch("OfflineJet_mass"            	         ,&OffJet_m                         );
    tree->Branch("OfflineJet_area"            	         ,&OffJet_area                      );
    // Not needed for Type-1 MET
    // tree->Branch("OfflineJet_chargedHadronEnergy"        ,&OffJet_chargedHadronEnergy       );
    // tree->Branch("OfflineJet_neutralHadronEnergy"        ,&OffJet_neutralHadronEnergy       );
    // tree->Branch("OfflineJet_photonEnergy"               ,&OffJet_photonEnergy 	            );
    // tree->Branch("OfflineJet_electronEnergy"             ,&OffJet_electronEnergy            );
    // tree->Branch("OfflineJet_muonEnergy"    	           ,&OffJet_muonEnergy                );
    // tree->Branch("OfflineJet_HFHadronEnergy"             ,&OffJet_HFHadronEnergy            );
    // tree->Branch("OfflineJet_HFEMEnergy"                 ,&OffJet_HFEMEnergy                );
    // tree->Branch("OfflineJet_HOEnergy"                   ,&OffJet_HOEnergy                  );
    //for jet id
    // tree->Branch("OfflineJet_chargedHadronEnergyFraction"        ,&OffJet_chargedHadronEnergyFraction       );
    // tree->Branch("OfflineJet_neutralHadronEnergyFraction"        ,&OffJet_neutralHadronEnergyFraction       );
    // tree->Branch("OfflineJet_chargedEmEnergyFraction"        ,&OffJet_chargedEmEnergyFraction       );
    // tree->Branch("OfflineJet_muonEnergyFraction"        ,&OffJet_muonEnergyFraction       );
    // tree->Branch("OfflineJet_neutralEmEnergyFraction"        ,&OffJet_neutralEmEnergyFraction       );
    // tree->Branch("OfflineJet_electronEnergyFraction"        ,&OffJet_electronEnergyFraction       );
    // tree->Branch("OfflineJet_photonEnergyFraction"        ,&OffJet_photonEnergyFraction       );
    tree->Branch("OfflineJet_nConstituents"              ,&OffJet_nConstituents       );
    // tree->Branch("OfflineJet_neutralMultiplicity"        ,&OffJet_neutralMultiplicity       );
    // tree->Branch("OfflineJet_chargedMultiplicity"        ,&OffJet_chargedMultiplicity       );

    // tree->Branch("OfflineJet_chargedHadronMultiplicity"  ,&OffJet_chargedHadronMultiplicity );
    // tree->Branch("OfflineJet_neutralHadronMultiplicity"  ,&OffJet_neutralHadronMultiplicity );
    // tree->Branch("OfflineJet_photonMultiplicity"         ,&OffJet_photonMultiplicity        );
    // tree->Branch("OfflineJet_electronMultiplicity"       ,&OffJet_electronMultiplicity      );
    // tree->Branch("OfflineJet_muonMultiplicity"           ,&OffJet_muonMultiplicity          );
    // tree->Branch("OfflineJet_HFHadronMultiplicity"       ,&OffJet_HFHadronMultiplicity      );
    // tree->Branch("OfflineJet_HFEMMultiplicity"           ,&OffJet_HFEMMultiplicity          );
    // tree->Branch("OfflineJet_passId"                     ,&OffJet_passId                    );

    //Offline AK4 Puppi Jets
    tree->Branch("nPuppiJet"            	        ,&n_jetoffpuppi                         ,"nOfflinePuppiJet/i");
    tree->Branch("nPuppiJetId"            	        ,&n_jetIdoffpuppi                       ,"nOfflinePuppiJetId/i");
    tree->Branch("OfflinePuppiJet_pt"            	           ,&OffPuppiJet_pt                        );
    tree->Branch("OfflinePuppiJet_eta"            	           ,&OffPuppiJet_eta                       );
    tree->Branch("OfflinePuppiJet_phi"            	           ,&OffPuppiJet_phi                       );
    tree->Branch("OfflinePuppiJet_mass"            	         ,&OffPuppiJet_m                         );
    tree->Branch("OfflinePuppiJet_area"            	         ,&OffPuppiJet_area                      );
    // Not needed for Type-1 MET
    // tree->Branch("OfflinePuppiJet_chargedHadronEnergy"        ,&OffPuppiJet_chargedHadronEnergy       );
    // tree->Branch("OfflinePuppiJet_neutralHadronEnergy"        ,&OffPuppiJet_neutralHadronEnergy       );
    // tree->Branch("OfflinePuppiJet_photonEnergy"               ,&OffPuppiJet_photonEnergy 	            );
    // tree->Branch("OfflinePuppiJet_electronEnergy"             ,&OffPuppiJet_electronEnergy            );
    // tree->Branch("OfflinePuppiJet_muonEnergy"    	           ,&OffPuppiJet_muonEnergy                );
    // tree->Branch("OfflinePuppiJet_HFHadronEnergy"             ,&OffPuppiJet_HFHadronEnergy            );
    // tree->Branch("OfflinePuppiJet_HFEMEnergy"                 ,&OffPuppiJet_HFEMEnergy                );
    // tree->Branch("OfflinePuppiJet_HOEnergy"                   ,&OffPuppiJet_HOEnergy                  );
    //for jet id
    // tree->Branch("OfflinePuppiJet_chargedHadronEnergyFraction"        ,&OffPuppiJet_chargedHadronEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_neutralHadronEnergyFraction"        ,&OffPuppiJet_neutralHadronEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_chargedEmEnergyFraction"        ,&OffPuppiJet_chargedEmEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_muonEnergyFraction"        ,&OffPuppiJet_muonEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_neutralEmEnergyFraction"        ,&OffPuppiJet_neutralEmEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_electronEnergyFraction"        ,&OffPuppiJet_electronEnergyFraction       );
    // tree->Branch("OfflinePuppiJet_photonEnergyFraction"        ,&OffPuppiJet_photonEnergyFraction       );
    tree->Branch("OfflinePuppiJet_nConstituents"              ,&OffPuppiJet_nConstituents       );
    // tree->Branch("OfflinePuppiJet_neutralMultiplicity"        ,&OffPuppiJet_neutralMultiplicity       );
    // tree->Branch("OfflinePuppiJet_chargedMultiplicity"        ,&OffPuppiJet_chargedMultiplicity       );

    // tree->Branch("OfflinePuppiJet_chargedHadronMultiplicity"  ,&OffPuppiJet_chargedHadronMultiplicity );
    // tree->Branch("OfflinePuppiJet_neutralHadronMultiplicity"  ,&OffPuppiJet_neutralHadronMultiplicity );
    // tree->Branch("OfflinePuppiJet_photonMultiplicity"         ,&OffPuppiJet_photonMultiplicity        );
    // tree->Branch("OfflinePuppiJet_electronMultiplicity"       ,&OffPuppiJet_electronMultiplicity      );
    // tree->Branch("OfflinePuppiJet_muonMultiplicity"           ,&OffPuppiJet_muonMultiplicity          );
    // tree->Branch("OfflinePuppiJet_HFHadronMultiplicity"       ,&OffPuppiJet_HFHadronMultiplicity      );
    // tree->Branch("OfflinePuppiJet_HFEMMultiplicity"           ,&OffPuppiJet_HFEMMultiplicity          );
    // tree->Branch("OfflinePuppiJet_passId"                     ,&OffPuppiJet_passId                    );

    
    //CZZ: added Offline AK8 PFJets (built AK8 from Offline PFCands using FastJet)
    tree->Branch("nOfflineFatJet"                       ,&n_fatjet                      ,"nOfflineFatJet/i");
    tree->Branch("OfflineFatJet_area"                    ,&OfflineFatJet_area                   );
    tree->Branch("OfflineFatJet_eta"                     ,&OfflineFatJet_eta                    );
    //tree->Branch("OfflineFatJet_n2b1"                    ,&OfflineFatJet_n2b1                   );
    //tree->Branch("OfflineFatJet_n3b1"                    ,&OfflineFatJet_n3b1                   );
    tree->Branch("OfflineFatJet_phi"                     ,&OfflineFatJet_phi                    );
    tree->Branch("OfflineFatJet_pt"                      ,&OfflineFatJet_pt                     );
    // tree->Branch("OfflineFatJet_tau1"                    ,&OfflineFatJet_tau1                   );
    // tree->Branch("OfflineFatJet_tau2"                    ,&OfflineFatJet_tau2                   );
    // tree->Branch("OfflineFatJet_tau3"                    ,&OfflineFatJet_tau3                   );
    // tree->Branch("OfflineFatJet_tau4"                    ,&OfflineFatJet_tau4                   );
    // tree->Branch("OfflineFatJet_tau21"                   ,&OfflineFatJet_tau21                  );
    // tree->Branch("OfflineFatJet_tau32"                   ,&OfflineFatJet_tau32                  );
    tree->Branch("OfflineFatJet_mass"                    ,&OfflineFatJet_mass                   );
    // tree->Branch("OfflineFatJet_msoftdrop"               ,&OfflineFatJet_msoftdrop              );
    // tree->Branch("OfflineFatJet_mtrim"                   ,&OfflineFatJet_mtrim                  );
    tree->Branch("OfflineFatJet_nConstituents"                  ,&OfflineFatJet_nConstituents                 );
    

    //Offline AK8 Puppi PFJets
    tree->Branch("nOfflinePuppiFatJet"            	             ,&n_fatjet_off_puppi                         ,"nOfflinePuppiFatJet/i");
    tree->Branch("nOfflinePuppiFatJetId"            	           ,&n_fatjetIdoffpuppi                         ,"nOfflinePuppiFatJetId/i");
    tree->Branch("OfflinePuppiFatJet_pt"            	           ,&OffPuppiFatJet_pt                        );
    tree->Branch("OfflinePuppiFatJet_eta"            	           ,&OffPuppiFatJet_eta                       );
    tree->Branch("OfflinePuppiFatJet_phi"            	           ,&OffPuppiFatJet_phi                       );
    tree->Branch("OfflinePuppiFatJet_mass"            	         ,&OffPuppiFatJet_m                         );
    tree->Branch("OfflinePuppiFatJet_area"            	         ,&OffPuppiFatJet_area                      );
    tree->Branch("OfflinePuppiFatJet_chargedHadronEnergy"        ,&OffPuppiFatJet_chargedHadronEnergy       );
    tree->Branch("OfflinePuppiFatJet_neutralHadronEnergy"        ,&OffPuppiFatJet_neutralHadronEnergy       );
    tree->Branch("OfflinePuppiFatJet_photonEnergy"               ,&OffPuppiFatJet_photonEnergy 	            );
    tree->Branch("OfflinePuppiFatJet_electronEnergy"             ,&OffPuppiFatJet_electronEnergy            );
    tree->Branch("OfflinePuppiFatJet_muonEnergy"    	           ,&OffPuppiFatJet_muonEnergy                );
    tree->Branch("OfflinePuppiFatJet_HFHadronEnergy"             ,&OffPuppiFatJet_HFHadronEnergy            );
    tree->Branch("OfflinePuppiFatJet_HFEMEnergy"                 ,&OffPuppiFatJet_HFEMEnergy                );
    tree->Branch("OfflinePuppiFatJet_HOEnergy"                   ,&OffPuppiFatJet_HOEnergy                  );
    tree->Branch("OfflinePuppiFatJet_chargedHadronMultiplicity"  ,&OffPuppiFatJet_chargedHadronMultiplicity );
    tree->Branch("OfflinePuppiFatJet_neutralHadronMultiplicity"  ,&OffPuppiFatJet_neutralHadronMultiplicity );
    tree->Branch("OfflinePuppiFatJet_photonMultiplicity"         ,&OffPuppiFatJet_photonMultiplicity        );
    tree->Branch("OfflinePuppiFatJet_electronMultiplicity"       ,&OffPuppiFatJet_electronMultiplicity      );
    tree->Branch("OfflinePuppiFatJet_muonMultiplicity"           ,&OffPuppiFatJet_muonMultiplicity          );
    tree->Branch("OfflinePuppiFatJet_HFHadronMultiplicity"       ,&OffPuppiFatJet_HFHadronMultiplicity      );
    tree->Branch("OfflinePuppiFatJet_HFEMMultiplicity"           ,&OffPuppiFatJet_HFEMMultiplicity          );
    tree->Branch("OfflinePuppiFatJet_passId"                     ,&OffPuppiFatJet_passId                    );

    //for jet id
    tree->Branch("OfflinePuppiFatJet_chargedHadronEnergyFraction"        ,&OffPuppiFatJet_chargedHadronEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_neutralHadronEnergyFraction"        ,&OffPuppiFatJet_neutralHadronEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_chargedEmEnergyFraction"        ,&OffPuppiFatJet_chargedEmEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_muonEnergyFraction"        ,&OffPuppiFatJet_muonEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_neutralEmEnergyFraction"        ,&OffPuppiFatJet_neutralEmEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_muonEnergyFraction"        ,&OffPuppiFatJet_muonEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_photonEnergyFraction"        ,&OffPuppiFatJet_photonEnergyFraction       );
    tree->Branch("OfflinePuppiFatJet_nConstituents"              ,&OffPuppiFatJet_nConstituents      );
    tree->Branch("OfflinePuppiFatJet_neutralMultiplicity"        ,&OffPuppiFatJet_neutralMultiplicity       );
    tree->Branch("OfflinePuppiFatJet_chargedMultiplicity"        ,&OffPuppiFatJet_chargedMultiplicity       );

    /*
    //offline PF Cands
    tree->Branch("nOfflinePFCands"            	        ,&n_offpfcand 		        ,"nOfflinePFCands/i");	
    tree->Branch("nOfflinePFMuons"            	        ,&n_offpfMu 		        ,"nOfflinePFMuons/i");	
    tree->Branch("nOfflinePFElectrons"            	        ,&n_offpfEl 		        ,"nOfflinePFElectrons/i");	
    tree->Branch("OfflinePFCands_pt"                 ,&OffPFcand_pt     );
    tree->Branch("OfflinePFCands_mass"                 ,&OffPFcand_m     );
    tree->Branch("OfflinePFCands_eta"                ,&OffPFcand_eta    );
    tree->Branch("OfflinePFCands_phi"                ,&OffPFcand_phi    );
    tree->Branch("OfflinePFCands_pdgId"                   ,&OffPFcand_pdgid                  );
    tree->Branch("OfflinePFCands_charge"                       ,&OffPFcand_q                      );

    //CZZ: added jet indices for Offline PF candidates
    tree->Branch("OfflineFatJetPFCands_jetIdx"                   ,&OfflineFatJetPFCands_jetIdx 	                );
    tree->Branch("OfflineFatJetPFCands_pFCandsIdx"                   ,&OfflineFatJetPFCands_pFCandsIdx 	                );
    */

    //offline Puppi PF Cands
    tree->Branch("nOfflinePFCands"            	        ,&n_offpuppipfcand 		      ,"nOfflinePuppiPFCands/i");	
    tree->Branch("OfflinePFCands_pt"                   ,&OffPuppiPFcand_pt                                  );
    tree->Branch("OfflinePFCands_mass"                 ,&OffPuppiPFcand_m                                   );
    tree->Branch("OfflinePFCands_eta"                  ,&OffPuppiPFcand_eta                                 );
    tree->Branch("OfflinePFCands_phi"                  ,&OffPuppiPFcand_phi                                 );
    tree->Branch("OfflinePFCands_pdgId"                ,&OffPuppiPFcand_pdgid                               );
    tree->Branch("OfflinePFCands_charge"               ,&OffPuppiPFcand_q                                   );

    //CZZ: added jet indices for Offline Puppi PF candidates
    tree->Branch("OfflineFatJetPFCands_jetIdx"                   ,&OfflinePuppiFatJetPFCands_jetIdx 	        );
    tree->Branch("OfflineFatJetPFCands_pFCandsIdx"                   ,&OfflinePuppiFatJetPFCands_pFCandsIdx 	);
  }

  //CZZ: added MET collections
  tree->Branch("ScoutMET_pt",&met_pt);
  tree->Branch("ScoutMET_phi",&met_phi);
  if(!onlyScouting){
    tree->Branch("OfflineMET_pt",&met_pt_reco);
    tree->Branch("OfflineMET_phi",&met_phi_reco);
  }
  
  //Gen MET branches
  tree->Branch("genMET_pt",&genMET_pt);
  tree->Branch("genMET_phi",&genMET_phi);



}


ScoutingNanoAOD_fromAOD::~ScoutingNanoAOD_fromAOD() {
}

void ScoutingNanoAOD_fromAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
    

  Handle<reco::VertexCollection> recoverticesH;
  Handle<std::vector<reco::PFJet>> ak4_pfjetsoffH;
  Handle<std::vector<reco::PFJet>> puppi_ak4_pfjetsoffH;
  Handle<std::vector<reco::PFJet>> puppi_ak8_pfjetsoffH;
  Handle<std::vector<reco::GsfElectron> > electronsoffH;
  Handle<edm::View<reco::Muon>> muonsoffH;

  Handle<vector<ScoutingElectron> > electronsH;
  Handle<vector<ScoutingMuon> > muonsH;
  Handle<vector<ScoutingPhoton> > photonsH;
  Handle<vector<ScoutingPFJet> > pfjetsH;
  Handle<vector<ScoutingParticle> > pfcandsH;
  Handle<vector<ScoutingVertex> > verticesH;
  Handle<std::vector<reco::PFCandidate>> pfcandsoffH;

  Handle<vector<reco::GenJet> > genak4jetsH;
  Handle<vector<reco::GenJet> > genak8jetsH;
  Handle<std::vector<reco::GenMET> > genMetH;

  Handle<std::vector<reco::PFMET>> metReco;
  Handle<double> metPt;
  Handle<double> metPhi;


  //define particles
  vector<double> neutralHadrons_ids = {111,130,310,2112};
  vector<double> chargedHadrons_ids = {211,321,999211,2212};
  vector<double> photons_ids = {22};
  vector<double> electrons_ids = {11};
  vector<double> muons_ids = {13};


  if(auto handle = iEvent.getHandle(pfcandsToken)){
    runScouting = true;
  }

  
  if(auto handle = iEvent.getHandle(recoPfCandidateToken)){
    runOffline = true;
  }
  
  if(auto handle = iEvent.getHandle(genak8jetsToken)){
    runGen = true;
  }

  if(runScouting){
    iEvent.getByToken(electronsToken, electronsH);
    iEvent.getByToken(muonsToken, muonsH);
    iEvent.getByToken(photonsToken, photonsH);
    iEvent.getByToken(pfjetsToken, pfjetsH);
    iEvent.getByToken(pfcandsToken, pfcandsH);
    iEvent.getByToken(verticesToken, verticesH);
    iEvent.getByToken(metPtToken, metPt);
    iEvent.getByToken(metPhiToken, metPhi);
  }

  if(runGen){
    iEvent.getByToken(genak4jetsToken, genak4jetsH);
    iEvent.getByToken(genak8jetsToken, genak8jetsH);
    iEvent.getByToken(genMetToken, genMetH);
  }

  if(runOffline){
    iEvent.getByToken(recoverticeToken  , recoverticesH  );
    iEvent.getByToken(recoPfCandidateToken, pfcandsoffH);
    iEvent.getByToken(recoElectronToken, electronsoffH);
    iEvent.getByToken(recoMuonToken, muonsoffH);
    iEvent.getByToken(recoak4PFJetToken, ak4_pfjetsoffH);
    iEvent.getByToken(recoak4PuppiJetToken, puppi_ak4_pfjetsoffH);
    iEvent.getByToken(recoak8PuppiJetToken, puppi_ak8_pfjetsoffH);
    iEvent.getByToken(recoMetToken, metReco);
  }


  Handle<vector<PileupSummaryInfo> > puInfo;
  if(auto handle = iEvent.getHandle(pileupInfoToken2)){
      iEvent.getByToken(pileupInfoToken2, puInfo);
  }
  else {
      iEvent.getByToken(pileupInfoToken, puInfo);
  }
  
  Handle<GenEventInfoProduct > genEvtInfo;
  iEvent.getByToken(genEvtInfoToken, genEvtInfo);

  run = iEvent.eventAuxiliary().run();
  event_ = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();  

  // Which triggers fired
  hltResult_.clear();
  hltResultName_.clear();

  iEvent.getByToken(triggerBits_, triggerBits);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  scouting_trig=0; 
  scouting_trig_prescaled=0;
  scouting_trig_zero_bias=0;
  offline_trig=0; 
  veto_trig=0; 
  for(size_t j = 0; j < hltSeeds_.size(); j++){
        TPRegexp pattern(hltSeeds_[j]);
        TPRegexp pattern1("DST_HT410_PFScouting_v");
        TPRegexp pattern2("HLT_IsoMu24_v*");
        TPRegexp pattern3("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        TPRegexp pattern4("HLT_Ele32_WPTight_Gsf_v*");
        TPRegexp pattern5("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        TPRegexp pattern6("HLT_PFHT1050_v*");
        TPRegexp pattern7("DST_L1HTT_CaloScouting_PFScouting_v");
        TPRegexp pattern8("DST_CaloJet40_CaloScouting_PFScouting_v*");
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      const std::string& hltbitName = names.triggerName(i);
      std::string hltpathName = hltbitName;
      bool hltpassFinal = triggerBits->accept(i);
        if( 
          TString(hltpathName).Contains(pattern1) and hltpassFinal)
          {
          scouting_trig=1;
          }
        //these are all events passing the prescaled trigger used to measure the reference trigger efficiency
        if( 
          TString(hltpathName).Contains(pattern7) and hltpassFinal)
          {
          scouting_trig_prescaled=1;
          }
        //these is a zero bias trigger to measure the L1 HT seed trigger efficiency
        if(
          TString(hltpathName).Contains(pattern8) and hltpassFinal)
          {
          scouting_trig_zero_bias=1;
          }
        if( 
          TString(hltpathName).Contains(pattern6) and hltpassFinal)
          {
          offline_trig=1;
          }
        if( hltpassFinal and (
          TString(hltpathName).Contains(pattern2)
          or TString(hltpathName).Contains(pattern3)
          or TString(hltpathName).Contains(pattern4)
          or TString(hltpathName).Contains(pattern5)
          ) 
          ){
          veto_trig=1;
        } 
        if( TString(hltpathName).Contains(pattern)){
          hltResult_.push_back(hltpassFinal);
          hltResultName_.push_back(hltbitName);
        }
      }  
  }


  // *
  // Electrons here
  // *

  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_d0.clear();
  Electron_dz.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_HoE.clear();
  Electron_ooEMOop.clear();
  Electron_mHits.clear();
  Electron_charge.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_trkiso.clear();
  Electron_combinediso.clear();
  Electron_ID.clear();
  Electron_chargedMiniIso.clear();
  Electron_combinedMiniIso.clear();
  n_ele = 0;

  vector<ScoutingParticle> PFcands;
  PFcands.clear();


  if(runScouting){
    for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) 
      {
        Electron_pt.push_back(electrons_iter->pt());
        Electron_eta.push_back(electrons_iter->eta());
        Electron_phi.push_back(electrons_iter->phi());	
        Electron_m.push_back(electrons_iter->m());
        Electron_detain.push_back(electrons_iter->dEtaIn());
        Electron_dphiin.push_back(electrons_iter->dPhiIn());
        Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
        Electron_HoE.push_back(electrons_iter->hOverE());	
        Electron_ooEMOop.push_back(electrons_iter->ooEMOop());
        Electron_mHits.push_back(electrons_iter->missingHits());
        Electron_charge.push_back(electrons_iter->charge());
        Electron_trkiso.push_back(electrons_iter->trackIso());
        Electron_ecaliso.push_back(electrons_iter->ecalIso());
        Electron_hcaliso.push_back(electrons_iter->hcalIso());
        Electron_d0.push_back(electrons_iter->d0());
        Electron_dz.push_back(electrons_iter->dz());
        n_ele++;

        //Notice: here pushing back the electrons in the PFcands vector (for scouting electrons are not stored as PFcands, but as a separate collection)
        ScoutingParticle tmp(electrons_iter->pt(),electrons_iter->eta(),electrons_iter->phi(),electrons_iter->m(),(-11)*electrons_iter->charge(),0);
        PFcands.push_back(tmp);
        TLorentzVector electron_p4 = TLorentzVector();
        electron_p4.SetPtEtaPhiM(electrons_iter->pt(), electrons_iter->eta(), electrons_iter->phi(), electrons_iter->m());
        float combinediso; 
        bool electronID = false;
        if(abs(electrons_iter->eta())<1.479){
        combinediso = (electrons_iter->trackIso() + std::max(0.f,electrons_iter->ecalIso() -1.f) + electrons_iter->hcalIso()) / electron_p4.Et();
          electronID = 
          (fabs(electrons_iter->dEtaIn()) < 0.007)
          & (fabs(electrons_iter->dPhiIn()) < 0.15)
          & (electrons_iter->sigmaIetaIeta() < 0.01)
          & (electrons_iter->hOverE() < 0.12)
          & (fabs(electrons_iter->d0()) < 0.02)
          & (fabs(electrons_iter->dz()) < 0.2)
          & (electrons_iter->ooEMOop() < 0.05)
          & (combinediso < 0.15);
        }else{
        combinediso = (electrons_iter->trackIso() + electrons_iter->ecalIso()  + electrons_iter->hcalIso()) / electron_p4.Et();
          electronID = 
          (fabs(electrons_iter->dEtaIn()) < 0.009)
          & (fabs(electrons_iter->dPhiIn()) < 0.10)
          & (electrons_iter->sigmaIetaIeta() < 0.03)
          & (electrons_iter->hOverE() < 0.10)
          & (fabs(electrons_iter->d0()) < 0.02)
          & (fabs(electrons_iter->dz()) < 0.2)
          & (electrons_iter->ooEMOop() < 0.05)
          & (combinediso < 0.15);
        }
        Electron_combinediso.push_back(combinediso);
        Electron_ID.push_back(electronID);
    }
  }


  OffElectron_pt.clear();
  OffElectron_eta.clear();
  OffElectron_phi.clear();
  OffElectron_m.clear();
  OffElectron_d0.clear();
  OffElectron_dz.clear();
  OffElectron_detain.clear();
  OffElectron_dphiin.clear();
  OffElectron_sigmaietaieta.clear();
  OffElectron_HoE.clear();
  OffElectron_ooEMOop.clear();
  OffElectron_mHits.clear();
  OffElectron_charge.clear();
  OffElectron_ecaliso.clear();
  OffElectron_hcaliso.clear();
  OffElectron_trkiso.clear();
  OffElectron_combinediso.clear();
  OffElectron_ID.clear();
  n_ele_off = 0;

  double ooEmooP_ = 1e30;
  double relIso = 0.0;

  if (runOffline){
    if (recoverticesH->empty()) return; // skip the event if no PV found
    //const reco::Vertex &pv = vertices->front();
    
    nPV_ = recoverticesH -> size();

    VertexCollection::const_iterator firstGoodVertex = recoverticesH->end();
    int firstGoodVertexIdx = 0;
    for (VertexCollection::const_iterator vtx = recoverticesH->begin(); 
        vtx != recoverticesH->end(); ++vtx, ++firstGoodVertexIdx) {
      // The "good vertex" selection is borrowed from Giovanni Zevi Della Porta
      // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
      // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
      if (  /*!vtx->isFake() &&*/ 
        !(vtx->chi2()==0 && vtx->ndof()==0) 
        &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
        && std::abs(vtx->position().Z())<=24.0) {
            firstGoodVertex = vtx;
            break;
      }
    }
    
    if ( firstGoodVertex==recoverticesH->end() )
      return; // skip event if there are no good PVs


    for (auto electronsoff_iter = electronsoffH->begin(); electronsoff_iter != electronsoffH->end(); ++electronsoff_iter) 
      {
        
        OffElectron_pt.push_back(electronsoff_iter->pt());
        OffElectron_eta.push_back(electronsoff_iter->eta());
        OffElectron_phi.push_back(electronsoff_iter->phi());	
        OffElectron_m.push_back(electronsoff_iter->mass());
        OffElectron_charge.push_back(electronsoff_iter->charge());
        OffElectron_detain.push_back(electronsoff_iter->deltaEtaSuperClusterTrackAtVtx());
        OffElectron_dphiin.push_back(electronsoff_iter->deltaPhiSuperClusterTrackAtVtx());
        OffElectron_sigmaietaieta.push_back(electronsoff_iter->full5x5_sigmaIetaIeta());
        OffElectron_HoE.push_back(electronsoff_iter->hcalOverEcal());	
      
        if( electronsoff_iter -> ecalEnergy() == 0 ){
          printf("Electron energy is zero!\n");
          ooEmooP_ = 1e30;
        }else if( !std::isfinite(electronsoff_iter -> ecalEnergy())){
          printf("Electron energy is not finite!\n");
          ooEmooP_ = 1e30;
        }else{
          ooEmooP_ = std::abs(1.0/electronsoff_iter -> ecalEnergy() - electronsoff_iter -> eSuperClusterOverP()/electronsoff_iter -> ecalEnergy() );
        }

        OffElectron_ooEMOop.push_back(ooEmooP_);
        OffElectron_mHits.push_back(electronsoff_iter->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));


        // Isolation
        relIso = ((electronsoff_iter->pfIsolationVariables()).sumChargedHadronPt +  (electronsoff_iter->pfIsolationVariables()).sumNeutralHadronEt + (electronsoff_iter->pfIsolationVariables()).sumPhotonEt)/electronsoff_iter->pt();

        OffElectron_d0.push_back((-1) * electronsoff_iter-> gsfTrack()->dxy(firstGoodVertex->position()));
        OffElectron_dz.push_back(electronsoff_iter -> gsfTrack()->dz( firstGoodVertex->position() ));
        n_ele_off++;

        //Notice, applying tghe same ID as for the scouting electrons for comparison
        TLorentzVector electronoff_p4 = TLorentzVector();
        electronoff_p4.SetPtEtaPhiM(electronsoff_iter->pt(), electronsoff_iter->eta(), electronsoff_iter->phi(), electronsoff_iter->mass());
        bool electronID_off = false;
        if(abs(electronsoff_iter->eta())<1.479){
          electronID_off = 
          (fabs(electronsoff_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.007)
          & (fabs(electronsoff_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.15)
          & (electronsoff_iter->full5x5_sigmaIetaIeta() < 0.01)
          & (electronsoff_iter->hcalOverEcal() < 0.12)
          & (fabs((-1) * electronsoff_iter-> gsfTrack()->dxy(firstGoodVertex->position())) < 0.02)
          & (fabs(electronsoff_iter -> gsfTrack()->dz( firstGoodVertex->position() )) < 0.2)
          & (ooEmooP_ < 0.05)
          & (relIso < 0.15);
        }else{
          electronID_off = 
          (fabs(electronsoff_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.009)
          & (fabs(electronsoff_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.10)
          & (electronsoff_iter->full5x5_sigmaIetaIeta() < 0.03)
          & (electronsoff_iter->hcalOverEcal() < 0.10)
          & (fabs((-1) * electronsoff_iter-> gsfTrack()->dxy(firstGoodVertex->position())) < 0.02)
          & (fabs(electronsoff_iter -> gsfTrack()->dz( firstGoodVertex->position() )) < 0.2)
          & (ooEmooP_ < 0.05)
          & (relIso < 0.15);
        }
        OffElectron_combinediso.push_back(relIso);
        OffElectron_ID.push_back(electronID_off);
    }

  }


  // *
  // Photons here
  // *
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_m.clear();
  Photon_sigmaietaieta.clear();
  Photon_HoE.clear();
  Photon_ecaliso.clear();
  Photon_hcaliso.clear();
  n_pho = 0;

  if(runScouting){
    for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
      Photon_pt.push_back(photons_iter->pt());
      Photon_eta.push_back(photons_iter->eta());
      Photon_phi.push_back(photons_iter->phi());
      Photon_m.push_back(photons_iter->m());
      Photon_sigmaietaieta.push_back(photons_iter->sigmaIetaIeta());
      Photon_HoE.push_back(photons_iter->hOverE());
      Photon_ecaliso.push_back(photons_iter->ecalIso());
      Photon_hcaliso.push_back(photons_iter->hcalIso());
      
      n_pho++;
    }
  }

  // *
  // Primary vertices
  // * 
  n_pvs = 0;
  Vertex_x.clear();
  Vertex_y.clear();
  Vertex_z.clear();
  Vertex_tracksSize.clear();
  Vertex_chi2.clear();
  Vertex_ndof.clear();
  Vertex_isValidVtx.clear();
  if(runScouting){
  for (auto vertices_iter = verticesH->begin(); vertices_iter != verticesH->end(); ++vertices_iter) {
        Vertex_x.push_back( vertices_iter->x() );
        Vertex_y.push_back( vertices_iter->y() );
        Vertex_z.push_back( vertices_iter->z() );
        Vertex_tracksSize.push_back( vertices_iter->tracksSize() );
        Vertex_chi2.push_back( vertices_iter->chi2() );
        Vertex_ndof.push_back( vertices_iter->ndof() );
        Vertex_isValidVtx.push_back( vertices_iter->isValidVtx() );
        n_pvs++;
    }
  }


  if (!doData) {
    for(auto PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI){
      int pu_bunchcrossing = PVI->getBunchCrossing();
      if(pu_bunchcrossing ==0){
        PU_num = PVI->getTrueNumInteractions();
      }
    }
  }


  // * 
  // Particle Flow candidates 
  // *

  PFcand_pt.clear();
  PFcand_eta.clear();
  PFcand_phi.clear();
  PFcand_m.clear();
  PFcand_pdgid.clear();
  PFcand_q.clear();
  PFcand_vertex.clear();
  FatJetPFCands_jetIdx.clear();
  FatJetPFCands_pFCandsIdx.clear();

  vector<PseudoJet> fj_part;
  n_pfcand = 0;
  n_pfMu =0;
  n_pfEl =0;


  if(runScouting){

    for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter) {
      ScoutingParticle tmp(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->pt())),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->eta())),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->phi())),pfcands_iter->m(),pfcands_iter->pdgId(),pfcands_iter->vertex());
      PFcands.push_back(tmp);
    }
    

    //sort PFcands according to pT
    struct {
      bool operator()(ScoutingParticle a, ScoutingParticle b) const { return a.pt() > b.pt(); }
    } custompT;

    std::sort(PFcands.begin(), PFcands.end(), custompT);


    for(auto & pfcands_iter : PFcands ){ //fills PFcand track info
      
      if (pfcands_iter.pt() < 0.5) continue;
      if (abs(pfcands_iter.eta()) >= 2.4 ) continue;

      PFcand_pt.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.pt())));
      PFcand_eta.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.eta())));
      PFcand_phi.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.phi())));
      if(abs(pfcands_iter.pdgId()) == 13){
        n_pfMu ++;
      }
      if(abs(pfcands_iter.pdgId()) == 11){
        n_pfEl ++;
      }
    
      PFcand_m.push_back(pfcands_iter.m());
      PFcand_pdgid.push_back(pfcands_iter.pdgId());
      PFcand_q.push_back(getCharge(pfcands_iter.pdgId()));
      PFcand_vertex.push_back(pfcands_iter.vertex());

      // For clustering fat jets
      PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);
      temp_jet.reset_PtYPhiM(pfcands_iter.pt(), pfcands_iter.eta(), pfcands_iter.phi(), pfcands_iter.m());
      temp_jet.set_user_index(n_pfcand);
      temp_jet.set_user_info(new PdgIdInfo(pfcands_iter.pdgId()));
      fj_part.push_back(temp_jet);
      
      n_pfcand++;
    }

  }


OffPFcand_pt.clear();
OffPFcand_eta.clear();
OffPFcand_phi.clear();
OffPFcand_m.clear();
OffPFcand_pdgid.clear();
OffPFcand_q.clear();
n_offpfcand = 0;
n_offpfMu =0;
n_offpfEl =0;
vector<PseudoJet> off_fj_part;


if(runOffline){

    for(auto pfcandsoff_iter = pfcandsoffH->begin(); pfcandsoff_iter != pfcandsoffH->end(); ++pfcandsoff_iter ){ //fills PFcand track info

      if (pfcandsoff_iter->pt() < 0.5) continue;
      if (abs(pfcandsoff_iter->eta()) >= 2.4 ) continue;

      OffPFcand_pt.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcandsoff_iter->pt())));
      OffPFcand_eta.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcandsoff_iter->eta())));
      OffPFcand_phi.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcandsoff_iter->phi())));

      if(abs(pfcandsoff_iter->pdgId()) == 13){
      n_offpfMu ++;
      }
      if(abs(pfcandsoff_iter->pdgId()) == 11){
      n_offpfEl ++;
      }
    
      OffPFcand_m.push_back(pfcandsoff_iter->mass());
      OffPFcand_pdgid.push_back(pfcandsoff_iter->pdgId());
      OffPFcand_q.push_back(getCharge(pfcandsoff_iter->pdgId()));

      // For clustering fat jets
      PseudoJet temp_offjet = PseudoJet(0, 0, 0, 0);
      temp_offjet.reset_PtYPhiM(pfcandsoff_iter->pt(), pfcandsoff_iter->eta(), pfcandsoff_iter->phi(), pfcandsoff_iter->mass());
      temp_offjet.set_user_index(n_offpfcand);
      off_fj_part.push_back(temp_offjet);

      n_offpfcand++;
    }

}


  // 
  // Scouting muons   
  //

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_chi2.clear();
  Muon_isGlobalMuon.clear();
  Muon_isTrackerMuon.clear();
  Muon_ndof.clear();
  Muon_charge.clear();
  Muon_dxy.clear();
  Muon_dz.clear();
  Muon_nvalidmuon_hits.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nmatchedstations.clear();
  Muon_type.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkpterror.clear();
  Muon_trkphierror.clear();
  Muon_trketaerror.clear();
  Muon_trkdszerror.clear();
  Muon_trkdsz.clear();
  Muon_chargedMiniIso.clear();
  Muon_combinedMiniIso.clear();
  n_mu=0;  
  
  if(runScouting){
  for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
 	  Muon_pt.push_back(muons_iter->pt()); 
   	Muon_eta.push_back(muons_iter->eta());
   	Muon_phi.push_back(muons_iter->phi());
   	Muon_m.push_back(muons_iter->m());
   	Muon_ecaliso.push_back(muons_iter->ecalIso());
   	Muon_hcaliso.push_back(muons_iter->hcalIso());
   	Muon_trkiso.push_back(muons_iter->trackIso());
   	Muon_chi2.push_back(muons_iter->chi2());
   	Muon_ndof.push_back(muons_iter->ndof());
   	Muon_charge.push_back(muons_iter->charge());
   	Muon_dxy.push_back(muons_iter->dxy());
   	Muon_dz.push_back(muons_iter->dz());
   	Muon_nvalidmuon_hits.push_back(muons_iter->nValidMuonHits());
   	Muon_nvalidpixelhits.push_back(muons_iter->nValidPixelHits());
   	Muon_nmatchedstations.push_back(muons_iter->nMatchedStations());
    Muon_type.push_back(muons_iter->type());
    Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
    Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
    Muon_trklambda.push_back(muons_iter->trk_lambda());
    Muon_trkpt.push_back(muons_iter->trk_pt());
    Muon_trkphi.push_back(muons_iter->trk_phi());
    Muon_trketa.push_back(muons_iter->trk_eta());
    Muon_trkqoverperror.push_back(muons_iter->dxyError());
    Muon_trklambdaerror.push_back(muons_iter->dzError());
    Muon_trkpterror.push_back(muons_iter->trk_qoverpError());
    Muon_trkphierror.push_back(muons_iter->trk_lambdaError());
    Muon_trketaerror.push_back(muons_iter->trk_phiError());
    Muon_trkdszerror.push_back(muons_iter->trk_dsz());
    Muon_trkdsz.push_back(muons_iter->trk_dszError());
    Muon_isGlobalMuon.push_back(muons_iter->isGlobalMuon());
    Muon_isTrackerMuon.push_back(muons_iter->isTrackerMuon());
    n_mu++;
  }
}

  // 
  // Offline muons   
  //

  OffMuon_pt.clear();
  OffMuon_eta.clear();
  OffMuon_phi.clear();
  OffMuon_m.clear();
  OffMuon_isGlobalMuon.clear();
  OffMuon_isTrackerMuon.clear();
  OffMuon_isStandAloneMuon.clear();
  OffMuon_isCaloMuon.clear();
  OffMuon_isPFMuon.clear();
  OffMuon_charge.clear();
  n_mu_off=0;  
  
  if(runOffline){
    for (auto muonsoff_iter = muonsoffH->begin(); muonsoff_iter != muonsoffH->end(); ++muonsoff_iter) {
      OffMuon_pt.push_back(muonsoff_iter->pt()); 
      OffMuon_eta.push_back(muonsoff_iter->eta());
      OffMuon_phi.push_back(muonsoff_iter->phi());
      OffMuon_m.push_back(muonsoff_iter->mass());
      OffMuon_charge.push_back(muonsoff_iter->charge());
      OffMuon_isGlobalMuon.push_back(muonsoff_iter->isGlobalMuon());
      OffMuon_isTrackerMuon.push_back(muonsoff_iter->isTrackerMuon());
      OffMuon_isStandAloneMuon.push_back(muonsoff_iter->isStandAloneMuon());
      OffMuon_isCaloMuon.push_back(muonsoff_iter->isCaloMuon());
      OffMuon_isPFMuon.push_back(muonsoff_iter->isPFMuon());
      n_mu_off++;
  }
}

  //* loop over electrons to compute to compute the mini-isolation
  if(runScouting){
    for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) 
      {
        TLorentzVector electron_p4 = TLorentzVector();
        electron_p4.SetPtEtaPhiM(electrons_iter->pt(), electrons_iter->eta(), electrons_iter->phi(), electrons_iter->m());
        reco::Candidate::PolarLorentzVector electron_p4_reco(electron_p4.Pt(), electron_p4.Eta(), electron_p4.Phi(), electron_p4.M());
        pat::PFIsolation ele_miniiso =  getMiniPFIsolationScout(PFcands, electron_p4_reco);
        auto chg = ele_miniiso.chargedHadronIso();
        auto neu = ele_miniiso.neutralHadronIso();
        auto pho = ele_miniiso.photonIso();
        float scale = 1.0 / electrons_iter->pt();
        Electron_chargedMiniIso.push_back(scale * chg);
        Electron_combinedMiniIso.push_back(scale * (chg + std::max(0.0, static_cast<double>(neu + pho) )));
        
    }
  }

  // * Loop over Muons and compute mini-isolaton
  if(runScouting){
    for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) 
      {
        TLorentzVector muon_p4 = TLorentzVector();
        muon_p4.SetPtEtaPhiM(muons_iter->pt(), muons_iter->eta(), muons_iter->phi(), muons_iter->m());
        reco::Candidate::PolarLorentzVector muon_p4_reco(muon_p4.Pt(), muon_p4.Eta(), muon_p4.Phi(), muon_p4.M());
        pat::PFIsolation mu_miniiso =  getMiniPFIsolationScout(PFcands, muon_p4_reco);
        auto chg = mu_miniiso.chargedHadronIso();
        auto neu = mu_miniiso.neutralHadronIso();
        auto pho = mu_miniiso.photonIso();
        float scale = 1.0 / muons_iter->pt();
        Muon_chargedMiniIso.push_back(scale * chg);
        Muon_combinedMiniIso.push_back(scale * (chg + std::max(0.0, static_cast<double>(neu + pho) )));
        
    }
  }


  // * 
  // AK4 Jets 
  // * 
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_m.clear();
  Jet_area.clear();
  // Not needed for Type-1 MET
  // Jet_chargedHadronEnergy.clear();
  // Jet_neutralHadronEnergy.clear();
  // Jet_photonEnergy.clear();
  // Jet_electronEnergy.clear();
  // Jet_muonEnergy.clear();
  // Jet_HFHadronEnergy.clear();
  // Jet_HFEMEnergy.clear();
  // Jet_HOEnergy.clear();
  //
  // // for jet id
  // Jet_chargedHadronEnergyFraction.clear();
  // Jet_neutralHadronEnergyFraction.clear();
  // Jet_electronEnergyFraction.clear();
  // Jet_photonEnergyFraction.clear();
  // Jet_muonEnergyFraction.clear();
  // Jet_chargedMultiplicity.clear();
  // Jet_neutralMultiplicity.clear();
  //
  // Jet_chargedHadronMultiplicity.clear();
  // Jet_neutralHadronMultiplicity.clear();
  // Jet_photonMultiplicity.clear();
  // Jet_electronMultiplicity.clear();
  // Jet_muonMultiplicity.clear();
  // Jet_HFHadronMultiplicity.clear();
  // Jet_HFEMMultiplicity.clear();
  // Jet_csv.clear();
  // Jet_mvaDiscriminator.clear();
  Jet_nConstituents.clear();
  // Jet_passId.clear();
  Jet_genJetIdx.clear();
  
  //OffJet (regular PF jets)
  OffJet_pt.clear();
  OffJet_eta.clear();
  OffJet_phi.clear();
  OffJet_m.clear();
  OffJet_area.clear();
  // Not needed for Type-1 MET
  // OffJet_chargedHadronEnergy.clear();
  // OffJet_neutralHadronEnergy.clear();
  // OffJet_photonEnergy.clear();
  // OffJet_electronEnergy.clear();
  // OffJet_muonEnergy.clear();
  // OffJet_HFHadronEnergy.clear();
  // OffJet_HFEMEnergy.clear();
  // OffJet_HOEnergy.clear();
  //
  // // for jet id
  // OffJet_chargedHadronEnergyFraction.clear();
  // OffJet_neutralHadronEnergyFraction.clear();
  // OffJet_chargedEmEnergyFraction.clear();
  // OffJet_neutralEmEnergyFraction.clear();
  // OffJet_muonEnergyFraction.clear();
  // OffJet_electronEnergyFraction.clear();
  // OffJet_photonEnergyFraction.clear();
  // OffJet_chargedMultiplicity.clear();
  // OffJet_neutralMultiplicity.clear();
  OffJet_nConstituents.clear();
  //
  // OffJet_chargedHadronMultiplicity.clear();
  // OffJet_neutralHadronMultiplicity.clear();
  // OffJet_photonMultiplicity.clear();
  // OffJet_electronMultiplicity.clear();
  // OffJet_muonMultiplicity.clear();
  // OffJet_HFHadronMultiplicity.clear();
  // OffJet_HFEMMultiplicity.clear();
  // OffJet_passId.clear();
  
  //OffPuppiJet (Puppi jets)
  OffPuppiJet_pt.clear();
  OffPuppiJet_eta.clear();
  OffPuppiJet_phi.clear();
  OffPuppiJet_m.clear();
  OffPuppiJet_area.clear();
  // Not needed for Type-1 MET
  // OffPuppiJet_chargedHadronEnergy.clear();
  // OffPuppiJet_neutralHadronEnergy.clear();
  // OffPuppiJet_photonEnergy.clear();
  // OffPuppiJet_electronEnergy.clear();
  // OffPuppiJet_muonEnergy.clear();
  // OffPuppiJet_HFHadronEnergy.clear();
  // OffPuppiJet_HFEMEnergy.clear();
  // OffPuppiJet_HOEnergy.clear();
  //
  // // for jet id
  // OffPuppiJet_chargedHadronEnergyFraction.clear();
  // OffPuppiJet_neutralHadronEnergyFraction.clear();
  // OffPuppiJet_chargedEmEnergyFraction.clear();
  // OffPuppiJet_neutralEmEnergyFraction.clear();
  // OffPuppiJet_muonEnergyFraction.clear();
  // OffPuppiJet_electronEnergyFraction.clear();
  // OffPuppiJet_photonEnergyFraction.clear();
  // OffPuppiJet_chargedMultiplicity.clear();
  // OffPuppiJet_neutralMultiplicity.clear();
  OffPuppiJet_nConstituents.clear();
  //
  // OffPuppiJet_chargedHadronMultiplicity.clear();
  // OffPuppiJet_neutralHadronMultiplicity.clear();
  // OffPuppiJet_photonMultiplicity.clear();
  // OffPuppiJet_electronMultiplicity.clear();
  // OffPuppiJet_muonMultiplicity.clear();
  // OffPuppiJet_HFHadronMultiplicity.clear();
  // OffPuppiJet_HFEMMultiplicity.clear();
  // OffPuppiJet_passId.clear();
  
  n_jet = 0;
  n_jetId = 0;
  n_jetIdoff = 0;
  n_jetoff = 0;
  n_jetIdoffpuppi = 0;
  n_jetoffpuppi = 0;
  ht = 0;
  htoff = 0;
  htoffpuppi = 0;
  passJetId = false;

  if(runScouting){
    for (auto pfjet = pfjetsH->begin(); pfjet != pfjetsH->end(); ++pfjet) {

      // Calculate HT at start of loop
      passJetId = jetID(*pfjet);
      if ( passJetId == true && pfjet->pt() >= 30) {
        ht += pfjet->pt();
        n_jetId++;
      }

      // Create LorentzVector for energy calculations
      math::PtEtaPhiMLorentzVector j(pfjet->pt(), pfjet->eta(), pfjet->phi(), pfjet->m());
      
      // Check EM energy fraction < 0.9 (approximated with electronEnergy + photonEnergy since chargedEmEnergy and neutralEmEnergy are not available in scouting)
      double emEnergyFraction = (pfjet->electronEnergy() + pfjet->photonEnergy()) / j.E();
      if (emEnergyFraction >= 0.9) continue;

      // Remove muon overlap (avoids double counting muon energy in MET calculation for Type-1 MET)
      math::PtEtaPhiMLorentzVector correctedJetP4 = j;
      const std::vector<int> & constituentIndices = pfjet->constituents();
      for (const auto & idx : constituentIndices) {
        if (idx >= 0 && idx < (int)pfcandsH->size()) {
          const auto & cand = (*pfcandsH)[idx];
          if (abs(cand.pdgId()) == 13) { // muon
            math::PtEtaPhiMLorentzVector muonP4(cand.pt(), cand.eta(), cand.phi(), cand.m());
            correctedJetP4 -= muonP4;
          }
        }
      }

      // Apply pt cut on corrected jet pt
      if (correctedJetP4.pt() < jetAK4ScoutPtMin) continue;

      // Save jet variables with corrected momentum
      Jet_pt .push_back( correctedJetP4.pt() );
      Jet_eta.push_back( correctedJetP4.eta());
      Jet_phi.push_back( correctedJetP4.phi());
      Jet_m  .push_back( correctedJetP4.mass()  );

      Jet_area.push_back( pfjet->jetArea());

      // Not needed for Type-1 MET
      // Jet_chargedHadronEnergy.push_back( pfjet->chargedHadronEnergy());
      // Jet_neutralHadronEnergy.push_back( pfjet->neutralHadronEnergy());
      // Jet_photonEnergy       .push_back( pfjet->photonEnergy()       );
      // Jet_electronEnergy     .push_back( pfjet->electronEnergy()     );
      // Jet_muonEnergy         .push_back( pfjet->muonEnergy()     );
      // Jet_HFHadronEnergy     .push_back( pfjet->HFHadronEnergy() );
      // Jet_HFEMEnergy         .push_back( pfjet->HFEMEnergy()     );
      // Jet_HOEnergy           .push_back( pfjet->HOEnergy()       );
      // 
      // Jet_chargedHadronMultiplicity.push_back( pfjet->chargedHadronMultiplicity());
      // Jet_neutralHadronMultiplicity.push_back( pfjet->neutralHadronMultiplicity());
      // Jet_photonMultiplicity       .push_back( pfjet->photonMultiplicity()       );
      // Jet_electronMultiplicity     .push_back( pfjet->electronMultiplicity()     );
      // Jet_muonMultiplicity         .push_back( pfjet->muonMultiplicity()         );
      // Jet_HFHadronMultiplicity     .push_back( pfjet->HFHadronMultiplicity()     );
      // Jet_HFEMMultiplicity         .push_back( pfjet->HFEMMultiplicity()         );
      //
      // Jet_chargedHadronEnergyFraction.push_back( pfjet->chargedHadronEnergy()/j.E());
      // Jet_neutralHadronEnergyFraction.push_back( pfjet->neutralHadronEnergy()/j.E());
      // Jet_electronEnergyFraction.push_back( pfjet->electronEnergy()/j.E());
      // Jet_photonEnergyFraction.push_back( pfjet->photonEnergy()/j.E());
      // Jet_muonEnergyFraction.push_back( pfjet->muonEnergy()/j.E());
      // Jet_chargedMultiplicity.push_back( pfjet->muonMultiplicity() + pfjet->electronMultiplicity() + pfjet->chargedHadronMultiplicity());
      // Jet_neutralMultiplicity.push_back( pfjet->photonMultiplicity() + pfjet->neutralHadronMultiplicity());
      //
      // Jet_csv             .push_back( pfjet->csv() );
      // Jet_mvaDiscriminator.push_back( pfjet->mvaDiscriminator()    );
      Jet_nConstituents   .push_back( pfjet->constituents().size() );
      // Jet_passId.push_back( passJetId );
      
      n_jet++; 

      }
  }

  // Fill OffJet (regular PF jets)
  if(runOffline){
    for (auto it = ak4_pfjetsoffH->begin(); it != ak4_pfjetsoffH->end(); ++it) {
      auto pfjet = *it;

      // Calculate HT at start of loop
      passJetId = jetIDoff(pfjet);
      if ( passJetId == true && pfjet.pt() >= 30) {
        htoff += pfjet.pt();
        n_jetIdoff++;
      }

      // Check EM energy fraction < 0.9
      double emEnergyFraction = (pfjet.chargedEmEnergy() + pfjet.neutralEmEnergy()) / pfjet.energy();
      if (emEnergyFraction >= 0.9) continue;

      // Remove muon overlap
      reco::Candidate::LorentzVector correctedJetP4 = pfjet.p4();
      const std::vector<reco::CandidatePtr> & cands = pfjet.getJetConstituents();
      for (const auto & cand : cands) {
        const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand.get());
        const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand.get());
        if ( mu != 0 && abs(mu->pdgId()) == 13 ) {
          reco::Candidate::LorentzVector muonP4 = cand->p4();
          correctedJetP4 -= muonP4;
        }
      }

      // Apply pt cut on corrected jet pt
      if (correctedJetP4.pt() < jetAK4PtMin) continue;

      // Save jet variables with corrected momentum
      OffJet_pt .push_back( correctedJetP4.pt() );
      OffJet_eta.push_back( correctedJetP4.eta());
      OffJet_phi.push_back( correctedJetP4.phi());
      OffJet_m  .push_back( correctedJetP4.mass() );

      OffJet_area.push_back( pfjet.jetArea());

      // Not needed for Type-1 MET
      // OffJet_chargedHadronEnergy.push_back( pfjet.chargedHadronEnergy());
      // OffJet_neutralHadronEnergy.push_back( pfjet.neutralHadronEnergy());
      // OffJet_photonEnergy       .push_back( pfjet.photonEnergy()       );
      // OffJet_electronEnergy     .push_back( pfjet.electronEnergy()     );
      // OffJet_muonEnergy         .push_back( pfjet.muonEnergy()     );
      // OffJet_HFHadronEnergy     .push_back( pfjet.HFHadronEnergy() );
      // OffJet_HFEMEnergy         .push_back( pfjet.HFEMEnergy()     );
      // OffJet_HOEnergy           .push_back( pfjet.hoEnergy()       );
      //
      // OffJet_chargedHadronEnergyFraction.push_back( pfjet.chargedHadronEnergy()/pfjet.energy());
      // OffJet_neutralHadronEnergyFraction.push_back( pfjet.neutralHadronEnergy()/pfjet.energy());
      // OffJet_chargedEmEnergyFraction.push_back( pfjet.chargedEmEnergy()/pfjet.energy());
      // OffJet_neutralEmEnergyFraction.push_back( pfjet.neutralEmEnergy()/pfjet.energy());
      // OffJet_muonEnergyFraction.push_back( pfjet.muonEnergy()/pfjet.energy());
      // OffJet_electronEnergyFraction.push_back( pfjet.electronEnergy()/pfjet.energy());
      // OffJet_photonEnergyFraction.push_back( pfjet.photonEnergy()/pfjet.energy());
      // OffJet_chargedMultiplicity.push_back( pfjet.chargedMultiplicity());
      // OffJet_neutralMultiplicity.push_back( pfjet.neutralMultiplicity());
      OffJet_nConstituents.push_back( pfjet.chargedMultiplicity() + pfjet.neutralMultiplicity());
      // 
      // OffJet_chargedHadronMultiplicity.push_back( pfjet.chargedHadronMultiplicity());
      // OffJet_neutralHadronMultiplicity.push_back( pfjet.neutralHadronMultiplicity());
      // OffJet_photonMultiplicity       .push_back( pfjet.photonMultiplicity()       );
      // OffJet_electronMultiplicity     .push_back( pfjet.electronMultiplicity()     );
      // OffJet_muonMultiplicity         .push_back( pfjet.muonMultiplicity()         );
      // OffJet_HFHadronMultiplicity     .push_back( pfjet.HFHadronMultiplicity()     );
      // OffJet_HFEMMultiplicity         .push_back( pfjet.HFEMMultiplicity()         );
      // OffJet_passId.push_back( passJetId );

      n_jetoff++;

    } 
  }

  // Fill OffPuppiJet (Puppi jets)
  if(runOffline){
    for (auto it = puppi_ak4_pfjetsoffH->begin(); it != puppi_ak4_pfjetsoffH->end(); ++it) {
      auto pfjet = *it;

      // Calculate HT at start of loop
      passJetId = jetIDoff(pfjet);
      if ( passJetId == true && pfjet.pt() >= 30) {
        htoffpuppi += pfjet.pt();
        n_jetIdoffpuppi++;
      }

      // Check EM energy fraction < 0.9
      double emEnergyFraction = (pfjet.chargedEmEnergy() + pfjet.neutralEmEnergy()) / pfjet.energy();
      if (emEnergyFraction >= 0.9) continue;

      // Remove muon overlap
      reco::Candidate::LorentzVector correctedJetP4 = pfjet.p4();
      const std::vector<reco::CandidatePtr> & cands = pfjet.getJetConstituents();
      for (const auto & cand : cands) {
        const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand.get());
        const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand.get());
        if ( mu != 0 && abs(mu->pdgId()) == 13 ) {
          reco::Candidate::LorentzVector muonP4 = cand->p4();
          correctedJetP4 -= muonP4;
        }
      }

      // Apply pt cut on corrected jet pt
      if (correctedJetP4.pt() < jetAK4PtMin) continue;

      // Save jet variables with corrected momentum
      OffPuppiJet_pt .push_back( correctedJetP4.pt() );
      OffPuppiJet_eta.push_back( correctedJetP4.eta());
      OffPuppiJet_phi.push_back( correctedJetP4.phi());
      OffPuppiJet_m  .push_back( correctedJetP4.mass() );

      OffPuppiJet_area.push_back( pfjet.jetArea());

      // Not needed for Type-1 MET
      // OffPuppiJet_chargedHadronEnergy.push_back( pfjet.chargedHadronEnergy());
      // OffPuppiJet_neutralHadronEnergy.push_back( pfjet.neutralHadronEnergy());
      // OffPuppiJet_photonEnergy       .push_back( pfjet.photonEnergy()       );
      // OffPuppiJet_electronEnergy     .push_back( pfjet.electronEnergy()     );
      // OffPuppiJet_muonEnergy         .push_back( pfjet.muonEnergy()     );
      // OffPuppiJet_HFHadronEnergy     .push_back( pfjet.HFHadronEnergy() );
      // OffPuppiJet_HFEMEnergy         .push_back( pfjet.HFEMEnergy()     );
      // OffPuppiJet_HOEnergy           .push_back( pfjet.hoEnergy()       );
      //
      // OffPuppiJet_chargedHadronEnergyFraction.push_back( pfjet.chargedHadronEnergy()/pfjet.energy());
      // OffPuppiJet_neutralHadronEnergyFraction.push_back( pfjet.neutralHadronEnergy()/pfjet.energy());
      // OffPuppiJet_chargedEmEnergyFraction.push_back( pfjet.chargedEmEnergy()/pfjet.energy());
      // OffPuppiJet_neutralEmEnergyFraction.push_back( pfjet.neutralEmEnergy()/pfjet.energy());
      // OffPuppiJet_muonEnergyFraction.push_back( pfjet.muonEnergy()/pfjet.energy());
      // OffPuppiJet_electronEnergyFraction.push_back( pfjet.electronEnergy()/pfjet.energy());
      // OffPuppiJet_photonEnergyFraction.push_back( pfjet.photonEnergy()/pfjet.energy());
      // OffPuppiJet_chargedMultiplicity.push_back( pfjet.chargedMultiplicity());
      // OffPuppiJet_neutralMultiplicity.push_back( pfjet.neutralMultiplicity());
      OffPuppiJet_nConstituents.push_back( pfjet.chargedMultiplicity() + pfjet.neutralMultiplicity());
      // 
      // OffPuppiJet_chargedHadronMultiplicity.push_back( pfjet.chargedHadronMultiplicity());
      // OffPuppiJet_neutralHadronMultiplicity.push_back( pfjet.neutralHadronMultiplicity());
      // OffPuppiJet_photonMultiplicity       .push_back( pfjet.photonMultiplicity()       );
      // OffPuppiJet_electronMultiplicity     .push_back( pfjet.electronMultiplicity()     );
      // OffPuppiJet_muonMultiplicity         .push_back( pfjet.muonMultiplicity()         );
      // OffPuppiJet_HFHadronMultiplicity     .push_back( pfjet.HFHadronMultiplicity()     );
      // OffPuppiJet_HFEMMultiplicity         .push_back( pfjet.HFEMMultiplicity()         );
      // OffPuppiJet_passId.push_back( passJetId );

      n_jetoffpuppi++;

    } 
  }

  // * 
  // FatJets 
  // *

  //here introduce AK definition
  JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);
  //double jet_pt_min = 10.0;
  //double fatjet_pt_min = 100.0;
  //double sd_z_cut = 0.10;
  //double sd_beta = 0;
  //SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
  //Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));

  //double beta = 1.0;
  //Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  //Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  //Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  //Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  //Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

  //EnergyCorrelatorN2 N2=EnergyCorrelatorN2(1.0);
  //EnergyCorrelatorN3 N3=EnergyCorrelatorN3(1.0);
  
  
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  //here add for scouting AK8 jets
  FatJet_area.clear();
  FatJet_eta .clear();
  FatJet_phi .clear();
  FatJet_pt  .clear();
  FatJet_mass.clear();
  //FatJet_n2b1.clear();
  //FatJet_n3b1.clear();
  //FatJet_tau1.clear();
  //FatJet_tau2.clear();
  //FatJet_tau3.clear();
  //FatJet_tau4.clear();
  //FatJet_tau21.clear();
  //FatJet_tau32.clear();
  //FatJet_msoftdrop.clear();
  //FatJet_mtrim.clear();
  FatJet_nConstituents.clear();

  // for jet id
  FatJet_chargedHadronEnergyFraction.clear();
  FatJet_neutralHadronEnergyFraction.clear();
  FatJet_electronEnergyFraction.clear();
  FatJet_photonEnergyFraction.clear();
  FatJet_muonEnergyFraction.clear();
  FatJet_chargedMultiplicity.clear();
  FatJet_neutralMultiplicity.clear();
  FatJet_passId.clear();
  FatJet_genJetAK8Idx.clear();

  n_fatjet = 0;


  if(runScouting){

    ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
    vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(jetAK8ScoutPtMin)); //pt min cut in config file
    for (auto &j: ak8_jets) {
      if (abs(j.eta()) > 2.4) continue; //eta cut
      FatJet_eta.push_back(j.eta());
      FatJet_phi.push_back(j.phi_std());
      FatJet_pt .push_back(j.pt());
      FatJet_mass.push_back(j.m());
      FatJet_area.push_back(j.area());
      FatJet_nConstituents.push_back(j.constituents().size());

      // needs to be done manually, not included in scouting 
      // charged and neutral EM fractions are approximated with just electrons/photons (not available in scouting)
      // loop over constituents of the jet and compute the fractions of photons, electrons, muons, charged hadrons, neutral hadrons 
      double sum_photon = 0;
      double sum_electron = 0;
      double sum_charged_hadron = 0;
      double sum_neutral_hadron = 0;
      double sum_muon = 0;
      int sum_charged_hadron_multiplicity = 0;
      int sum_neutral_hadron_multiplicity = 0;
      int sum_photon_multiplicity = 0;
      int sum_electron_multiplicity = 0;
      int sum_muon_multiplicity = 0;


      for (auto &k: j.constituents()) {
        //check if the particle is a photon using the pdgid in the user info and comparing it with photons_ids
        auto found_photon = std::find(photons_ids.begin(), photons_ids.end(), abs(k.user_info<PdgIdInfo>().pdg_id()));
        if (found_photon!= photons_ids.end()) {
            sum_photon += k.E();
            sum_photon_multiplicity++;
        }

        //check if the particle is an electron using the pdgid in the user info and comparing it with electrons_ids
        auto found_electron = std::find(electrons_ids.begin(), electrons_ids.end(), abs(k.user_info<PdgIdInfo>().pdg_id()));
        if (found_electron!= electrons_ids.end()) {
            sum_electron += k.E();
            sum_electron_multiplicity++;
        }
        
        
        //check if the particle is a muon using the pdgid in the user info and comparing it with muons_ids
        auto found_muon = std::find(muons_ids.begin(), muons_ids.end(), abs(k.user_info<PdgIdInfo>().pdg_id()));
        if (found_muon!= muons_ids.end()) {
            sum_muon += k.E();
            sum_muon_multiplicity++;
        }


        //check if the particle is a charged hadron using the pdgid in the user info and comparing it with charged_hadrons_ids
        auto found_charged_hadron = std::find(chargedHadrons_ids.begin(), chargedHadrons_ids.end(), abs(k.user_info<PdgIdInfo>().pdg_id()));
        if (found_charged_hadron!= chargedHadrons_ids.end()) {
            sum_charged_hadron += k.E();
            sum_charged_hadron_multiplicity++;
        }


        //check if the particle is a neutral hadron using the pdgid in the user info and comparing it with neutral_hadrons_ids
        auto found_neutral_hadron = std::find(neutralHadrons_ids.begin(), neutralHadrons_ids.end(), abs(k.user_info<PdgIdInfo>().pdg_id()));
        if (found_neutral_hadron!= neutralHadrons_ids.end()) {
            sum_neutral_hadron += k.E();
            sum_neutral_hadron_multiplicity++;
        }
      
      }

      //define variables for jetID
      //define the fraction of electrons
      double jet_charged_em_fraction = (sum_electron) / j.E();
      //define the fraction of photons
      double jet_neutral_em_fraction = (sum_photon) / j.E();
      //define the fraction of muons
      double jet_muon_fraction = sum_muon / j.E();
      //define the fraction of charged hadrons
      double jet_charged_hadron_fraction = sum_charged_hadron / j.E();
      //define the fraction of neutral hadrons
      double jet_neutral_hadron_fraction = sum_neutral_hadron / j.E();
      //define multiplicity of charged hadrons, neutral hadrons, photons, electrons, muons
      double NumConst = sum_charged_hadron_multiplicity + sum_neutral_hadron_multiplicity + sum_photon_multiplicity + sum_electron_multiplicity + sum_muon_multiplicity;
      double CHM =  sum_muon_multiplicity + sum_electron_multiplicity + sum_charged_hadron_multiplicity;

      bool passID = (abs(j.eta())<=2.6 && jet_charged_em_fraction < 0.8 && jet_neutral_em_fraction < 0.9 && jet_muon_fraction < 0.8 && jet_neutral_hadron_fraction < 0.9 && jet_charged_hadron_fraction > 0 && NumConst > 1 && CHM > 0);

      FatJet_chargedHadronEnergyFraction.push_back(jet_charged_hadron_fraction);
      FatJet_neutralHadronEnergyFraction.push_back(jet_neutral_hadron_fraction);
      FatJet_electronEnergyFraction.push_back(jet_charged_em_fraction);
      FatJet_photonEnergyFraction.push_back(jet_neutral_em_fraction);
      FatJet_muonEnergyFraction.push_back(jet_muon_fraction);
      FatJet_chargedMultiplicity.push_back(sum_muon_multiplicity + sum_electron_multiplicity + sum_charged_hadron_multiplicity);
      FatJet_neutralMultiplicity.push_back(sum_photon_multiplicity + sum_neutral_hadron_multiplicity);

      FatJet_passId.push_back(passID) ;
      n_fatjet++;
    } 

    unsigned int n_pfcand_tot = 0;
    for (auto & pfcands_iter : PFcands ) {
      if (pfcands_iter.pt() < 1.) continue;
      if (abs(pfcands_iter.eta()) >= 2.4 ) continue;    
      int tmpidx = -1;
      int ak8count = 0;
      for (auto &j: ak8_jets) {
        for (auto &k: j.constituents()){
          if ((UInt_t)k.user_index() == n_pfcand_tot){
            tmpidx = ak8count;
            ak8count++;
            break;
          }
        }
        if (tmpidx>-1){
          FatJetPFCands_jetIdx.push_back(tmpidx);
          FatJetPFCands_pFCandsIdx.push_back(n_pfcand_tot);
          break;
        }else{
          ak8count++;
        }
      }
      n_pfcand_tot++;
    }
  }

  
  //here add for offline AK8 jets
  OfflineFatJet_area.clear();
  OfflineFatJet_eta .clear();
  OfflineFatJet_phi .clear();
  OfflineFatJet_pt  .clear();
  OfflineFatJet_mass.clear();
  // OfflineFatJet_n2b1.clear();
  // OfflineFatJet_n3b1.clear();
  // OfflineFatJet_tau1.clear();
  // OfflineFatJet_tau2.clear();
  // OfflineFatJet_tau3.clear();
  // OfflineFatJet_tau4.clear();
  // OfflineFatJet_tau21.clear();
  // OfflineFatJet_tau32.clear();
  // OfflineFatJet_msoftdrop.clear();
  // OfflineFatJet_mtrim.clear();
  OfflineFatJet_nConstituents.clear();
  
  if (runOffline){

    ClusterSequenceArea ak8_cs_offline(off_fj_part, ak8_def, area_def);
    vector<PseudoJet> ak8_jets_offline = sorted_by_pt(ak8_cs_offline.inclusive_jets(jetAK8PtMin));

    n_fatjet_off = 0;
    
    for(auto &j: ak8_jets_offline) {
      if (abs(j.pseudorapidity()) > 2.4) continue; //eta cut
      OfflineFatJet_area.push_back(j.area());
      OfflineFatJet_eta.push_back(j.pseudorapidity());
      OfflineFatJet_mass.push_back(j.m());

      OfflineFatJet_nConstituents.push_back(j.constituents().size());
      
      // PseudoJet sd_ak8_off = sd_groomer(j);
      // OfflineFatJet_msoftdrop.push_back(sd_ak8_off.m());
      
      // PseudoJet trimmed_ak8_off = trimmer(j);
      // OfflineFatJet_mtrim.push_back(trimmed_ak8_off.m());
      
      // OfflineFatJet_n2b1.push_back(N2(sd_ak8_off));
      // OfflineFatJet_n3b1.push_back(N3(sd_ak8_off));
      OfflineFatJet_phi.push_back(j.phi_std());
      OfflineFatJet_pt.push_back(j.pt());
      // OfflineFatJet_tau1.push_back(nSub1.result(j));
      // OfflineFatJet_tau2.push_back(nSub2.result(j));
      // OfflineFatJet_tau3.push_back(nSub3.result(j));
      // OfflineFatJet_tau4.push_back(nSub4.result(j));

      n_fatjet_off++;
    }

    // unsigned int n_pfcand_off_tot = 0;
    // for(auto pfcandsoff_iter = pfcandsoffH->begin(); pfcandsoff_iter != pfcandsoffH->end(); ++pfcandsoff_iter ){

    //   if (pfcandsoff_iter->pt() < 1.) continue;
    //   if (abs(pfcandsoff_iter->eta()) >= 2.4 ) continue;    
    //   int tmpidx_off = -1;
    //   int ak8count_off = 0;
    //   for (auto &j: ak8_jets_offline) {
    //     for (auto &k: j.constituents()){
    //       if ((UInt_t)k.user_index() == n_pfcand_off_tot){
    //         tmpidx_off = ak8count_off;
    //         ak8count_off++;
    //         break;
    //       }
    //     }
    //     if (tmpidx_off>-1){
    //       OfflineFatJetPFCands_jetIdx.push_back(tmpidx_off);
    //       OfflineFatJetPFCands_pFCandsIdx.push_back(n_pfcand_off_tot);
    //       break;
    //     }else{
    //       ak8count_off++;
    //     }
    //   }
    //   n_pfcand_off_tot++;
    // }

}


OffPuppiFatJet_pt.clear();
OffPuppiFatJet_eta.clear();
OffPuppiFatJet_phi.clear();
OffPuppiFatJet_m.clear();
OffPuppiFatJet_area.clear();
OffPuppiFatJet_chargedHadronEnergy.clear();
OffPuppiFatJet_neutralHadronEnergy.clear();
OffPuppiFatJet_photonEnergy.clear();
OffPuppiFatJet_electronEnergy.clear();
OffPuppiFatJet_muonEnergy.clear();
OffPuppiFatJet_HFHadronEnergy.clear();
OffPuppiFatJet_HFEMEnergy.clear();
OffPuppiFatJet_HOEnergy.clear();
OffPuppiFatJet_chargedHadronMultiplicity.clear();
OffPuppiFatJet_neutralHadronMultiplicity.clear();
OffPuppiFatJet_photonMultiplicity.clear();
OffPuppiFatJet_electronMultiplicity.clear();
OffPuppiFatJet_muonMultiplicity.clear();
OffPuppiFatJet_HFHadronMultiplicity.clear();
OffPuppiFatJet_HFEMMultiplicity.clear();
OffPuppiFatJet_passId.clear();
n_fatjet_off_puppi = 0;
n_fatjetIdoffpuppi = 0;
passOffPuppiFatJetId = false;

// for jet id
OffPuppiFatJet_chargedHadronEnergyFraction.clear();
OffPuppiFatJet_neutralHadronEnergyFraction.clear();
OffPuppiFatJet_chargedEmEnergyFraction.clear();
OffPuppiFatJet_neutralEmEnergyFraction.clear();
OffPuppiFatJet_muonEnergyFraction.clear();
OffPuppiFatJet_electronEnergyFraction.clear();
OffPuppiFatJet_photonEnergyFraction.clear();
OffPuppiFatJet_chargedMultiplicity.clear();
OffPuppiFatJet_neutralMultiplicity.clear();
OffPuppiFatJet_nConstituents.clear();


vector<vector<reco::CandidatePtr>> OffPFcandsAK8Puppi;

if(runOffline){

  for (auto it = puppi_ak8_pfjetsoffH->begin(); it != puppi_ak8_pfjetsoffH->end(); ++it) {
    auto pfjet = *it;

    if (pfjet.pt() < jetAK8PtMin) continue; //pt cut
    if (abs(pfjet.eta()) > 2.4) continue; //eta cut

    OffPuppiFatJet_pt .push_back( pfjet.pt());
    OffPuppiFatJet_eta.push_back( pfjet.eta());
    OffPuppiFatJet_phi.push_back( pfjet.phi());
    OffPuppiFatJet_m  .push_back( pfjet.mass());

    OffPuppiFatJet_area.push_back( pfjet.jetArea());

    OffPuppiFatJet_chargedHadronEnergy.push_back( pfjet.chargedHadronEnergy());
    OffPuppiFatJet_neutralHadronEnergy.push_back( pfjet.neutralHadronEnergy());
    OffPuppiFatJet_photonEnergy       .push_back( pfjet.photonEnergy()       );
    OffPuppiFatJet_electronEnergy     .push_back( pfjet.electronEnergy()     );
    OffPuppiFatJet_muonEnergy         .push_back( pfjet.muonEnergy()     );
    OffPuppiFatJet_HFHadronEnergy     .push_back( pfjet.HFHadronEnergy() );
    OffPuppiFatJet_HFEMEnergy         .push_back( pfjet.HFEMEnergy()     );
    OffPuppiFatJet_HOEnergy           .push_back( pfjet.hoEnergy()       );

    OffPuppiFatJet_chargedHadronMultiplicity.push_back( pfjet.chargedHadronMultiplicity());
    OffPuppiFatJet_neutralHadronMultiplicity.push_back( pfjet.neutralHadronMultiplicity());
    OffPuppiFatJet_photonMultiplicity       .push_back( pfjet.photonMultiplicity()       );
    OffPuppiFatJet_electronMultiplicity     .push_back( pfjet.electronMultiplicity()     );
    OffPuppiFatJet_muonMultiplicity         .push_back( pfjet.muonMultiplicity()         );
    OffPuppiFatJet_HFHadronMultiplicity     .push_back( pfjet.HFHadronMultiplicity()     );
    OffPuppiFatJet_HFEMMultiplicity         .push_back( pfjet.HFEMMultiplicity()         );

    OffPuppiFatJet_chargedHadronEnergyFraction.push_back( pfjet.chargedHadronEnergy()/pfjet.energy());
    OffPuppiFatJet_neutralHadronEnergyFraction.push_back( pfjet.neutralHadronEnergy()/pfjet.energy());
    OffPuppiFatJet_chargedEmEnergyFraction.push_back( pfjet.chargedEmEnergy()/pfjet.energy());
    OffPuppiFatJet_neutralEmEnergyFraction.push_back( pfjet.neutralEmEnergy()/pfjet.energy());
    OffPuppiFatJet_muonEnergyFraction.push_back( pfjet.muonEnergy()/pfjet.energy());
    OffPuppiFatJet_electronEnergyFraction.push_back( pfjet.electronEnergy()/pfjet.energy());
    OffPuppiFatJet_photonEnergyFraction.push_back( pfjet.photonEnergy()/pfjet.energy());
    OffPuppiFatJet_chargedMultiplicity.push_back( pfjet.chargedMultiplicity());
    OffPuppiFatJet_neutralMultiplicity.push_back( pfjet.neutralMultiplicity());
    OffPuppiFatJet_nConstituents.push_back( pfjet.chargedMultiplicity() + pfjet.neutralMultiplicity());

    n_fatjet_off_puppi++;

    passOffPuppiFatJetId = jetIDoff(pfjet);
    OffPuppiFatJet_passId.push_back( passOffPuppiFatJetId );
    OffPFcandsAK8Puppi.push_back(pfjet.getJetConstituents());

    // apply jet ID 
    if ( passOffPuppiFatJetId == false ) continue; 
    n_fatjetIdoffpuppi++ ; 
  }

}
//offline PF Cands
n_offpuppipfcand = 0;
int n_offpuppipfcand_tot = 0;
OffPuppiPFcand_pt.clear();
OffPuppiPFcand_m.clear();
OffPuppiPFcand_eta.clear();
OffPuppiPFcand_phi.clear();
OffPuppiPFcand_pdgid.clear();
OffPuppiPFcand_q.clear();
//CZZ: added jet indices for Offline PF candidates
OfflinePuppiFatJetPFCands_jetIdx.clear();
OfflinePuppiFatJetPFCands_pFCandsIdx.clear();

if(runOffline){
  for (unsigned j=0;j<OffPFcandsAK8Puppi.size();++j){
  //Adding PFcands to AK8 jets
  //iter over jets
    for (unsigned ic=0;ic<OffPFcandsAK8Puppi[j].size();++ic){
      const reco::Candidate* pfc = dynamic_cast <const reco::Candidate*> (OffPFcandsAK8Puppi[j][ic].get());
      if (pfc->pt() < 1.) continue;
      if (abs(pfc->eta()) >= 2.4 ) continue; 
      OffPuppiPFcand_pt.push_back(pfc->pt()); 
      OffPuppiPFcand_eta.push_back(pfc->eta());
      OffPuppiPFcand_phi.push_back(pfc->phi());
      OffPuppiPFcand_m.push_back(pfc->mass());
      OffPuppiPFcand_pdgid.push_back(pfc->pdgId());
      OffPuppiPFcand_q.push_back(getCharge(pfc->pdgId()));
      // save indices 
      OfflinePuppiFatJetPFCands_pFCandsIdx.push_back(n_offpuppipfcand_tot);
      OfflinePuppiFatJetPFCands_jetIdx.push_back(j);
      n_offpuppipfcand_tot++;
    }
  } 

  n_offpuppipfcand = OffPuppiPFcand_pt.size(); 

}

//  Handle<double> rhoH;
  Handle<double> rhoH2;
  if(runScouting){
    iEvent.getByToken(rhoToken2, rhoH2);
    rho2 = *rhoH2;
  }else{// rho=0;
    rho2=0;}

  if(doSignal or isMC){
    PSweights = genEvtInfo->weights();
    Handle<double> prefirewgt;
    iEvent.getByToken(prefireToken, prefirewgt);
    prefire = *prefirewgt;
    Handle<double> prefirewgtup;
    iEvent.getByToken(prefireTokenup, prefirewgtup);
    prefireup = *prefirewgtup;
    Handle<double> prefirewgtdown;
    iEvent.getByToken(prefireTokendown, prefirewgtdown);
    prefiredown = *prefirewgtdown;
  }


  // * Gen jets Ak4 // *

  GenJet_pt.clear();
  GenJet_eta.clear();
  GenJet_phi.clear();
  GenJet_mass.clear();
  
  n_genjet = 0;
  if(runGen){
    for (auto genjet = genak4jetsH->begin(); genjet != genak4jetsH->end(); ++genjet) {
      if (genjet->pt() > jetAK4PtMin){
        if(abs(genjet->eta()) > 2.4) continue;
        GenJet_pt.push_back( genjet->pt() );
        GenJet_eta.push_back( genjet->eta());
        GenJet_phi.push_back( genjet->phi());
        GenJet_mass.push_back( genjet->mass()  );
        n_genjet++;
      }
    }
  }

  // * Match AK4 reco jets to gen jets // *
  if(runGen && runScouting){
    const double maxDeltaR_AK4 = 0.4;  // Use jet radius for matching
    std::vector<bool> genJetMatched(GenJet_pt.size(), false);  // Track which gen jets are already matched
    
    for(size_t ijet = 0; ijet < Jet_pt.size(); ++ijet) {
      int bestMatch = -1;
      double minDeltaR = maxDeltaR_AK4;
      
      for(size_t igen = 0; igen < GenJet_pt.size(); ++igen) {
        if(genJetMatched[igen]) continue;  // Skip already matched gen jets (resolveAmbiguities)
        
        // Use CMSSW's standard deltaR function
        double dR = reco::deltaR(Jet_eta[ijet], Jet_phi[ijet], GenJet_eta[igen], GenJet_phi[igen]);
        
        if(dR < minDeltaR) {
          minDeltaR = dR;
          bestMatch = igen;
        }
      }
      
      if(bestMatch >= 0) {
        genJetMatched[bestMatch] = true;  // Mark this gen jet as matched
      }
      Jet_genJetIdx.push_back(bestMatch);
    }
  } else {
    // Fill with -1 if no gen info
    for(size_t ijet = 0; ijet < Jet_pt.size(); ++ijet) {
      Jet_genJetIdx.push_back(-1);
    }
  }


  // * Gen jets AK8 // *

  GenFatJet_pt.clear();
  GenFatJet_eta.clear();
  GenFatJet_phi.clear();
  GenFatJet_mass.clear();
    
  n_genfatjet = 0;

  if(runGen){
    for (auto genjet = genak8jetsH->begin(); genjet != genak8jetsH->end(); ++genjet) {
      if (genjet->pt() > jetAK8PtMin){
        if(abs(genjet->eta()) > 2.4) continue;
        GenFatJet_pt.push_back( genjet->pt() );
        GenFatJet_eta.push_back( genjet->eta());
        GenFatJet_phi.push_back( genjet->phi());
        GenFatJet_mass.push_back( genjet->mass()  );
        n_genfatjet++;
      }
    }
  }

  // * Match AK8 reco jets to gen jets // *
  if(runGen && runScouting){
    const double maxDeltaR_AK8 = 0.8;  // Use jet radius for matching
    std::vector<bool> genFatJetMatched(GenFatJet_pt.size(), false);  // Track which gen jets are already matched
    
    for(size_t ijet = 0; ijet < FatJet_pt.size(); ++ijet) {
      int bestMatch = -1;
      double minDeltaR = maxDeltaR_AK8;
      
      for(size_t igen = 0; igen < GenFatJet_pt.size(); ++igen) {
        if(genFatJetMatched[igen]) continue;  // Skip already matched gen jets (resolveAmbiguities)
        
        // Use CMSSW's standard deltaR function
        double dR = reco::deltaR(FatJet_eta[ijet], FatJet_phi[ijet], GenFatJet_eta[igen], GenFatJet_phi[igen]);
        
        if(dR < minDeltaR) {
          minDeltaR = dR;
          bestMatch = igen;
        }
      }
      
      if(bestMatch >= 0) {
        genFatJetMatched[bestMatch] = true;  // Mark this gen jet as matched
      }
      FatJet_genJetAK8Idx.push_back(bestMatch);
    }
  } else {
    // Fill with -1 if no gen info
    for(size_t ijet = 0; ijet < FatJet_pt.size(); ++ijet) {
      FatJet_genJetAK8Idx.push_back(-1);
    }
  }


  // * 
  // L1 info
  // *
  l1Result_.clear();
  l1Prescale_.clear();

  if (doL1) {

    //I seem to recall this function being slow so perhaps cache for a given lumi 
    //(it only changes on lumi boundaries)  
    //note to the reader, what I'm doing is extremely dangerious (a const cast), never do this!           
    //however in this narrow case, it fixes a bug in l1t::L1TGlobalUtil (the method should be const)          
    //and it is safe for this specific instance                                                                                     
    l1t::L1TGlobalUtil& l1GtUtils = const_cast<l1t::L1TGlobalUtil&> (hltPSProv_.l1tGlobalUtil());

    // For debugging: from https://github.com/Sam-Harper/usercode/blob/09e2252601da473ba02de966930863df57512438/TrigTools/plugins/L1MenuExample.cc
    std::cout <<"l1 menu: name decisions prescale "<<std::endl;

    for(size_t bitNr=0;bitNr<l1GtUtils.decisionsFinal().size();bitNr++){
        const std::string& bitName = l1GtUtils.decisionsFinal()[bitNr].first; // l1GtUtils.decisionsFinal() is of type std::vector<std::pair<std::string,bool> >
        bool passInitial = l1GtUtils.decisionsInitial()[bitNr].second; //before masks and prescales, so if we have a 15 GeV electron passing L1_SingleEG10, it will show up as true but will likely not cause a L1 acccept due to the seeds high prescale
        bool passInterm = l1GtUtils.decisionsInterm()[bitNr].second; //after mask (?, unsure what this is)
        bool passFinal = l1GtUtils.decisionsFinal()[bitNr].second; //after masks & prescales, true means it gives a L1 accept to the HLT
        int prescale = l1GtUtils.prescales()[bitNr].second;
        std::cout <<"   "<<bitNr<<" "<<bitName<<" "<<passInitial<<" "<<passInterm<<" "<<passFinal<<" "<<prescale<<std::endl;
        for(size_t i = 0; i < l1Seeds_.size(); i++){
          std::string l1Name = l1Seeds_[i];
          std::string pathName = bitName;
          if(bitName.find(l1Name) != std::string::npos ){
            l1Result_  .push_back(passFinal);
            l1Prescale_.push_back(prescale);
          }
        }
    }


  }

  //Adding met info
  met_pt = -1;
  met_phi = -1;
  met_pt_reco = -1;
  met_phi_reco = -1;
  genMET_pt = -1;
  genMET_phi = -1;

  if (runScouting){
    met_pt = *metPt;
    met_phi = *metPhi;
  }

  if (runOffline){
    met_pt_reco = metReco->front().pt();
    met_phi_reco = metReco->front().phi();
  }
  
  if (runGen){
    genMET_pt = genMetH->front().pt();
    genMET_phi = genMetH->front().phi();
  }
  
  tree->Fill();	
	
}


void ScoutingNanoAOD_fromAOD::beginJob() {
  
}

void ScoutingNanoAOD_fromAOD::endJob() {
}

void ScoutingNanoAOD_fromAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  //we need to initalise the menu each run (menu can and will change on run boundaries)           
  //for L1
  bool changed=false;
  hltPSProv_.init(iRun,iSetup,hltProcess_,changed);
  const l1t::L1TGlobalUtil& l1GtUtils = hltPSProv_.l1tGlobalUtil();
  std::cout <<"l1 menu "<<l1GtUtils.gtTriggerMenuName()<<" version "<<l1GtUtils.gtTriggerMenuVersion()<<" comment "<<std::endl;
  std::cout <<"hlt name "<<hltPSProv_.hltConfigProvider().tableName()<<std::endl;

}

void ScoutingNanoAOD_fromAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD_fromAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
edm::Handle<GenLumiInfoHeader> genLumiInfoHead;
    iLumi.getByToken(genLumiInfoHeadTag_, genLumiInfoHead);
    if (!genLumiInfoHead.isValid())
      edm::LogWarning("LHETablesProducer")
          << "No GenLumiInfoHeader product found, will not fill generator model string.\n";

    //std::string label;
    if (genLumiInfoHead.isValid()) {
      label = genLumiInfoHead->configDescription();
      //printf("label: %s\n",label.c_str());
      boost::replace_all(label, "-", "_");
      boost::replace_all(label, "/", "_");
      label = std::string("GenModel_") + label;
    }
}

void ScoutingNanoAOD_fromAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

int ScoutingNanoAOD_fromAOD::getCharge(int pdgId) {
  // following workbook 
  if      (abs(pdgId) == 11 ) return 1; // electron
  else if (abs(pdgId) == 13 ) return 1; // muon
  else if (abs(pdgId) == 211) return 1; // pion
  return 0;
  // 130 = KLong - neutral hadron 
  // 22 = photon 
  // 1 = HF hadron, where HF means forward calo
  // 2 = HF em particle, where HF means forward calo
}
//bool ScoutingNanoAOD_fromAOD::jetIDoff(const reco::PFJet &pfjet){
bool ScoutingNanoAOD_fromAOD::jetIDoff(const reco::PFJet &pfjet){
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
// same definition is used in SUEP Scouting search EXO-23-001
    //moved HT cut
    TLorentzVector jet; 
    jet.SetPtEtaPhiM(pfjet.pt(), pfjet.eta(), pfjet.phi(), pfjet.mass() );
    
    float NHF  = pfjet.neutralHadronEnergy()/jet.E();
    float NEMF = pfjet.neutralEmEnergy()/jet.E();
    float CHF  = pfjet.chargedHadronEnergy()/jet.E();
    float MUF  = pfjet.muonEnergy()/jet.E();
    float CEMF = pfjet.chargedEmEnergy()/jet.E();
    float NumConst = pfjet.chargedMultiplicity()+pfjet.neutralMultiplicity();
    float CHM      = pfjet.chargedMultiplicity(); 
    bool passID = (abs(pfjet.eta())<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );

    return passID;
}
bool ScoutingNanoAOD_fromAOD::jetID(const ScoutingPFJet &pfjet){
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
// same definition is used in SUEP Scouting search EXO-23-001
// scouting version uses photon/electron energy instead of neutral/charged EM energy
    //moved HT cut
    TLorentzVector jet; 
    jet.SetPtEtaPhiM(pfjet.pt(), pfjet.eta(), pfjet.phi(), pfjet.m() );
    
    float NHF  = pfjet.neutralHadronEnergy()/jet.E();
    float NEMF = pfjet.photonEnergy()/jet.E();
    float CHF  = pfjet.chargedHadronEnergy()/jet.E();
    float MUF  = pfjet.muonEnergy()/jet.E();
    float CEMF = pfjet.electronEnergy()/jet.E();
    float NumConst = pfjet.chargedHadronMultiplicity()+pfjet.neutralHadronMultiplicity()+pfjet.photonMultiplicity() + pfjet.electronMultiplicity() + pfjet.muonMultiplicity() + pfjet.HFHadronMultiplicity() + pfjet.HFEMMultiplicity();
    float CHM      = pfjet.chargedHadronMultiplicity() +pfjet.electronMultiplicity() + pfjet.muonMultiplicity(); 
    
    bool passID = (abs(pfjet.eta())<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );

    return passID;
}

void ScoutingNanoAOD_fromAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD_fromAOD);
