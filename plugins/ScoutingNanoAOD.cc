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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

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



using namespace std;


class ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingNanoAOD(const edm::ParameterSet&);
  ~ScoutingNanoAOD();
		
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

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;
  const edm::EDGetTokenT<std::vector<ScoutingMuon> >            muonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingElectron> >  	electronsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPhoton> >  	photonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingParticle> >  	pfcandsToken;
  const edm::EDGetTokenT<std::vector<ScoutingPFJet> >  		pfjetsToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet> >  		pfjetsoffToken;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> >       pileupInfoToken;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> >       pileupInfoToken2;
  const edm::EDGetTokenT<GenEventInfoProduct>                  genEvtInfoToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >  	gensToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >  	gensToken2;
  const edm::EDGetTokenT<double>  	rhoToken2;
  const edm::EDGetTokenT<double>  	prefireToken;
  const edm::EDGetTokenT<double>  	prefireTokenup;
  const edm::EDGetTokenT<double>  	prefireTokendown;
  const edm::EDGetTokenT<GenLumiInfoHeader>  	genLumiInfoHeadTag_;

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
  bool era_16;
  bool runScouting = false;
  bool runOffline =false;
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
  UInt_t offline_trig; 
  UInt_t veto_trig;
  //Photon
  UInt_t n_pho;
  vector<Float16_t> 	       Photon_pt;
  vector<Float16_t>            Photon_eta;
  vector<Float16_t>            Photon_phi;
  vector<Float16_t>	       Photon_m;
  vector<Float16_t>	       Photon_sigmaietaieta;
  vector<Float16_t>	       Photon_HoE;
  vector<Float16_t>            Photon_ecaliso;
  vector<Float16_t>	       Photon_hcaliso;

  //Electron
  UInt_t n_ele;
  vector<Float16_t> 	       Electron_pt;
  vector<Float16_t>            Electron_eta;
  vector<Float16_t>            Electron_phi;
  vector<Float16_t>	       Electron_m;
  vector<Float16_t>            Electron_d0;
  vector<Float16_t>	       Electron_dz;
  vector<Float16_t>	       Electron_detain;
  vector<Float16_t>	       Electron_dphiin;
  vector<Float16_t>	       Electron_sigmaietaieta;
  vector<Float16_t>	       Electron_HoE;
  vector<Float16_t>	       Electron_ooEMOop;
  vector<Float16_t>	       Electron_mHits;
  vector<Float16_t>            Electron_charge;
  vector<Float16_t>            Electron_ecaliso;
  vector<Float16_t>	       Electron_hcaliso;
  vector<Float16_t>            Electron_trkiso;
  vector<Float16_t>            Electron_combinediso;
  vector<bool>            Electron_ID;

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

  UInt_t                       PU_num;
  //PFJets
  UInt_t                       n_jet;
  UInt_t                       n_jetId;
  float                        ht;
  float                        htoff;
  float                        Muon_totPt;
  float                        Electron_totPt;
  bool                         passJetId;
  vector<Float16_t> 	       Jet_pt;
  vector<Float16_t>            Jet_eta;
  vector<Float16_t>            Jet_phi;
  vector<Float16_t>	       Jet_m;
  vector<Float16_t>	       Jet_area;
  vector<Float16_t>	       Jet_chargedHadronEnergy;
  vector<Float16_t>            Jet_neutralHadronEnergy;
  vector<Float16_t>	       Jet_photonEnergy;
  vector<Float16_t>	       Jet_electronEnergy;
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
  vector<bool>                 Jet_passId;

  
 
  //PFCand
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
  vector<Float16_t>            PFcand_fjidx;
  vector<Float16_t>            PFcand_dR;
  vector<Float16_t>            PFcand_alldR;
  
  
  // SUEP decay products
  float                        scalar_pt;
  float                        scalar_eta;
  float                        scalar_phi;
  float                        scalar_m;
  vector<Float16_t>	       truth_pts;
  vector<Float16_t>	       truth_etas;
  vector<Float16_t>	       truth_phis;
  vector<Float16_t>	       truth_dR;
  vector<Float16_t>	       truth_mass;
  vector<UInt_t>	       truth_PV;
  vector<Float16_t>	       truth_PVdZ;

  // Fatjets 
  UInt_t                       n_fatjet;
  vector<Float16_t>            FatJet_area;
  vector<Float16_t>            FatJet_eta;
  vector<Float16_t>            FatJet_n2b1;
  vector<Float16_t>            FatJet_n3b1;
  vector<Float16_t>            FatJet_phi;
  vector<Float16_t>            FatJet_pt;
  vector<Float16_t>            FatJet_tau1;
  vector<Float16_t>            FatJet_tau2;
  vector<Float16_t>            FatJet_tau3;
  vector<Float16_t>            FatJet_tau4;
  vector<Float16_t>            FatJet_tau21;
  vector<Float16_t>            FatJet_tau32;
  vector<Float16_t>            FatJet_mass;
  vector<Float16_t>            FatJet_msoftdrop;
  vector<Float16_t>            FatJet_mtrim;
  vector<Float16_t>            FatJet_nconst;

  // Fatjets Cambridge Aachen 
  UInt_t                       n_fatjet_CA;
  vector<Float16_t>            FatJet_area_CA;
  vector<Float16_t>            FatJet_eta_CA;
  vector<Float16_t>            FatJet_n2b1_CA;
  vector<Float16_t>            FatJet_n3b1_CA;
  vector<Float16_t>            FatJet_phi_CA;
  vector<Float16_t>            FatJet_pt_CA;
  vector<Float16_t>            FatJet_tau1_CA;
  vector<Float16_t>            FatJet_tau2_CA;
  vector<Float16_t>            FatJet_tau3_CA;
  vector<Float16_t>            FatJet_tau4_CA;
  vector<Float16_t>            FatJet_tau21_CA;
  vector<Float16_t>            FatJet_tau32_CA;
  vector<Float16_t>            FatJet_mass_CA;
  vector<Float16_t>            FatJet_msoftdrop_CA;
  vector<Float16_t>            FatJet_mtrim_CA;
  vector<Float16_t>            FatJet_nconst_CA;

 
//  float                        rho;
  float                        rho2;
  float                        prefire;
  float                        prefireup;
  float                        prefiredown;

  // Event shape variables
  float                        event_isotropy;
  float                        event_circularity;
  float                        event_sphericity;
  float                        event_thrust; // need to save actual reco objects for thrust
  
        
  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;
  int event_;

};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  muonsToken               (consumes<std::vector<ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken           (consumes<std::vector<ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken             (consumes<std::vector<ScoutingPhoton> >           (iConfig.getParameter<edm::InputTag>("photons"))), 
  pfcandsToken             (consumes<std::vector<ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
  pfjetsToken              (consumes<std::vector<ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))), 
  pfjetsoffToken           (consumes<std::vector<reco::PFJet> >              (iConfig.getParameter<edm::InputTag>("pfjetsoff"))), 
  pileupInfoToken          (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo"))),
  pileupInfoToken2         (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo_sig"))),
  genEvtInfoToken          (consumes<GenEventInfoProduct>                    (iConfig.getParameter<edm::InputTag>("geneventinfo"))),    
  gensToken                (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("gens"))),
  gensToken2               (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("gens_sig"))),
  rhoToken2                (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho2"))),
  prefireToken             (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProb"))),
  prefireTokenup           (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
  prefireTokendown         (consumes<double>                                 (edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
  genLumiInfoHeadTag_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),
  doL1                     (iConfig.existsAs<bool>("doL1")              ?    iConfig.getParameter<bool>  ("doL1")            : false),
  doData                   (iConfig.existsAs<bool>("doData")            ?    iConfig.getParameter<bool>  ("doData")            : false),
  doSignal                 (iConfig.existsAs<bool>("doSignal")          ?    iConfig.getParameter<bool>  ("doSignal")            : false),
  isMC                     (iConfig.existsAs<bool>("isMC")              ?    iConfig.getParameter<bool>  ("isMC")            : true),
  era_16                   (iConfig.existsAs<bool>("era_16")            ?    iConfig.getParameter<bool>  ("era_16")            : false),


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
  tree = fs->make<TTree>("tree"       , "tree");

  // Event weights
    
  tree->Branch("lumSec"		                  ,&lumSec            ,"lumSec/i");
  tree->Branch("run"		                    ,&run                  ,"run/i");
  tree->Branch("event"		                    ,&event_                  ,"event/i");
  tree->Branch("PSweights"            	    ,&PSweights 	                 );
  tree->Branch("prefire"		                ,&prefire                      );
  tree->Branch("prefireup"		              ,&prefireup                    );
  tree->Branch("prefiredown"		            ,&prefiredown                  );
    
  // Triggers
  tree->Branch("hltResult"                      ,&hltResult_                    );              
  tree->Branch("hltResultName"                  ,&hltResultName_                );              
  tree->Branch("l1Result"		                    ,&l1Result_	                );		
  tree->Branch("l1Prescale"		                  ,&l1Prescale_                   );		
  //Electrons
  tree->Branch("n_ele"               	         ,&n_ele                        ,"n_ele/i");
  tree->Branch("Electron_pt"                    ,&Electron_pt                   );
  tree->Branch("Electron_eta"                   ,&Electron_eta 	                );
  tree->Branch("Electron_phi"                   ,&Electron_phi                  );
  tree->Branch("Electron_charge"                ,&Electron_charge               );
  tree->Branch("Electron_m"            	        ,&Electron_m                    );
  tree->Branch("Electron_trkiso"                 ,&Electron_trkiso 	        );
  tree->Branch("Electron_HoE"                   ,&Electron_HoE                  );
  tree->Branch("Electron_sigmaietaieta"         ,&Electron_sigmaietaieta        );
  tree->Branch("Electron_dphiin"                ,&Electron_dphiin 	        );
  tree->Branch("Electron_detain"                ,&Electron_detain 	        );
  tree->Branch("Electron_mHits"                 ,&Electron_mHits 	        );
  tree->Branch("Electron_ooEMOop"               ,&Electron_ooEMOop              );
  tree->Branch("Electron_ecaliso"               ,&Electron_ecaliso              );
  tree->Branch("Electron_hcaliso"               ,&Electron_hcaliso              );
  tree->Branch("Electron_combinediso"               ,&Electron_combinediso   );
  tree->Branch("Electron_ID"               ,&Electron_ID   );
  tree->Branch("Electron_d0"               ,&Electron_d0              );
  tree->Branch("Electron_dz"               ,&Electron_dz              );

  tree->Branch("scouting_trig"            	        ,&scouting_trig 			,"scounting_trig/i");
  tree->Branch("offline_trig"            	        ,&offline_trig 			,"offline_trig/i");
  tree->Branch("veto_trig"            	        ,&veto_trig 			,"veto_trig/i");
  tree->Branch("genModel"            	        ,&label 			);
  //Photons
  tree->Branch("n_pho"            	        ,&n_pho 			,"n_pho/i");
  tree->Branch("Photon_pt"            	        ,&Photon_pt                     );
  tree->Branch("Photon_eta"            	        ,&Photon_eta                    );
  tree->Branch("Photon_phi"            	        ,&Photon_phi                    );	
  tree->Branch("Photon_m"            	        ,&Photon_m 	                );
  tree->Branch("Photon_hcaliso"                 ,&Photon_hcaliso 		);
  tree->Branch("Photon_ecaliso"                 ,&Photon_ecaliso 		);
  tree->Branch("Photon_HoE"            	        ,&Photon_HoE                    );
  tree->Branch("Photon_sigmaietaieta"           ,&Photon_sigmaietaieta	        );


  tree->Branch("n_pfcand"            	        ,&n_pfcand 		        ,"n_pfcand/i");	
  tree->Branch("n_pfMu"            	        ,&n_pfMu 		        ,"n_pfMu/i");	
  tree->Branch("n_pfEl"            	        ,&n_pfEl 		        ,"n_pfEl/i");	
  tree->Branch("PFcand_pt"        	        ,&PFcand_pt 		        );
  tree->Branch("PFcand_eta"            	        ,&PFcand_eta 	                );
  tree->Branch("PFcand_phi"            	        ,&PFcand_phi		        );
  tree->Branch("PFcand_m"            	        ,&PFcand_m 		        );
  tree->Branch("PFcand_pdgid"                   ,&PFcand_pdgid                  );
  tree->Branch("PFcand_q"                       ,&PFcand_q                      );
  tree->Branch("PFcand_vertex"                  ,&PFcand_vertex                 );
  tree->Branch("PFcand_fjidx"                   ,&PFcand_fjidx 	                );
  tree->Branch("PFcand_dR"        	        ,&PFcand_dR 	                );
  tree->Branch("PFcand_alldR"        	        ,&PFcand_alldR 	                );

  tree->Branch("scalar_pt"                      ,&scalar_pt                     );
  tree->Branch("scalar_eta"                     ,&scalar_eta                    );
  tree->Branch("scalar_phi"                     ,&scalar_phi                    );
  tree->Branch("scalar_m"                       ,&scalar_m                      );
  tree->Branch("gen_pt"                         ,&truth_pts                     );
  tree->Branch("gen_eta"                        ,&truth_etas                    );
  tree->Branch("gen_phi"                        ,&truth_phis                    );
  tree->Branch("gen_mass"                       ,&truth_mass                    );
  tree->Branch("gen_dR"                         ,&truth_dR                      );
  tree->Branch("gen_PV"                         ,&truth_PV                      );
  tree->Branch("gen_PVdZ"                       ,&truth_PVdZ                    );

  tree->Branch("n_mu"            	        ,&n_mu 	                        ,"n_mu/i");
  tree->Branch("Muon_pt"                        ,&Muon_pt                       );
  tree->Branch("Muon_eta"                       ,&Muon_eta                      );
  tree->Branch("Muon_phi"                       ,&Muon_phi                      );
  tree->Branch("Muon_m"                         ,&Muon_m                        );
  tree->Branch("Muon_ecaliso"                   ,&Muon_ecaliso                  );
  tree->Branch("Muon_hcaliso"                   ,&Muon_hcaliso                  );
  tree->Branch("Muon_trkiso"                    ,&Muon_trkiso                   );
  tree->Branch("Muon_chi2"                      ,&Muon_chi2                     );
  tree->Branch("Muon_isGlobalMuon"              ,&Muon_isGlobalMuon             );
  tree->Branch("Muon_isTrackerMuon"             ,&Muon_isTrackerMuon            );
  tree->Branch("Muon_ndof"                      ,&Muon_ndof                     );
  tree->Branch("Muon_charge"                    ,&Muon_charge	                  );
  tree->Branch("Muon_dxy"                       ,&Muon_dxy                      );
  tree->Branch("Muon_dz"                        ,&Muon_dz                       );
  tree->Branch("Muon_nvalidmuon_hits"           ,&Muon_nvalidmuon_hits          );
  tree->Branch("Muon_validpixelhits"            ,&Muon_nvalidpixelhits          );
  tree->Branch("Muon_nmatchedstations"          ,&Muon_nmatchedstations         );
  tree->Branch("Muon_type"                      ,&Muon_type                     );
  tree->Branch("Muon_nvalidstriphits"           ,&Muon_nvalidstriphits          );
  tree->Branch("Muon_trkqoverp"                 ,&Muon_trkqoverp                );
  tree->Branch("Muon_trklambda"                 ,&Muon_trklambda                );
  tree->Branch("Muon_trkpt"                     ,&Muon_trkpt                    );
  tree->Branch("Muon_trkphi"                    ,&Muon_trkphi                   );
  tree->Branch("Muon_trketa"                    ,&Muon_trketa                   );
  tree->Branch("Muon_trkqoverperror"            ,&Muon_trkqoverperror           );
  tree->Branch("Muon_trklambdaerror"            ,&Muon_trklambdaerror           );
  tree->Branch("Muon_trkpterror"                ,&Muon_trkpterror               );
  tree->Branch("Muon_trkphierror"               ,&Muon_trkphierror              );
  tree->Branch("Muon_trketaerror"               ,&Muon_trketaerror              );
  tree->Branch("Muon_trkdszerror"               ,&Muon_trkdszerror              );
  tree->Branch("Muon_trkdsz"                    ,&Muon_trkdsz                   );

  tree->Branch("ht"                             ,&ht                            );
  tree->Branch("htoff"                             ,&htoff                            );
  tree->Branch("Muon_totPt"                     ,&Muon_totPt                    );
  tree->Branch("Electron_totPt"                 ,&Electron_totPt                );
  tree->Branch("PU_num"            	        ,&PU_num                        ,"PU_num/i");
  tree->Branch("n_jet"            	        ,&n_jet                         ,"n_jet/i");
  tree->Branch("n_jetId"            	        ,&n_jetId                       ,"n_jetId/i");
  tree->Branch("Jet_pt"            	        ,&Jet_pt                        );
  tree->Branch("Jet_eta"            	        ,&Jet_eta                       );
  tree->Branch("Jet_phi"            	        ,&Jet_phi                       );
  tree->Branch("Jet_m"            	        ,&Jet_m                         );
  tree->Branch("Jet_area"            	        ,&Jet_area                      );
  tree->Branch("Jet_chargedHadronEnergy"        ,&Jet_chargedHadronEnergy       );
  tree->Branch("Jet_neutralHadronEnergy"        ,&Jet_neutralHadronEnergy       );
  tree->Branch("Jet_photonEnergy"               ,&Jet_photonEnergy 	        );
  tree->Branch("Jet_electronEnergy"             ,&Jet_electronEnergy            );
  tree->Branch("Jet_muonEnergy"    	        ,&Jet_muonEnergy                );
  tree->Branch("Jet_HFHadronEnergy"             ,&Jet_HFHadronEnergy            );
  tree->Branch("Jet_HFEMEnergy"                 ,&Jet_HFEMEnergy                );
  tree->Branch("Jet_HOEnergy"                   ,&Jet_HOEnergy                  );
  tree->Branch("Jet_chargedHadronMultiplicity"  ,&Jet_chargedHadronMultiplicity );
  tree->Branch("Jet_neutralHadronMultiplicity"  ,&Jet_neutralHadronMultiplicity );
  tree->Branch("Jet_photonMultiplicity"         ,&Jet_photonMultiplicity        );
  tree->Branch("Jet_electronMultiplicity"       ,&Jet_electronMultiplicity      );
  tree->Branch("Jet_muonMultiplicity"           ,&Jet_muonMultiplicity          );
  tree->Branch("Jet_HFHadronMultiplicity"       ,&Jet_HFHadronMultiplicity      );
  tree->Branch("Jet_HFEMMultiplicity"           ,&Jet_HFEMMultiplicity          );
  tree->Branch("Jet_csv"            	        ,&Jet_csv                       );
  tree->Branch("Jet_mvaDiscriminator"           ,&Jet_mvaDiscriminator          );
  tree->Branch("Jet_nConstituents"              ,&Jet_nConstituents             );
  tree->Branch("Jet_passId"                     ,&Jet_passId                    );
  
  
  tree->Branch("n_fatjet"                       ,&n_fatjet                      ,"n_fatjet/i");
  tree->Branch("FatJet_area"                    ,&FatJet_area                   );
  tree->Branch("FatJet_eta"                     ,&FatJet_eta                    );
  tree->Branch("FatJet_n2b1"                    ,&FatJet_n2b1                   );
  tree->Branch("FatJet_n3b1"                    ,&FatJet_n3b1                   );
  tree->Branch("FatJet_phi"                     ,&FatJet_phi                    );
  tree->Branch("FatJet_pt"                      ,&FatJet_pt                     );
  tree->Branch("FatJet_tau1"                    ,&FatJet_tau1                   );
  tree->Branch("FatJet_tau2"                    ,&FatJet_tau2                   );
  tree->Branch("FatJet_tau3"                    ,&FatJet_tau3                   );
  tree->Branch("FatJet_tau4"                    ,&FatJet_tau4                   );
  tree->Branch("FatJet_tau21"                   ,&FatJet_tau21                  );
  tree->Branch("FatJet_tau32"                   ,&FatJet_tau32                  );
  tree->Branch("FatJet_mass"                    ,&FatJet_mass                   );
  tree->Branch("FatJet_msoftdrop"               ,&FatJet_msoftdrop              );
  tree->Branch("FatJet_mtrim"                   ,&FatJet_mtrim                  );
  tree->Branch("FatJet_nconst"                  ,&FatJet_nconst                 );
  
  tree->Branch("n_fatjet_CA"                       ,&n_fatjet_CA                      ,"n_fatjet_CA/i");
  tree->Branch("FatJet_area_CA"                    ,&FatJet_area_CA                   );
  tree->Branch("FatJet_eta_CA"                     ,&FatJet_eta_CA                    );
  tree->Branch("FatJet_n2b1_CA"                    ,&FatJet_n2b1_CA                   );
  tree->Branch("FatJet_n3b1_CA"                    ,&FatJet_n3b1_CA                   );
  tree->Branch("FatJet_phi_CA"                     ,&FatJet_phi_CA                    );
  tree->Branch("FatJet_pt_CA"                      ,&FatJet_pt_CA                     );
  tree->Branch("FatJet_tau1_CA"                    ,&FatJet_tau1_CA                   );
  tree->Branch("FatJet_tau2_CA"                    ,&FatJet_tau2_CA                   );
  tree->Branch("FatJet_tau3_CA"                    ,&FatJet_tau3_CA                   );
  tree->Branch("FatJet_tau4_CA"                    ,&FatJet_tau4_CA                   );
  tree->Branch("FatJet_tau21_CA"                   ,&FatJet_tau21_CA                  );
  tree->Branch("FatJet_tau32_CA"                   ,&FatJet_tau32_CA                  );
  tree->Branch("FatJet_mass_CA"                    ,&FatJet_mass_CA                   );
  tree->Branch("FatJet_msoftdrop_CA"               ,&FatJet_msoftdrop_CA              );
  tree->Branch("FatJet_mtrim_CA"                   ,&FatJet_mtrim_CA                  );
  tree->Branch("FatJet_nconst_CA"                  ,&FatJet_nconst_CA                 );

  tree->Branch("rho"                            ,&rho2                           );

  tree->Branch("event_isotropy"                 ,&event_isotropy                );
  tree->Branch("event_circularity"              ,&event_circularity             );
  tree->Branch("event_sphericity"               ,&event_sphericity              );
  tree->Branch("event_thrust"                   ,&event_thrust                  );
  
  
}


ScoutingNanoAOD::~ScoutingNanoAOD() {
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;


  Handle<vector<reco::PFJet> > pfjetsoffH;

  Handle<vector<ScoutingElectron> > electronsH;
  Handle<vector<ScoutingMuon> > muonsH;
  Handle<vector<ScoutingPhoton> > photonsH;
  Handle<vector<ScoutingPFJet> > pfjetsH;
  Handle<vector<ScoutingParticle> > pfcandsH;
  Handle<vector<reco::PFCandidate> > tracksH1;
  Handle<vector<pat::PackedCandidate> > tracksH2;
   
  if(auto handle = iEvent.getHandle(pfcandsToken)){
    runScouting = true;
  }
 
  if(runScouting){
    iEvent.getByToken(electronsToken, electronsH);

    iEvent.getByToken(muonsToken, muonsH);

    iEvent.getByToken(photonsToken, photonsH);

    iEvent.getByToken(pfjetsToken, pfjetsH);
      
    iEvent.getByToken(pfcandsToken, pfcandsH);

 
   
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
       
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {                                                          
      const std::string& hltbitName = names.triggerName(i);
      std::string hltpathName = hltbitName;
      bool hltpassFinal = triggerBits->accept(i);
       
        if( 
          TString(hltpathName).Contains(pattern1) and hltpassFinal)
          {
          scouting_trig=1;
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
  // Electrons here, also electrons are not contained in pf candidate collection. need to merge them explicitly
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
        & (combinediso/electrons_iter->pt() < 0.15);
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
        & (combinediso/electrons_iter->pt() < 0.15);
      }
      Electron_combinediso.push_back(combinediso);
      Electron_ID.push_back(electronID);
    }}

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
  }}

  
  
  if (!doData) {
    for(auto PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI){
      int pu_bunchcrossing = PVI->getBunchCrossing();
      if(pu_bunchcrossing ==0){
        PU_num = PVI->getTrueNumInteractions();
      }
    }
  }


if(runScouting and isMC){  //do not run for data
  Handle<vector<reco::GenParticle> > genP;
  iEvent.getByToken(gensToken2, genP);

  truth_pts.clear();
  truth_etas.clear();
  truth_phis.clear();
  truth_mass.clear();
  truth_dR.clear();
  truth_PV.clear();
  truth_PVdZ.clear();

  for (auto genp_iter = genP->begin(); genp_iter != genP->end(); ++genp_iter ) {
   
    if (genp_iter->status()!=1){continue;}
    if (genp_iter->charge()==0){continue;}
      truth_pts.push_back(genp_iter->pt());
      truth_etas.push_back(genp_iter->eta());
      truth_phis.push_back(genp_iter->phi());
      truth_mass.push_back(genp_iter->mass());

      
  }
}
  // * 
  // Particle Flow candidates 
  // *

  
  if(runScouting){
     for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter) {
       ScoutingParticle tmp(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->pt())),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->eta())),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->phi())),pfcands_iter->m(),pfcands_iter->pdgId(),pfcands_iter->vertex());  
    
      PFcands.push_back(tmp);
    }
    }

    //sort PFcands according to pT
    struct {
      bool operator()(ScoutingParticle a, ScoutingParticle b) const { return a.pt() > b.pt(); }
    } custompT;

    std::sort(PFcands.begin(), PFcands.end(), custompT);


    PFcand_pt.clear();
    PFcand_eta.clear();
    PFcand_phi.clear();
    PFcand_m.clear();
    PFcand_pdgid.clear();
    PFcand_q.clear();
    PFcand_vertex.clear();
    PFcand_fjidx.clear();
    PFcand_dR.clear();
    PFcand_alldR.clear();

    
    vector<unsigned int> daughters_used; // for 1-to-1 reco-truth matching
    daughters_used.clear();
//////////////////////////////////////
//new code 
  vector<vector<float> > dr_vector;
  vector<PseudoJet> fj_part;
  vector<math::XYZVector> event_tracks; // all event tracks
  math::XYZVector trk = math::XYZVector(0,0,0); 
  n_pfcand = 0;
  Muon_totPt =0; 
  Electron_totPt =0; 
  n_pfMu =0;
  n_pfEl =0;
  for(auto & pfcands_iter : PFcands ){ //fills PFcand track info
    vector<float> dr_vector_row; //sets all dR values between pFcands and gen tracks
    if (pfcands_iter.pt() < 0.5) continue;
    if (abs(pfcands_iter.eta()) >= 2.4 ) continue;

    PFcand_pt.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.pt())));
    PFcand_eta.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.eta())));
    PFcand_phi.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.phi())));
    if(abs(pfcands_iter.pdgId()) == 13){
     Muon_totPt += pfcands_iter.pt(); 
     n_pfMu ++;
    }
    if(abs(pfcands_iter.pdgId()) == 11){
     Electron_totPt += pfcands_iter.pt(); 
     n_pfEl ++;
    }
    if(doSignal){
      Handle<vector<reco::GenParticle> > genP;
      iEvent.getByToken(gensToken2, genP);
      for (auto genp_iter = genP->begin(); genp_iter != genP->end(); ++genp_iter ) {
        if (genp_iter->status()!=1){continue;}
        if (genp_iter->charge()==0){continue;}
        
        auto dR = deltaR2(genp_iter->eta(),genp_iter->phi(),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.eta())),MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter.phi())));

        if(pfcands_iter.vertex() ==0){ ///////////////////EXTRA SELECTION
        dr_vector_row.push_back(dR);//fills are dR values for this pFcand and all gen
        if(dR < 0.3 && abs(genp_iter->eta()) < 2.4){
        PFcand_alldR.push_back(dR);
        }
        }
        
      }
      dr_vector.push_back(dr_vector_row);// makes matrix of all pFcand-gen dR values
    }
    PFcand_m.push_back(pfcands_iter.m());
    PFcand_pdgid.push_back(pfcands_iter.pdgId());
    PFcand_q.push_back(getCharge(pfcands_iter.pdgId()));
    PFcand_vertex.push_back(pfcands_iter.vertex());
   
    // Cluster charged PF candidates into fat jets
    if (pfcands_iter.vertex() != 0) continue;
    if (abs(pfcands_iter.eta()) >= 2.4 ) continue;
    if (pfcands_iter.pt() < 0.5) continue; 
    if (getCharge(pfcands_iter.pdgId()) == 0 ) continue;

    // For clustering fat jets
    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);
    temp_jet.reset_PtYPhiM(pfcands_iter.pt(), pfcands_iter.eta(), pfcands_iter.phi(), pfcands_iter.m());
    temp_jet.set_user_index(n_pfcand);
    if (pfcands_iter.vertex() == 0 && getCharge(pfcands_iter.pdgId()) != 0 ){
      fj_part.push_back(temp_jet);
    
      // Event shape variables on whole event
      trk = math::XYZVector(0,0,0);
      trk.SetXYZ(temp_jet.px(), temp_jet.py(), temp_jet.pz() );
      event_tracks.push_back(trk);
    }

    n_pfcand++;
  }

if(doSignal){  //do not run for other samples to save time
  // 1 to 1 gen matching
  std::vector<int> used_pf;
  std::vector<int> used_gen;
  std::vector<float> dr_matched;
  std::vector<int> gen_pv;
  float min;
  do{
    min = std::numeric_limits<float>::max();
    int row =0;
    int minrow=-1;
    int mincol=-1;
    for( auto & v : dr_vector){ //loops over dR values for each PFcand
      if(std::find(used_pf.begin(), used_pf.end(), row) != used_pf.end()){row++; continue;}// skip if this row is already a matched pf Candidates and increase the row counter
      int col=0;
      for( auto & e : v){ // loops over dR values for each gen associated with this PFCand
        if(std::find(used_gen.begin(), used_gen.end(), col) != used_gen.end()){col++; continue;}// skip if this col is already a matched gen and increase the col counter
        if(e <= min){ //finds min dR value within matrix. sets the used col and row position for the minimum as well as the new min dR value.
           mincol = col;
           minrow = row;
           min = e;
        }
        col++; // gen counter position increment
      }
      row++; // pfcand counter positsion increment
    }// all dR values have been looped over and the min dR has been found with position in gen and pfcand
    if(minrow != -1 && mincol != -1){ // if there is a dR match
        used_pf.push_back(minrow);// index of Pf cand with match
        used_gen.push_back(mincol);// index of gen cand with match
        dr_matched.push_back(min);// dR between pF cand and gen
        gen_pv.push_back(PFcand_vertex.at(minrow));
    }
  }while(min < 0.3); //cut off value for min dR

for(int e = 0; e < static_cast<int>(PFcand_pt.size()); e++){//loop over pf cands again to set dR values in proper positions
    auto it = find(used_pf.begin(), used_pf.end(), e); // see if PF cand has a match
    if(it != used_pf.end()){
      float dR = dr_matched.at(it-used_pf.begin());// get dR associated with this PF cand
      PFcand_dR.push_back(dR); //push back dR at that match positon into proper position.
    }
    else{ 
      PFcand_dR.push_back(0.3);//no match found set as fake value
    }
  }
for(int e = 0; e < static_cast<int>(truth_pts.size()); e++){//loop over pf cands again to set dR values in proper positions
    auto it = find(used_gen.begin(), used_gen.end(), e); // see if PF cand has a match
    if(it != used_gen.end()){
      float dR = dr_matched.at(it-used_gen.begin());// get dR associated with this PF cand      
      truth_dR.push_back(dR); //push back dR at that match positon into proper position.
      int genPV = gen_pv.at(it-used_gen.begin());// get dR associated with this PF cand
      if(genPV ==0){
      truth_PV.push_back(0);
      }
      else if(genPV ==-1){
      truth_PV.push_back(2);
      truth_PVdZ.push_back(-1);
      }
      else{
      truth_PV.push_back(1);     
      }     
    }
    else{ 
      truth_dR.push_back(0.3);//no match found set as fake value
      truth_PV.push_back(3);
      truth_PVdZ.push_back(-1);     
    }
  }

}




///////////////////////////////////////

  // 
  // Muons   
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
  n_mu=0;  
  
  if(runScouting){
   for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
 	Muon_pt.push_back(muons_iter->pt()); 
   	Muon_eta.push_back(muons_iter->eta());
   	Muon_phi.push_back(muons_iter->phi());
   	Muon_m.push_back(muons_iter->m());
   	Muon_ecaliso.push_back(muons_iter->ecalIso());
   	Muon_hcaliso.push_back(muons_iter->hcalIso());
   	Muon_trkiso.push_back(muons_iter->chi2());
   	Muon_chi2.push_back(muons_iter->ndof());
   	Muon_ndof.push_back(muons_iter->charge());
   	Muon_charge.push_back(muons_iter->dxy());
   	Muon_dxy.push_back(muons_iter->dz());
   	Muon_dz.push_back(muons_iter->nValidMuonHits());
   	Muon_nvalidmuon_hits.push_back(muons_iter->nValidPixelHits());
   	Muon_nvalidpixelhits.push_back(muons_iter->nMatchedStations());
   	Muon_nmatchedstations.push_back(muons_iter->nTrackerLayersWithMeasurement());
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
  // * 
  // Jets 
  // * 
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_m.clear();
  Jet_area.clear();
  Jet_chargedHadronEnergy.clear();
  Jet_neutralHadronEnergy.clear();
  Jet_photonEnergy.clear();
  Jet_electronEnergy.clear();
  Jet_muonEnergy.clear();
  Jet_HFHadronEnergy.clear();
  Jet_HFEMEnergy.clear();
  Jet_HOEnergy.clear();
  Jet_chargedHadronMultiplicity.clear();
  Jet_neutralHadronMultiplicity.clear();
  Jet_photonMultiplicity.clear();
  Jet_electronMultiplicity.clear();
  Jet_muonMultiplicity.clear();
  Jet_HFHadronMultiplicity.clear();
  Jet_HFEMMultiplicity.clear();
  Jet_csv.clear();
  Jet_mvaDiscriminator.clear();
  Jet_nConstituents.clear();
  Jet_passId.clear();
  
  n_jet = 0;
  n_jetId = 0;
  ht = 0;
  passJetId = false;

  if(runScouting){
   for (auto pfjet = pfjetsH->begin(); pfjet != pfjetsH->end(); ++pfjet) {

    Jet_pt .push_back( pfjet->pt() );
    Jet_eta.push_back( pfjet->eta());
    Jet_phi.push_back( pfjet->phi());
    Jet_m  .push_back( pfjet->m()  );

    Jet_area.push_back( pfjet->jetArea());

    Jet_chargedHadronEnergy.push_back( pfjet->chargedHadronEnergy());
    Jet_neutralHadronEnergy.push_back( pfjet->neutralHadronEnergy());
    Jet_photonEnergy       .push_back( pfjet->photonEnergy()       );
    Jet_electronEnergy     .push_back( pfjet->electronEnergy()     );
    Jet_muonEnergy         .push_back( pfjet->muonEnergy()     );
    Jet_HFHadronEnergy     .push_back( pfjet->HFHadronEnergy() );
    Jet_HFEMEnergy         .push_back( pfjet->HFEMEnergy()     );
    Jet_HOEnergy           .push_back( pfjet->HOEnergy()       );
    
    Jet_chargedHadronMultiplicity.push_back( pfjet->chargedHadronMultiplicity());
    Jet_neutralHadronMultiplicity.push_back( pfjet->neutralHadronMultiplicity());
    Jet_photonMultiplicity       .push_back( pfjet->photonMultiplicity()       );
    Jet_electronMultiplicity     .push_back( pfjet->electronMultiplicity()     );
    Jet_muonMultiplicity         .push_back( pfjet->muonMultiplicity()         );
    Jet_HFHadronMultiplicity     .push_back( pfjet->HFHadronMultiplicity()     );
    Jet_HFEMMultiplicity         .push_back( pfjet->HFEMMultiplicity()         );

    Jet_csv             .push_back( pfjet->csv() );
    Jet_mvaDiscriminator.push_back( pfjet->mvaDiscriminator()    );
    Jet_nConstituents   .push_back( pfjet->constituents().size() );
    
    n_jet++;

    passJetId = jetID(*pfjet);
    Jet_passId.push_back( passJetId );

    // apply jet ID 
    if ( passJetId == false ) continue; 
    if (pfjet->pt() < 30){continue;}//raise pt threshold for HT calculation 
    ht += pfjet->pt() ; 
    n_jetId++ ; 

  }
  }
 
  // loop through constituents & save

  // * 
  // FatJets 
  // *
  FatJet_area.clear();
  FatJet_eta .clear();
  FatJet_phi .clear();
  FatJet_pt  .clear();
  FatJet_mass.clear();
  FatJet_n2b1.clear();
  FatJet_n3b1.clear();
  FatJet_tau1.clear();
  FatJet_tau2.clear();
  FatJet_tau3.clear();
  FatJet_tau4.clear();
  FatJet_tau21.clear();
  FatJet_tau32.clear();
  FatJet_msoftdrop.clear();
  FatJet_mtrim.clear();
  FatJet_nconst.clear();
  
  // * 
  // FatJets_CA 
  // *
  FatJet_area_CA.clear();
  FatJet_eta_CA .clear();
  FatJet_phi_CA .clear();
  FatJet_pt_CA  .clear();
  FatJet_mass_CA.clear();
  FatJet_n2b1_CA.clear();
  FatJet_n3b1_CA.clear();
  FatJet_tau1_CA.clear();
  FatJet_tau2_CA.clear();
  FatJet_tau3_CA.clear();
  FatJet_tau4_CA.clear();
  FatJet_tau21_CA.clear();
  FatJet_tau32_CA.clear();
  FatJet_msoftdrop_CA.clear();
  FatJet_mtrim_CA.clear();
  FatJet_nconst_CA.clear();

  JetDefinition ak08_def = JetDefinition(antikt_algorithm, 0.8);
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
  Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));

  JetDefinition CA08_def = JetDefinition(cambridge_algorithm, 0.8);  

 
  double beta = 1.0;                                                                                  
  Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

  Nsubjettiness nSub1_CA = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta)); 
  Nsubjettiness nSub2_CA = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub3_CA = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub4_CA = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub5_CA = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

  EnergyCorrelatorN2 N2=EnergyCorrelatorN2(1.0);
  EnergyCorrelatorN3 N3=EnergyCorrelatorN3(1.0);

  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  ClusterSequenceArea ak08_cs(fj_part, ak08_def, area_def);
  vector<PseudoJet> ak08_jets = sorted_by_pt(ak08_cs.inclusive_jets(30.0)); //pt min

  ClusterSequenceArea CA08_cs(fj_part, CA08_def, area_def);//CA 
  vector<PseudoJet> CA08_jets = sorted_by_pt(CA08_cs.inclusive_jets(30.0)); //pt min

  
  n_fatjet = 0;
  for(auto &j: ak08_jets) {
    FatJet_area.push_back(j.area());
    FatJet_eta .push_back(j.pseudorapidity());
    FatJet_phi .push_back(j.phi_std());
    FatJet_pt  .push_back(j.pt());
    FatJet_mass.push_back(j.m());

    FatJet_nconst.push_back(j.constituents().size());

    PseudoJet sd_ak8 = sd_groomer(j);
    FatJet_msoftdrop.push_back(sd_ak8.m());
    
    PseudoJet trimmed_ak8 = trimmer(j);
    FatJet_mtrim.push_back(trimmed_ak8.m());
    
    // Energy correlation
    FatJet_n2b1.push_back(N2(sd_ak8));
    FatJet_n3b1.push_back(N3(sd_ak8));
    
    // Nsubjettiness, tau 
    FatJet_tau1.push_back(nSub1.result(j));
    FatJet_tau2.push_back(nSub2.result(j));
    FatJet_tau3.push_back(nSub3.result(j));
    FatJet_tau4.push_back(nSub4.result(j));
    FatJet_tau21.push_back(nSub2.result(j)/nSub1.result(j));
    FatJet_tau32.push_back(nSub3.result(j)/nSub2.result(j));

    
    n_fatjet++;
  }


    //Cambridge Aachen
  n_fatjet_CA = 0;
  for(auto &j: CA08_jets) {
    FatJet_area_CA.push_back(j.area());
    FatJet_eta_CA .push_back(j.pseudorapidity());
    FatJet_phi_CA .push_back(j.phi_std());
    FatJet_pt_CA  .push_back(j.pt());
    FatJet_mass_CA.push_back(j.m());

    FatJet_nconst_CA.push_back(j.constituents().size());

    PseudoJet sd_CA8 = sd_groomer(j);
    FatJet_msoftdrop_CA.push_back(sd_CA8.m());
    
    PseudoJet trimmed_CA8 = trimmer(j);
    FatJet_mtrim_CA.push_back(trimmed_CA8.m());
    
    // Energy correlation
    FatJet_n2b1_CA.push_back(N2(sd_CA8));
    FatJet_n3b1_CA.push_back(N3(sd_CA8));
    
    // Nsubjettiness, tau 
    FatJet_tau1_CA.push_back(nSub1_CA.result(j));  
    FatJet_tau2_CA.push_back(nSub2_CA.result(j));
    FatJet_tau3_CA.push_back(nSub3_CA.result(j));
    FatJet_tau4_CA.push_back(nSub4_CA.result(j));
    FatJet_tau21_CA.push_back(nSub2_CA.result(j)/nSub1_CA.result(j));
    FatJet_tau32_CA.push_back(nSub3_CA.result(j)/nSub2_CA.result(j));
    
    
    n_fatjet_CA++;
  }



  unsigned int n_pfcand_tot = 0;
  for (auto & pfcands_iter : PFcands ) {
    if (pfcands_iter.pt() < 1.) continue;
    if (abs(pfcands_iter.eta()) >= 2.4 ) continue;    
    int tmpidx = -1;
    int ak08count = 0;
    for (auto &j: ak08_jets) {
      for (auto &k: j.constituents()){
        if ((UInt_t)k.user_index() == n_pfcand_tot){
          tmpidx = ak08count;
          ak08count++;
          break;
        }
      }
      if (tmpidx>-1)
        break;
      else
        ak08count++;
    }
    PFcand_fjidx.push_back(tmpidx);
    n_pfcand_tot++;
  }

//  Handle<double> rhoH;
  Handle<double> rhoH2;
  if(runScouting){
  iEvent.getByToken(rhoToken2, rhoH2);
  rho2 = *rhoH2;
  }else{// rho=0;
    rho2=0;}

  if(doSignal or (isMC and not era_16)){
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

 // done for all events, no need to reset?
 EventShapeVariables event_algo(event_tracks);
 event_isotropy    = event_algo.isotropy();
 event_sphericity  = event_algo.sphericity();
 event_circularity = event_algo.circularity();

 
  
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


 tree->Fill();	
	
}


void ScoutingNanoAOD::beginJob() {
  
}

void ScoutingNanoAOD::endJob() {
}

void ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  

  //we need to initalise the menu each run (menu can and will change on run boundaries)           
  //for L1
  bool changed=false;
  hltPSProv_.init(iRun,iSetup,hltProcess_,changed);
  const l1t::L1TGlobalUtil& l1GtUtils = hltPSProv_.l1tGlobalUtil();
  std::cout <<"l1 menu "<<l1GtUtils.gtTriggerMenuName()<<" version "<<l1GtUtils.gtTriggerMenuVersion()<<" comment "<<std::endl;
  std::cout <<"hlt name "<<hltPSProv_.hltConfigProvider().tableName()<<std::endl;

}

void ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
edm::Handle<GenLumiInfoHeader> genLumiInfoHead;
    iLumi.getByToken(genLumiInfoHeadTag_, genLumiInfoHead);
    if (!genLumiInfoHead.isValid())
      edm::LogWarning("LHETablesProducer")
          << "No GenLumiInfoHeader product found, will not fill generator model string.\n";
    
    if (genLumiInfoHead.isValid()) {
      label = genLumiInfoHead->configDescription();
      label = std::string("GenModel_") + label;
    }
}

void ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

int ScoutingNanoAOD::getCharge(int pdgId) {
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
bool ScoutingNanoAOD::jetIDoff(const reco::PFJet &pfjet){
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    //moved HT cut
    TLorentzVector jet; 
    jet.SetPtEtaPhiM(pfjet.pt(), pfjet.eta(), pfjet.phi(), pfjet.mass() );
    
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
bool ScoutingNanoAOD::jetID(const ScoutingPFJet &pfjet){
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
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

void ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD);
