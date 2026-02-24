import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
from PhysicsTools.SVJScouting.pdfrecalculator_cfi import PDFRecalculator

# Set parameters externally 
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterMuons', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to contain at least two muons'
)

params.register(
    'reducedInfo', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to store just the reduced information'
)

params.register(
    'trigProcess', 
    'HLT', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData', 
    #'101X_dataRun2_HLT_v7',
    '101X_dataRun2_Prompt_v11', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagMC', 
    '102X_upgrade2018_realistic_v15', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)# check this

params.register(
    'xsec', 
    0.001, 
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)#fix this

params.register(
    'fileList', 
    'none', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'input list of files'
)

params.setDefault(
    'maxEvents', 
    -1
)

params.setDefault(
    'outputFile', 
    'test.root' 
)

params.register(
  "era",
  "2018",
  VarParsing.multiplicity.singleton,VarParsing.varType.string,
  "era"
)

params.register(
    'signal', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not signal is run'
)

params.register(
    'MatrixElementInfo', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not MatrixElementInfo is available'
)
params.register(
    'UL2016preVFP', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample belongs to the 2016 preVFP era'
)


# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100 #1000 
#debug
# process.MessageLogger.categories.append('ScoutingNanoAOD_fromMiniAOD')  # Add your category
# process.MessageLogger.categories.append('PDFRecalculator')  # Add another category if needed
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(0)  # Disable general INFO messages
# )
# process.MessageLogger.cerr.ScoutingNanoAOD_fromMiniAOD = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)  # Enable all messages for this category
# )
# process.MessageLogger.cerr.PDFRecalculator = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)  # Enable all messages for this category
# )
#debug

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    #Rethrow = cms.untracked.vstring()
    FailPath = cms.untracked.vstring("ProductNotFound")
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(params.maxEvents) )

# Input EDM files
#list = FileUtils.loadListFromFile(options.inputFiles)
#readFiles = cms.untracked.vstring(*list)

if params.fileList == "none" : readFiles = params.inputFiles
else : 
    readFiles = cms.untracked.vstring( FileUtils.loadListFromFile (os.environ['CMSSW_BASE']+'/src/PhysicsTools/ScoutingNanoAOD/test/'+params.fileList) )
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(readFiles) 
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC : 
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(params.outputFile)
)

# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheInfo = cms.InputTag("externalLHEProducer"),
    genInfo = cms.InputTag("generator"),
    useLHEWeights = cms.bool(params.useWeights)
)

from PhysicsTools.SVJScouting.addSVJPDFsAndScales_cff import addSVJPDFsAndScales
# check if these settings are correct
addSVJPDFsAndScales(process, debug_flag=True, recalculatePDFs_flag=True, recalculateScales_flag=True, nEM_flag=2, nQCD_flag=0, pythiaSettings_flag='', pdfSetName_flag='NNPDF31_nnlo_as_0118_mc_hessian_pdfas')


# # Load the PDFRecalculator
# process.load("PhysicsTools.SVJScouting.pdfrecalculator_cfi")
# #process.PDFRecalculator = cms.EDProducer("PDFRecalculator")
# print(process.PDFRecalculator)



HLTInfo = [
    "DST_DoubleMu1_noVtx_CaloScouting_v*",
    "DST_DoubleMu3_noVtx_CaloScouting_v*",
    "DST_DoubleMu3_noVtx_Mass10_PFScouting_v*",
    "DST_L1HTT_CaloScouting_PFScouting_v*",
    "DST_CaloJet40_CaloScouting_PFScouting_v*",
    "DST_HT250_CaloScouting_v*",
    "DST_HT410_PFScouting_v*",
    "DST_HT450_PFScouting_v*"]
L1Info = [
    'L1_HTT200er',
    'L1_HTT255er',
    'L1_HTT280er',
    'L1_HTT320er',
    'L1_HTT360er',
    'L1_HTT400er',
    'L1_HTT450er',
    'L1_SingleJet180',
    'L1_SingleJet200',
    'L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5',
    'L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5',
    'L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5',
    'L1_ETT2000']
runSig = False
if params.signal:
  runSig = True


# Load Puppi
#
process.load('CommonTools/PileupAlgos/Puppi_cff')
# The following modification ensures that tune V15 is setup for puppi
# https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_26/CommonTools/PileupAlgos/python/Puppi_cff.py#L123-L130
process.puppi.EtaMinUseDeltaZ = 2.4
process.puppi.PtMaxCharged = 20.
process.puppi.PtMaxNeutralsStartSlope = 20.
process.puppi.NumOfPUVtxsForCharged = 2
process.puppi.algos[0].etaMin = [-0.01]

process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.clonePackedCands   = cms.bool(True)
process.puppi.useExistingWeights = cms.bool(False)

#define a new jet collection
from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJetsPuppi
process.ak8PFJetsPuppi  = ak8PFJetsPuppi.clone (jetPtMin = 100.)


#and add the jet collection with
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(process,labelName = 'AK8PFPUPPI', jetSource = cms.InputTag('ak8PFJetsPuppi'), algo = 'AK', rParam=0.8, genJetCollection=cms.InputTag('slimmedGenJetsAK8'),pfCandidates = cms.InputTag('packedPFCandidates'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    muSource =cms.InputTag( 'slimmedMuons'),
    elSource = cms.InputTag('slimmedElectrons')
)


process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD_fromMiniAOD',
    doL1              = cms.bool(False),
    doData            = cms.bool(not params.isMC and not params.signal),
    doSignal          = cms.bool(runSig),
    addMatrixElementInfo = cms.bool(params.MatrixElementInfo),
    isMC              = cms.bool(params.isMC),
    era = cms.string(params.era),
    stageL1Trigger    = cms.uint32(2),

    hltProcess=cms.string("HLT"),
    bits              = cms.InputTag("TriggerResults", "", "HLT"),
    
    triggerresults   = cms.InputTag("TriggerResults", "", params.trigProcess),
    triggerConfiguration = cms.PSet(
    	hltResults               = cms.InputTag('TriggerResults','','HLT'),
    	l1tResults               = cms.InputTag(''),
    	daqPartitions            = cms.uint32(1),
    	l1tIgnoreMaskAndPrescale = cms.bool(False),
    	throw                    = cms.bool(False)
  	),
    ReadPrescalesFromFile = cms.bool( False ),
    AlgInputTag = cms.InputTag("gtStage2Digis"),
    l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
    l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
    l1Seeds           = cms.vstring(L1Info),
    hltSeeds          = cms.vstring(HLTInfo),
    GetLumiInfoHeader=cms.InputTag("generator"),

    #PDF weights
    PDFweights = cms.InputTag("PDFRecalculator", "LHEPdfWeight"),
    ScaleWeights = cms.InputTag("PDFRecalculator", "LHEScaleWeight"),

    #scouting objects
    muons             = cms.InputTag("hltScoutingMuonPackerCalo"),
    electrons         = cms.InputTag("hltScoutingEgammaPacker"),
    photons           = cms.InputTag("hltScoutingEgammaPacker"),
    pfcands           = cms.InputTag("hltScoutingPFPacker"),
    pfjets            = cms.InputTag("hltScoutingPFPacker"),
    vertices          = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
    metPt             = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),
    metPhi            = cms.InputTag("hltScoutingPFPacker", "pfMetPhi"),

    #HLT AK4 PF jets
    jetAK4ScoutPtMin=cms.double(15),
    
    #HLT AK8 PF jets
    jetAK8ScoutPtMin=cms.double(100),

    #offline objects
    pfcandsReco=cms.InputTag("packedPFCandidates"),
    ak4pfjetsReco=cms.InputTag("slimmedJetsPuppi"),
    jetAK4PtMin=cms.double(15),

    verticesReco=cms.InputTag('offlineSlimmedPrimaryVertices'),
    electronsReco=cms.InputTag("slimmedElectrons"),
    muonsReco=cms.InputTag("slimmedMuons"),
    metReco=cms.InputTag("slimmedMETs"),
    PuppimetReco=cms.InputTag("slimmedMETsPuppi"),

    #Puppi AK8 PF
    ak8pfjetsReco=cms.InputTag("ak8PFJetsPuppi"),
    jetAK8PtMin=cms.double(100),

    #gen info and pileup
    genak4jets        = cms.InputTag("slimmedGenJets"),
    genak8jets        = cms.InputTag("slimmedGenJetsAK8"),
    pileupinfo        = cms.InputTag("addPileupInfo"),
    pileupinfo_sig    = cms.InputTag("slimmedAddPileupInfo"),
    geneventinfo     = cms.InputTag("generator"),
    gens             = cms.InputTag("prunedGenParticles"),
    rho               = cms.InputTag("fixedGridRhoFastjetAllScouting"),
    rho2              = cms.InputTag("hltScoutingPFPacker","rho"),
    genMet            = cms.InputTag("genMetTrue"),

)


if(params.isMC):
    # taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe

    if (params.era == "2016") and (params.UL2016preVFP):
        prefiring_era = "2016preVFP"
    elif (params.era == "2016") and (not params.UL2016preVFP):
        prefiring_era = "2016postVFP"
    else:
        prefiring_era = params.era


    prefiring_dict = {
        '2016ULpreVFP' : {'ECAL':'UL2016preVFP', 'Muon':'2016preVFP'},
        '2016ULpostVFP' : {'ECAL':'UL2016postVFP', 'Muon':'2016postVFP'},
        '2017' : {'ECAL':'UL2017BtoF', 'Muon':'20172018'},
        '2018' : {'ECAL':'None', 'Muon':'20172018'}
    }

    from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
    process.prefiringweight = l1PrefiringWeightProducer.clone(
    ThePhotons           = cms.InputTag("hltScoutingEgammaPacker"),
    TheMuons             = cms.InputTag("hltScoutingMuonPacker"),
    TheJets            = cms.InputTag("hltScoutingPFPacker"),
    DataEraECAL = cms.string(prefiring_dict[prefiring_era]["ECAL"]),
    DataEraMuon = cms.string(prefiring_dict[prefiring_era]["Muon"]),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
    )
    process.p = cms.Path(process.puppi   * process.ak8PFJetsPuppi * process.prefiringweight * process.theoryWeightsMC * process.mmtree)
else:

    process.p = cms.Path(process.puppi  * process.ak8PFJetsPuppi * process.theoryWeightsMC * process.mmtree)



