# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: ScoutingAOD --filein file:AOD_QCD_HT300to500.root --fileout file:ScoutingAOD.root --eventcontent AOD --datatier AOD --conditions 106X_upgrade2018_realistic_v11_L1v1-v3 --mc --step RAW2DIGI,RECO --processName=ScoutingAOD --no_exec --python_filename=ScoutingAOD_cfg.py
import FWCore.ParameterSet.Config as cms

# automatically sets up inputFiles, maxEvents, etc.
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')
params.parseArguments()

process = cms.Process('ScoutingAOD')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(params.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(params.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('ScoutingAOD nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string(params.outputFile),
    #outputCommands = process.AODEventContent.outputCommands
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_GEN",
        "keep *_*_*_SIM",
        "keep *_*_*_DIGI2RAW",
        "keep *_*_*_HLT",
        #"keep *_offlinePrimaryVertices__RECO",
        #"keep *_offlinePrimaryVerticesWithBS__RECO",
        #"keep *_TriggerResults__RECO",
        #"keep *_particleFlow_electrons_RECO",
        #"keep *_particleFlow_muons_RECO",
        #"keep *_particleBasedIsolation_gedGsfElectrons_RECO",
        #"keep *_gedGsElectrons__RECO",
        #"keep *_muons__RECO",
        #"keep *_particleFlow__RECO",
        #"keep *_pfMet__RECO",
        #"keep *_inclusiveSecondaryVertices__RECO",
        )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v11_L1v1', '')

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
process.schedule = cms.Schedule(process.endjob_step,process.AODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
