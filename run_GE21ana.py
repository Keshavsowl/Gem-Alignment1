import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('analyzer',Phase2C17I13M9)
# from Configuration.Eras.Era_Phase2_cff import Phase2

# process = cms.Process('analyzer',Phase2)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag

### This is the misalignment part
misalign = False
do_GEM = False
do_CSC = True

if misalign:
  #db_file = 'sqlite_file:dummy_dx1.db'
  # gem_db_file = 'sqlite_file:/afs/cern.ch/work/t/toakhter/public/GEM_Alignment/2024H_prompt_reco_v0.db' #for GEM
  # csc_db_file = 'sqlite_file:csc.db' #for csc alignment only in this case
  # gpr_db_file = 'sqlite_file:gpr.db' #for gpr only in this case
  process.GlobalTag.toGet = cms.VPSet(
    #GE11 rec/tag
    # cms.PSet(
    #     connect = cms.string(gem_db_file),
    #     record = cms.string('GEMAlignmentRcd'),
    #     tag = cms.string('GEMAlignmentRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(gem_db_file),
    #     record = cms.string('GEMAlignmentErrorExtendedRcd'),
    #     tag = cms.string('GEMAlignmentErrorExtendedRcd')
    # )
    #ME11 rec/tag
    # cms.PSet(
    #     connect = cms.string(csc_db_file),
    #     record = cms.string('CSCAlignmentRcd'),
    #     tag = cms.string('CSCAlignmentRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(csc_db_file),
    #     record = cms.string('CSCAlignmentErrorExtendedRcd'),
    #     tag = cms.string('CSCAlignmentErrorExtendedRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(gpr_db_file), 
    #     record = cms.string('GlobalPositionRcd'), 
    #     tag = cms.string('GlobalPositionRcd') #cms.string('IdealGeometry')
    # )
  )

process.gemGeometry.applyAlignment = cms.bool(do_GEM)
process.CSCGeometryESModule.applyAlignment = cms.bool(do_CSC)
################################


#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_prompt', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

# process.maxEvents = cms.untracked.PSet(
#   input = cms.untracked.int32(options.nEvents)
# )
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
			fileNames = cms.untracked.vstring(options.inputFiles), 
			inputCommands = cms.untracked.vstring(
			  "keep *",
        'keep *_*muonGEMDigis*_*_*',
			  "drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			  "drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			)
      # SelectEvents = cms.untracked.PSet(
      #   SelectEvents = cms.vstring('rechit_step')
      # )
		)

#testfile = "/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_83.root"
#process.source.fileNames.append('file:'+testfile)
outfile = "output_test.root"

#process.source.fileNames.append('root://cms-xrd-global.cern.ch/')
#process.source.fileNames.append('file:/eos/cms/store/relval/CMSSW_15_0_0/RelValSingleMuPt10/GEN-SIM-RECO/141X_mcRun4_realistic_v3_STD_RegeneratedGS_Run4D110_noPU-v2/2580000/0b0d313e-56e1-4e64-aa58-8a2e61767bf5.root')

process.source.fileNames.append(
    'root://cms-xrd-global.cern.ch//store/mc/Phase2Spring24DIGIRECOMiniAOD/DYTo2L_Bin-M-50_TuneCP5_14TeV_pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU140_Trk1GeV_140X_mcRun4_realistic_v5-v1/2810000/00021d9f-9ae0-4cf9-96b9-574212bb82dc.root')

#process.source.fileNames.append(
#    'file:/eos/cms/store/relval/CMSSW_15_0_0/RelValZMM_14/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v3/2580000/ff884c68-8c5d-463c-9335-ab6b7ef48f79.root')


process.options = cms.untracked.PSet(
                        TryToContinue = cms.untracked.vstring('ProductNotFound') #SkipEvent parameter does not work for CMSSW_13_3_X and above
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) 

from RecoLocalMuon.CSCSegment.cscSegments_cfi import *
process.cscSegments = cscSegments.clone()

process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')
process.gemRecHits = cms.EDProducer("GEMRecHitProducer",
    recAlgoConfig = cms.PSet(),
    recAlgo = cms.string('GEMRecHitStandardAlgo'),
    gemDigiLabel = cms.InputTag("muonGEMDigis"),
    ge21Off = cms.bool(False),
)
# process.gemGeometry.fromDD4hep = cms.bool(True)
# process.muonGeometryConstants.fromDD4hep = cms.bool(True)


process.analyzer = cms.EDAnalyzer('ge21analyzer', 
	      process.MuonServiceProxy,
      # cscSegments = cms.InputTag("cscSegments"),
       cscSegments = cms.InputTag("slimmedDisplacedMuons"),
	      gemRecHits = cms.InputTag("gemRecHits"),
        # gemRecHits = cms.InputTag("gemRecHits", "", "GEMLocalRECO"), #uncomment this if reconstructing GE21 hits
	      gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("slimmedDisplacedMuons"),#("ALCARECOMuAlCalIsolatedMu:SelectedMuons"),
        # ref_track = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted"),
	      vertexCollection = cms.InputTag("offlinePrimaryVertices"),
        tracker_prop = cms.bool(True),
        CSC_prop = cms.bool(False),
        Segment_prop = cms.bool(True),
        trackerRefit_prop = cms.bool(False),
        SegmentReco_prop = cms.bool(False),
        debug = cms.bool(True),
        isCosmic = cms.bool(False)
)

process.p = cms.Path(process.gemRecHits * process.analyzer)
# process.p = cms.Path(process.gemRecHits * process.analyzer)
