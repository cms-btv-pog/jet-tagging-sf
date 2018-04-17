import FWCore.ParameterSet.Config as cms

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#### caution: use the correct global tag for MC or Data 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8' ## '80X_mcRun2_asymptotic_2016_miniAODv2_v1'  ##MC

# Load the producer for MVA IDs. Make sure it is also added to the sequence!
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.ak4PFCHSL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK4PFchs'),
    srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )

process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
	'ak4PFCHSL1Fastjet', 
        'ak4PFchsL2Relative', 
        'ak4PFchsL3Absolute')
)

######## update JetCollection with the HIP mitigation
##process.load("Configuration.StandardSequences.MagneticField_cff")
##process.load("Configuration.Geometry.GeometryRecoDB_cff")
##
##jetSource = 'slimmedJets'
##jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
##bTagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags', 'pfCombinedMVAV2BJetTags']
##
##from PhysicsTools.PatAlgos.tools.jetTools import *
##updateJetCollection(
##        process,
##        jetSource = cms.InputTag(jetSource),
##        jetCorrections = jetCorrectionsAK4,
####        btagInfos = bTagInfos, #list of bTagInfos
##        btagDiscriminators = bTagDiscriminators, #list of taggers
##        runIVF = True,
##        hipMitigation = True,
##        explicitJTA = False,
####        postfix = postfix
##)
##########

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#        '/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0030B9D6-72C1-E611-AE49-02163E00E602.root'
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/0EA60289-18C4-E611-8A8F-008CFA110AB4.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/101D622A-85C4-E611-A7C2-C4346BC80410.root'
#        '/store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/60000/0005201C-D41B-E611-8E37-002481E0D398.root'
#        '/store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/08F114C5-7A01-E611-B29D-0002C94D54FA.root'
#        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/0251DBB7-201B-E611-8653-0CC47A4F1C2E.root'

#        '/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/00087FEB-236E-E511-9ACB-003048FF86CA.root'
#        '/store/mc/RunIISpring15DR74/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/141B9915-1F08-E511-B9FF-001E675A6AB3.root',
#        '/store/user/puigh/TTHSync/ttjets_phys14_20bx25_withfatjets_v2.root'
            )
)



process.ttHTreeMaker = cms.EDAnalyzer('csvTreeMaker',
#    inSample = cms.int32(2500),##
#    sampleName = cms.string("TTJets"),##
    inSample = cms.int32(2300),##
    sampleName = cms.string("ZJets"),##
#    inSample = cms.int32(2310),##
#    sampleName = cms.string("LowMassZJets"),##
#    inSample = cms.int32(2514),##
#    sampleName = cms.string("singletW"),##
#    inSample = cms.int32(2515),##
#    sampleName = cms.string("singletbarW"),##
#    inSample = cms.int32(2600),##
#    sampleName = cms.string("WW"),##
    XS = cms.double(1.),
    nGen = cms.double(1.),
    lumi = cms.double(10000),

    )

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('csv_treeMaker.root')
)

#process.p = cms.Path(process.electronMVAValueMapProducer * process.ttHTreeMaker)
process.p = cms.Path(process.ttHTreeMaker)
