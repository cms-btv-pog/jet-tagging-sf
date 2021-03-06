# -*- coding: utf-8 -*-

"""
TreeMaker module initialization file.
"""


__all__ = ["treeMaker"]


import FWCore.ParameterSet.Config as cms


treeMaker = cms.EDAnalyzer("TreeMaker",
    verbose=cms.untracked.bool(False),
    outputFile=cms.string("output.root"),
    campaign=cms.string(""),
    metaDataFile=cms.string("meta.root"),
    isData=cms.bool(False),
    leptonChannel=cms.string(""),
    eeTriggers=cms.vstring(),
    emuTriggers=cms.vstring(),
    mumuTriggers=cms.vstring(),
    eTriggers=cms.vstring(),
    muTriggers=cms.vstring(),
    metFilters=cms.vstring(),
    jesFiles=cms.vstring(),
    jesRanges=cms.vint32(),
    jesUncFiles=cms.vstring(),
    jesUncSrcFile=cms.string(""),
    jesUncSources=cms.vstring(),
    jerPtResolutionFile=cms.string(""),
    jerScaleFactorFile=cms.string(""),
    deepCSVWP=cms.double(0.),
    deepJetWP=cms.double(0.),
    applyHEMFilter=cms.bool(False),
    applyL1Weights=cms.bool(False),
    genInfoCollection=cms.InputTag("generator"),
    triggerBitsCollection=cms.InputTag("TriggerResults", "", "HLT"),
    metFilterBitsCollection=cms.InputTag("TriggerResults", "", "PAT"),
    pileupInfoCollection=cms.InputTag("slimmedAddPileupInfo"),
    beamSpotCollection=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    electronCollection=cms.InputTag("slimmedElectrons"),
    muonCollection=cms.InputTag("slimmedMuons"),
    metCollection=cms.InputTag("slimmedMETs"),
    jetCollection=cms.InputTag("slimmedJets"),
    genJetCollection=cms.InputTag("slimmedGenJets"),
    rhoCollection=cms.InputTag("fixedGridRhoFastjetAll"),
)
