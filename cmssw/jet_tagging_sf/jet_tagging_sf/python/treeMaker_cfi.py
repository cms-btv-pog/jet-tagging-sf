# -*- coding: utf-8 -*-

"""
TreeMaker module initialization file.
"""


__all__ = ["treeMaker"]


import FWCore.ParameterSet.Config as cms


treeMaker = cms.EDAnalyzer("TreeMaker",
    verbose=cms.untracked.bool(False),
    isData=cms.bool(False),
    leptonChannel=cms.string(""),
    eeTriggers=cms.vstring(),
    emuTriggers=cms.vstring(),
    mumuTriggers=cms.vstring(),
    metFilters=cms.vstring(),
    jesFiles=cms.vstring(),
    jesRanges=cms.vint32(),
    jesUncFiles=cms.vstring(),
    jesUncSrcFile=cms.string(""),
    jesUncSources=cms.vstring(),
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
    rhoCollection=cms.InputTag("fixedGridRhoFastjetAll"),
    eleVIDCollection=cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),
)
