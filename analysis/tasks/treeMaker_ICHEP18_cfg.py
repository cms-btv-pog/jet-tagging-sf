# -*- coding: utf-8 -*-

"""
Config file to create jet tagging scale factor trees.
"""


import os

import FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing


try:
    # create options
    options = VarParsing("python")

    # set defaults of common options
    options.setDefault("inputFiles", "root://xrootd-cms.infn.it//store/data/Run2017B/DoubleEG/MINIAOD/17Nov2017-v1/20000/065312BE-A3D5-E711-A0C7-0CC47A1E0DCC.root")
    options.setDefault("outputFile", "output.root")
    options.setDefault("maxEvents", -1)

    # add custom options
    options.register(
        "metaDataFile",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "path to the meta data file to write",
    )
    options.register(
        "globalTag",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "the global tag to use",
    )
    options.register(
        "lumiFile",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "file for selecting runs and lumis",
    )
    options.register(
        "isData",
        False,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.bool,
        "input dataset contains real data",
    )
    options.register(
        "leptonChannel",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "the lepton channel name when running on real data",
    )
    options.register(
        "eeTriggers",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "ee triggers to use",
    )
    options.register(
        "emuTriggers",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "emu triggers to use",
    )
    options.register(
        "mumuTriggers",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "mumu triggers to use",
    )
    options.register(
        "metFilters",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "MET filters to use",
    )
    options.register(
        "jesFiles",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "txt files containing jes infos",
    )
    options.register(
        "jesRanges",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.int,
        "a flat list of range pairs",
    )
    options.register(
        "jesUncFiles",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "txt files containing the combined jes uncertainty infos",
    )
    options.register(
        "jesUncSrcFile",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "txt file containing the per-source jes uncertainty infos",
    )
    options.register(
        "jesUncSources",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "jes uncertainty sources to consider",
    )
    options.register(
        "reportEvery",
        1000,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.int,
        "number of events after which a report message is written",
    )
    options.register(
        "summary",
        False,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.bool,
        "print a summary at the end?",
    )
    options.parseArguments()

    # sanity checks
    if options.isData and not options.leptonChannel:
        raise Exception("a lepton channel is required when running on real data")

    # create the process and a sequence for additional modules
    process = cms.Process("JTSF")
    seq = cms.Sequence()

    # some default collections
    electronCollection = cms.InputTag("slimmedElectrons")
    muonCollection = cms.InputTag("slimmedMuons")
    metCollection = cms.InputTag("slimmedMETs")
    jetCollection = cms.InputTag("slimmedJets")
    metFilterBitsCollection = cms.InputTag("TriggerResults", "", "RECO")

    # message logger
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

    # source defintion
    process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))

    # good run and lumi selection
    if options.isData and options.lumiFile:
        lumiList = LumiList(filename=options.lumiFile)
        process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()

    # standard seuquences with global tag
    if options.globalTag:
        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
        process.GlobalTag.globaltag = options.globalTag

    # standard and geometry sequences
    process.load("Configuration.StandardSequences.GeometryDB_cff")
    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

    # electron ID on uncorrected electrons
    # no option to configure the electron collection available here
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    setupEgammaPostRecoSeq(
        process,
        isMiniAOD=True,
        applyEnergyCorrections=False,
        applyVIDOnCorrectedEgamma=False,
        era="2017-Nov17ReReco",
    )
    seq += process.egammaScaleSmearSeq
    seq += process.egammaPostRecoSeq
    electronCollection = cms.InputTag("slimmedElectrons", "", process.name_())

    # electron energy calibration
    from RecoEgamma.EgammaTools.calibratedEgammas_cff import calibratedPatElectrons
    process.correctedElectrons = calibratedPatElectrons.clone(
        src=electronCollection,
        produceCalibratedObjs=cms.bool(True),
        semiDeterministic=cms.bool(True),
    )
    seq += process.correctedElectrons
    electronCollection = cms.InputTag("correctedElectrons", "", process.name_())

    # load and configure the tree maker
    process.load("jet_tagging_sf.jet_tagging_sf.treeMaker_cfi")
    process.treeMaker.verbose = cms.untracked.bool(False)
    process.treeMaker.outputFile = cms.string(options.__getattr__("outputFile", noTags=True))
    process.treeMaker.metaDataFile = cms.string(options.metaDataFile)
    process.treeMaker.isData = cms.bool(options.isData)
    process.treeMaker.leptonChannel = cms.string(options.leptonChannel)
    process.treeMaker.eeTriggers = cms.vstring(options.eeTriggers)
    process.treeMaker.emuTriggers = cms.vstring(options.emuTriggers)
    process.treeMaker.mumuTriggers = cms.vstring(options.mumuTriggers)
    process.treeMaker.metFilters = cms.vstring(options.metFilters)
    process.treeMaker.jesFiles = cms.vstring(options.jesFiles)
    process.treeMaker.jesRanges = cms.vint32(options.jesRanges)
    process.treeMaker.jesUncFiles = cms.vstring(options.jesUncFiles)
    process.treeMaker.jesUncSrcFile = cms.string(options.jesUncSrcFile)
    process.treeMaker.jesUncSources = cms.vstring(options.jesUncSources)
    process.treeMaker.metFilterBitsCollection = metFilterBitsCollection
    process.treeMaker.electronCollection = electronCollection
    process.treeMaker.muonCollection = muonCollection
    process.treeMaker.metCollection = metCollection
    process.treeMaker.jetCollection = jetCollection

    # additional configuration
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

    # process options
    process.options = cms.untracked.PSet(
        allowUnscheduled=cms.untracked.bool(True),
        wantSummary=cms.untracked.bool(options.summary),
    )

    # tell the process what to run
    process.p = cms.Path(seq + process.treeMaker)

except:
    import traceback
    traceback.print_exc()
    raise
