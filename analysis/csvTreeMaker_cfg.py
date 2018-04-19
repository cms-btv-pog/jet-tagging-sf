# -*- coding: utf-8 -*-

"""
Config file to create CSV SF tuples.
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
    options.setDefault("maxEvents", 100)

    # add custom options
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
    process = cms.Process("CSVSF")
    seq = cms.Sequence()

    # some default collections
    defaultProcess = "RECO" if options.isData else "PAT"
    electronCollection = cms.InputTag("slimmedElectrons", "", defaultProcess)
    muonCollection = cms.InputTag("slimmedMuons", "", defaultProcess)
    metCollection = cms.InputTag("slimmedMETs", "", defaultProcess)
    jetCollection = cms.InputTag("slimmedJets", "", defaultProcess)

    # message logger
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

    # source defintion
    process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))

    # good run and lumi selection
    if options.isData and options.lumiFile:
        lumi_list = LumiList(filename=options.lumiFile)
        process.source.lumisToProcess = lumi_list.getVLuminosityBlockRange()

    # standard seuquences with global tag
    if options.globalTag:
        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
        process.GlobalTag.globaltag = options.globalTag

    # geometry sequences
    process.load("Configuration.StandardSequences.GeometryDB_cff")

    # particle data table
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

    # configure the tfile service
    output_file = options.__getattr__("outputFile", noTags=True)
    process.TFileService = cms.Service("TFileService", fileName=cms.string(output_file))

    # load and configure the csv tree maker
    process.load("jet_tagging_sf.jet_tagging_sf.csvTreeMaker_cfi")
    process.csvTreeMaker.verbose = cms.untracked.bool(True)
    process.csvTreeMaker.isData = cms.bool(options.isData)
    process.csvTreeMaker.leptonChannel = cms.string(options.leptonChannel)
    process.csvTreeMaker.eeTriggers = cms.vstring(options.eeTriggers)
    process.csvTreeMaker.emuTriggers = cms.vstring(options.emuTriggers)
    process.csvTreeMaker.mumuTriggers = cms.vstring(options.mumuTriggers)
    process.csvTreeMaker.metFilters = cms.vstring(options.metFilters)
    process.csvTreeMaker.jesFiles = cms.vstring(options.jesFiles)
    process.csvTreeMaker.jesRanges = cms.vint32(options.jesRanges)
    process.csvTreeMaker.jesUncFiles = cms.vstring(options.jesUncFiles)
    process.csvTreeMaker.jesUncSrcFile = cms.string(options.jesUncSrcFile)
    process.csvTreeMaker.jesUncSources = cms.vstring(options.jesUncSources)
    process.csvTreeMaker.electronCollection = electronCollection
    process.csvTreeMaker.muonCollection = muonCollection
    process.csvTreeMaker.metCollection = metCollection
    process.csvTreeMaker.jetCollection = jetCollection
    process.csvTreeMaker.rhoCollection = cms.InputTag("fixedGridRhoFastjetAll")
    process.csvTreeMaker.eleVIDCollection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight")

    # additional configuration
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

    process.options = cms.untracked.PSet(
        allowUnscheduled=cms.untracked.bool(True),
        wantSummary=cms.untracked.bool(options.summary),
    )

    # tell the process what to run
    process.p = cms.Path(seq + process.csvTreeMaker)

except:
    import traceback
    traceback.print_exc()
    raise
