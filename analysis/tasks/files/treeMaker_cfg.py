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
        "campaign",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "campaign which the dataset to process belongs to",
    )
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
        "eTriggers",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "e triggers to use",
    )
    options.register(
        "muTriggers",
        [],
        VarParsing.multiplicity.list,
        VarParsing.varType.string,
        "mu triggers to use",
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
        "jerPtResolutionFile",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "JER pt resolution file",
    )
    options.register(
        "jerScaleFactorFile",
        "",
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "JER scale factor file",
    )
    options.register(
        "deepCSVWP",
        0.,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.float,
        "Working point to count number of deepcsv tagged jets",
    )
    options.register(
        "deepJetWP",
        0.,
        VarParsing.multiplicity.singleton,
        VarParsing.varType.float,
        "Working point to count number of deepjet tagged jets",
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

    miniAODProcess = "RECO" if options.isData else "PAT"

    # some default collections
    electronCollection = cms.InputTag("slimmedElectrons")
    muonCollection = cms.InputTag("slimmedMuons")
    metCollection = cms.InputTag("slimmedMETs")
    jetCollection = cms.InputTag("slimmedJets")
    metFilterBitsCollection = cms.InputTag("TriggerResults", "", miniAODProcess)

    # message logger
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

    # source defintion
    process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))

    # good run and lumi selection
    if options.isData and options.lumiFile:
        lumiList = LumiList(filename=options.lumiFile)
        process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()

    # standard sequences with global tag
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
    # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    params = {
            "isMiniAOD": True,
            "applyEnergyCorrections": False,
            "applyVIDOnCorrectedEgamma": False,
    }
    if options.campaign == "Run2_pp_13TeV_Legacy18":
        params["era"] = "2018-Prompt"
    elif options.campaign == "Run2_pp_13TeV_Legacy17":
        params["era"] = "2017-Nov17ReReco"
    elif options.campaign == "Run2_pp_13TeV_Legacy16":
        params["runEnergyCorrections"] = False
        params["era"] = "2016-Legacy"
    elif options.campaign == "Run2_pp_13TeV_UltraLegacy17":
        params["era"] = "2017-UL"
    else:
        raise ValueError("Unknown campaign {}".format(options.campaign))
    setupEgammaPostRecoSeq(process, **params)

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

    # updated MET Filter:
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    if options.campaign in ["Run2_pp_13TeV_Legacy18", "Run2_pp_13TeV_Legacy17"]:
        process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

        baddetEcallist = cms.vuint32(
            [872439604,872422825,872420274,872423218,
            872423215,872416066,872435036,872439336,
            872420273,872436907,872420147,872439731,
            872436657,872420397,872439732,872439339,
            872439603,872422436,872439861,872437051,
            872437052,872420649,872422436,872421950,
            872437185,872422564,872421566,872421695,
            872421955,872421567,872437184,872421951,
            872421694,872437056,872437057,872437313]
        )

        process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
            "EcalBadCalibFilter",
            EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
            ecalMinEt        = cms.double(50.),
            baddetEcal    = baddetEcallist,
            taggingMode = cms.bool(True),
            debug = cms.bool(False)
        )
        seq += process.ecalBadCalibReducedMINIAODFilter

    # MET correction
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    params = {
        "isData" : options.isData,
        "jecUncFile" : os.path.basename(options.jesUncFiles[0]),
        "electronColl" : electronCollection.value(),
        "muonColl" : muonCollection.value(),
        "jetCollUnskimmed" : jetCollection.value(),
    }
    if options.campaign == "Run2_pp_13TeV_Legacy17":
        params["fixEE2017"] = True
        params["fixEE2017Params"] = {"userawPt": True, "ptThreshold": 50.0, "minEtaThreshold": 2.65, "maxEtaThreshold": 3.139}

    runMetCorAndUncFromMiniAOD(process, **params)
    seq += process.fullPatMetSequence
    metCollection = cms.InputTag("slimmedMETs", "", process.name_())

    # add DeepJet discriminators
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

    if options.campaign != "Run2_pp_13TeV_Legacy18":
        updateJetCollection(
           process,
           jetSource = jetCollection,
           pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
           svSource = cms.InputTag('slimmedSecondaryVertices'),
           # Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
           jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
           btagDiscriminators = [
              'pfDeepFlavourJetTags:probb',
              'pfDeepFlavourJetTags:probbb',
              'pfDeepFlavourJetTags:problepb',
              'pfDeepFlavourJetTags:probc',
              'pfDeepFlavourJetTags:probuds',
              'pfDeepFlavourJetTags:probg'
              ],
           postfix='NewDFTraining'
        )
        process.deepFlavour = cms.Task(
                             process.patJetCorrFactorsNewDFTraining,
                             process.updatedPatJetsNewDFTraining,
                             process.patJetCorrFactorsTransientCorrectedNewDFTraining,
                             process.updatedPatJetsTransientCorrectedNewDFTraining,
                             process.pfDeepFlavourJetTagsNewDFTraining,
                             process.pfDeepFlavourTagInfosNewDFTraining,
                             process.pfDeepCSVTagInfosNewDFTraining,
                             process.selectedUpdatedPatJetsNewDFTraining,
                             process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining,
                             process.pfImpactParameterTagInfosNewDFTraining
                             )
        seq.associate(process.deepFlavour)
        jetCollection = cms.InputTag("selectedUpdatedPatJetsNewDFTraining", "", process.name_())

    # L1 prefiring weight
    if options.campaign.endswith(("16", "17")) and not options.isData:
        applyL1Weights = True
        if options.campaign.endswith("17"):
            data_era = "2017BtoF"
        elif options.campaign.endswith("16"):
            data_era = "2016BtoH"
        else:
            raise ValueError("campaign {} should not have l1 prefiring weights applied".format(options.campaign))

        from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
        process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
            DataEra = cms.string(data_era),
            UseJetEMPt = cms.bool(False),
            PrefiringRateSystematicUncty = cms.double(0.2),
            SkipWarnings = False)
        seq += process.prefiringweight
    else:
        applyL1Weights = False

    # deterministic seeds
    process.load("PhysicsTools.PatUtils.deterministicSeeds_cfi")
    process.deterministicSeeds.produceCollections = cms.bool(True)
    process.deterministicSeeds.produceValueMaps = cms.bool(False)
    process.deterministicSeeds.electronCollection = electronCollection
    process.deterministicSeeds.muonCollection = muonCollection
    #process.deterministicSeeds.tauCollection = tauCollection
    #process.deterministicSeeds.photonCollection = photonCollection
    process.deterministicSeeds.jetCollection = jetCollection
    process.deterministicSeeds.METCollection = metCollection
    seq += process.deterministicSeeds

    # overwrite output collections
    muonCollection = cms.InputTag("deterministicSeeds", "muonsWithSeed", process.name_())
    jetCollection = cms.InputTag("deterministicSeeds", "jetsWithSeed", process.name_())
    metCollection = cms.InputTag("deterministicSeeds", "METsWithSeed", process.name_())
    electronCollection = cms.InputTag("deterministicSeeds", "electronsWithSeed", process.name_())

    # load and configure the tree maker
    process.load("JetTaggingSF.JetTaggingSF.treeMaker_cfi")
    process.treeMaker.verbose = cms.untracked.bool(False)
    process.treeMaker.outputFile = cms.string(options.__getattr__("outputFile", noTags=True))
    process.treeMaker.campaign = cms.string(options.campaign)
    process.treeMaker.metaDataFile = cms.string(options.metaDataFile)
    process.treeMaker.isData = cms.bool(options.isData)
    process.treeMaker.leptonChannel = cms.string(options.leptonChannel)
    process.treeMaker.eeTriggers = cms.vstring(options.eeTriggers)
    process.treeMaker.emuTriggers = cms.vstring(options.emuTriggers)
    process.treeMaker.mumuTriggers = cms.vstring(options.mumuTriggers)
    process.treeMaker.eTriggers = cms.vstring(options.eTriggers)
    process.treeMaker.muTriggers = cms.vstring(options.muTriggers)
    process.treeMaker.metFilters = cms.vstring(options.metFilters)
    process.treeMaker.jesFiles = cms.vstring(options.jesFiles)
    process.treeMaker.jesRanges = cms.vint32(options.jesRanges)
    process.treeMaker.jesUncFiles = cms.vstring(options.jesUncFiles)
    process.treeMaker.jesUncSrcFile = cms.string(options.jesUncSrcFile)
    process.treeMaker.jesUncSources = cms.vstring(options.jesUncSources)
    process.treeMaker.jerPtResolutionFile = cms.string(options.jerPtResolutionFile)
    process.treeMaker.jerScaleFactorFile = cms.string(options.jerScaleFactorFile)
    process.treeMaker.deepJetWP = cms.double(options.deepJetWP)
    process.treeMaker.deepCSVWP = cms.double(options.deepCSVWP)
    process.treeMaker.metFilterBitsCollection = metFilterBitsCollection
    process.treeMaker.electronCollection = electronCollection
    process.treeMaker.muonCollection = muonCollection
    process.treeMaker.metCollection = metCollection
    process.treeMaker.jetCollection = jetCollection
    process.treeMaker.applyHEMFilter = cms.bool(True) if options.campaign == "Run2_pp_13TeV_Legacy18" else cms.bool(False)
    process.treeMaker.applyL1Weights = applyL1Weights

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
