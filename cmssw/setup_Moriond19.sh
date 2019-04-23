#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="slc${JTSF_DIST_VERSION}_amd64_gcc700"
    export CMSSW_VERSION="CMSSW_10_2_11"
    [ -z "$CMSSW_BASE" ] && export CMSSW_BASE="$JTSF_DATA/cmssw/$SCRAM_ARCH/$CMSSW_VERSION"

    source "/cvmfs/cms.cern.ch/cmsset_default.sh"

    if [ ! -d "$CMSSW_BASE" ]; then
        mkdir -p "$( dirname "$CMSSW_BASE" )"
        cd "$( dirname "$CMSSW_BASE" )"
        scramv1 project CMSSW "$CMSSW_VERSION"
        cd "$CMSSW_VERSION/src"
        eval `scramv1 runtime -sh`
        scram b

        #
        # custom topics
        #

        git cms-merge-topic cms-egamma:EgammaPostRecoTools
        git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029 #optional but speeds up the photon ID value module so things fun faster
        git cms-merge-topic cms-egamma:slava77-btvDictFix_10210 #fixes the Run2018D dictionary issue, see https://github.com/cms-sw/cmssw/issues/26182

        # E-gamma
        git cms-addpkg EgammaAnalysis/ElectronTools  #check out the package otherwise code accessing it will crash
        rm EgammaAnalysis/ElectronTools/data -rf   #delete the data directory so we can populate it ourselves
        git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
        cd EgammaAnalysis/ElectronTools/data
        git checkout ScalesSmearing2018_Dev
        cd -
        git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev

        # b-tagging
        git cms-addpkg RecoBTag
        git cms-addpkg PhysicsTools/PatAlgos
        git cms-merge-topic rauser:PrunedTraining_NoPuppi_10_2_11
        git clone -b PrunedTraining_NoPuppi https://github.com/emilbols/RecoBTag-Combined RecoBTag/Combined/data

        git clone -b 10_2_X_v1.06 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

        # MET
        git cms-addpkg RecoMET/METFilters

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
    fi

    cd "$origin"

    # set default campaign
    export JTSF_CAMPAIGN="2018_Run2_pp_13TeV_MORIOND19"

    export PYTHONPATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/lcg/root/6.12.07-gnimlf5/lib:$PYTHONPATH
}
action "$@"
