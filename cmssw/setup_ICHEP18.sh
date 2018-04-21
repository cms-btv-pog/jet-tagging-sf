#!/usr/bin/env 

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="slc6_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_9_4_6_patch1"
    export CMSSW_BASE="$JTSF_DATA/cmssw/$CMSSW_VERSION"

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

        # CutBased ID
        git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP

        # ECAL scale and resolution corrections
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations
        git cms-merge-topic cms-egamma:EgammaPostRecoTools_940

        # electron ID and tool files
        git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git \
            -b CMSSW_9_4_0_pre3_TnP \
            "$CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data"
        git clone https://github.com/Sam-Harper/EgammaAnalysis-ElectronTools.git \
            -b ReReco17NovScaleAndSmearing \
            "$CMSSW_BASE/external/$SCRAM_ARCH/data/EgammaAnalysis/ElectronTools/data"

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
    fi

    cd "$origin"
}
action "$@"
