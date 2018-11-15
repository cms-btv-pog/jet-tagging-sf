#!/usr/bin/env 

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="slc6_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_9_4_8"
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

        # ECAL scale and resolution corrections
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations
        git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
        git cms-merge-topic cms-egamma:Egamma80XMiniAODV2_946

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
    fi

    cd "$origin"
}
action "$@"
