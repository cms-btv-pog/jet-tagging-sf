#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="${JTSF_DIST_VERSION}_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_9_4_12"
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

        # ECAL scale and resolution corrections
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations
        git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
        # Electron VID
        git cms-merge-topic cms-egamma:EgammaID_949

        # new DeepJet Training
        git cms-addpkg RecoBTag/TensorFlow
        git cherry-pick 94ceae257f846998c357fcad408986cc8a039152

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
    fi

    # set default campaign
    export JTSF_CAMPAIGN="2018_Run2_pp_13TeV_MORIOND19legacy"

    cd "$origin"
}
action "$@"
