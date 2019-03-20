#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="slc${JTSF_DIST_VERSION}_amd64_gcc700"
    export CMSSW_VERSION="CMSSW_10_2_11"
    export CMSSW_BASE="$JTSF_DATA/cmssw/$SCRAM_ARCH/$CMSSW_VERSION"

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
        git cms-addpkg RecoBTag/TensorFlow
        git cherry-pick 94ceae257f846998c357fcad408986cc8a039152

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
