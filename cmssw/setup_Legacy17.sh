#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="${JTSF_DIST_VERSION}_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_9_4_9"
    [ -z "$CMSSW_BASE" ] && export CMSSW_BASE="$JTSF_DATA/cmssw/$JTSF_CAMPAIGN/$SCRAM_ARCH/$CMSSW_VERSION"

    source "/cvmfs/cms.cern.ch/cmsset_default.sh"

    if [ ! -d "$CMSSW_BASE" ]; then
        mkdir -p "$( dirname "$CMSSW_BASE" )"
        cd "$( dirname "$CMSSW_BASE" )"
        scramv1 project CMSSW "$CMSSW_VERSION"
        cd "$CMSSW_VERSION/src"
        eval `scramv1 runtime -sh`
        scram b

        if [ "$JTSF_ON_GRID" == "1" ]; then # unpack custom installation from .tgz
            cd "$( dirname "$CMSSW_BASE" )"
            if [ -f "cmssw.tgz" ]; then
                tar -xzvf "cmssw.tgz" --directory $CMSSW_BASE
                rm "cmssw.tgz"
            fi
            cd "$CMSSW_BASE/src"
        else
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

            # Updated MET filter
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond%202018
            git cms-addpkg RecoMET/METFilters

            # MET EE fix
            git cms-merge-topic cms-met:METFixEE2017_949_v2

            # deterministics seeds
            git cms-merge-topic yrath:deterministicSeeds
        fi

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
        scram build
    fi

    cd "$origin"
}
action "$@"
