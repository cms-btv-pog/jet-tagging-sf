#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="${JTSF_DIST_VERSION}_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_10_6_8_patch1"
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
        git clone git@github.com:cms-egamma/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
        cd  EgammaUser/EgammaPostRecoTools
        git checkout master
        cd -

        # Updated MET filter
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond%202018
        git cms-addpkg RecoMET/METFilters

        # deterministics seeds
        git cms-merge-topic yrath:deterministicSeeds

        scram b -j "$scram_cores"

    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
    fi

    cd "$origin"
}
action "$@"
