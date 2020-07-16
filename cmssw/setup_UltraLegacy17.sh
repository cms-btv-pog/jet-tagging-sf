#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    if [ $JTSF_DIST_VERSION == "slc6" ]; then # dummy cmssw environment, only for sandbox launching
        echo "Warning: Cannot setup CMSSW_10_6_X on slc6, creating dummy cmssw environment"
        export CMSSW_VERSION="CMSSW_10_2_17"
    else
        export CMSSW_VERSION="CMSSW_10_6_8_patch1"
    fi

    export SCRAM_ARCH="${JTSF_DIST_VERSION}_amd64_gcc700"
    [ -z "$CMSSW_BASE" ] && export CMSSW_BASE="$JTSF_DATA/cmssw/$SCRAM_ARCH/$CMSSW_VERSION"

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
            if [ $JTSF_DIST_VERSION != "slc6" ]; then
                #
                # custom topics
                #

                # Updated MET filter
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond%202018
                git cms-addpkg RecoMET/METFilters

                # deterministics seeds
                git cms-addpkg PhysicsTools/PatUtils
                git fetch git@github.com:yrath/cmssw.git deterministicSeeds_106X && git cherry-pick aa1ff1709e66e3c44570d562e77ad125559cc6f2

                # ECAL scale and resolution corrections
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations
                git clone git@github.com:cms-egamma/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
                cd  EgammaUser/EgammaPostRecoTools
                git checkout master
            fi
            cd -
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
