#!/usr/bin/env

action() {
    local origin="$( pwd )"
    local scram_cores="$SCRAM_CORES"
    [ -z "$scram_cores" ] && scram_cores="1"

    export SCRAM_ARCH="${JTSF_DIST_VERSION}_amd64_gcc700"
    export CMSSW_VERSION="CMSSW_10_2_18"
    [ -z "$CMSSW_BASE" ] && export CMSSW_BASE="$JTSF_DATA/cmssw/$SCRAM_ARCH/$CMSSW_VERSION"

    source "/cvmfs/cms.cern.ch/cmsset_default.sh"

    export PYTHONPATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/lcg/root/6.12.07-gnimlf5/lib:$PYTHONPATH

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

            # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes
            git cms-merge-topic cms-egamma:EgammaPostRecoTools
            git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029 #optional but speeds up the photon ID value module so things fun faster
            git cms-merge-topic cms-egamma:slava77-btvDictFix_10210 #fixes the Run2018D dictionary issue, see https://github.com/cms-sw/cmssw/issues/26182
            scram b -j "$scram_cores"

            # E-gamma
            git cms-addpkg EgammaAnalysis/ElectronTools  #check out the package otherwise code accessing it will crash
            rm EgammaAnalysis/ElectronTools/data -rf   #delete the data directory so we can populate it ourselves
            git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data

            # fix for new JER version
            # git cms-merge-topic ahinzmann:resolutionSmearingFix102

            # deterministic seed producer
            #git cms-merge-topic yrath:deterministicSeeds_102X
            git fetch git@github.com:yrath/cmssw.git deterministicSeeds_106X && git cherry-pick aa1ff1709e66e3c44570d562e77ad125559cc6f2

            # MET
            git cms-addpkg RecoMET/METFilters
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
