#!/usr/bin/env bash

action() {
    #
    # global variables
    #

    export JTSF_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && /bin/pwd )"
    [ -z "$JTSF_DATA" ] && export JTSF_DATA="/user/public/jet-tagging-sf"
    export JTSF_SOFTWARE="$JTSF_DATA/software"
    export JTSF_STORE="$JTSF_DATA/store"
    export JTSF_LOCAL_CACHE="$JTSF_DATA/cache"


    #
    # helper functions
    #

    _install_pip() {
        pip install --ignore-installed --prefix "$JTSF_SOFTWARE" "$1"
    }

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:$PYTHONPATH"
    }

    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:$PATH"
    }


    #
    # CMSSW setup
    #

    source "/cvmfs/cms.cern.ch/cmsset_default.sh"
    export SCRAM_ARCH="slc6_amd64_gcc630"
    export CMSSW_VERSION="CMSSW_9_4_6_patch1"
    export CMSSW_BASE="$JTSF_DATA/cmssw/$CMSSW_VERSION"

    if [ ! -d "$CMSSW_BASE" ]; then
        mkdir -p "$( dirname "$CMSSW_BASE" )"
        cd "$( dirname "$CMSSW_BASE" )"
        scramv1 project CMSSW "$CMSSW_VERSION"
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
        scram b
        cd "$JTSF_BASE"
    else
        cd "$CMSSW_BASE/src"
        eval `scramv1 runtime -sh`
        cd "$JTSF_BASE"
    fi


    #
    # install minimal software stack once
    #

    # software paths
    _addbin "$JTSF_SOFTWARE/bin"
    _addpy "$JTSF_SOFTWARE/lib/python2.7/site-packages"

    # software that is used in this project
    if [ ! -d "$JTSF_SOFTWARE" ]; then
        echo "installing development software in $JTSF_SOFTWARE"
        mkdir -p "$JTSF_SOFTWARE"

        _install_pip luigi
        _install_pip six
        _install_pip scinum
        _install_pip order
        LAW_INSTALL_CUSTOM_SCRIPT=1 _install_pip git+https://github.com/riga/law.git

        # gfal2
        cd "$JTSF_SOFTWARE"
        wget https://www.dropbox.com/s/3nylghi0xtqaiyy/gfal2.tgz
        tar -xzf gfal2.tgz
        rm gfal2.tgz
        cd "$JTSF_BASE"
    fi

    # source gfal2
    source "$JTSF_SOFTWARE/gfal2/setup.sh"


    #
    # env setup
    #

    # add _this_ repo
    _addpy "$JTSF_BASE"

    # law setup
    export LAW_HOME="$JTSF_BASE/.law"
    export LAW_CONFIG_FILE="$JTSF_BASE/law.cfg"
    export LUIGI_CONFIG_PATH="$JTSF_BASE/luigi.cfg"
    source "$( law completion )"
}
action "$@"
