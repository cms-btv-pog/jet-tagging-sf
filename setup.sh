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
    [ -z "$JTSF_CMSSW_SETUP" ] && export JTSF_CMSSW_SETUP="ICHEP18"


    #
    # CMSSW setup
    #

    if [ "$JTSF_CMSSW_SETUP" = "ICHEP18" ]; then
        source "$JTSF_BASE/cmssw/setup_ICHEP18.sh"
    else
        2>&1 echo "unknown JTSF_CMSSW_SETUP '$JTSF_CMSSW_SETUP'"
        return "1"
    fi


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
    # minimal software stack
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
