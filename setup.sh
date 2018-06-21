#!/usr/bin/env bash

action() {
    #
    # global variables
    #

    export JTSF_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && /bin/pwd )"

    # check if we're on lxplus
    [[ "$( hostname )" = lxplus*.cern.ch ]] && JTSF_ON_LXPLUS="1" || JTSF_ON_LXPLUS="0"

    # default data directory
    if [ -z "$JTSF_DATA" ]; then
        if [ "$JTSF_ON_LXPLUS" = "1" ]; then
            JTSF_DATA="$JTSF_BASE/.data"
        else
            JTSF_DATA="/user/public/jet-tagging-sf"
        fi
    fi

    # default grid user
    if [ -z "$JTSF_GRID_USER" ]; then
        if [ "$JTSF_ON_LXPLUS" = "1" ]; then
            JTSF_GRID_USER="$( whoami )"
            echo "setting JTSF_GRID_USER to $JTSF_GRID_USER"
        else
            2>&1 echo "please set the JTSF_GRID_USER to your grid user name and try again"
            return "1"
        fi
    fi

    # other defaults
    [ -z "$JTSF_SOFTWARE" ] && JTSF_SOFTWARE="$JTSF_DATA/software"
    [ -z "$JTSF_STORE" ] && JTSF_STORE="$JTSF_DATA/store"
    [ -z "$JTSF_LOCAL_CACHE" ] && JTSF_LOCAL_CACHE="$JTSF_DATA/cache"
    [ -z "$JTSF_CMSSW_SETUP" ] && JTSF_CMSSW_SETUP="ICHEP18"

    # export variables
    export JTSF_DATA
    export JTSF_ON_LXPLUS
    export JTSF_GRID_USER
    export JTSF_SOFTWARE
    export JTSF_STORE
    export JTSF_LOCAL_CACHE
    export JTSF_CMSSW_SETUP


    #
    # CMSSW setup
    #

    if [ "$JTSF_CMSSW_SETUP" = "ICHEP18" ]; then
        source "$JTSF_BASE/cmssw/setup_ICHEP18.sh" || return "$?"
    else
        2>&1 echo "unknown JTSF_CMSSW_SETUP '$JTSF_CMSSW_SETUP'"
        return "1"
    fi


    #
    # helper functions
    #

    jtsf_install_pip() {
        pip install --ignore-installed --no-cache-dir --prefix "$JTSF_SOFTWARE" "$@"
    }
    export -f jtsf_install_pip

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
    jtsf_install_software() {
        local origin="$( pwd )"
        local mode="$1"

        if [ -d "$JTSF_SOFTWARE" ]; then
            if [ "$mode" = "force" ]; then
                echo "remove software in $JTSF_SOFTWARE"
                rm -rf "$JTSF_SOFTWARE"
            else
                if [ "$mode" != "silent" ]; then
                    echo "software already installed in $JTSF_SOFTWARE"
                fi
                return "0"
            fi
        fi

        echo "installing development software in $JTSF_SOFTWARE"
        mkdir -p "$JTSF_SOFTWARE"

        jtsf_install_pip --no-dependencies uproot
        jtsf_install_pip slackclient
        jtsf_install_pip order
        jtsf_install_pip git+https://github.com/riga/luigi.git@fix/dynamicParamKwargs
        LAW_INSTALL_CUSTOM_SCRIPT="1" jtsf_install_pip --no-dependencies git+https://github.com/riga/law.git

        # gfal2
        cd "$JTSF_SOFTWARE"
        wget https://www.dropbox.com/s/3nylghi0xtqaiyy/gfal2.tgz
        tar -xzf gfal2.tgz
        rm gfal2.tgz
        cd "$origin"
    }
    export -f jtsf_install_software
    jtsf_install_software silent

    # setup gfal2 separately
    source "$JTSF_SOFTWARE/gfal2/setup.sh" || return "$?"


    #
    # env setup
    #

    # add _this_ repo
    _addpy "$JTSF_BASE"

    # law and luigi setup
    export LAW_HOME="$JTSF_BASE/.law"
    export LAW_CONFIG_FILE="$JTSF_BASE/law.cfg"
    source "$( law completion )"

    if [ -z "$JTSF_SCHEDULER_HOST" ]; then
        2>&1 echo "NOTE: \$JTSF_SCHEDULER_HOST is not set, use '--local-scheduler' in your tasks!"
    fi
}
action "$@"
