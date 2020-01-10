#!/usr/bin/env bash

action() {
    #
    # global variables
    #

    export JTSF_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && /bin/pwd )"

    # check if we're on lxplus
    [[ "$( hostname )" = lxplus*.cern.ch ]] && JTSF_ON_LXPLUS="1" || JTSF_ON_LXPLUS="0"
    export JTSF_ON_LXPLUS

    # check if we're on VISPA
    [ ! -z "$( hostname | grep vispa )" ] && JTSF_ON_VISPA="1" || JTSF_ON_VISPA="0"
    export JTSF_ON_VISPA

    # figure out distribution version
    export JTSF_DIST_VERSION="slc$( lsb_release -rs | head -c 1 )"

    # default data directory
    if [ -z "$JTSF_DATA" ]; then
        if [ "$JTSF_ON_LXPLUS" = "1" ]; then
            export JTSF_DATA="$JTSF_BASE/.data"
        elif [ "$JTSF_ON_VISPA" = "1" ]; then
            export JTSF_DATA="/net/scratch/cms/jet-tagging-sf"
        else
            export JTSF_DATA="/user/$( whoami )/jet-tagging-sf"
        fi
    fi

    # default grid user
    if [ -z "$JTSF_GRID_USER" ]; then
        if [ "$JTSF_ON_LXPLUS" = "1" ]; then
            export JTSF_GRID_USER="$( whoami )"
            echo "NOTE: lxplus detected, setting JTSF_GRID_USER to $JTSF_GRID_USER"
        else
            2>&1 echo "please set the JTSF_GRID_USER to your grid user name and try again"
            return "1"
        fi
    fi

    # default CMSSW setup when on VISPA
    [ "$JTSF_ON_VISPA" = "1" ] && export JTSF_CMSSW_SETUP="NONE"

    # other defaults
    [ -z "$JTSF_SOFTWARE" ] && export JTSF_SOFTWARE="$JTSF_DATA/$JTSF_DIST_VERSION/software/$( whoami )"
    [ -z "$JTSF_STORE" ] && export JTSF_STORE="$JTSF_DATA/store"
    [ -z "$JTSF_LOCAL_CACHE" ] && export JTSF_LOCAL_CACHE="$JTSF_DATA/cache"
    [ -z "$JTSF_CMSSW_SETUP" ] && export JTSF_CMSSW_SETUP="Moriond19"

    # law and luigi setup
    export LAW_HOME="$JTSF_BASE/.law"
    export LAW_CONFIG_FILE="$JTSF_BASE/law.cfg"
    [ "$JTSF_ON_LXPLUS" == "0" ] && export LAW_TARGET_TMP_DIR="$JTSF_DATA/tmp"

    if [ "$JTSF_ON_GRID" == "1" ]; then
        export JTSF_LUIGI_WORKER_KEEP_ALIVE="False"
        export JTSF_LUIGI_WORKER_FORCE_MULTIPROCESSING="True"
    else
        export JTSF_LUIGI_WORKER_KEEP_ALIVE="True"
        export JTSF_LUIGI_WORKER_FORCE_MULTIPROCESSING="False"
    fi

    if [ -z "$JTSF_SCHEDULER_HOST" ]; then
        2>&1 echo "NOTE: \$JTSF_SCHEDULER_HOST is not set, use '--local-scheduler' in your tasks!"
    fi


    #
    # CMSSW setup
    #

    if [ "$JTSF_CMSSW_SETUP" = "NONE" ]; then
        echo "NOTE: skipping CMSSW setup"
    else
        source "$JTSF_BASE/cmssw/setup_$JTSF_CMSSW_SETUP.sh" || return "$?"
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

        jtsf_install_pip slackclient
        jtsf_install_pip docutils
        jtsf_install_pip git+https://github.com/riga/order.git@4b78ad6c06caee65f42e470a3c88fb61bba2d8f8
        jtsf_install_pip git+https://github.com/spotify/luigi.git
        LAW_INSTALL_CUSTOM_SCRIPT="1" jtsf_install_pip --no-dependencies git+https://github.com/riga/law.git

        # gfal2
        if [ $JTSF_DIST_VERSION == "slc6" ]; then
            cd "$JTSF_SOFTWARE"
            wget https://www.dropbox.com/s/3nylghi0xtqaiyy/gfal2.tgz
            tar -xzf gfal2.tgz
            rm gfal2.tgz
            cd "$origin"
        fi
    }
    export -f jtsf_install_software
    jtsf_install_software silent

    # setup gfal2 separately
    if [ $JTSF_DIST_VERSION == "slc6" ]; then
        source "$JTSF_SOFTWARE/gfal2/setup.sh" || return "$?"
    else
        export GLOBUS_THREAD_MODEL="none"
        export LD_LIBRARY_PATH="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib64:/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib:$LD_LIBRARY_PATH"
        export PATH="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/bin:/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/sbin:$PATH"
        export PYTHONPATH="$PYTHONPATH:/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib64/python2.7/site-packages/"
    fi

    # add _this_ repo
    _addpy "$JTSF_BASE"

    # source law's bash completion scipt
    source "$( law completion )"
}
action "$@"
