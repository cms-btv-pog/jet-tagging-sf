#!/usr/bin/env bash

action() {
    #
    # global variables
    #

    export JTSF_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && /bin/pwd )"
    export JTSF_SOFTWARE="$JTSF_BASE/software"
    export JTSF_DATA="/user/djschmidt/btaggingSF/jet-tagging-sf-output"
    export JTSF_LOCAL_CACHE="/user/djschmidt/btaggingSF/jet-tagging-sf-cache"


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
    # install minimal software stack once
    #

    if [ ! -d "$JTSF_SOFTWARE" ]; then
        echo "installing development software in $JTSF_SOFTWARE"

        _install_pip luigi
        _install_pip six
        _install_pip scinum
        _install_pip order
        # _install_pip law
        _install_pip git+https://github.com/riga/law.git
    fi


    #
    # env setup
    #

    # software
    _addbin "$JTSF_SOFTWARE/bin"
    _addpy "$JTSF_SOFTWARE/lib/python2.7/site-packages"

    # add _this_ repo
    _addpy "$JTSF_BASE"

    # law setup
    export LAW_HOME="$JTSF_BASE/.law"
    export LAW_CONFIG_FILE="$JTSF_BASE/law.cfg"
    source "$( law completion )"
}
action "$@"
