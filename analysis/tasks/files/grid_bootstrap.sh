#!/usr/bin/env bash

load_replica() {
    local remote_base="$1"
    local bundle_re="$2"
    local arc_path="$3"

    local arc="$( gfal-ls "$remote_base" | grep -Po "$bundle_re" | shuf -n 1 )"
    if [ -z "$arc" ]; then
        >&2 echo "could not determine archive to load from $remote_base"
        return "1"
    fi

    gfal-copy "$remote_base/$arc" "$arc_path"
    if [ "$?" != "0" ]; then
        >&2 echo "could not load archive $arc from $remote_base"
        return "1"
    fi
}

action() {
    # figure out distribution version
    export JTSF_DIST_VERSION="slc$( lsb_release -rs | head -c 1 )"

    #
    # set env variables
    #

    export PATH_ORIG="$PATH"
    export PYTHONPATH_ORIG="$PYTHONPATH"
    export LD_LIBRARY_PATH_ORIG="$LD_LIBRARY_PATH"
    export GFAL_PLUGIN_DIR_ORIG="$GFAL_PLUGIN_DIR"

    export JTSF_DATA="$TMP/jtsf_data"
    export JTSF_SOFTWARE="$JTSF_DATA/$JTSF_DIST_VERSION/software"
    export SANDBOX_JTSF_SOFTWARE="$JTSF_DATA/{{sandbox_jtsf_dist_version}}/software"

    export JTSF_STORE="$JTSF_DATA/store"
    export JTSF_LOCAL_CACHE="$JTSF_DATA/cache"

    export JTSF_GRID_USER="{{jtsf_grid_user}}"

    export JTSF_CMSSW_SETUP="{{jtsf_cmssw_setup}}"
    export SANDBOX_CMSSW_VERSION="{{sandbox_cmssw_version}}"
    export SANDBOX_CMSSW_BASE="$JTSF_DATA/cmssw/{{sandbox_scram_arch}}/$SANDBOX_CMSSW_VERSION"

    export JTSF_ON_GRID="1"

    mkdir -p "$JTSF_DATA"

    if [ $JTSF_DIST_VERSION == "slc7" ]; then # need to fix gfal before anything is downloaded
        export PYTHONPATH="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib64/python2.7/site-packages:/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib/python2.7/site-packages:$PYTHONPATH"
        export PATH="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/bin/:$PATH"
        export GLOBUS_THREAD_MODEL="none"
        export LD_LIBRARY_PATH="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib64:/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib:$LD_LIBRARY_PATH"
        export GFAL_PLUGIN_DIR="/cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/usr/lib64/gfal2-plugins"
    fi

    #
    # load CMSSW
    #

    mkdir -p "$( dirname "$SANDBOX_CMSSW_BASE" )"
    cd "$( dirname "$SANDBOX_CMSSW_BASE" )"
    load_replica "{{cmssw_base_url}}" "$SANDBOX_CMSSW_VERSION\.\d+\.tgz" "cmssw.tgz"
    cd "$TMP"

    #
    # load the software bundle
    #

    mkdir -p "$JTSF_SOFTWARE"
    cd "$JTSF_SOFTWARE"
    load_replica "{{software_base_url}}" "software\.\d+\.tgz" "software.tgz"
    tar -xzf "software.tgz"
    rm "software.tgz"
    cd "$TMP"

    if [[ "$JTSF_SOFTWARE" != "$SANDBOX_JTSF_SOFTWARE" ]]; then
        mkdir -p "$SANDBOX_JTSF_SOFTWARE"
        cd "$SANDBOX_JTSF_SOFTWARE"
        load_replica "{{sandbox_software_base_url}}" "software\.\d+\.tgz" "software.tgz"
        tar -xzf "software.tgz"
        rm "software.tgz"
        cd "$TMP"
    fi

    #
    # load the repo bundle
    #

    load_replica "{{repo_base}}" "jet-tagging-sf\.{{repo_checksum}}\.\d+\.tgz" "repo.tgz"
    tar -xzf "repo.tgz"
    rm "repo.tgz"

    # copy user proxy
    cp $X509_USER_PROXY "grid.proxy"
    export SANDBOX_X509_USER_PROXY="$TMP/grid.proxy"

    # source the repo setup
    source "jet-tagging-sf/setup.sh"
    cd $LAW_JOB_HOME
}
action "$@"
