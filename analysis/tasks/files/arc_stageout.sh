#!/usr/bin/env bash

action() {
    echo "ARC stageout"

    # check if output_uri is set
    local output_uri="{{output_uri}}"
    if [ -z "$output_uri" ]; then
        2>&1 echo "output_uri empty, abort"
        return "1"
    fi
    echo "output_uri is $output_uri"

    #
    # helpers
    #

    stageout_file() {
        local src_dir="$1"
        local file_name="$2"
        local msg_name="$3"

        if [ -f "$src_dir/$file_name" ]; then
            echo "stageout $msg_name $src_dir/$file_name"
            python -c "import gfal2;\
                ctx = gfal2.creat_context();\
                params = ctx.transfer_parameters();\
                params.overwrite = True;\
                ctx.filecopy(params, 'file://$src_dir/$file_name', '$output_uri/$file_name')"
        else
            2>&1 echo "cannot stageout missing $msg_name $src_dir/$file_name"
        fi
    }


    #
    # actual stageout
    #

    # log file
    stageout_file "$LAW_JOB_INIT_DIR" "{{log_file}}" "log file"
}
action "$@"
