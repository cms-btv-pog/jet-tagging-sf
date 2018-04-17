# -*- coding: utf-8 -*-

__all__ = ["calc_checksum", "wget"]


import os
import subprocess

import law


def calc_checksum(*paths, **kwargs):
    exclude = law.util.make_list(kwargs.get("exclude", ["*.pyc", "*.git*"]))
    exclude = " ".join("! -path '%s'" % p for p in exclude)

    sums = []
    for path in paths:
        path = os.path.expandvars(os.path.expanduser(path))
        if os.path.isfile(path):
            cmd = "sha1sum \"%s\"" % path
        elif os.path.isdir(path):
            cmd = "files=\"$( find \"%s\" -type f %s -print | sort -z )\"; "\
                "(for f in $files; do sha1sum $f; done) | sha1sum" % (path, exclude)
        else:
            raise IOError("file or directory '%s' does not exist" % path)

        code, out, _ = law.util.interruptable_popen(cmd, stdout=subprocess.PIPE, shell=True,
            executable="/bin/bash")
        if code != 0:
            raise Exception("checksum calculation failed")

        sums.append(out.strip().split(" ")[0])

    if len(sums) == 1:
        return sums[0]
    else:
        cmd = "echo \"%s\" | sha1sum" % ",".join(sums)
        code, out, _ = law.util.interruptable_popen(cmd, stdout=subprocess.PIPE, shell=True,
            executable="/bin/bash")
        if code != 0:
            raise Exception("checksum combination failed")

        return out.strip().split(" ")[0]


def wget(url, dst, verbose=False):
    """
    Better version of `urllib.urlretrieve` using wget. Downloads a file from *url* and saves to the
    path *dst*. Returns *True* on success, *False* otherwise.
    """
    dirname, filename = os.path.split(dst)
    cwd = dirname or None
    std = None if verbose else subprocess.PIPE
    cmd = "wget -O {} {}".format(filename, url)
    code = law.interruptable_popen(cmd, shell=True, executable="/bin/bash", cwd=cwd, stdout=std,
        stderr=std)[0]
    return code == 0
