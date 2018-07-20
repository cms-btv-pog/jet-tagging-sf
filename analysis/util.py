# -*- coding: utf-8 -*-


__all__ = [
    "calc_checksum", "wget", "call_thread", "determine_xrd_redirector", "parse_leaf_list",
    "get_tree_names", "get_trees", "copy_trees", "TreeExtender", "TreeInFileExtender",
]


import os
import array
import collections
import subprocess
import threading
import Queue
import re

import law


def calc_checksum(*paths, **kwargs):
    exclude = law.util.make_list(kwargs.get("exclude", ["*.pyc", "*.git*"]))
    exclude = " ".join("! -path '{}'".format(p) for p in exclude)

    sums = []
    for path in paths:
        path = os.path.expandvars(os.path.expanduser(path))
        if os.path.isfile(path):
            cmd = "sha1sum \"{}\"".format(path)
        elif os.path.isdir(path):
            cmd = "files=\"$( find \"{}\" -type f {} -print | sort -z )\"; "\
                "(for f in $files; do sha1sum $f; done) | sha1sum".format(path, exclude)
        else:
            raise IOError("file or directory '{}' does not exist".format(path))

        code, out, _ = law.util.interruptable_popen(cmd, stdout=subprocess.PIPE, shell=True,
            executable="/bin/bash")
        if code != 0:
            raise Exception("checksum calculation failed")

        sums.append(out.strip().split(" ")[0])

    if len(sums) == 1:
        return sums[0]
    else:
        cmd = "echo \"{}\" | sha1sum".format(",".join(sums))
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
    code = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash", cwd=cwd,
        stdout=std, stderr=std)[0]
    return code == 0


def call_thread(func, args=(), kwargs=None, timeout=None):
    """
    Execute a function *func* in a thread and aborts the call when *timeout* is reached. *args* and
    *kwargs* are forwarded to the function. The return value is a 3-tuple:
    ``(func(), err, timeout satisified)``.
    """
    def wrapper(q, *args, **kwargs):
        try:
            ret = func(*args, **kwargs), None
        except Exception as e:
            ret = None, str(e)
        q.put(ret)

    q = Queue.Queue(1)

    thread = threading.Thread(target=wrapper, args=(q,) + args, kwargs=kwargs)
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        return None, None, True
    else:
        return q.get() + (False,)


def determine_xrd_redirector(lfn, timeout=30):
    import ROOT
    ROOT.gROOT.SetBatch()

    redirectors = ["xrootd-cms.infn.it", "cms-xrd-global.cern.ch", "cmsxrootd.fnal.gov"]
    pfn = lambda rdr: "root://{}/{}".format(rdr, lfn)

    def check(pfn):
        t = ROOT.TFile.Open(pfn)
        t.Close()

    for rdr in 2 * redirectors:
        _, err, timedout = call_thread(check, (pfn(rdr),), timeout=timeout)
        if not timedout and err is None:
            redirector = rdr
            break
    else:
        raise Exception("could not determine redirector to load %s" % lfn)

    return redirector


root_array_types = {
    "D": "d",
    "F": "f",
    "I": "i",
}


def parse_leaf_list(leaf_list, default_type="D"):
    """
    Parses a string *leaf_list* (e.g. ``"leafA/F:leafB/D:leafC"``) and returns a list of tuples
    containing information in the leaves: ``(name, root_type, python_type)``. Supported ROOT types
    are ``F``, ``D`` and ``I``.
    """
    data = []
    for leaf in leaf_list.split(":"):
        parts = leaf.rsplit("/", 1)
        if len(parts) == 1:
            parts.append(default_type)
        name, root_type = parts
        if root_type not in root_array_types:
            raise Exception("unknown root type: " + root_type)
        data.append((name, root_type, root_array_types[root_type]))
    return data


def get_tree_names(tfile, name_pattern="*"):
    """
    Returns the names of all trees found in a *tfile* that pass *name_pattern*.
    """
    import ROOT

    names = []
    for tkey in tfile.GetListOfKeys():
        name = tkey.GetName()

        if not law.util.multi_match(name, name_pattern):
            continue

        tobj = tfile.Get(name)
        if not isinstance(tobj, ROOT.TTree):
            continue

        names.append(name)

    return names


def get_trees(tfile, name_pattern="*"):
    """
    Returns all trees found in a *tfile* that pass *name_pattern*.

    .. code-block:: python

        tfile = ROOT.TFile("/tmp/file.root")
        get_trees(tfile)
        # -> [<ROOT.TTree object ("myTreeA") at 0x2cb6400>,
        #     <ROOT.TTree object ("myTreeB") at 0x2df6420>,
        #     <ROOT.TTree object ("fooTree") at 0x2aa6480>,
        #     ... ]

        get_trees(tfile, "myTree*")
        # -> [<ROOT.TTree object ("myTreeA") at 0x2cb6400>,
        #     <ROOT.TTree object ("myTreeB") at 0x2df6420>]
    """
    names = get_tree_names(tfile, name_pattern=name_pattern)
    return [tfile.Get(name) for name in names]


def copy_trees(src, dst, name_pattern="*", force=False):
    """
    Copies all trees from a *src* file that match *name_pattern* into an other file *dst*. When the
    target file exists and *force* is *False*, an *IOError* is raised.
    """
    import ROOT

    src = os.path.expandvars(os.path.expanduser(src))
    dst = os.path.expandvars(os.path.expanduser(dst))

    if os.path.exists(dst):
        if not force:
            raise IOError("destination file '{}' exists, force is False".format(dst))
        else:
            os.remove(dst)

    # simple cp when all trees should be copied
    tfile_src = ROOT.TFile.Open(src, "READ")
    all_names = get_tree_names(tfile_src)
    names = get_tree_names(tfile_src, name_pattern=name_pattern)
    if len(names) == len(all_names):
        tfile_src.Close()
        ROOT.TFile.Cp(src, dst, ROOT.kFALSE)
        return dst

    # create the dst directory
    dst_dir = os.path.dirname(dst)
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    # open the dst file
    tfile_dst = ROOT.TFile.Open(dst, "RECREATE")
    tfile_dst.cd()

    # copy all input trees
    for name in names:
        tree = tfile_src.Get(name)
        tfile_dst.cd()
        copy = tree.CloneTree(-1, "fast")
        copy.SetName(tree.GetName())
        copy.Write()

    tfile_dst.Close()
    tfile_src.Close()

    return dst


class TreeExtender(object):

    DEFAULT_VALUE = {
        "D": -1e5,
        "F": -1e5,
        "I": int(-1e5),
    }

    class Entry(object):
        pass

    def __init__(self, trees, default_type=None, write_every=-1):
        super(TreeExtender, self).__init__()
        import ROOT

        self.trees = law.util.make_list(trees)
        self.write_every = write_every

        # check default type
        if isinstance(self.trees[0], ROOT.TNtuple):
            default_type = "F"
        elif default_type is None:
            default_type = "D"
        self.default_type = default_type

        # red existing branches
        self.existing_branches = []
        self.existing_aliases = []
        for tree in self.trees:
            self.existing_branches += [branch.GetName() for branch in tree.GetListOfBranches()]
            if tree.GetListOfAliases():
                self.existing_aliases += [alias.GetName() for alias in tree.GetListOfAliases()]
            # do not read unneeded branches
            tree.SetBranchStatus("*", 0)

        self.existing_branches = list(set(self.existing_branches))
        self.existing_aliases = list(set(self.existing_aliases))

        # store data for new branches (name -> leaf data), and branches + aliases to unpack (name)
        self.new_branches = collections.OrderedDict()
        self.unpack_branches = []
        self.unpack_aliases = []

    def __enter__(self):
        return self

    def __exit__(self, err_type, err_value, traceback):
        for tree in self.trees:
            tree.SetBranchStatus("*", 1)
        return

    def add_branch(self, name, leaf_list=None, unpack=None):
        """ add_branch(name, leaf_list=name, unpack=None)
        Adds a new branch *name* with leaves defined by *leaf_list* which defaults to *name*.
        *unpack* is forwarded to :py:meth:`unpack_branch`.
        """
        if leaf_list is None:
            leaf_list = name
        name = name.split("/", 1)[0]

        if name in self.existing_branches:
            raise ValueError("cannot add branch '{}', already present".format(name))
        if name in self.new_branches:
            raise ValueError("cannot add branch '{}', will be already added".format(name))

        self.new_branches[name] = parse_leaf_list(leaf_list, self.default_type)

        if unpack is not None:
            self.unpack_branch(unpack)

    def unpack_branch(self, name):
        """
        Adds a branch *name* to the list of branches that will be unpacked for reading during
        iteration. It can also be a list of names, a pattern or a list of patterns.
        If *name* is an alias, all contained branch names will be unpacked.
        """
        names = law.util.make_list(name)
        def add_name(existing, unpack):
            for name in existing:
                if name in unpack:
                    continue
                if law.util.multi_match(name, names):
                    unpack.append(name)

        add_name(self.existing_branches, self.unpack_branches)
        add_name(self.existing_aliases, self.unpack_aliases)

    def __iter__(self):
        import ROOT

        for i, tree in enumerate(self.trees):
            # tree data containing new branches, arrays and their defaults in tuples
            tree_data = []

            # wrapper to the current entry for easy array and tree access during iteration
            entry = self.Entry()
            setattr(entry, "_tree_", tree)
            setattr(entry, "_i_", i)

            # new branches
            for name, leaf_data in self.new_branches.items():
                root_type, py_type = leaf_data[0][1:]
                leaf_string = leaf_data[0][0] + "/" + root_type
                if len(leaf_data) > 1:
                    leaf_string += ":" + ":".join(name for name, _, _ in leaf_data[1:])

                arr = array.array(py_type, len(leaf_data) * [self.DEFAULT_VALUE[root_type]])
                branch = tree.Branch(name, arr, leaf_string)
                defaults = array.array(py_type, len(leaf_data) * [self.DEFAULT_VALUE[root_type]])

                tree_data.append((branch, arr, defaults))
                setattr(entry, name, arr)

            # unpack aliases
            add_unpack_branches = []
            for name in self.unpack_aliases:
                alias = tree.GetAlias(name)
                # add branches used in alias to unpacked branches
                used_branches = parse_branch_names(alias, tree)
                add_unpack_branches.extend(used_branches)

            # unpack branches
            for name in self.unpack_branches + add_unpack_branches:
                branch = tree.GetBranch(name)
                first_leaf = branch.GetListOfLeaves()[0]
                if isinstance(first_leaf, ROOT.TLeafI):
                    root_type = "I"
                elif isinstance(first_leaf, ROOT.TLeafF):
                    root_type = "F"
                else:
                    root_type = "D"
                py_type = root_array_types[root_type]
                arr = array.array(py_type, branch.GetNleaves() * [self.DEFAULT_VALUE[root_type]])

                tree.SetBranchStatus(name, 1)
                tree.SetBranchAddress(name, arr)
                setattr(entry, name, arr)

            # define formulas to calculate alias values
            aliases = {}
            for name in self.unpack_aliases:
                alias = tree.GetAlias(name)
                formula = ROOT.TTreeFormula(alias, name, tree)
                aliases[name] = formula

            # start iterating
            for j in xrange(tree.GetEntries()):
                # read the next entry
                tree.GetEntry(j)

                # evaluate aliases
                for name, formula in aliases.items():
                    setattr(entry, name, [formula.EvalInstance(0)])

                # reset values of new branches
                for _, arr, defaults in tree_data:
                    arr[:] = defaults[:]

                # actual yield
                yield entry

                # fill branches
                for branch, _, _ in tree_data:
                    branch.Fill()

                # write the tree
                if self.write_every > 0 and not (j + 1) % self.write_every:
                    tree.Write()

            tree.Write()

def parse_branch_names(expression, tree, expandAliases=True):
    """
    Parses an *expression* string and returns the contained branch names in a list. The *tree*
    is required to validate the found branch names. When multiple
    expressions are passed, a joined list with unique branch names is returned. When
    *expandAliases* is *True*, all expressions are tested for aliases which get expanded.
    """
    expressions = list(expression) if isinstance(expression, (list, tuple, set)) \
        else [expression]

    allBranches = [b.GetName() for b in tree.GetListOfBranches()]
    branches = []

    while expressions:
        expression = expressions.pop(0)
        spaced = re.sub("(\+|\-|\*|\/|\(|\)|\:|\,)", " ", expression)
        parts = [s.strip() for s in spaced.split(" ") if s.strip()]
        for b in parts:
            if expandAliases and tree.GetAlias(b):
                expressions.insert(0, tree.GetAlias(b))
                continue
            if b in allBranches and b not in branches:
                branches.append(b)

    return branches


class TreeInFileExtender(TreeExtender):

    def __init__(self, path, name=None, copy_path=None, force=False, **kwargs):
        import ROOT

        if copy_path is not None:
            path = copy_trees(path, copy_path, name_pattern=name, force=force)

        self.path = os.path.expandvars(os.path.expanduser(path))
        self.file = ROOT.TFile.Open(self.path, "UPDATE")
        trees = get_trees(self.file, name_pattern=name)

        super(TreeExtender, self).__init__(trees, **kwargs)

    def __del__(self):
        self.close()

    def __exit__(self, err_type, err_value, traceback):
        self.close()

    def close(self):
        if self.file and self.file.IsOpen():
            self.file.Close()
