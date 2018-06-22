# -*- coding: utf-8 -*-

"""
Helper classes and functions.
"""


import os
import array
import collections

import law


class TreeExtender(object):

    DEFAULT_VALUE = {
        "D": -1e5,
        "F": -1e5,
        "I": int(-1e5),
    }

    class Entry(object):
        pass

    def __init__(self, path, name=None, copy_path=None, force=False, default_type=None,
            write_every=-1):
        super(TreeExtender, self).__init__()
        import ROOT

        if copy_path is not None:
            path = copy_trees(path, copy_path, name_pattern=name, force=force)

        self.path = os.path.expandvars(os.path.expanduser(path))
        self.file = ROOT.TFile.Open(self.path, "UPDATE")
        self.trees = get_trees(self.file, name_pattern=name)

        # check default type
        if isinstance(self.trees[0], ROOT.TNtuple):
            default_type = "F"
        elif default_type is None:
            default_type = "D"
        self.default_type = default_type

        # red existing branches
        self.existing_branches = []
        for tree in self.trees:
            self.existing_branches += [branch.GetName() for branch in tree.GetListOfBranches()]

        # store data for new branches (name -> leaf data), and branches to unpack (name)
        self.new_branches = collections.OrderedDict()
        self.unpack_branches = []

        # other attributes
        self.write_every = write_every

    def __enter__(self):
        return self

    def __exit__(self, cls, err, traceback):
        if self.file and self.file.IsOpen():
            self.file.Close()

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
        """
        if not isinstance(name, (list, tuple, set)):
            names = [name]
        else:
            names = name

        for name in self.existing_branches:
            if name in self.unpack_branches:
                continue
            if law.util.multi_match(name, names):
                self.unpack_branches.append(name)
                break

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
                root_type, python_type = leaf_data[0][1:]
                leaf_string = leaf_data[0][0] + "/" + root_type
                if len(leaf_data) > 1:
                    leaf_string += ":" + ":".join(name for name, _, _ in leaf_data[1:])

                arr = array.array(python_type, len(leaf_data) * [self.DEFAULT_VALUE[root_type]])
                branch = tree.Branch(name, arr, leaf_string)
                defaults = array(python_type, len(leaf_data) * [self.DEFAULT_VALUE[root_type]])

                tree_data.append((branch, arr, defaults))

                setattr(entry, name, arr)

            # unpack branches
            for name in self.unpack_branches:
                branch = tree.GetBranch(name)
                first_leaf = branch.GetListOfLeaves()[0]
                if isinstance(first_leaf, ROOT.TLeafI):
                    root_type = "I"
                elif isinstance(first_leaf, ROOT.TLeafF):
                    root_type = "F"
                else:
                    root_type = "D"
                python_type = root_array_types[root_type]
                arr = array(python_type, branch.GetNleaves() * [self.DEFAULT_VALUE[root_type]])
                tree.SetBranchAddress(name, arr)
                setattr(entry, name, arr)

            # start iterating
            for j in xrange(tree.GetEntries()):
                # read the next entry
                tree.GetEntry(j)

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
    Returns the names of all trees found in a *tfile* that pass *name_pattern*. Files opened by this
    method the first time are not cached.
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
    Returns all trees in a *tfile*. *name_pattern* is used as a matching pattern for all contained
    tree names.

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

    if os.path.exists(dst) and not force:
        raise IOError("destination file '{}' exists, force is False".format(dst))

    # simple cp when all trees should be copied
    tfile_src = ROOT.TFile.Open(src, "READ")
    all_names = get_tree_names(tfile_src)
    names = get_tree_names(tfile_src, name_pattern=name_pattern)
    if len(names) == len(all_names):
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
