# Jet Tagging Scale Factors Measurement

### Implementation Info

The *Luigi Analysis Workflow* (**law**) package is used to structure the scale factor measurement.
It is based on the *luigi* pipelining tool, formerly developed by Spotify.
Logical units of the measurement code are structured in so-called `Task` classes.
Tasks produce output `Target`'s and can require each other.
This defines the actual *workflow* to run.
If a task is run from the command-line, all required tasks whose outputs do not exist yet are triggered as well.
See the links below for more information on **luigi**.

**law** is an extension ontop of **luigi**.
It adds useful classes  to work with infrastructure available to the HEP community.
E.g. it supports submission to the WLCG as well HTCondor and LSF batch systems, and also storage of files on remote resources using a bunch of protocols (EOS, xRootD, dCache, GridFTP, WebDav, Dropbox, CERNBox, etc.).

Useful links:

- **luigi**: [Repo](https://github.com/spotify/luigi), [Docs](https://luigi.readthedocs.io), [hello_world.py](https://github.com/spotify/luigi/blob/master/examples/hello_world.py)
- **law**: [Repo](https://github.com/riga/law), [Docs](https://law.readthedocs.io/en/latest/) (in progress), [Examples](https://github.com/riga/law/tree/master/examples)


### Setup

```shell
# source the main setup file
scram_cores=4 source setup.sh

# create the law task index file
# (required for bash/zsh autocompletion)
law index

# install CMSSW modules once
# (only needed when creating new trees)
law run InstallCMSSWCode
```

To let the tasks talk to a central luigi scheduler, you need to set `$JRST_SCHEDULER_HOST`.

### Running Tasks

##### `WriteTrees`

```shell
# example
law run WriteTrees --dataset data_B_ee --version prod1 --grid-ce CNAF --poll-interval 1 --transfer-logs
```

### Available Grid CEs

- `RWTH`
- `RWTH_short`
- `DESY`
- `CNAF`
- `IRFU`
- `IIHE`
- `CIEMAT`
