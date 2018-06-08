# Jet Tagging Scale Factors Measurement

### Setup

```shell
# source the main setup file
# (ICHEP18 is currently the only available CMSSW setup)
JTSF_CMSSW_SETUP=ICHEP18 scram_cores=4 source setup.sh

# create the law task index file
law db

# install CMSSW module(s)
law run InstallCMSSWCode
```

### Create Trees

```shell
# example
law run CreateTrees --dataset data_B_ee --version prod1 --grid-ce CNAF --poll-interval 1 --transfer-logs
```

### Available Grid CEs

- `RWTH`
- `RWTH\_short`
- `DESY`
- `CNAF`
- `IRFU`
- `IIHE`
- `CIEMAT`
