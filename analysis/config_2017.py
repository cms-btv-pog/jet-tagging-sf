# -*- coding: utf-8 -*-


__all__ = ["campaign", "analysis"]


import order as od
import scinum as sn


# campaign
campaign = od.Campaign("2017_13Tev_25ns", 1, ecm=13, bx=25)

# processes
process_data_ee = od.Process(
    "data_ee", 1,
    is_data=True,
    label="data",
)

# datasets
dataset_data_ee = od.Dataset(
    "data_ee", 1,
    campaign=campaign,
    is_data=True,
    n_files=759,
    n_events=58088760,
    keys=["/DoubleEG/Run2017B-17Nov2017-v1/MINIAOD"],
)

# link processes to datasets
dataset_data_ee.add_process(process_data_ee)

# add the analysis and a config for the 2017 campaign
analysis = od.Analysis("jet-tagging-sf", 1)
cfg = analysis.add_config(campaign=campaign)

# link processes
cfg.add_process(process_data_ee)

# add datasets
cfg.add_dataset(dataset_data_ee)
