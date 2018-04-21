import enum


class Tagger(enum.IntEnum):
    cMVA = 0
    csv = 1


class LeptonType(enum.IntEnum):
    DoubleEG = -100
    DoubleMuon = -200
    MuonEG = -300
    All = 0


class RunEra(enum.Enum):
    B = 1
    C = 2
    D = 3
    E = 4
    F = 5


class FinalState(enum.IntEnum):
    zjets = 2300
    lowMasszjets = 2310
    # WJetsToLNu = 2400
    ttjets = 2500
    # TToLepton_s = 2510
    # TBarToLepton_s = 2511
    # TToLeptons_t = 2512
    # TBarToLeptons_t = 2513
    singletW = 2514
    singletbarW = 2515
    # TTZJets = 2523
    # TTWJets = 2524
    WW = 2600
