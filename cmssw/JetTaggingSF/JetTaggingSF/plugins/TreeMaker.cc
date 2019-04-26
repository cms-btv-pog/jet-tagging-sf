/*
 * CMSSW module that does some event selection and outputs a tree that can be used to measure
 * jet tagging scale factors.
 *
 * Authors:
 *   - Marcel Rieger
 *   - Yannik Rath
 */

#include <memory>
#include <cmath>
#include <fnmatch.h>

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom3.h"

#include "JetTaggingSF/JetTaggingSF/interface/VarMap.h"

typedef std::string string;
typedef std::vector<string> vstring;
typedef std::vector<int32_t> vint;
typedef std::pair<string, string> stringPair;
typedef math::XYZTLorentzVector LorentzVector;

double EMPTY_VALUE = -1e5;

enum VertexID
{
    V_INVALID = 0,
    V_VALID
};

enum ElectronID
{
    E_INVALID = 0,
    E_LOOSE,
    E_TIGHT
};

enum MuonID
{
    M_INVALID = 0,
    M_LOOSE,
    M_TIGHT
};

enum JetID
{
    J_INVALID = 0,
    J_VALID,
};

enum LeptonChannel
{
    C_INVALID = 0,
    C_EE,
    C_EMU,
    C_MUMU,
//    C_E,
//    C_MU
};

bool comparePt(const pat::Jet& jet1, const pat::Jet& jet2)
{
    return jet1.pt() > jet2.pt();
}

float electronEffectiveArea(const pat::Electron& electron)
{
    // numbers from https://github.com/lsoffi/cmssw/blob/CMSSW_9_4_0_pre3_TnP/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
    double absSCEta = fabs(electron.superCluster()->eta());

    if (absSCEta <= 1.0)
    {
        return 0.1566;
    }
    else if (absSCEta <= 1.479)
    {
        return 0.1626;
    }
    else if (absSCEta <= 2.0)
    {
        return 0.1073;
    }
    else if (absSCEta <= 2.2)
    {
        return 0.0854;
    }
    else if (absSCEta <= 2.3)
    {
        return 0.1051;
    }
    else if (absSCEta <= 2.4)
    {
        return 0.1204;
    }
    else if (absSCEta <= 5.0)
    {
        return 0.1524;
    }
    else
    {
        return 0.;
    }
}

class TreeMaker : public edm::EDAnalyzer
{
public:
    explicit TreeMaker(const edm::ParameterSet&);
    ~TreeMaker();

private:
    // interface methods
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // methods for handling variables and output objects
    void setupJetCorrectionObjects();
    void setupVariables();
    void setVariable(const string& name, double value);

    // selection methods
    bool metFilterSelection(const edm::Event&);
    bool triggerSelection(const edm::Event&, LeptonChannel&);
    bool electronSelection(const edm::Event&, reco::Vertex&, double, std::vector<pat::Electron>&,
        std::vector<pat::Electron>&);
    bool muonSelection(const edm::Event&, reco::Vertex&, std::vector<pat::Muon>&,
        std::vector<pat::Muon>&);
    bool leptonSelection(std::vector<pat::Electron>&, std::vector<pat::Electron>&,
        std::vector<pat::Muon>&, std::vector<pat::Muon>&, reco::RecoCandidate*&,
        reco::RecoCandidate*&, LeptonChannel&, bool&);
    bool jetMETSelection(const edm::Event&, double, reco::RecoCandidate*, reco::RecoCandidate*,
        const pat::MET&, const string&, const string&, std::vector<pat::Jet>&, pat::MET&, bool);
    VertexID vertexID(reco::Vertex&);
    ElectronID electronID(pat::Electron&, reco::Vertex&, double);
    MuonID muonID(pat::Muon&, reco::Vertex&);
    JetID jetID(pat::Jet&, reco::RecoCandidate*, reco::RecoCandidate*, bool, bool);
    bool tightJetID_2016(pat::Jet&);
    bool tightJetID_2017(pat::Jet&);
    bool tightJetID_2018(pat::Jet&);

    // random helpers
    double readGenWeight(const edm::Event&);
    float readPU(const edm::Event&);
    double readRho(const edm::Event&);
    void applyJES(pat::Jet&, const string&, const string&, int64_t, double);
    void applyJER(pat::Jet&, const std::vector<reco::GenJet>*, const string&, const string&,
        double);

    // options
    bool verbose_;
    string outputFile_;
    string campaign_;
    string metaDataFile_;
    bool isData_;
    string leptonChannel_;
    vstring eeTriggers_;
    vstring emuTriggers_;
    vstring mumuTriggers_;
    //vstring eTriggers_;
    //vstring muTriggers_;
    vstring metFilters_;
    vstring jesFiles_;
    vint jesRanges_;
    vstring jesUncFiles_;
    string jesUncSrcFile_;
    vstring jesUncSources_;
    string jerPtResolutionFile_;
    string jerScaleFactorFile_;
    bool (TreeMaker::*tightJetID_)(pat::Jet&);
    double maxJetEta_;

    // tokens
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
    edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > vertexToken_;
    edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
    edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
    edm::EDGetTokenT<std::vector<pat::MET> > metToken_;
    edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet> > genJetToken_;
    edm::EDGetTokenT<double> rhoToken_;

    // additional members
    VarMap varMap_;
    TFile* tfile_;
    TFile* tfileMeta_;
    TTree* tree_;
    TH1F* eventHist_;
    TH1F* selectedEventHist_;
    TH1F* weightHist_;
    TH1F* selectedWeightHist_;
    TH1F* cutflowHist_;
    TH1F* pileupHist_;
    size_t nJESRanges_;
    size_t nJESFilesPerRange_;
    std::vector<std::pair<string, string> > jetVariations_;
    std::vector<FactorizedJetCorrector*> jetCorrectors_; // per jes range
    std::vector<JetCorrectionUncertainty*> jetCorrectorUncs_; // per jes range
    std::vector<JetCorrectionUncertainty*> jetCorrectorUncSources_; // per source, one for all ranges
    JME::JetResolution* jerResolution_;
    JME::JetResolutionScaleFactor* jerScaleFactor_;
    TRandom3* rnd_;
};

TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)
    : verbose_(iConfig.getUntrackedParameter<bool>("verbose"))
    , outputFile_(iConfig.getParameter<string>("outputFile"))
    , campaign_(iConfig.getParameter<string>("campaign"))
    , metaDataFile_(iConfig.getParameter<string>("metaDataFile"))
    , isData_(iConfig.getParameter<bool>("isData"))
    , leptonChannel_(iConfig.getParameter<string>("leptonChannel"))
    , eeTriggers_(iConfig.getParameter<vstring>("eeTriggers"))
    , emuTriggers_(iConfig.getParameter<vstring>("emuTriggers"))
    , mumuTriggers_(iConfig.getParameter<vstring>("mumuTriggers"))
    //, eTriggers_(iConfig.getParameter<vstring>("eTriggers"))
    //, muTriggers_(iConfig.getParameter<vstring>("muTriggers"))
    , metFilters_(iConfig.getParameter<vstring>("metFilters"))
    , jesFiles_(iConfig.getParameter<vstring>("jesFiles"))
    , jesRanges_(iConfig.getParameter<vint>("jesRanges"))
    , jesUncFiles_(iConfig.getParameter<vstring>("jesUncFiles"))
    , jesUncSrcFile_(iConfig.getParameter<string>("jesUncSrcFile"))
    , jesUncSources_(iConfig.getParameter<vstring>("jesUncSources"))
    , jerPtResolutionFile_(iConfig.getParameter<string>("jerPtResolutionFile"))
    , jerScaleFactorFile_(iConfig.getParameter<string>("jerScaleFactorFile"))
    , genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoCollection")))
    , triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBitsCollection")))
    , metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBitsCollection")))
    , pileupInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupInfoCollection")))
    , beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotCollection")))
    , vertexToken_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexCollection")))
    , electronToken_(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection")))
    , muonToken_(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonCollection")))
    , metToken_(consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metCollection")))
    , jetToken_(consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetCollection")))
    , genJetToken_(consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetCollection")))
    , rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection")))
    , tfile_(nullptr)
    , tfileMeta_(nullptr)
    , tree_(nullptr)
    , eventHist_(nullptr)
    , selectedEventHist_(nullptr)
    , weightHist_(nullptr)
    , selectedWeightHist_(nullptr)
    , cutflowHist_(nullptr)
    , pileupHist_(nullptr)
    , jerResolution_(nullptr)
    , jerScaleFactor_(nullptr)
    , rnd_(0)
{
    if (verbose_)
    {
        std::cout << "running TreeMaker in verbose mode" << std::endl;
    }

    setupJetCorrectionObjects();
    setupVariables();

    rnd_ = new TRandom3(0);
}

TreeMaker::~TreeMaker()
{
}

void TreeMaker::setupJetCorrectionObjects()
{
    // build jet variations, i.e. a vector a string pairs (variation, direction)
    // begin with the nominal one
    jetVariations_.push_back(stringPair("", ""));

    if (!isData_ && !jerPtResolutionFile_.empty() && !jerScaleFactorFile_.empty())
    {
        jerResolution_ = new JME::JetResolution(jerPtResolutionFile_);
        jerScaleFactor_ = new JME::JetResolutionScaleFactor(jerScaleFactorFile_);

        // add jet variations
        jetVariations_.push_back(stringPair("jer", "up"));
        jetVariations_.push_back(stringPair("jer", "down"));
    }

    // stop here, when no jes files are given
    if (jesFiles_.size() == 0)
    {
        return;
    }

    // sanity checks
    if (jesRanges_.size() % 2 != 0)
    {
        throw std::runtime_error("please provide an even number of jesRanges");
    }
    nJESRanges_ = jesRanges_.size() / 2;
    if (jesFiles_.size() % nJESRanges_ != 0)
    {
        throw std::runtime_error("the number of JESFiles does not match the number of JESRanges");
    }
    nJESFilesPerRange_ = jesFiles_.size() / nJESRanges_;
    if (jesUncFiles_.size() > 0 && jesUncFiles_.size() != nJESRanges_)
    {
        throw std::runtime_error("the number of JESUncFiles does not match the number of JESRanges");
    }

    // initialize the jet correctors
    // example: when there are 8 jes files and 2 jes ranges, each JETCorrector has 4 files which
    // should correspond to the number of correction levels applied
    for (size_t i = 0; i < nJESRanges_; i++)
    {
        std::vector<JetCorrectorParameters> corrParams;
        for (size_t j = 0; j < nJESFilesPerRange_; j++)
        {
            JetCorrectorParameters params(jesFiles_[i * nJESFilesPerRange_ + j]);
            corrParams.push_back(params);
        }
        jetCorrectors_.push_back(new FactorizedJetCorrector(corrParams));
    }

    // also process uncertaintes, but - of course - only for MC
    if (!isData_)
    {
        // initialize the total jet corrector uncertainties
        if (jesUncFiles_.size() > 0)
        {
            for (size_t i = 0; i < jesUncFiles_.size(); i++)
            {
                JetCorrectorParameters params(jesUncFiles_[i]);
                jetCorrectorUncs_.push_back(new JetCorrectionUncertainty(params));
            }
        }

        // initialize the factorized jet corrector uncertainties
        for (size_t i = 0; i < jesUncSources_.size(); ++i)
        {
            if (!jesUncSrcFile_.empty())
            {
                JetCorrectorParameters params(jesUncSrcFile_, jesUncSources_[i]);
                jetCorrectorUncSources_.push_back(new JetCorrectionUncertainty(params));
            }
            else
            {
                jetCorrectorUncSources_.push_back(0);
            }
        }
    }

    // add jet variations for each jes source
    for (size_t i = 0; i < jesUncSources_.size(); i++)
    {
        jetVariations_.push_back(stringPair("jes" + jesUncSources_[i], "up"));
        jetVariations_.push_back(stringPair("jes" + jesUncSources_[i], "down"));
    }
}

void TreeMaker::setupVariables()
{
    // event variables
    varMap_.addInt32("is_data");
    varMap_.addInt64("event");
    varMap_.addInt32("run");
    varMap_.addInt32("lumi");
    varMap_.addDouble("gen_weight");
    varMap_.addFloat("pu");
    varMap_.addDouble("rho");
    varMap_.addInt32("channel");

    // leptons
    for (size_t i = 1; i <= 2; i++)
    {
        varMap_.addDouble("lep" + std::to_string(i) + "_E");
        varMap_.addDouble("lep" + std::to_string(i) + "_px");
        varMap_.addDouble("lep" + std::to_string(i) + "_py");
        varMap_.addDouble("lep" + std::to_string(i) + "_pz");
        varMap_.addDouble("lep" + std::to_string(i) + "_eta");
        varMap_.addInt32("lep" + std::to_string(i) + "_charge");
        varMap_.addInt32("lep" + std::to_string(i) + "_pdg");
        varMap_.addFloat("lep" + std::to_string(i) + "_iso");
        varMap_.addInt32("lep" + std::to_string(i) + "_tight");
        varMap_.addDouble("lep" + std::to_string(i) + "_eta_sc");
    }
    varMap_.addDouble("mll");
    varMap_.addDouble("dr_ll");

    // jet and MET and other JES/JER dependent variables
    for (size_t i = 0; i < (isData_ ? 1 : jetVariations_.size()); i++)
    {
        string postfix = "";
        if (!jetVariations_[i].first.empty())
        {
            postfix += "_" + jetVariations_[i].first + "_" + jetVariations_[i].second;
        }

        varMap_.addInt32("jetmet_pass" + postfix);
        varMap_.addInt32("n_jets" + postfix);
        varMap_.addInt32("pass_z_mask" + postfix);
        varMap_.addDouble("met_px" + postfix);
        varMap_.addDouble("met_py" + postfix);
        varMap_.addDouble("mht" + postfix);

        // jets
        for (size_t j = 1; j <= 4; j++)
        {
            varMap_.addDouble("jet" + std::to_string(j) + "_E" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_px" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_py" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_pz" + postfix);

            varMap_.addDouble("jet" + std::to_string(j) + "_eta" + postfix);
            varMap_.addInt32("jet" + std::to_string(j) + "_flavor" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_csvv2" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepcsv_b" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepcsv_bb" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepcsv_c" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepcsv_udsg" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_b" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_bb" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_lepb" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_c" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_uds" + postfix);
            varMap_.addDouble("jet" + std::to_string(j) + "_deepjet_udsg" + postfix);
        }
    }
}

void TreeMaker::beginJob()
{
    // create the tree
    tfile_ = new TFile(outputFile_.c_str(), "RECREATE");
    tfile_->cd();
    tree_ = new TTree("tree", "tree");

    // create the counting histograms, 2 bins each to separate information for events with positive
    // or negative weights
    tfileMeta_ = new TFile(metaDataFile_.c_str(), "RECREATE");
    tfileMeta_->cd();
    eventHist_ = new TH1F("events", "", 2, -1., 1.);
    selectedEventHist_ = new TH1F("selected_events", "", 2, -1., 1.);
    weightHist_ = new TH1F("event_weights", "", 2, -1., 1.);
    selectedWeightHist_ = new TH1F("selected_event_weights", "", 2, -1., 1.);
    cutflowHist_ = new TH1F("cutflow", "", 6, 0., 6.);
    pileupHist_ = new TH1F("pileup", "", 100, 0., 100.);

    // add branches based on added variables
    for (size_t i = 0; i < varMap_.size(); i++)
    {
        string name = varMap_.getName(i);
        string typeFlag = varMap_.getFlag(name);

        if (verbose_)
        {
            std::cout << "adding variable: " << name << std::endl;
        }

        if (typeFlag == "F")
        {
            tree_->Branch(name.c_str(), &varMap_.getFloat(name), (name + "/F").c_str());
        }
        else if (typeFlag == "D")
        {
            tree_->Branch(name.c_str(), &varMap_.getDouble(name), (name + "/D").c_str());
        }
        else if (typeFlag == "I")
        {
            tree_->Branch(name.c_str(), &varMap_.getInt32(name), (name + "/I").c_str());
        }
        else if (typeFlag == "L")
        {
            tree_->Branch(name.c_str(), &varMap_.getInt64(name), (name + "/L").c_str());
        }
    }
    std::cout << "total variables: " << varMap_.size() << std::endl;

    // set campaign specific information
    if (campaign_ == "2018_Run2_pp_13TeV_MORIOND19legacy")
    {
        tightJetID_ = &TreeMaker::tightJetID_2016;
        maxJetEta_ = 2.4;
    }
    else if (campaign_ == "2017_Run2_pp_13TeV_ICHEP18")
    {
        tightJetID_ = &TreeMaker::tightJetID_2017;
        maxJetEta_ = 2.5;
    }
    else if (campaign_ == "2018_Run2_pp_13TeV_MORIOND19")
    {
        tightJetID_ = &TreeMaker::tightJetID_2018;
        maxJetEta_ = 2.5;
    }
    else
    {
        throw std::runtime_error("Unknown campaign " + campaign_);
    }
}

void TreeMaker::endJob()
{
    // log some stats
    std::cout << "total events        : " << eventHist_->Integral() << std::endl;
    std::cout << "selected events     : " << selectedEventHist_->Integral() << std::endl;
    std::cout << "sum weights         : " << weightHist_->Integral() << std::endl;
    std::cout << "sum selected weights: " << selectedWeightHist_->Integral() << std::endl;

    tfile_->cd();
    tree_->Write();
    tfile_->Close();

    tfileMeta_->cd();
    eventHist_->Write();
    selectedEventHist_->Write();
    weightHist_->Write();
    selectedWeightHist_->Write();
    cutflowHist_->Write();
    pileupHist_->Write();
    tfileMeta_->Close();
}

void TreeMaker::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
    varMap_.reset();

    // start cutflow
    double cutflowBin = 0.5;
    cutflowHist_->Fill(cutflowBin++);

    // read the generator generator infos
    double genWeight = isData_ ? 1. : readGenWeight(event);

    // fill hists
    double histPos = genWeight < 0 ? -0.5 : 0.5;
    eventHist_->Fill(histPos, 1.);
    weightHist_->Fill(histPos, genWeight);

    // get vertex
    edm::Handle<std::vector<reco::Vertex> > verticesHandle;
    event.getByToken(vertexToken_, verticesHandle);
    std::vector<reco::Vertex> vertices(*verticesHandle);
    reco::Vertex vertex = vertices.front();

    // read pu infos
    float pu = isData_ ? vertices.size() : readPU(event);

    // fill pileup hist
    pileupHist_->Fill(pu, genWeight);

    // read the rho value
    double rho = readRho(event);

    // check MET filters
    if (!metFilterSelection(event))
    {
        return;
    }
    cutflowHist_->Fill(cutflowBin++);

    // select the vertex

    if (vertexID(vertex) != V_VALID)
    {
        return;
    }
    cutflowHist_->Fill(cutflowBin++);

    // read and select electrons
    std::vector<pat::Electron> electrons;
    std::vector<pat::Electron> tightElectrons;
    if (!electronSelection(event, vertex, rho, electrons, tightElectrons))
    {
        return;
    }

    // read and select muons
    std::vector<pat::Muon> muons;
    std::vector<pat::Muon> tightMuons;
    if (!muonSelection(event, vertex, muons, tightMuons))
    {
        return;
    }

    // combined lepton selection
    reco::RecoCandidate* lep1 = nullptr;
    reco::RecoCandidate* lep2 = nullptr;
    LeptonChannel channel = C_INVALID;
    bool is_sl = false;
    bool passLeptonSelection = leptonSelection(
        electrons, tightElectrons, muons, tightMuons, lep1, lep2, channel, is_sl);
    if (!passLeptonSelection)
    {
        return;
    }
    cutflowHist_->Fill(cutflowBin++);

    // trigger selection
    if (!triggerSelection(event, channel))
    {
        return;
    }
    cutflowHist_->Fill(cutflowBin++);

    // read the MET
    edm::Handle<std::vector<pat::MET> > metHandle;
    event.getByToken(metToken_, metHandle);
    pat::MET metOrig(metHandle->at(0));
    metOrig.setP4(metOrig.shiftedP4(pat::MET::NoShift, pat::MET::Type1XY));

    // read and select jets plus MET
    std::vector<std::vector<pat::Jet> > jets;
    std::vector<pat::MET> mets;
    std::vector<bool> passJetMETSelection;
    bool passOneJetMETSelection = false;
    size_t rndSeed = rnd_->GetSeed();
    for (size_t i = 0; i < (isData_ ? 1 : jetVariations_.size()); i++)
    {
        string variation = jetVariations_[i].first;
        string direction = jetVariations_[i].second;
        std::vector<pat::Jet> jets2;
        pat::MET met(metOrig);

        rnd_->SetSeed(rndSeed);

        bool pass = jetMETSelection(
            event, rho, lep1, lep2, metOrig, variation, direction, jets2, met, is_sl);
        jets.push_back(jets2);
        mets.push_back(met);
        passJetMETSelection.push_back(pass);

        passOneJetMETSelection |= pass;
    }
    if (!passOneJetMETSelection)
    {
        return;
    }
    // count only the nominal jet met decision
    if (passJetMETSelection[0])
    {
        cutflowHist_->Fill(cutflowBin++);
    }

    // fill hists
    selectedEventHist_->Fill(histPos, 1.);
    selectedWeightHist_->Fill(histPos, genWeight);

    // event variables
    varMap_.setInt32("is_data", isData_);
    varMap_.setInt64("event", event.id().event());
    varMap_.setInt32("run", event.id().run());
    varMap_.setInt32("lumi", event.luminosityBlock());
    varMap_.setDouble("gen_weight", genWeight);
    varMap_.setFloat("pu", pu);
    varMap_.setDouble("rho", rho);
    varMap_.setInt32("channel", int32_t(channel));

    // lepton variables
    size_t n_leps = is_sl ? 1 : 2;
    for (size_t i = 1; i <= n_leps; i++)
    {
        reco::RecoCandidate* lep = i == 1 ? lep1 : lep2;
        varMap_.setDouble("lep" + std::to_string(i) + "_E", lep->energy());
        varMap_.setDouble("lep" + std::to_string(i) + "_px", lep->px());
        varMap_.setDouble("lep" + std::to_string(i) + "_py", lep->py());
        varMap_.setDouble("lep" + std::to_string(i) + "_pz", lep->pz());
        varMap_.setDouble("lep" + std::to_string(i) + "_eta", lep->eta());
        varMap_.setInt32("lep" + std::to_string(i) + "_charge", lep->charge());
        varMap_.setInt32("lep" + std::to_string(i) + "_pdg", lep->pdgId());
        double absPdgId = abs(lep->pdgId());
        if (absPdgId == 11)
        {
            pat::Electron* e = dynamic_cast<pat::Electron*>(lep);
            varMap_.setDouble("lep" + std::to_string(i) + "_eta_sc", e->superCluster()->eta());
            varMap_.setFloat("lep" + std::to_string(i) + "_iso", e->userFloat("iso"));
            varMap_.setInt32("lep" + std::to_string(i) + "_tight", e->userInt("tight"));
        }
        else if (absPdgId == 13)
        {
            pat::Muon* mu = dynamic_cast<pat::Muon*>(lep);
            varMap_.setDouble("lep" + std::to_string(i) + "_eta_sc", mu->eta());
            varMap_.setFloat("lep" + std::to_string(i) + "_iso", mu->userFloat("iso"));
            varMap_.setInt32("lep" + std::to_string(i) + "_tight", mu->userInt("tight"));
        }
        else
        {
            throw std::runtime_error("cannot handle lepton pdg id " + std::to_string(absPdgId));
        }
    }
    double mll = 0.;
    if (!is_sl)
    {
        mll = (lep1->p4() + lep2->p4()).M();
        double dr_ll = deltaR(lep1->p4(), lep2->p4());
        varMap_.setDouble("mll", mll);
        varMap_.setDouble("dr_ll", dr_ll);
    }

    // jet and MET variables
    for (size_t i = 0; i < (isData_ ? 1 : jetVariations_.size()); i++)
    {
        string postfix = "";
        if (!jetVariations_[i].first.empty())
        {
            postfix += "_" + jetVariations_[i].first + "_" + jetVariations_[i].second;
        }

        varMap_.setInt32("jetmet_pass" + postfix, passJetMETSelection[i]);
        varMap_.setInt32("n_jets" + postfix, jets[i].size());
        varMap_.setDouble("met_px" + postfix, mets[i].px());
        varMap_.setDouble("met_py" + postfix, mets[i].py());

        // jets
        double mht_x = 0.;
        double mht_y = 0.;
        if (!is_sl)
        {
            mht_x += lep1->px() + lep2->px();
            mht_y += lep1->py() + lep2->py();
            for (size_t j = 1; j <= jets[i].size(); j++)
            {
                mht_x += jets[i][j - 1].px();
                mht_y += jets[i][j - 1].py();
            }

            double mht = sqrt(mht_x * mht_x + mht_y * mht_y);
            varMap_.setDouble("mht" + postfix, mht);

            bool pass_z_mask = (mll < (65.5 + 3 * mht / 8)) ||
                               (mll > (108. - mht / 4)) ||
                               (mll < (79. - 3 * mht / 4)) ||
                               (mll > (99. + mht / 2));
            varMap_.setInt32("pass_z_mask" + postfix, (int)pass_z_mask);
        }

        for (size_t j = 1; j <= 4; j++)
        {
            if (jets[i].size() < j)
            {
                break;
            }

            pat::Jet* jet = &jets[i][j - 1];
            varMap_.setDouble("jet" + std::to_string(j) + "_E" + postfix, jet->energy());
            varMap_.setDouble("jet" + std::to_string(j) + "_px" + postfix, jet->px());
            varMap_.setDouble("jet" + std::to_string(j) + "_py" + postfix, jet->py());
            varMap_.setDouble("jet" + std::to_string(j) + "_pz" + postfix, jet->pz());
            varMap_.setDouble("jet" + std::to_string(j) + "_eta" + postfix, jet->eta());
            varMap_.setInt32("jet" + std::to_string(j) + "_flavor" + postfix,
                jet->hadronFlavour());
            varMap_.setDouble("jet" + std::to_string(j) + "_csvv2" + postfix,
                jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepcsv_b" + postfix,
                jet->bDiscriminator("pfDeepCSVJetTags:probb"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepcsv_bb" + postfix,
                jet->bDiscriminator("pfDeepCSVJetTags:probbb"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepcsv_c" + postfix,
                jet->bDiscriminator("pfDeepCSVJetTags:probc"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepcsv_udsg" + postfix,
                jet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_b" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:probb"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_bb" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_lepb" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_c" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:probc"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_uds" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
            varMap_.setDouble("jet" + std::to_string(j) + "_deepjet_udsg" + postfix,
                jet->bDiscriminator("pfDeepFlavourJetTags:probg"));
        }
    }

    // finally, fill the tree
    tree_->Fill();
}

bool TreeMaker::metFilterSelection(const edm::Event& event)
{
    // check MET filters
    // unlike triggers all of them have to be accepted
    edm::Handle<edm::TriggerResults> metFilterBitsHandle;
    event.getByToken(metFilterBitsToken_, metFilterBitsHandle);

    const edm::TriggerNames& metFilterNames = event.triggerNames(*metFilterBitsHandle);
    for (size_t i = 0; i < metFilterBitsHandle->size(); i++)
    {
        string name = metFilterNames.triggerName(i);
        if (std::find(metFilters_.begin(), metFilters_.end(), name) != metFilters_.end())
        {
            if (!metFilterBitsHandle->accept(i))
            {
                return false;
            }
        }
    }

    return true;
}

bool TreeMaker::triggerSelection(const edm::Event& event, LeptonChannel& channel)
{
    edm::Handle<edm::TriggerResults> triggerBitsHandle;
    event.getByToken(triggerBitsToken_, triggerBitsHandle);
    if (!triggerBitsHandle.isValid())
    {
        return false;
    }

    std::vector<string> triggersToPass;
    if (channel == C_EE)
    {
        triggersToPass = eeTriggers_;
    }
    else if (channel == C_EMU)
    {
        triggersToPass = emuTriggers_;
    }
    else if (channel == C_MUMU)
    {
        triggersToPass = mumuTriggers_;
    }
    //else if (channel == C_E)
    //{
    //    triggersToPass = eTriggers_;
    //}
    //else if (channel == C_MU)
    //{
    //    triggersToPass = muTriggers_;
    //}

    const edm::TriggerNames& triggerNames = event.triggerNames(*triggerBitsHandle);
    for (size_t i = 0; i < triggerBitsHandle->size(); ++i)
    {
        if (!triggerBitsHandle->accept(i))
        {
            continue;
        }
        string name = triggerNames.triggerName(i);
        for (const string& triggerToPass : triggersToPass)
        {
            if (fnmatch(triggerToPass.c_str(), name.c_str(), 0) == 0)
            {
                return true;
            }
        }
    }

    return false;
}

bool TreeMaker::electronSelection(const edm::Event& event, reco::Vertex& vertex, double rho,
    std::vector<pat::Electron>& electrons, std::vector<pat::Electron>& tightElectrons)
{
    edm::Handle<std::vector<pat::Electron> > electronsHandle;
    event.getByToken(electronToken_, electronsHandle);

    std::vector<pat::Electron> electronsCopy(*electronsHandle);
    for (pat::Electron& electron : electronsCopy)
    {
        ElectronID type = electronID(electron, vertex, rho);
        if (type >= E_LOOSE)
        {
            electrons.push_back(electron);
            if (type == E_TIGHT)
            {
                tightElectrons.push_back(electron);
            }
        }
    }

    return true;
}

bool TreeMaker::muonSelection(const edm::Event& event, reco::Vertex& vertex,
    std::vector<pat::Muon>& muons, std::vector<pat::Muon>& tightMuons)
{
    edm::Handle<std::vector<pat::Muon> > muonsHandle;
    event.getByToken(muonToken_, muonsHandle);

    std::vector<pat::Muon> muonsCopy(*muonsHandle);
    for (pat::Muon& muon : muonsCopy)
    {
        MuonID type = muonID(muon, vertex);
        if (type >= M_LOOSE)
        {
            muons.push_back(muon);
            if (type == M_TIGHT)
            {
                tightMuons.push_back(muon);
            }
        }
    }

    return true;
}

bool TreeMaker::leptonSelection(std::vector<pat::Electron>& electrons,
    std::vector<pat::Electron>& tightElectrons, std::vector<pat::Muon>& muons,
    std::vector<pat::Muon>& tightMuons, reco::RecoCandidate*& lep1,
    reco::RecoCandidate*& lep2, LeptonChannel& channel, bool& is_sl)
{
    // require exactly 2 leptons, at least 1 tight, or exactly 1 tight lepton for the SL control region
    size_t nElectrons = electrons.size();
    size_t nTightElectrons = tightElectrons.size();

    size_t nMuons = muons.size();
    size_t nTightMuons = tightMuons.size();

    size_t nLeptons = nElectrons + nMuons;
    size_t nTightLeptons = nTightElectrons + nTightMuons;

    if (nLeptons != 2 || nTightLeptons < 1)
    {
        return false;
    }

    // extract leptons into RecoCandidate pointers and define the lepton channel
    if (nElectrons == 2 && nMuons == 0)
    {
        channel = C_EE;
        lep1 = &electrons[0];
        lep2 = &electrons[1];
    }
    else if (nElectrons == 1 && nMuons == 1)
    {
        channel = C_EMU;
        lep1 = &electrons[0];
        lep2 = &muons[0];
    }
    else if (nElectrons == 0 && nMuons == 2)
    {
        channel = C_MUMU;
        lep1 = &muons[0];
        lep2 = &muons[1];
    }
    //else if (nElectrons == 1 && nMuons == 0)
    //{
    //    is_sl = true;
    //    channel = C_E;
    //    lep1 = &electrons[0];
    //}
    //else if (nElectrons == 0 && nMuons == 1)
    //{
    //    is_sl = true;
    //    channel = C_MU;
    //    lep1 = &muons[0];
    //}
    else
    {
        throw std::runtime_error("could not determine lepton channel!");
    }

    // check for opposite charges (only DL channels)
    if (!is_sl && (lep1->charge() == lep2->charge()))
    {
        return false;
    }

    // a particular lepton channel might be required (by name), e.g. for data
    if (!leptonChannel_.empty())
    {
        if ((leptonChannel_ == "ee" && channel != C_EE) ||
            (leptonChannel_ == "emu" && channel != C_EMU) ||
            (leptonChannel_ == "mumu" && channel != C_MUMU))
            //(leptonChannel_ == "e" && channel != C_E) ||
            //(leptonChannel_ == "mu" && channel != C_MU))
        {
            return false;
        }
    }

    // re-order py pt (only DL channels)
    if (!is_sl && (lep2->pt() > lep1->pt()))
    {
        std::swap(lep1, lep2);
    }

    return true;
}

bool TreeMaker::jetMETSelection(const edm::Event& event, double rho,
    reco::RecoCandidate* lep1, reco::RecoCandidate* lep2, const pat::MET& metOrig,
    const string& variation, const string& direction, std::vector<pat::Jet>& jets, pat::MET& met,
    bool is_sl)
{
    // read jets
    edm::Handle<std::vector<pat::Jet> > jetsHandle;
    event.getByToken(jetToken_, jetsHandle);

    // read gen jets for matching within JER smearing
    const std::vector<reco::GenJet>* genJets = nullptr;
    if (!isData_)
    {
        edm::Handle<std::vector<reco::GenJet> > genJetsHandle;
        event.getByToken(genJetToken_, genJetsHandle);
        genJets = genJetsHandle.product();
    }

    // keep track of changes to MET due to jet corrections
    LorentzVector correctedMetP4(metOrig.p4());

    for (const pat::Jet& jetOrig : *jetsHandle)
    {
        pat::Jet jet(jetOrig);

        // make sure that jetID is applied at correct point because functionality has changed in the past
        double totalEnergyFraction = jet.neutralHadronEnergyFraction() + jet.neutralEmEnergyFraction()
            + jet.chargedHadronEnergyFraction() + jet.chargedEmEnergyFraction() + jet.muonEnergyFraction();
        if (std::fabs(totalEnergyFraction - 1. > 0.0001))
        {
            throw std::runtime_error("Jet energy fractions do not add up to 1.");
        }

        bool passPOGID;
        passPOGID = (this->*tightJetID_)(jet);

        // when variation is empty, apply nominal JES and JER,
        // when variation is JER, apply nominal JES and JER variation,
        // otherwise, apply JES variation and nominal JER
        if (variation.empty())
        {
            applyJES(jet, "", "", event.run(), rho);
            applyJER(jet, genJets, "", "", rho);
        }
        else if (variation == "jer")
        {
            applyJES(jet, "", "", event.run(), rho);
            applyJER(jet, genJets, variation, direction, rho);
        }
        else
        {
            applyJES(jet, variation, direction, event.run(), rho);
            applyJER(jet, genJets, "", "", rho);
        }

        JetID type = jetID(jet, lep1, lep2, is_sl, passPOGID);
        if (type == J_VALID)
        {
            jets.push_back(jet);

            // propagate to met
            correctedMetP4.SetPx(correctedMetP4.Px() - (jet.px() - jetOrig.px()));
            correctedMetP4.SetPy(correctedMetP4.Py() - (jet.py() - jetOrig.py()));
        }
    }

    // propagate the change in px and py to met
    met.setP4(correctedMetP4);

    size_t min_njets = is_sl ? 4 : 2;
    if (jets.size() < min_njets)
    {
        return false;
    }

    // sort jets
    std::sort(jets.begin(), jets.end(), comparePt);

    return true;
}

VertexID TreeMaker::vertexID(reco::Vertex& vertex)
{
    bool isValid = !vertex.isFake() &&
        vertex.ndof() >= 4.0 &&
        abs(vertex.z()) <= 24.0 &&
        abs(vertex.position().Rho()) <= 2.0;
    return isValid ? V_VALID : V_INVALID;
}

ElectronID TreeMaker::electronID(pat::Electron& electron, reco::Vertex& vertex, double rho)
{
    // loose pt cut
    if (electron.pt() <= 15.)
    {
        return E_INVALID;
    }

    // eta cut
    double absEta = fabs(electron.eta());
    if (absEta >= 2.4)
    {
        return E_INVALID;
    }

    // barrel-endcap crack filter
    double absSCEta = fabs(electron.superCluster()->eta());
    if (absSCEta >= 1.4442 && absSCEta <= 1.5660)
    {
        return E_INVALID;
    }

    // IP cuts
    const auto& gsfTrackRef = electron.gsfTrack();
    if (gsfTrackRef.isNull())
    {
        return E_INVALID;
    }
    bool isBarrel = absSCEta <= 1.479;
    if (fabs(gsfTrackRef.get()->dz(vertex.position())) >= (isBarrel ? 0.1 : 0.2))
    {
        return E_INVALID;
    }
    if (fabs(gsfTrackRef.get()->dxy(vertex.position())) >= (isBarrel ? 0.05 : 0.1))
    {
        return E_INVALID;
    }

    // electron VID
    if (electron.electronID("cutBasedElectronID-Fall17-94X-V2-tight") != 1.)
    {
        return E_INVALID;
    }

    // calculate the relative isolation for later use
    // no need to cut on since VID already did
    const auto& isoVars = electron.pfIsolationVariables();
    float effArea = electronEffectiveArea(electron);
    float neutralSum = isoVars.sumNeutralHadronEt + isoVars.sumPhotonEt - rho * effArea;
    float iso = (isoVars.sumChargedHadronPt + fmax(0.0, neutralSum)) / electron.pt();
    electron.addUserFloat("iso", iso, true);

    // tight or loose decision only depends on pt
    bool isTight = electron.pt() > 25;
    electron.addUserInt("tight", isTight, true);

    return isTight ? E_TIGHT : E_LOOSE;
}

MuonID TreeMaker::muonID(pat::Muon& muon, reco::Vertex& vertex)
{
    // loose pt cut
    if (muon.pt() <= 15.)
    {
        return M_INVALID;
    }

    // eta cut
    double absEta = fabs(muon.eta());
    if (absEta >= 2.4)
    {
        return M_INVALID;
    }

    // isolation cut
    const auto& isoVars = muon.pfIsolationR04();
    float neutralSum = isoVars.sumNeutralHadronEt + isoVars.sumPhotonEt - 0.5 * isoVars.sumPUPt;
    float iso = (isoVars.sumChargedHadronPt + fmax(0.0, neutralSum)) / muon.pt();
    muon.addUserFloat("iso", iso, true);
    if (iso >= 0.15)
    {
        return M_INVALID;
    }

    // tight ID criteria from muon POG
    const auto& globalTrackRef = muon.globalTrack();
    const auto& bestTrackRef = muon.muonBestTrack();
    const auto& innerTrackRef = muon.innerTrack();
    bool passPOGID = !globalTrackRef.isNull() &&
        !bestTrackRef.isNull() &&
        muon.isGlobalMuon() &&
        muon.isPFMuon() &&
        globalTrackRef.get()->normalizedChi2() < 10. &&
        globalTrackRef.get()->hitPattern().numberOfValidMuonHits() > 0 &&
        muon.numberOfMatchedStations() > 1 &&
        std::min(fabs(bestTrackRef.get()->dxy(vertex.position())), fabs(muon.dB())) < 0.2 &&
        fabs(bestTrackRef.get()->dz(vertex.position())) < 0.5 &&
        innerTrackRef.get()->hitPattern().numberOfValidPixelHits() > 0 &&
        innerTrackRef.get()->hitPattern().trackerLayersWithMeasurement() > 5;
    if (!passPOGID)
    {
        return M_INVALID;
    }

    // tight or loose decision depends on pt
    bool isTight = muon.pt() > 25.;
    muon.addUserInt("tight", isTight, true);

    return isTight ? M_TIGHT : M_LOOSE;
}

JetID TreeMaker::jetID(pat::Jet& jet, reco::RecoCandidate* lep1, reco::RecoCandidate* lep2,
    bool is_sl, bool passPOGID)
{
    // pt cut
    if (jet.pt() <= 20.)
    {
        return J_INVALID;
    }

    // eta cut
    double absEta = fabs(jet.eta());
    if (absEta >= maxJetEta_)
    {
        return J_INVALID;
    }

    // loose PU jet ID
    if (jet.userInt("pileupJetId:fullId") < 4)
    {
        return J_INVALID;
    }

    // check distance to selected leptons
    if (deltaR(jet.p4(), lep1->p4()) < 0.4)
    {
        return J_INVALID;
    }
    if (!is_sl && (deltaR(jet.p4(), lep2->p4()) < 0.4))
    {
        return J_INVALID;
    }

    if (!passPOGID)
    {
        return J_INVALID;
    }

    return J_VALID;
}

double TreeMaker::readGenWeight(const edm::Event& event)
{
    edm::Handle<GenEventInfoProduct> genInfoHandle;
    event.getByToken(genInfoToken_, genInfoHandle);
    return genInfoHandle->weight();
}

float TreeMaker::readPU(const edm::Event& event)
{
    edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoHandle;
    event.getByToken(pileupInfoToken_, pileupInfoHandle);

    for (const PileupSummaryInfo& info : *pileupInfoHandle)
    {
        if (info.getBunchCrossing() == 0)
        {
            return info.getTrueNumInteractions();
        }
    }
    throw std::runtime_error("no valid bunch crossing found to read PU infos");
}

double TreeMaker::readRho(const edm::Event& event)
{
    edm::Handle<double> rhoHandle;
    event.getByToken(rhoToken_, rhoHandle);
    return *rhoHandle;
}

void TreeMaker::applyJES(pat::Jet& jet, const string& variation, const string& direction,
    int64_t run, double rho)
{
    // before we start, lookup the proper objects, as all values depend on the run number
    FactorizedJetCorrector* jetCorrector = nullptr;
    JetCorrectionUncertainty* jetCorrectorUnc = nullptr;
    for (size_t i = 0; i < nJESRanges_; i++)
    {
        if (run >= jesRanges_[2 * i] && run <= jesRanges_[2 * i + 1])
        {
            jetCorrector = jetCorrectors_[i];
            if (jetCorrectorUncs_.size() > 0)
            {
                jetCorrectorUnc = jetCorrectorUncs_[i];
            }
        }
    }
    if (!jetCorrector)
    {
        throw std::runtime_error("no jetCorrector found for run " + std::to_string(run));
    }

    // first, uncorrect the jet
    // note: uncorrectFactor would actually lead to the original p4, i.e., before applying any JES
    // so "p4 * uncorrectFactor = p4_orig"
    // or "p4 = uncorrectFactor / p4_orig"
    double uncorrectFactor = jet.jecFactor("Uncorrected");

    jet.setP4(jet.p4() * uncorrectFactor);

    // then, calculate the new correction factor
    // so "p4_new = p4 * recorrectFactor"
    // rho is from fixedGridRhoFastJetAll
    jetCorrector->setJetPt(jet.pt());
    jetCorrector->setJetEta(jet.eta());
    jetCorrector->setJetA(jet.jetArea());
    jetCorrector->setRho(rho);
    double recorrectFactor = jetCorrector->getCorrection();

    // apply a variation / uncertainty?
    if (!variation.empty())
    {
        // variation is "jes<source>"
        string source = variation.substr(3);
        JetCorrectionUncertainty* corr = nullptr;
        if (source == "total")
        {
            // total variation
            corr = jetCorrectorUnc;
        }
        else
        {
            // factorized variation
            for (size_t i = 0; i < jesUncSources_.size(); i++)
            {
                if (jesUncSources_[i] == source)
                {
                    corr = jetCorrectorUncSources_[i];
                    break;
                }
            }
        }
        if (!corr)
        {
            throw std::runtime_error("could not find jet corrector for source " + source);
        }
        // update the recorrectFactor accordingly (use corrected jet)
        corr->setJetPt(jet.pt() * recorrectFactor);
        corr->setJetEta(jet.eta());
        if (direction == "up")
        {
            recorrectFactor *= 1. + corr->getUncertainty(true);
        }
        else
        {
            recorrectFactor *= 1. - corr->getUncertainty(false);
        }
    }

    // finally, scale the four-vector
    jet.setP4(jet.p4() * recorrectFactor);
}

void TreeMaker::applyJER(pat::Jet& jet, const std::vector<reco::GenJet>* genJets,
    const string& variation, const string& direction, double rho)
{
    if (isData_ || jerResolution_ == nullptr || jerScaleFactor_ == nullptr || jet.pt() == 0)
    {
        return;
    }

    Variation v = Variation::NOMINAL;
    if (variation == "jer")
    {
        v = direction == "up" ? Variation::UP : Variation::DOWN;
    }

    double res = jerResolution_->getResolution({ { JME::Binning::JetPt, jet.pt() },
        { JME::Binning::JetEta, jet.eta() }, { JME::Binning::Rho, rho } });
    double sf = jerScaleFactor_->getScaleFactor({ { JME::Binning::JetEta, jet.eta() } }, v);

    // try to find a matched gen jet
    double minDR = 0.2;
    double maxDPt = 3 * res * jet.pt();
    const reco::GenJet* matchedGenJet = nullptr;
    for (size_t i = 0; i < genJets->size(); i++)
    {
        double dR = deltaR(genJets->at(i), jet);
        if (dR >= minDR)
        {
            continue;
        }

        double dPt = std::fabs(genJets->at(i).pt() - jet.pt());
        if (dPt > maxDPt)
        {
            continue;
        }

        minDR = dR;
        matchedGenJet = &(genJets->at(i));
    }

    // scaling when a gen jet was found, random smearing otherwise (and sf > 1)
    double jerFactor = 1.;
    if (matchedGenJet)
    {
        jerFactor = 1 + (sf - 1) * (jet.pt() - matchedGenJet->pt()) / jet.pt();
        rnd_->Gaus(0.0, 1.);
    }
    else if (sf > 1)
    {
        double width = res * std::sqrt(sf * sf - 1);
        jerFactor = 1 + rnd_->Gaus(0.0, width);
    }
    else {
        rnd_->Gaus(0.0, 1.);
    }

    // truncate the jer factor when the energy is too small
    jerFactor = std::max(jerFactor, 1e-2 / jet.energy());

    // finally, update the jet
    jet.setP4(jet.p4() * jerFactor);
}

// campaign specific functions
// tight ID criteria from jetmet POG
bool TreeMaker::tightJetID_2016(pat::Jet& jet) {
    bool passPOGID;
    double absEta = fabs(jet.eta());

    if (absEta <= 2.4)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.9 &&
            jet.nConstituents() > 1 &&
            jet.chargedHadronEnergyFraction() > 0. &&
            jet.chargedMultiplicity() > 0 &&
            jet.chargedEmEnergyFraction() < 0.99;
    }
    else if (absEta <= 2.7)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.9 &&
            jet.nConstituents() > 1;
    }
    else if (absEta <= 3.)
    {
        passPOGID = jet.neutralEmEnergyFraction() > 0.01 &&
            jet.neutralHadronEnergyFraction() < 0.98 &&
            jet.neutralMultiplicity() > 2;
    }
    else
    {
        passPOGID = jet.neutralEmEnergyFraction() < 0.9 &&
            jet.neutralMultiplicity() > 10;
    }
    return passPOGID;
}

bool TreeMaker::tightJetID_2017(pat::Jet& jet) {
    bool passPOGID;
    double absEta = fabs(jet.eta());

    if (absEta <= 2.4)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.9 &&
            jet.nConstituents() > 1 &&
            jet.chargedHadronEnergyFraction() > 0. &&
            jet.chargedMultiplicity() > 0;
    }
    else if (absEta <= 2.7)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.9 &&
            jet.nConstituents() > 1;
    }
    else if (absEta <= 3.)
    {
        passPOGID = jet.neutralEmEnergyFraction() > 0.02 &&
            jet.neutralEmEnergyFraction() < 0.99 &&
            jet.neutralMultiplicity() > 2;
    }
    else
    {
        passPOGID = jet.neutralEmEnergyFraction() < 0.9 &&
            jet.neutralHadronEnergyFraction() > 0.02 &&
            jet.neutralMultiplicity() > 10;
    }
    return passPOGID;
}

bool TreeMaker::tightJetID_2018(pat::Jet& jet) {
    bool passPOGID;
    double absEta = fabs(jet.eta());
    if (absEta <= 2.6)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.9 &&
            jet.nConstituents() > 1 &&
            jet.chargedHadronEnergyFraction() > 0. &&
            jet.chargedMultiplicity() > 0;
    }
    else if (absEta <= 2.7)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.9 &&
            jet.neutralEmEnergyFraction() < 0.99 &&
            jet.chargedMultiplicity() > 0;
    }
    else if (absEta <= 3.)
    {
        passPOGID = jet.neutralEmEnergyFraction() > 0.02 &&
            jet.neutralEmEnergyFraction() < 0.99 &&
            jet.neutralMultiplicity() > 2;
    }
    else
    {
        passPOGID = jet.neutralEmEnergyFraction() < 0.9 &&
            jet.neutralHadronEnergyFraction() > 0.2 &&
            jet.neutralMultiplicity() > 10;
    }
    return passPOGID;
}

DEFINE_FWK_MODULE(TreeMaker);
