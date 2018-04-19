/*
 * CMSSW module that does some event selection and outputs a tree that can be used to measure
 * b-tagging scale factors.
 *
 * Authors:
 *   - Marcel Rieger
 *   - David Schmidt
 */

#include <memory>
#include <cmath>
#include <fnmatch.h>

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "TH1F.h"
#include "TTree.h"

typedef std::string string;
typedef std::vector<string> vstring;
typedef std::vector<int32_t> vint;
typedef math::XYZTLorentzVector LorentzVector;

double EMPTY_VALUE = -1e5;

enum VertexType
{
    V_INVALID = 0,
    V_VALID
};

enum ElectronType
{
    E_INVALID = 0,
    E_LOOSE,
    E_TIGHT
};

enum MuonType
{
    M_INVALID = 0,
    M_LOOSE,
    M_TIGHT
};

enum JetType
{
    J_INVALID = 0,
    J_LOOSE,
    J_TIGHT
};

enum LeptonChannel
{
    C_INVALID = 0,
    C_EE,
    C_EMU,
    C_MUMU
};

vstring treeVariables = {
    "event", "run", "lumi", "channel", "rho", "weight",
    "lep1_E", "lep1_px", "lep1_py", "lep1_pz", "lep1_iso", "lep1_pdg", "lep1_tight",
    "lep2_E", "lep2_px", "lep2_py", "lep2_pz", "lep2_iso", "lep2_pdg", "lep2_tight"
};

class CSVTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit CSVTreeMaker(const edm::ParameterSet&);
    ~CSVTreeMaker();

private:
    // interface methods
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // methods for handling variables and output objects
    void resetVariables();
    void setVariable(const string& name, double value);
    void prepareOutput();

    // selection methods
    bool metFilterSelection(const edm::Event&);
    bool leptonSelection(std::vector<pat::Electron>&, std::vector<pat::Electron>&,
        std::vector<pat::Muon>&, std::vector<pat::Muon>&, reco::RecoCandidate*&,
        reco::RecoCandidate*&, LeptonChannel&);
    bool triggerSelection(const edm::Event&, LeptonChannel&);
    VertexType vertexID(reco::Vertex&);
    ElectronType electronID(pat::Electron&, reco::Vertex&);
    MuonType muonID(pat::Muon&, reco::Vertex&);
    JetType jetID(pat::Jet&, reco::RecoCandidate*, reco::RecoCandidate*);

    // options
    bool verbose_;
    bool isData_;
    string leptonChannel_;
    vstring eeTriggers_;
    vstring emuTriggers_;
    vstring mumuTriggers_;
    vstring metFilters_;
    vstring jesFiles_;
    vint jesRanges_;
    vstring jesUncFiles_;
    string jesUncSrcFile_;
    vstring jesUncSources_;

    // tokes
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
    edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > vertexToken_;
    edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
    edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
    edm::EDGetTokenT<std::vector<pat::MET> > metToken_;
    edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVIDToken_;

    // additional members
    std::map<string, double> vars_;
    size_t n_events_;
    size_t n_selected_;
    TTree* tree_;
    TH1F* weight_hist_;
    size_t nJESRanges_;
    size_t nJESFilesPerRange_;
    std::vector<FactorizedJetCorrector*> jetCorrectors_; // per jes range
    std::vector<JetCorrectionUncertainty*> jetCorrectorUncs_; // per jes range
    std::vector<JetCorrectionUncertainty*> jetCorrectorUncSources_; // per source, one for all ranges

    // edm::EDGetTokenT<reco::GenJetCollection> genJetsToken;

    // edm::EDGetTokenT <edm::TriggerResults> triggerResultsToken;
    // edm::EDGetTokenT <edm::TriggerResults> filterResultsToken;

    // // new MVAelectron
    // // edm::EDGetTokenT< edm::View<pat::Electron> > EDMElectronsToken;
    // // MVA values and categories
    // // edm::EDGetTokenT<edm::ValueMap<float> > EDMeleMVAvaluesToken;
    // // edm::EDGetTokenT<edm::ValueMap<int> > EDMeleMVAcategoriesToken;

    // edm::EDGetTokenT <reco::VertexCollection> vertexToken;
    // edm::EDGetTokenT <pat::ElectronCollection> electronToken;
    // edm::EDGetTokenT <pat::MuonCollection> muonToken;
    // edm::EDGetTokenT <pat::JetCollection> jetToken;
    // // edm::EDGetTokenT <edm::View<pat::Jet > > jetToken;

    // edm::EDGetTokenT <pat::METCollection> metToken;
    // // edm::EDGetTokenT <pat::METCollection> metNoHFToken;

    // edm::EDGetTokenT <pat::PackedCandidateCollection> packedpfToken;

    // edm::EDGetTokenT <reco::BeamSpot> beamspotToken;
    // // edm::EDGetTokenT <reco::ConversionCollection> EDMConversionCollectionToken;
    // edm::EDGetTokenT <double> rhoToken;
    // edm::EDGetTokenT <reco::GenParticleCollection> mcparicleToken;
    // edm::EDGetTokenT <std::vector< PileupSummaryInfo > > puInfoToken;

    // edm::EDGetTokenT <GenEventInfoProduct> genInfoProductToken;

    // HLTConfigProvider hlt_config_;

    // bool verbose_;
    // bool isData;
    // const int insample_;
    // const std::string sampleName_;
    // const double xSec_;
    // const double nGen_;
    // const double intLumi_;

    // std::string hltTag;
    // std::string filterTag;

    // int nevents;
    // double nevents_wgt;

    // int nevents_clean;
    // double nevents_clean_wgt;

    // EventVars *eve;
    // TH1F* h_nEvents;
    // TH1F* h_numPV;

    // // EGammaMvaEleEstimatorCSA14* myMVATrig;
};

CSVTreeMaker::CSVTreeMaker(const edm::ParameterSet& iConfig)
    : verbose_(iConfig.getUntrackedParameter<bool>("verbose"))
    , isData_(iConfig.getParameter<bool>("isData"))
    , leptonChannel_(iConfig.getParameter<string>("leptonChannel"))
    , eeTriggers_(iConfig.getParameter<vstring>("eeTriggers"))
    , emuTriggers_(iConfig.getParameter<vstring>("emuTriggers"))
    , mumuTriggers_(iConfig.getParameter<vstring>("mumuTriggers"))
    , metFilters_(iConfig.getParameter<vstring>("metFilters"))
    , jesFiles_(iConfig.getParameter<vstring>("jesFiles"))
    , jesRanges_(iConfig.getParameter<vint>("jesRanges"))
    , jesUncFiles_(iConfig.getParameter<vstring>("jesUncFiles"))
    , jesUncSrcFile_(iConfig.getParameter<string>("jesUncSrcFile"))
    , jesUncSources_(iConfig.getParameter<vstring>("jesUncSources"))
    , genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoCollection")))
    , lheEventToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventCollection")))
    , triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBitsCollection")))
    , metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBitsCollection")))
    , pileupInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupInfoCollection")))
    , beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotCollection")))
    , vertexToken_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexCollection")))
    , electronToken_(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection")))
    , muonToken_(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonCollection")))
    , metToken_(consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metCollection")))
    , jetToken_(consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetCollection")))
    , rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection")))
    , eleVIDToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVIDCollection")))
    , n_events_(0)
    , n_selected_(0)
    , tree_(0)
    , weight_hist_(0)
{
    if (verbose_)
    {
        std::cout << "running CSVTreeMaker in verbose mode" << std::endl;
    }

    usesResource("TFileService");

    // eve = 0;
    // worldTree->Branch("eve.", "EventVars", &eve, 8000, 1);
    // h_nEvents = fs_->make<TH1F>("nEvents", "nEvents", 3, 0, 3);
    // h_nEvents->GetXaxis()->SetBinLabel(1, "All");
    // h_nEvents->GetXaxis()->SetBinLabel(2, "Pos");
    // h_nEvents->GetXaxis()->SetBinLabel(3, "Neg");

    // h_numPV = fs_->make<TH1F>("numPVs", "numPVs", 50, 0, 50);

    // nevents = 0;
    // nevents_wgt = 0;

    // nevents_clean = 0;
    // nevents_clean_wgt = 0;
}

CSVTreeMaker::~CSVTreeMaker()
{
    if (tree_)
    {
        delete tree_;
        tree_ = 0;
    }
    if (weight_hist_)
    {
        delete weight_hist_;
        weight_hist_ = 0;
    }
}

void CSVTreeMaker::resetVariables()
{
    std::map<string, double>::iterator it;
    for (it = vars_.begin(); it != vars_.end(); it++)
    {
        it->second = EMPTY_VALUE;
    }
}

void CSVTreeMaker::setVariable(const string& name, double value)
{
    vars_[name] = value;
}

void CSVTreeMaker::beginJob()
{
    edm::Service<TFileService> fs;

    // create the histogram containing the sum of event weights in its only bin
    weight_hist_ = fs->make<TH1F>("event_weights", "", 1, 0., 1.);

    // create the tree
    tree_ = fs->make<TTree>("events", "events");

    // define variables and add branches based on added variables
    for (size_t i = 0; i < treeVariables.size(); i++)
    {
        string name = treeVariables[i];

        if (verbose_)
        {
            std::cout << "adding variable: " << name << std::endl;
        }

        vars_[name] = EMPTY_VALUE;
        tree_->Branch(name.c_str(), &vars_[name], (name + "/D").c_str());
    }
}

void CSVTreeMaker::endJob()
{
    tree_->Write();
    weight_hist_->Write();

    // log some stats
    if (verbose_)
    {
        std::cout << "processed events: " << n_events_ << std::endl;
        std::cout << "selected events: " << n_selected_ << std::endl;
    }
}

void CSVTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    resetVariables();
    n_events_ += 1;

    // check MET filters
    if (!metFilterSelection(iEvent))
    {
        return;
    }

    // read the select the vertex
    edm::Handle<std::vector<reco::Vertex> > verticesHandle;
    iEvent.getByToken(vertexToken_, verticesHandle);
    if (!verticesHandle.isValid())
    {
        return;
    }
    std::vector<reco::Vertex> vertices(*verticesHandle);
    reco::Vertex vertex = vertices.front();
    if (vertexID(vertex) != V_VALID)
    {
        return;
    }

    // read and select electrons
    edm::Handle<std::vector<pat::Electron> > electronsHandle;
    iEvent.getByToken(electronToken_, electronsHandle);
    std::vector<pat::Electron> electrons;
    std::vector<pat::Electron> tightElectrons;
    if (electronsHandle.isValid())
    {
        std::vector<pat::Electron> electronsCopy(*electronsHandle);
        for (pat::Electron& electron : electronsCopy)
        {
            ElectronType type = electronID(electron, vertex);
            if (type == E_LOOSE)
            {
                electrons.push_back(electron);
            }
            else if (type == E_TIGHT)
            {
                electrons.push_back(electron);
                tightElectrons.push_back(electron);
            }
        }
    }

    // read and select muons
    edm::Handle<std::vector<pat::Muon> > muonsHandle;
    iEvent.getByToken(muonToken_, muonsHandle);
    std::vector<pat::Muon> muons;
    std::vector<pat::Muon> tightMuons;
    if (muonsHandle.isValid())
    {
        std::vector<pat::Muon> muonsCopy(*muonsHandle);
        for (pat::Muon& muon : muonsCopy)
        {
            MuonType type = muonID(muon, vertex);
            if (type == M_LOOSE)
            {
                muons.push_back(muon);
            }
            else if (type == M_TIGHT)
            {
                muons.push_back(muon);
                tightMuons.push_back(muon);
            }
        }
    }

    // lepton selection
    reco::RecoCandidate* lep1 = nullptr;
    reco::RecoCandidate* lep2 = nullptr;
    LeptonChannel channel = C_INVALID;
    bool passedLeptonSelection = leptonSelection(
        electrons, tightElectrons, muons, tightMuons, lep1, lep2, channel);
    if (!passedLeptonSelection)
    {
        return;
    }

    // trigger selection
    if (!triggerSelection(iEvent, channel))
    {
        return;
    }

    if (0)
    {
        //     edm::Handle<pat::METCollection> pfmet;
        //     iEvent.getByToken(metToken, pfmet);

        //     // edm::Handle<pat::METCollection> pfmetNoHF;
        //     // iEvent.getByToken(metNoHFToken,pfmetNoHF);

        //     edm::Handle<double> rhoHandle;
        //     iEvent.getByToken(rhoToken, rhoHandle);
        //     ////------- set up rho for lepton effArea Isolation correction
        //     double rho_event = ((rhoHandle.isValid())) ? *rhoHandle : -99;
        //     miniAODhelper.SetRho(rho_event);

        //     edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
        //     if (!isData)
        //         iEvent.getByToken(puInfoToken, PupInfo);

        //     edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
        //     if (!isData)
        //         iEvent.getByToken(genInfoProductToken, GenEventInfoHandle);

        //     double GenEventInfoWeight = 1;
        //     if (!isData)
        //         GenEventInfoWeight = GenEventInfoHandle.product()->weight();
        //     h_nEvents->Fill(0.5);
        //     if (GenEventInfoWeight > 0)
        //         h_nEvents->Fill(1.5);
        //     else
        //         h_nEvents->Fill(2.5);

        //     eve->numPVs_ = numpv;

        //     if (GenEventInfoWeight > 0)
        //         h_numPV->Fill(numpv, 1);
        //     else
        //         h_numPV->Fill(numpv, -1);

        //     double numTruePV = -1;
        //     double numGenPV = -1;
        //     if ((PupInfo.isValid()))
        //     {
        //         for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
        //         {
        //             int BX = PVI->getBunchCrossing();
        //             if (BX == 0)
        //             {
        //                 numTruePV = PVI->getTrueNumInteractions();
        //                 numGenPV = PVI->getPU_NumInteractions();
        //             }
        //         }
        //     }

        //     eve->numTruePV_ = numTruePV;
        //     eve->numGenPV_ = numGenPV;

        //     double wgt_lumi = 1.;
        //     if (insample_ >= 0)
        //         wgt_lumi = xSec_ * intLumi_ * 1. / nGen_;

        //     // Increment event counter
        //     nevents++;
        //     nevents_wgt += wgt_lumi;

        //     // Number of events after filtering, triggers, etc. (now nothing)
        //     nevents_clean++;
        //     nevents_clean_wgt += wgt_lumi;

        //     int run = iEvent.id().run();
        //     int lumi = iEvent.luminosityBlock();
        //     long evt = iEvent.id().event();

        //     eve->run_ = run;
        //     eve->lumi_ = lumi;
        //     eve->evt_ = evt;

        //     const JetCorrector* corrector = JetCorrector::getJetCorrector("ak4PFchsL1L2L3", iSetup); //Get the jet corrector from the event setup

        //     miniAODhelper.SetJetCorrector(corrector);

        //     /////////
        //     ///
        //     /// Electrons
        //     ///
        //     ////////
        //     // miniAODhelper.SetElectronMVAinfo(h_conversioncollection, bsHandle);
        //     // /// Ele MVA ID
        //     // std::vector<pat::Electron> selectedElectrons_tight = miniAODhelper.GetSelectedElectrons( electrons, minTightLeptonPt, electronID::electronEndOf15MVA80iso0p15, 2.4 );
        //     // std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons( electrons, looseLeptonPt, electronID::electronEndOf15MVA80iso0p15, 2.4 );

        //     /// Ele Cut Based ID
        //     std::vector<pat::Electron> selectedElectrons_tight = miniAODhelper.GetSelectedElectrons(*electrons, minTightLeptonPt, electronID::electron80XCutBasedM, 2.4);
        //     std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons(*electrons, looseLeptonPt, electronID::electron80XCutBasedM, 2.4);

        //     int numTightElectrons = int(selectedElectrons_tight.size());
        //     int numLooseElectrons = int(selectedElectrons_loose.size()); // - numTightElectrons;

        //     /////////
        //     ///
        //     /// Muons
        //     ///
        //     ////////
        //     std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons(*muons, minTightLeptonPt, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.4);
        //     std::vector<pat::Muon> selectedMuons_loose = miniAODhelper.GetSelectedMuons(*muons, looseLeptonPt, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.4);
        //     // std::vector<pat::Muon> selectedMuons_loose = miniAODhelper.GetSelectedMuons( *muons, looseLeptonPt, muonID::muonLoose );

        //     int numTightMuons = int(selectedMuons_tight.size());
        //     int numLooseMuons = int(selectedMuons_loose.size()); // - numTightMuons;

        //     eve->numTightMuons_ = numTightMuons;
        //     eve->numLooseMuons_ = numLooseMuons;
        //     eve->numTightElectrons_ = numTightElectrons;
        //     eve->numLooseElectrons_ = numLooseElectrons;

        //     /////////////
        //     // --- Pass exactly two leptons cuts
        //     // ---------- the lepton selections might change in the future -------
        //     //
        //     /////////////

        //     bool passDIL = ((numTightMuons + numTightElectrons) >= 1 && (numLooseMuons + numLooseElectrons) == 2);
        //     // bool passDIL = (numTightMuons + numTightElectrons) == 2; //old selections

        //     // Event must be DIL, as requested
        //     if (!passDIL)
        //         return;

        //     eve->PassDIL_ = (passDIL) ? 1 : 0;

        //     //// DIL subcategories: ee, mumu, emu
        //     int TwoMuon = 0, TwoElectron = 0, MuonElectron = 0;
        //     if (numTightMuons >= 1 && numLooseMuons == 2 && numLooseElectrons == 0)
        //         TwoMuon = 1;
        //     if (numTightElectrons >= 1 && numLooseElectrons == 2 && numLooseMuons == 0)
        //         TwoElectron = 1;
        //     if ((numTightMuons + numTightElectrons) >= 1 && numLooseMuons == 1 && numLooseElectrons == 1)
        //         MuonElectron = 1;
        //     // if ( numTightMuons==2 && numTightElectrons==0 ) TwoMuon = 1; //old selections
        //     // if ( numTightElectrons==2 && numTightMuons==0 ) TwoElectron = 1;
        //     // if ( numTightMuons==1 && numTightElectrons==1 ) MuonElectron = 1;

        //     eve->TwoMuon_ = TwoMuon;
        //     eve->TwoElectron_ = TwoElectron;
        //     eve->MuonElectron_ = MuonElectron;

        //     ///////////
        //     /////  extra lepton variables: opposite charge, dR_leplep and mass_leplep
        //     //////////
        //     vint lepton_trkCharge, lepton_isMuon; //, lepton_isTight, lepton_isLoose;

        //     std::vector<TLorentzVector> vec_TLV_lep;
        //     TLorentzVector sum_lepton_vect;
        //     sum_lepton_vect.SetPxPyPzE(0., 0., 0., 0.);

        //     ///// --- loop over muons
        //     for (std::vector<pat::Muon>::const_iterator iMu = selectedMuons_loose.begin(); iMu != selectedMuons_loose.end(); iMu++)
        //     {
        //         // Get muon 4Vector and add to vecTLorentzVector for muons
        //         TLorentzVector leptonP4;
        //         leptonP4.SetPxPyPzE(iMu->px(), iMu->py(), iMu->pz(), iMu->energy());
        //         vec_TLV_lep.push_back(leptonP4);

        //         sum_lepton_vect += leptonP4;

        //         lepton_isMuon.push_back(1);

        //         int trkCharge = -99;
        //         if (iMu->muonBestTrack().isAvailable())
        //             trkCharge = iMu->muonBestTrack()->charge();
        //         lepton_trkCharge.push_back(trkCharge);
        //     }

        //     ///// --- loop over electrons
        //     for (std::vector<pat::Electron>::const_iterator iEle = selectedElectrons_loose.begin(); iEle != selectedElectrons_loose.end(); iEle++)
        //     {
        //         // Get electron 4Vector and add to vecTLorentzVector for electrons
        //         TLorentzVector leptonP4;
        //         leptonP4.SetPxPyPzE(iEle->px(), iEle->py(), iEle->pz(), iEle->energy());
        //         vec_TLV_lep.push_back(leptonP4);

        //         sum_lepton_vect += leptonP4;

        //         lepton_isMuon.push_back(0);

        //         int trkCharge = -99;
        //         if (iEle->gsfTrack().isAvailable())
        //             trkCharge = iEle->gsfTrack()->charge();
        //         lepton_trkCharge.push_back(trkCharge);
        //     }

        //     eve->lepton_trkCharge_ = lepton_trkCharge;
        //     eve->lepton_vect_TLV_ = vec_TLV_lep; /// non existing yet
        //     eve->lepton_isMuon_ = lepton_isMuon;

        //     /// DIL specific: opposite charged lepton pair
        //     int oppositeLepCharge = -9;
        //     if (lepton_trkCharge.size() == 2)
        //     {
        //         int chg0 = lepton_trkCharge[0];
        //         int chg1 = lepton_trkCharge[1];

        //         if ((chg0 * chg1) == -1)
        //             oppositeLepCharge = 1;
        //         else if ((chg0 * chg1) == 1)
        //             oppositeLepCharge = 0;
        //         else if (chg0 == -99)
        //             oppositeLepCharge = -1;
        //         else if (chg1 == -99)
        //             oppositeLepCharge = -2;
        //         else
        //             oppositeLepCharge = -3;
        //     }
        //     eve->oppositeLepCharge_ = oppositeLepCharge;

        //     /// DIL specific: dR and mass of the dilepton pair
        //     double mass_leplep = -99;
        //     if (vec_TLV_lep.size() == 2)
        //     {
        //         mass_leplep = sum_lepton_vect.M();
        //         eve->mass_leplep_ = mass_leplep;
        //         eve->dR_leplep_ = vec_TLV_lep[0].DeltaR(vec_TLV_lep[1]);
        //     }

        //     ///////////
        //     ///
        //     /// various weights: PU, lepton SF, trigger SF, etc
        //     ///
        //     //////////
        //     // PU scale factor
        //     double wgtPU = 1;
        //     eve->wgt_pu_ = wgtPU;

        //     // Lepton scale factor
        //     float leptonSF = 1.0;
        //     eve->wgt_lepSF_ = leptonSF;

        //     eve->wgt_lumi_ = intLumi_;
        //     eve->wgt_xs_ = xSec_; //mySample.xSec;
        //     eve->wgt_nGen_ = nGen_; //mySample.nGen;

        //     eve->wgt_generator_ = GenEventInfoWeight;

        //     /////////
        //     ///
        //     /// Jets
        //     ///
        //     ////////
        //     // Do jets stuff

        //     /// checking PU jet ID
        //     // for ( unsigned int i=0; i<pfjets->size(); ++i ) {
        //     //   const pat::Jet & patjet = pfjets->at(i);

        //     //   if(!patjet.hasUserFloat("pileupJetId:fullDiscriminant") || !patjet.hasUserInt("pileupJetId:fullId")) {
        //     // 	std::cout << "not mvaValue or idflag info for the jet ---" << std::endl;
        //     // 	continue;
        //     //   }

        //     //   float mvaValue   = patjet.userFloat("pileupJetId:fullDiscriminant");
        //     //   int    idflag = patjet.userInt("pileupJetId:fullId");

        //     //   cout << "jet " << i << " pt " << patjet.pt() << " eta " << patjet.eta() << " PU JetID MVA " << mvaValue << " PU JetID flag " << idflag;
        //     //   if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose ) ){
        //     //   	cout << " pass loose wp";
        //     //   }
        //     //   if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) ){
        //     //   	  cout << " pass medium wp";
        //     //   }
        //     //   if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight ) ){
        //     //     cout << " pass tight wp";
        //     //   }
        //     //   cout << endl;

        //     // }
        //     // cout << endl;

        //     //-------
        //     std::vector<pat::Jet> pfJets_ID = miniAODhelper.GetSelectedJets(*pfjets, 0., 999, jetID::jetLoose, '-');
        //     /// lepton-jet overlap cleanning, might change the input lepton collections in the future
        //     std::vector<pat::Jet> pfJets_ID_clean = miniAODhelper.GetDeltaRCleanedJets(pfJets_ID, selectedMuons_loose, selectedElectrons_loose, 0.4);

        //     std::vector<pat::Jet> unCorrectedJets = miniAODhelper.GetUncorrectedJets(pfJets_ID_clean);
        //     std::vector<pat::Jet> correctedJets_noSys = miniAODhelper.GetCorrectedJets(unCorrectedJets, iEvent, iSetup, genjetCollection);
        //     /// jet pt threshold, 20GeV for now
        //     std::vector<pat::Jet> selectedJets_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, 20., 2.4, jetID::none, '-');
        //     std::vector<pat::Jet> selectedJets_tag_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, 20., 2.4, jetID::none, 'M');

        //     // std::vector<pat::Jet> selectedJets_loose_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, 20., 3.0, jetID::none, '-' );
        //     // std::vector<pat::Jet> selectedJets_loose_tag_noSys_unsorted = miniAODhelper.GetSelectedJets( correctedJets_noSys, 20., 3.0, jetID::none, 'M' );

        //     // bool for passing jet/tag cuts, first is for tree filling, second is for event counting
        //     // bool hasNumJetNumTag = false;

        //     double totalWgt = leptonSF * wgtPU; //wgt_lumi;

        //     ///////
        //     // Loop over systematics: nominal, JESUp, JESDown
        //     ///////
        //     bool passingTwoJet = false;
        //     for (int iSys = 0; iSys < rNumSys; iSys++)
        //     {

        //         //no sys for data,
        //         if (isData)
        //         {
        //             if (iSys != 0)
        //                 continue;
        //         }

        //         eve->wgt_[iSys] = totalWgt;

        //         sysType::sysType iSysType = sysType::NA;
        //         switch (iSys)
        //         {
        //         case 1:
        //             iSysType = sysType::JESup;
        //             break;
        //         case 2:
        //             iSysType = sysType::JESdown;
        //             break;

        //         case 3:
        //             iSysType = sysType::JESFlavorQCDup;
        //             break;
        //         case 4:
        //             iSysType = sysType::JESFlavorQCDdown;
        //             break;

        //         case 5:
        //             iSysType = sysType::JESSinglePionHCALup;
        //             break;
        //         case 6:
        //             iSysType = sysType::JESSinglePionHCALdown;
        //             break;

        //         case 7:
        //             iSysType = sysType::JESAbsoluteScaleup;
        //             break;
        //         case 8:
        //             iSysType = sysType::JESAbsoluteScaledown;
        //             break;

        //         case 9:
        //             iSysType = sysType::JESAbsoluteMPFBiasup;
        //             break;
        //         case 10:
        //             iSysType = sysType::JESAbsoluteMPFBiasdown;
        //             break;

        //         case 11:
        //             iSysType = sysType::JESPileUpPtRefup;
        //             break;
        //         case 12:
        //             iSysType = sysType::JESPileUpPtRefdown;
        //             break;
        //         ////////
        //         case 13:
        //             iSysType = sysType::JESSinglePionECALup;
        //             break;
        //         case 14:
        //             iSysType = sysType::JESSinglePionECALdown;
        //             break;

        //         case 15:
        //             iSysType = sysType::JESPileUpPtBBup;
        //             break;
        //         case 16:
        //             iSysType = sysType::JESPileUpPtBBdown;
        //             break;

        //         case 17:
        //             iSysType = sysType::JESPileUpPtEC1up;
        //             break;
        //         case 18:
        //             iSysType = sysType::JESPileUpPtEC1down;
        //             break;

        //         case 19:
        //             iSysType = sysType::JESPileUpDataMCup;
        //             break;
        //         case 20:
        //             iSysType = sysType::JESPileUpDataMCdown;
        //             break;

        //         case 21:
        //             iSysType = sysType::JESRelativeFSRup;
        //             break;
        //         case 22:
        //             iSysType = sysType::JESRelativeFSRdown;
        //             break;

        //         case 23:
        //             iSysType = sysType::JESTimePtEtaup;
        //             break;
        //         case 24:
        //             iSysType = sysType::JESTimePtEtadown;
        //             break;

        //         default:
        //             iSysType = sysType::NA;
        //             break;
        //         }

        //         // nJets/Tags
        //         // THESE MUST BE INITIALIZED TO -1
        //         eve->numJets_float_[iSys] = -1;
        //         eve->numTags_float_[iSys] = -1;

        //         /////////
        //         ///
        //         /// Pfjets
        //         ///
        //         ////////

        //         std::vector<pat::Jet> correctedJets = (iSys == 0) ? correctedJets_noSys : miniAODhelper.GetCorrectedJets(unCorrectedJets, iEvent, iSetup, genjetCollection, iSysType);
        //         std::vector<pat::Jet> selectedJets_unsorted = (iSys == 0) ? selectedJets_noSys_unsorted : miniAODhelper.GetSelectedJets(correctedJets, 20., 2.4, jetID::none, '-');
        //         std::vector<pat::Jet> selectedJets_tag_unsorted = (iSys == 0) ? selectedJets_tag_noSys_unsorted : miniAODhelper.GetSelectedJets(correctedJets, 20., 2.4, jetID::none, 'M');

        //         // Sort jet collections by pT
        //         std::vector<pat::Jet> selectedJets = miniAODhelper.GetSortedByPt(selectedJets_unsorted);
        //         std::vector<pat::Jet> selectedJets_tag = miniAODhelper.GetSortedByPt(selectedJets_tag_unsorted);

        //         // Get numJets, numTags
        //         int numJet = int(selectedJets.size());
        //         int numTag = int(selectedJets_tag.size());

        //         if (numJet < 2)
        //             continue; //return;
        //         passingTwoJet = true;
        //         // Get Corrected MET (propagating JEC and JER)
        //         std::vector<pat::Jet> oldJetsForMET = miniAODhelper.GetSelectedJets(*pfjets, 0., 999, jetID::jetMETcorrection, '-');
        //         std::vector<pat::Jet> oldJetsForMET_uncorr = miniAODhelper.GetUncorrectedJets(oldJetsForMET);
        //         std::vector<pat::Jet> newJetsForMET = miniAODhelper.GetCorrectedJets(oldJetsForMET_uncorr, iEvent, iSetup, genjetCollection, iSysType);
        //         std::vector<pat::MET> newMETs = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *pfmet);
        //         // std::vector<pat::MET> newMETs_NoHF = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *pfmetNoHF);

        //         pat::MET correctedMET = newMETs.at(0);
        //         // TLorentzVector metV(correctedMET.px(),correctedMET.py(),0.0,correctedMET.pt());
        //         // pat::MET correctedMET_NoHF = newMETs_NoHF.at(0);

        //         // Initialize event vars, selected jets
        //         double mht_px = 0.;
        //         double mht_py = 0.;

        //         vecTLorentzVector jetV;
        //         std::vector<double> csvV;
        //         std::vector<double> cMVAV;
        //         std::vector<double> jet_ptV; /// non existing
        //         std::vector<double> jet_etaV;
        //         vint jet_flavour_vect;
        //         vint jet_partonflavour_vect;

        //         vdouble jet_PUID_mvaV;
        //         vint jet_PUID_flagV;
        //         vint jet_PUID_passWPLooseV;

        //         // int jeti = 0;
        //         // Loop over selected jets
        //         for (std::vector<pat::Jet>::const_iterator iJet = selectedJets.begin(); iJet != selectedJets.end(); iJet++)
        //         {

        //             jet_ptV.push_back(iJet->pt());
        //             jet_etaV.push_back(iJet->eta());
        //             jet_partonflavour_vect.push_back(iJet->partonFlavour());
        //             jet_flavour_vect.push_back(iJet->hadronFlavour());

        //             // Get jet 4Vector and add to vecTLorentzVector for jets
        //             TLorentzVector jet0p4;
        //             jet0p4.SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
        //             jetV.push_back(jet0p4);

        //             // MHT
        //             mht_px += -iJet->px();
        //             mht_py += -iJet->py();

        //             //      double myCSV = miniAODhelper.GetJetCSV(*iJet, "pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //             double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //             csvV.push_back(myCSV);

        //             //      double mycMVA = miniAODhelper.GetJetCSV(*iJet, "pfCombinedMVABJetTags");
        //             double mycMVA = iJet->bDiscriminator("pfCombinedMVAV2BJetTags");
        //             cMVAV.push_back(mycMVA);

        //             //PU jet ID
        //             if (!iJet->hasUserFloat("pileupJetId:fullDiscriminant") || !iJet->hasUserInt("pileupJetId:fullId"))
        //             {
        //                 std::cout << "no mvaValue or idflag info for the jet ---" << std::endl;
        //             }

        //             float mvaValue = iJet->userFloat("pileupJetId:fullDiscriminant");
        //             int idflag = iJet->userInt("pileupJetId:fullId");

        //             // cout << "Sys " << iSys << ": jet " << jeti << " pt " << iJet->pt() << " eta " << iJet->eta() << " PU JetID MVA " << mvaValue << " PU JetID flag " << idflag << endl;;
        //             // jeti++;
        //             bool passPUIDLoose = PileupJetIdentifier::passJetId(idflag, PileupJetIdentifier::kLoose);

        //             jet_PUID_mvaV.push_back(mvaValue);
        //             jet_PUID_flagV.push_back(idflag);
        //             jet_PUID_passWPLooseV.push_back(passPUIDLoose);

        //         } // end loop over iJet

        //         ////
        //         eve->jet_vect_TLV_[iSys] = jetV;
        //         eve->jet_CSV_[iSys] = csvV;
        //         eve->jet_cMVA_[iSys] = cMVAV;
        //         eve->jet_pt_[iSys] = jet_ptV;
        //         eve->jet_eta_[iSys] = jet_etaV;
        //         eve->jet_flavour_[iSys] = jet_flavour_vect;
        //         eve->jet_partonflavour_[iSys] = jet_partonflavour_vect;

        //         //PU jet ID
        //         eve->jet_PUID_mva_[iSys] = jet_PUID_mvaV;
        //         eve->jet_PUID_flag_[iSys] = jet_PUID_flagV;
        //         eve->jet_PUID_passWPLoose_[iSys] = jet_PUID_passWPLooseV;

        //         // Add lepton 4Vector quantities to MHT
        //         mht_px += -sum_lepton_vect.Px();
        //         mht_py += -sum_lepton_vect.Py();
        //         double MHT = sqrt(mht_px * mht_px + mht_py * mht_py);

        //         bool PassBigDiamondZmask = (MuonElectron || (mass_leplep < (65.5 + 3 * MHT / 8)) || (mass_leplep > (108 - MHT / 4)) || (mass_leplep < (79 - 3 * MHT / 4)) || (mass_leplep > (99 + MHT / 2)));
        //         eve->PassZmask_ = (PassBigDiamondZmask) ? 1 : 0;

        //         double MET = correctedMET.pt();
        //         bool PassBigDiamondZmaskMET = (MuonElectron || (mass_leplep < (65.5 + 3 * MET / 8)) || (mass_leplep > (108 - MET / 4)) || (mass_leplep < (79 - 3 * MET / 4)) || (mass_leplep > (99 + MET / 2)));
        //         eve->PassZmaskMET_ = (PassBigDiamondZmaskMET) ? 1 : 0;

        //         // Compute MHT_
        //         eve->MHT_[iSys] = MHT;
        //         // MET
        //         eve->MET_[iSys] = correctedMET.pt();
        //         eve->MET_phi_[iSys] = correctedMET.phi();
        //         // MET_NoHF
        //         // eve->METNoHF_[iSys]      = correctedMET_NoHF.pt();
        //         // eve->METNoHF_phi_[iSys]  = correctedMET_NoHF.phi();

        //         // nJets/Tags
        //         eve->numJets_float_[iSys] = numJet;
        //         eve->numTags_float_[iSys] = numTag;

        //         // // jet1-4 pT
        //         // if( selectedJets.size()>0 ) eve->first_jet_pt_[iSys]  = selectedJets.at(0).pt();
        //         // if( selectedJets.size()>1 ) eve->second_jet_pt_[iSys] = selectedJets.at(1).pt();

        //     } // end loop over systematics

        //     // Fill tree if pass full selection
        //     if (passingTwoJet)
    }

    tree_->Fill();
}

bool CSVTreeMaker::metFilterSelection(const edm::Event& iEvent)
{
    // check MET filters
    // unlike triggers all of them have to be accepted
    edm::Handle<edm::TriggerResults> metFilterBitsHandle;
    iEvent.getByToken(metFilterBitsToken_, metFilterBitsHandle);
    if (!metFilterBitsHandle.isValid())
    {
        return false;
    }

    const edm::TriggerNames& metFilterNames = iEvent.triggerNames(*metFilterBitsHandle);
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

bool CSVTreeMaker::triggerSelection(const edm::Event& iEvent, LeptonChannel& channel)
{
    edm::Handle<edm::TriggerResults> triggerBitsHandle;
    iEvent.getByToken(triggerBitsToken_, triggerBitsHandle);
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

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBitsHandle);
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

bool CSVTreeMaker::leptonSelection(std::vector<pat::Electron>& electrons,
    std::vector<pat::Electron>& tightElectrons, std::vector<pat::Muon>& muons,
    std::vector<pat::Muon>& tightMuons, reco::RecoCandidate*& lep1,
    reco::RecoCandidate*& lep2, LeptonChannel& channel)
{
    // require exactly 2 leptons, at least 1 tight
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
    else
    {
        throw std::runtime_error("could not determine lepton channel!");
    }

    // check for opposite charges
    if (lep1->charge() == lep2->charge())
    {
        return false;
    }

    // a particular lepton channel might be required (by name), e.g. for data
    if (!leptonChannel_.empty())
    {
        if ((leptonChannel_ == "ee" && channel != C_EE) ||
            (leptonChannel_ == "emu" && channel != C_EMU) ||
            (leptonChannel_ == "mumu" && channel != C_MUMU))
        {
            return false;
        }
    }

    // re-order py pt
    if (lep2->pt() > lep1->pt())
    {
        std::swap(lep1, lep2);
    }

    return true;
}

VertexType CSVTreeMaker::vertexID(reco::Vertex& vertex)
{
    bool isValid = !vertex.isFake() &&
        vertex.ndof() >= 4.0 &&
        abs(vertex.z()) <= 24.0 &&
        abs(vertex.position().Rho()) <= 2.0;
    return isValid ? V_VALID : V_INVALID;
}

ElectronType CSVTreeMaker::electronID(pat::Electron& electron, reco::Vertex& vertex)
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
    if (fabs(gsfTrackRef.get()->dz(vertex.position())) >= (isBarrel ? 0.05 : 0.1))
    {
        return E_INVALID;
    }
    if (fabs(gsfTrackRef.get()->dxy(vertex.position())) >= (isBarrel ? 0.1 : 0.2))
    {
        return E_INVALID;
    }

    // TODO: VID

    // tight or loose decision only depends on pt
    return electron.pt() > 25. ? E_TIGHT : E_LOOSE;
}

MuonType CSVTreeMaker::muonID(pat::Muon& muon, reco::Vertex& vertex)
{
    // loose pt cut
    if (muon.pt() <= 15.)
    {
        return M_INVALID;
    }

    // loose eta cut
    double absEta = fabs(muon.eta());
    if (absEta >= 2.4)
    {
        return M_INVALID;
    }

    // loose isolation cut
    const auto& isoVars = muon.pfIsolationR04();
    float neutralSumPt = isoVars.sumNeutralHadronEt + isoVars.sumPhotonEt - 0.5 * isoVars.sumPUPt;
    float iso = (isoVars.sumChargedHadronPt + fmax(0.0, neutralSumPt)) / muon.pt();
    muon.addUserFloat("iso", iso, true);
    if (iso >= 0.25)
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


    // tight or loose decision depends on pt, eta and iso
    return (muon.pt() > 26. && absEta < 2.1 && iso < 0.15)? M_TIGHT : M_LOOSE;
}

JetType CSVTreeMaker::jetID(pat::Jet& jet, reco::RecoCandidate* lep1, reco::RecoCandidate* lep2)
{
    // loose pt cut
    if (jet.pt() <= 20.)
    {
        return J_INVALID;
    }

    // eta cut
    double absEta = fabs(jet.eta());
    if (absEta >= 2.4)
    {
        return J_INVALID;
    }

    // PU jet ID
    // TODO: update to different WP?
    if (jet.userInt("pileupJetIdUpdated:fullId") < 4)
    {
        return J_INVALID;
    }

    // check distance to selected leptons
    if (deltaR(jet.p4(), lep1->p4()) < 0.4)
    {
        return J_INVALID;
    }

    // loose ID criteria from jetmet POG
    double passPOGID;
    if (absEta <= 2.4)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.99 &&
            jet.neutralEmEnergyFraction() < 0.99 &&
            jet.nConstituents() > 1 &&
            jet.chargedHadronEnergyFraction() > 0. &&
            jet.chargedMultiplicity() > 0 &&
            jet.chargedEmEnergyFraction() < 0.99;
    }
    else if (absEta <= 2.7)
    {
        passPOGID = jet.neutralHadronEnergyFraction() < 0.99 &&
            jet.neutralEmEnergyFraction() < 0.99 &&
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
    if (!passPOGID)
    {
        return J_INVALID;
    }

    // tight or loose decision depends ony on pt
    return jet.pt() > 30. ? J_TIGHT : J_LOOSE;
}

double deltaR(const LorentzVector& v1, const LorentzVector& v2)
{
    double deltaEta = v1.eta() - v2.eta();
    double deltaPhi = fabs(v1.phi() - v2.phi());
    if (deltaPhi > M_PI)
    {
        deltaPhi = 2. * M_PI - deltaPhi;
    }
    return sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
}

DEFINE_FWK_MODULE(CSVTreeMaker);
