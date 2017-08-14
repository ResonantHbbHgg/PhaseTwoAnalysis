// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/BasicPatDistrib
// Class:      BasicPatDistrib
// 
/**\class BasicPatDistrib BasicPatDistrib.cc PhaseTwoAnalysis/BasicPatDistrib/plugins/BasicPatDistrib.cc

Description: produces histograms of basic quantities from PAT collections

Implementation:
   - lepton isolation might need to be refined
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/UPGTrackerTDRStudies#Muon_identification
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
      /!\ no ID is implemented for forward electrons as:
      - PFClusterProducer does not run on miniAOD
      - jurassic isolation needs tracks
   - PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
   - no JEC applied
   - b-tagging WP come from Run-2 https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Supported_Algorithms_and_Operati 
      - for pfCombinedInclusiveSecondaryVertexV2BJetTags: L = 0.5426, M = 0.8484, T = 0.9535)
      - for deepCSV: L = 0.2219, M = 0.6324, T = 0.8958
*/
/*
    Modified: Maxime Gouzevitch
    Date: August 2017 for Barrel TDR upgrade

    Modifier: Sam Higginbotham 
    Date: July 6th,2017. 
    - Photon ID:Selection based on https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details 
    - Starting template taken DIRECTLY from RecoEgamma/Examples/plugins/PatPhotonSimpleAnalyzer.cc  (Author:  M.B. Anderson based on simple photon analyzer by:  J. Stilley, A. Askew)


*/
//
// Original Author:  Elvire Bouvier
//         Created:  Wed, 14 Jun 2017 14:16:22 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//for photons
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "TMath.h"
#include "TTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BasicPatDistrib : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit BasicPatDistrib(const edm::ParameterSet&);
    ~BasicPatDistrib();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    bool useDeepCSV_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    PFJetIDSelectionFunctor jetIDTight_;
    edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
    edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> allGenPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<pat::Photon>> photonsToken_;

    // MC truth in fiducial phase space
    TH1D* h_genJets_n_;
    TH1D* h_genJets_pt_;
    TH1D* h_genJets_phi_;
    TH1D* h_genJets_eta_;

    // Vertices
    TH1D* h_allVertices_n_;
    // ... that pass ID
    TH1D* h_goodVertices_n_;

    // Jets
    TH1D* h_allJets_n_;
    TH1D* h_allJets_pt_;
    TH1D* h_allJets_phi_;
    TH1D* h_allJets_eta_;
    TH1D* h_allJets_csv_;
    TH1D* h_allJets_id_;
    // ... that pass kin cuts, loose ID
    TH1D* h_goodJets_n_;
    TH1D* h_goodJets_nb_;
    TH1D* h_goodJets_pt_;
    TH1D* h_goodJets_phi_;
    TH1D* h_goodJets_eta_;
    TH1D* h_goodJets_csv_;
    TH1D* h_goodLJets_n_;
    TH1D* h_goodLJets_nb_;
    TH1D* h_goodLJets_pt_;
    TH1D* h_goodLJets_phi_;
    TH1D* h_goodLJets_eta_;
    TH1D* h_goodLJets_csv_;
    TH1D* h_goodBJets_n_;
    TH1D* h_goodBJets_nb_;
    TH1D* h_goodBJets_pt_;
    TH1D* h_goodBJets_phi_;
    TH1D* h_goodBJets_eta_;
    TH1D* h_goodBJets_csv_;
    //Sam added recoJet histos
    TH1D* h_goodrecoJetBJets_Higgs_n_;
    TH1F* h_recobjetHiggsMass_;
    
    TH1D* h_goodMET_pt_;
    TH1D* h_goodMET_phi_;      

    //Photons! Incorporated from PatPhotonSimpleAnalyzer.h (see .cc file in the intro)    
    //std::string outputFile_;   // output file
    double minPhotonEt_;       // minimum photon Et
    double minPhotonAbsEta_;   // min and
    double maxPhotonAbsEta_;   // max abs(eta)
    double minPhotonR9_;       // minimum R9 = E(3x3)/E(SuperCluster)
    double maxPhotonHoverE_;   // maximum HCAL / ECAL
    double maxIetaIeta_;       // min sigma Ieta Ieta... ECAL crystal location
    bool   createPhotonTTree_; // Create a TTree of photon variables

    // Will be used for creating TTree of photons.
    // These names did not have to match those from a phtn->...
    // but do match for clarity.
    struct struct_recPhoton {
        float isolationEcalRecHit;
        float isolationHcalRecHit;
        float isolationSolidTrkCone;
        float isolationHollowTrkCone;
        float nTrkSolidCone;
        float nTrkHollowCone;
        float isEBGap;
        float isEEGap;
        float isEBEEGap;
        float r9;
        float pt;
        float et;
        float eta;
        float phi;
        float hadronicOverEm;
        float ecalIso;
        float hcalIso;
        float trackIso;
    } ;
    struct_recPhoton recPhoton;

    float genPhoton1_pt, genPhoton1_eta, genPhoton1_phi;
    float genPhoton2_pt, genPhoton2_eta, genPhoton2_phi;
    float genPhotonDble_mass;
    float genBJetDble_Higgs_mass;
    float recoPhoton1_pt, recoPhoton1_eta, recoPhoton1_phi;
    float recoPhoton2_pt, recoPhoton2_eta, recoPhoton2_phi;
    float recoPhotonDble_mass;
    float recoBJetDble_Higgs_mass;

    float patGenBJetDble_Higgs_mass;
    float recoJetGenBJetDble_Higgs_mass;
    size_t nGenPhotons, nRecoPhotons, nGenB, nPatGenB, nrecoJetGenB;
    
  TLorentzVector genPhoton1, genPhoton2, genHiggs1, recoPhoton1, recoPhoton2, recoPhoton1_raw, recoPhoton2_raw, recoPhoton, recoHiggs, recoHiggs_raw, genB1, genB2, genHiggs2, genHH, patGenB1,patGenB2, recoJetGenB1, recoJetGenB2, patGenHiggs2, recoJetGenHiggs1, recoJetGenHiggs2;
   
    std::vector<TLorentzVector> genPho1, genPho2;

    // root file to store histograms
    //TFile*  rootFile_;

    // data members for histograms to be filled

    // PhotonID Histograms
    TH1F* h_isoEcalRecHit_;
    TH1F* h_isoHcalRecHit_;
    TH1F* h_trk_pt_solid_;
    TH1F* h_trk_pt_hollow_;
    TH1F* h_ntrk_solid_;
    TH1F* h_ntrk_hollow_;
    TH1F* h_ebgap_;
    TH1F* h_eeGap_;
    TH1F* h_ebeeGap_;
    TH1F* h_r9_;

    // Photon Histograms
    TH1F* h_photonPt_;
  //    TH1F* h_photonEt_;
    TH1F* h_photonEta_;
  //    TH1F* h_photonPhi_;
    TH1F* h_hadoverem_;
    TH1F* h_photonIetaIeta_;
    TH1F* h_phoIsoNeuHad_;
    TH1F* h_phoIsoCharHad_;
    TH1F* h_photonIso_;
  TH1F* h_puppiPhoIsoNeuHad_;
  TH1F* h_puppiPhoIsoCharHad_;
    TH1F* h_puppiPhotonIso_;


    TH1F* h_recoPhotonHiggsMass;
    TH1F* h_recoPhotonHiggsMass_raw;
    TH1F* h_recoPhotonHiggsMass_HM;
    TH1F* h_recoPhotonHiggsMass_HM_raw;
    TH1F* h_recoPhotonHiggsMass_LM;
    TH1F* h_recoPhotonHiggsMass_LM_raw;

    TH1F* h_genPhotonHiggsMass;

    TH1F* h_genHHMass;

  //    TH1F* h_matchPhotonPt;
  // TH1F* h_matchPhotonPtGen;
  //    TH1F* h_genDeltaR;
  //  TH1F* h_recoDeltaR;
  //  TProfile* tp_RecoGenEfficiencyPt;

  /*
    // Photon's SuperCluster Histograms
    TH1F* h_photonScEt_;
    TH1F* h_photonScEta_;
    TH1F* h_photonScPhi_;
    TH1F* h_photonScEtaWidth_;
  */
    // Composite or Other Histograms
    TH1F* h_photonInAnyGap_;
    TH1F* h_nPassingPho_;
    TH1F* h_nPho_;
    TH1F* h_bjetHiggsMass_;
    TH1F* h_patGenbjetHiggsMass_;
    TH1D* h_goodBJets_Higgs_n_;
    TH1D* h_goodPatBJets_Higgs_n_;

    // TProfile
    TProfile* tp_photonIso_nVtx_;
    TProfile* h_PhotonPtr;
    TProfile* h_PhotonPtr_raw;

    // TTree
    TTree* tree_PhotonAll_;
    TTree* tree_genPhotonAll_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BasicPatDistrib::BasicPatDistrib(const edm::ParameterSet& iConfig):
  useDeepCSV_(iConfig.getParameter<bool>("useDeepCSV")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),  
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
  jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
  metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
  genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  allGenPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("allGenParts"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  photonsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  minPhotonEt_(iConfig.getParameter<double>("minPhotonEt")),
  minPhotonAbsEta_(iConfig.getParameter<double>("minPhotonAbsEta")),
  maxPhotonAbsEta_(iConfig.getParameter<double>("maxPhotonAbsEta")),
  minPhotonR9_(iConfig.getParameter<double>("minPhotonR9")),
  maxPhotonHoverE_(iConfig.getParameter<double>("maxPhotonHoverE")),
  maxIetaIeta_(iConfig.getParameter<double>("maxIetaIeta"))
  //outputFile_(iConfig.getParameter<std::string>("outputFile"))
{
  //Photon PATS 
  // Read Parameters from configuration file

  // output filename
  // Read variables that must be passed to allow a
  //  supercluster to be placed in histograms as a photon.
  
  // Read variable to that decidedes whether
  // a TTree of photons is created or not
  createPhotonTTree_ = iConfig.getParameter<bool>("createPhotonTTree");
  std::cout<<"Hello User I am running!"<<std::endl;
  // open output file to store histograms
  //rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");



  //now do what ever initialization is needed
  usesResource("TFileService");

  // MC truth in fiducial phase space
  h_genJets_n_ = fs_->make<TH1D>("GenJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_genJets_pt_ = fs_->make<TH1D>("GenJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_genJets_phi_ = fs_->make<TH1D>("GenJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_genJets_eta_ = fs_->make<TH1D>("GenJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);

  // Vertices
  h_allVertices_n_ = fs_->make<TH1D>("AllVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);
  // ... that pass ID
  h_goodVertices_n_ = fs_->make<TH1D>("GoodVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);

  // Jets
  h_allJets_n_ = fs_->make<TH1D>("AllJetsN",";Jet multiplicity;Events / 1", 15, 0., 15.);
  h_allJets_pt_ = fs_->make<TH1D>("AllJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 100, 0., 200.);
  h_allJets_phi_ = fs_->make<TH1D>("AllJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_allJets_eta_ = fs_->make<TH1D>("AllJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_allJets_csv_ = fs_->make<TH1D>("AllJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_allJets_id_ = fs_->make<TH1D>("AllJetsID",";;Jets / 1", 3, 0., 3.);
  h_allJets_id_->SetOption("bar");
  h_allJets_id_->SetBarWidth(0.75);
  h_allJets_id_->SetBarOffset(0.125);
  h_allJets_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allJets_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allJets_id_->GetXaxis()->SetBinLabel(3,"Tight");
  // ... that pass kin cuts, loose ID
  h_goodJets_n_ = fs_->make<TH1D>("GoodJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_goodJets_nb_ = fs_->make<TH1D>("GoodJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodJets_pt_ = fs_->make<TH1D>("GoodJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodJets_phi_ = fs_->make<TH1D>("GoodJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodJets_eta_ = fs_->make<TH1D>("GoodJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodJets_csv_ = fs_->make<TH1D>("GoodJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_goodLJets_n_ = fs_->make<TH1D>("GoodLightJetsN",";Jet multiplicity;Events / 1", 12, 0., 12.);
  h_goodLJets_nb_ = fs_->make<TH1D>("GoodLightJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodLJets_pt_ = fs_->make<TH1D>("GoodLightJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodLJets_phi_ = fs_->make<TH1D>("GoodLightJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodLJets_eta_ = fs_->make<TH1D>("GoodLightJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodLJets_csv_ = fs_->make<TH1D>("GoodLightJetsCSV",";CSV discriminant;Events / 0.02", 50, 0., 1.);
  h_goodBJets_n_ = fs_->make<TH1D>("GoodBtaggedJetsN",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_nb_ = fs_->make<TH1D>("GoodBtaggedJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_pt_ = fs_->make<TH1D>("GoodBtaggedJetsPt",";p_{T}(jet) (GeV);Events / (5 GeV)", 36, 20., 200.);
  h_goodBJets_phi_ = fs_->make<TH1D>("GoodBtaggedJetsPhi",";#phi(jet);Events / 0.2", 30, -3., 3.);
  h_goodBJets_eta_ = fs_->make<TH1D>("GoodBtaggedJetsEta",";#eta(jet);Events / 0.2", 50, -5., 5.);
  h_goodBJets_csv_ = fs_->make<TH1D>("GoodBtaggedJetsCSV",";CSV discriminant;Events / 0.01", 20, 0.8, 1.);
  //Sam added reco jets 
  h_goodrecoJetBJets_Higgs_n_ = fs_->make<TH1D>("GoodRecoBtaggedJetsN_R2_Selection",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_recobjetHiggsMass_ = fs_->make<TH1F>("recoBJetDble_Higgs_Mass", "genBJetDble_Higgs_mass",10000,1e9,500e9);

  // MET
  h_goodMET_pt_ = fs_->make<TH1D>("GoodMETPt",";p_{T}(MET) (GeV);Events / (5 GeV)", 60, 0., 300.);
  h_goodMET_phi_ = fs_->make<TH1D>("GoodMETPhi",";#phi(MET);Events / 0.2", 30, -3., 3.);

  // PhotonID Histograms
  h_isoEcalRecHit_ = fs_->make<TH1F>("photonEcalIso",          "Ecal Rec Hit Isolation", 100, 0, 100);
  h_isoHcalRecHit_ = fs_->make<TH1F>("photonHcalIso",          "Hcal Rec Hit Isolation", 100, 0, 100);
  h_trk_pt_solid_  = fs_->make<TH1F>("photonTrackSolidIso",    "Sum of track pT in a cone of #DeltaR" , 100, 0, 100);
  h_trk_pt_hollow_ = fs_->make<TH1F>("photonTrackHollowIso",   "Sum of track pT in a hollow cone" ,     100, 0, 100);
  h_ntrk_solid_    = fs_->make<TH1F>("photonTrackCountSolid",  "Number of tracks in a cone of #DeltaR", 100, 0, 100);
  h_ntrk_hollow_   = fs_->make<TH1F>("photonTrackCountHollow", "Number of tracks in a hollow cone",     100, 0, 100);
  h_ebgap_         = fs_->make<TH1F>("photonInEBgap",          "Ecal Barrel gap flag",  2, -0.5, 1.5);
  h_eeGap_         = fs_->make<TH1F>("photonInEEgap",          "Ecal Endcap gap flag",  2, -0.5, 1.5);
  h_ebeeGap_       = fs_->make<TH1F>("photonInEEgap",          "Ecal Barrel/Endcap gap flag",  2, -0.5, 1.5);
  h_r9_            = fs_->make<TH1F>("photonR9",               "R9 = E(3x3) / E(SuperCluster)", 300, 0, 3);

  // Photon Histograms
  h_photonPt_      = fs_->make<TH1F>("photonPt",     "Photon P_{T}", 50,0., 1000.);
  h_photonEta_     = fs_->make<TH1F>("photonEta",    "Photon #eta",   200, -4,   4);

  h_hadoverem_     = fs_->make<TH1F>("photonHoverE", "Hadronic over EM", 200, 0, 1);
  h_photonIetaIeta_ = fs_->make<TH1F>("photonSigmaIetaIeta","Photon #sigma_{i#etai#eta}",1500,0.0,0.3);
  h_phoIsoNeuHad_ = fs_->make<TH1F>("photonIsolatedNeuHadron","Isolated Photon by Neutral Hadron",100,0.0,140.0);
  h_phoIsoCharHad_ = fs_->make<TH1F>("photonIsolatedCharHadron","Isolated Photon by Charged Hadron",100,0.0,250.0);
  h_photonIso_ = fs_->make<TH1F>("photonIsolated","Isolated Photon",100,0.0,200.0);
  h_puppiPhoIsoNeuHad_ = fs_->make<TH1F>("puppiPhotonIsolatedNeuHadron","Isolated Photon by Neutral Hadron PUPPI",100,0.0,140.0);
  h_puppiPhoIsoCharHad_ = fs_->make<TH1F>("puppiPhotonIsolatedCharHadron","Isolated Photon by Charged Hadron PUPPI",100,0.0,250.0);
  h_puppiPhotonIso_ = fs_->make<TH1F>("puppiPhotonIsolated","Isolated Photon PUPPI",100,0.0,200.0);

  h_recoPhotonHiggsMass = fs_->make<TH1F>("recoPhotonHiggsMass", "recoPhoton_Higgs_mass",60,110.,140.);
  h_recoPhotonHiggsMass_HM = fs_->make<TH1F>("recoPhotonHiggsMass_HM", "recoPhoton_Higgs_mass_HM",60,110.,140.);
  h_recoPhotonHiggsMass_LM = fs_->make<TH1F>("recoPhotonHiggsMass_LM", "recoPhoton_Higgs_mass_LM",60,110.,140.);

  h_recoPhotonHiggsMass_raw = fs_->make<TH1F>("recoPhotonHiggsMass_raw", "recoPhoton_Higgs_mass_raw",60,110.,140.);
  h_recoPhotonHiggsMass_HM_raw = fs_->make<TH1F>("recoPhotonHiggsMass_HM_raw", "recoPhoton_Higgs_mass_HM_raw",60,110.,140.);
  h_recoPhotonHiggsMass_LM_raw = fs_->make<TH1F>("recoPhotonHiggsMass_LM_raw", "recoPhoton_Higgs_mass_LM_raw",60,110.,140.);


  h_genPhotonHiggsMass = fs_->make<TH1F>("genPhotonHiggsMass", "genPhoton_Higgs_mass",60,110.,140.);

  h_genHHMass = fs_->make<TH1F>("genHHMass", "gen HH Mass", 75, 250., 1000.);

  h_PhotonPtr = fs_->make<TProfile>("PhotonPtr", "pTRatio",50,0.,250., 0., 2);
  h_PhotonPtr_raw = fs_->make<TProfile>("PhotonPtr_raw", "pTRatio_raw",50,0.,250., 0., 2);

  //  h_recoDeltaR= fs_->make<TH1F>("RecoDeltaR", "RecoDeltaR",1000,0.001,5.);
  //  h_genDeltaR= fs_->make<TH1F>("GenDeltaR", "GenDeltaR",1000,0.001,5.);
  //  tp_RecoGenEfficiencyPt = fs_->make<TProfile>("Efficiency", "RecoGenEfficienyPt",1000,1.,250.,0.0,3.);


  // Photon's SuperCluster Histograms
  //  h_photonScEt_       = fs_->make<TH1F>("photonScEt",  "Photon SuperCluster E_{T}", 200,  0, 200);
  //  h_photonScEta_      = fs_->make<TH1F>("photonScEta", "Photon #eta",               200, -4,   4);
  //  h_photonScPhi_      = fs_->make<TH1F>("photonScPhi", "Photon #phi", 200, -1.*TMath::Pi(), TMath::Pi());
  // h_photonScEtaWidth_ = fs_->make<TH1F>("photonScEtaWidth","#eta-width",            100,  0,  .1);

  // Composite or Other Histograms
  //  h_photonInAnyGap_   = fs_->make<TH1F>("photonInAnyGap",     "Photon in any gap flag",  2, -0.5, 1.5);
  h_nPassingPho_      = fs_->make<TH1F>("photonPassingCount", "Total number photons (0=NotPassing, 1=Passing)", 2, -0.5, 1.5);
  h_nPho_             = fs_->make<TH1F>("photonCount",        "Number of photons passing cuts in event",  10,  0,  10);
  // tp_photonIso_nVtx_ = fs_->make<TProfile>("photonsNVtx","photon and Number of Vertices",1000,0.,220.,0.,100.);



  //Bjets invariant mass distribution. 
  h_goodBJets_Higgs_n_ = fs_->make<TH1D>("GoodBtaggedJetsN_R2_Selection",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodPatBJets_Higgs_n_ = fs_->make<TH1D>("GoodBtaggedJetsN_R2_Selection",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_bjetHiggsMass_ = fs_->make<TH1F>("genBJetDble_Higgs_Mass", "genBJetDble_Higgs_mass",10000,1e9,500e9);
  h_patGenbjetHiggsMass_ = fs_->make<TH1F>("patgenBJetDble_Higgs_Mass", "patgenBJetDble_Higgs_mass",10000,1e9,500e9);



  // Create a TTree of photons if set to 'True' in config file
  if ( createPhotonTTree_ ) {
    tree_PhotonAll_     = fs_->make<TTree>("TreePhotonAll", "Reconstructed Photon");
    tree_PhotonAll_->Branch("recPhoton", &recPhoton.isolationEcalRecHit, "isolationEcalRecHit/F:isolationHcalRecHit:isolationSolidTrkCone:isolationHollowTrkCone:nTrkSolidCone:nTrkHollowCone:isEBGap:isEEGap:isEBEEGap:r9:pt:et:eta:phi:hadronicOverEm:ecalIso:hcalIso:trackIso");


    tree_genPhotonAll_     = fs_->make<TTree>("TreeGenPhotonAll", "Generated Photon from Higgs");
    tree_genPhotonAll_->Branch("genPhoton1_pt", &genPhoton1_pt, "genPhoton1_pt/F");
    tree_genPhotonAll_->Branch("genPhoton1_eta", &genPhoton1_eta, "genPhoton1_eta/F");
    tree_genPhotonAll_->Branch("genPhoton1_phi", &genPhoton1_phi, "genPhoton1_phi/F");

    tree_genPhotonAll_->Branch("genPhoton2_pt", &genPhoton2_pt, "genPhoton2_pt/F");
    tree_genPhotonAll_->Branch("genPhoton2_eta", &genPhoton2_eta, "genPhoton2_eta/F");
    tree_genPhotonAll_->Branch("genPhoton2_phi", &genPhoton2_phi, "genPhoton2_phi/F");
 

    tree_genPhotonAll_->Branch("nGenPhotons", &nGenPhotons, "nGenPhotons/I");
    tree_genPhotonAll_->Branch("nGenB", &nGenB, "nGenB/I");
    
    tree_genPhotonAll_->Branch("genPhotonDouble_Mass", &genPhotonDble_mass, "genPhotonDble_mass/I");
  }
}


BasicPatDistrib::~BasicPatDistrib()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //delete rootFile_;

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
BasicPatDistrib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //scope   
  using namespace edm;
  using namespace std;
  //Gathering Physics Objects 
  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);

  Handle<std::vector<pat::Electron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  //  Handle<reco::ConversionCollection> conversions;
  //  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();  

  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<std::vector<pat::MET>> mets;
  iEvent.getByToken(metsToken_, mets);

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<pat::PackedGenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenParticle>> allGenParts;//Photon H mass is computed
  iEvent.getByToken(allGenPartsToken_, allGenParts);


  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  //Handle< View<pat::Photon> >  photonHandle;
  //iEvent.getByLabel("slimmedPhotons", photonHandle);
  //View<pat::Photon> photons = *photonHandle; 
  
  Handle<std::vector<pat::Photon>>  photonHandle;
  iEvent.getByToken(photonsToken_, photonHandle);
  std::vector<pat::Photon> photons = *photonHandle; 

  genPhoton1.SetPtEtaPhiM(0,0,0,0);
  genPhoton2.SetPtEtaPhiM(0,0,0,0);

  genB1.SetPtEtaPhiM(0,0,0,0);
  genB2.SetPtEtaPhiM(0,0,0,0);

  genHiggs1.SetPtEtaPhiM(0,0,0,0);
  genHiggs2.SetPtEtaPhiM(0,0,0,0);

  genHH.SetPtEtaPhiM(0,0,0,0);

  genPhoton1_pt=0, genPhoton1_eta=0, genPhoton1_phi=0;
  genPhoton2_pt=0, genPhoton2_eta=0, genPhoton2_phi=0;
  //genHiggs_pt=0,genHiggs_eta=0,genHiggs_phi=0;
  genPhotonDble_mass=0;
  genBJetDble_Higgs_mass=0;
  patGenBJetDble_Higgs_mass=0;
  nGenPhotons = 0, nGenB = 0;
  nRecoPhotons= 0;

  // Vertices
  int prVtx = -1;
  size_t nVtx = 0;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
    ++ nVtx;
  }
  if (prVtx < 0) return;
  h_goodVertices_n_->Fill(nVtx);
  h_allVertices_n_->Fill(vertices->size());
   
  // MC truth in fiducial phase space
  std::vector<size_t> jGenJets;
  size_t nGenJets = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < genParts->size(); j++) {
      if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
      if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    jGenJets.push_back(i);

    if (genJets->at(i).pt() < 30.) continue;
    if (fabs(genJets->at(i).eta()) > 4.7) continue;
    h_genJets_pt_->Fill(genJets->at(i).pt());
    h_genJets_phi_->Fill(genJets->at(i).phi());
    h_genJets_eta_->Fill(genJets->at(i).eta());
    ++nGenJets;
  }
  h_genJets_n_->Fill(nGenJets);
  

  for (size_t i = 0; i < allGenParts->size(); i++) {

    /*
    if (allGenParts->at(i).numberOfMothers() == 0) cout << i << " \t \t" << allGenParts->at(i).pdgId() << " st " << allGenParts->at(i).status()<< endl;
    if (allGenParts->at(i).numberOfMothers() == 1) cout << i << " " <<  allGenParts->at(i).mother(0)->pdgId() << " \t " << allGenParts->at(i).pdgId() << " st " << allGenParts->at(i).status() << endl;
    if (allGenParts->at(i).numberOfMothers() == 2) cout << i << " " << allGenParts->at(i).mother(0)->pdgId() << " " << allGenParts->at(i).mother(1)->pdgId() << " " << allGenParts->at(i).pdgId() << " st " << allGenParts->at(i).status() << endl;
    */    
    //    if (fabs(allGenParts->at(i).eta()) >  1.479) continue;
    size_t nMothers = allGenParts->at(i).numberOfMothers();
    /* 
   if (nMothers > 1) {
      cout << "=================================== PDG ID = " << allGenParts->at(i).mother(0)->pdgId() << endl;
      cout << "=================================== PDG ID = " << allGenParts->at(i).mother(1)->pdgId() << endl;
      cout << "" << endl << endl;
    }
    */



    //GenPhotons.... 

    if (nMothers==1){
      if (allGenParts->at(i).mother(0)->pdgId() != 25) continue;
      //cout << "PDG ID = " << allGenParts->at(i).mother(0)->pdgId() << endl;
      if (allGenParts->at(i).pdgId() == 22) {
	nGenPhotons++;
	if (nGenPhotons == 1){
	  genPhoton1_pt = allGenParts->at(i).pt();
	  genPhoton1_eta = allGenParts->at(i).eta();
	  genPhoton1_phi = allGenParts->at(i).phi();
	  genPhoton1.SetPtEtaPhiM(genPhoton1_pt, genPhoton1_eta, genPhoton1_phi, 0);
	  //genPho1.push_back(genPhoton1);
	} else if (nGenPhotons == 2){
	  genPhoton2_pt = allGenParts->at(i).pt();
	  genPhoton2_eta = allGenParts->at(i).eta();
	  genPhoton2_phi = allGenParts->at(i).phi();
	  genPhoton2.SetPtEtaPhiM(genPhoton2_pt, genPhoton2_eta, genPhoton2_phi, 0);
	  genHiggs1 = genPhoton1+genPhoton2;
	  //genPho2.push_back(genPhoton2);
	  genPhotonDble_mass = genHiggs1.M();
	  h_genPhotonHiggsMass->Fill(genPhotonDble_mass);
	  //std::cout<<"The Gen Higgs Mass:  "<<genPhotonDble_mass<<std::endl;
      //genPhotonDble_mass = TMath::Sqrt(2*genPhoton1_pt*genPhoton2_pt*(TMath::CosH(genPhoton1_eta-genPhoton2_eta)-TMath::Cos(genPhoton1_phi-genPhoton2_phi)));
	}
      } else if (abs(allGenParts->at(i).pdgId()) == 5) {
      //B-Jets
	  nGenB++;
      //double btagDisc = -1.;
      if (useDeepCSV_){
          //btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
              //jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
        }
      else{
           //btagDisc = jets->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        }
	  if (nGenB == 1){
      if (useDeepCSV_){ //&& btagDisc > 0.131)
            //|| (!useDeepCSV_ && btagDisc > 0.8484))  
        genB1.SetPtEtaPhiM(allGenParts->at(i).pt(),allGenParts->at(i).eta(), allGenParts->at(i).phi(), allGenParts->at(i).mass());
            }
      }
	  else if (nGenB == 2){
        genB2.SetPtEtaPhiM(allGenParts->at(i).pt(),allGenParts->at(i).eta(), allGenParts->at(i).phi(), allGenParts->at(i).mass());
        genHiggs2 = genB1+genB2;
        genBJetDble_Higgs_mass = genHiggs2.M();
        h_bjetHiggsMass_->Fill(genBJetDble_Higgs_mass);
	//        std::cout<<"Found a Higgs from Bjet!"<<std::endl;
        }
      }      
    //End 1Mother
    }
  //End loop over allGen particles 
  }

  genHH = genHiggs1+genHiggs2;

  h_genHHMass->Fill(genHH.M());
  
  h_goodBJets_Higgs_n_->Fill(nGenB);

  //Gen Particles Histograms 
  if (createPhotonTTree_)
    tree_genPhotonAll_->Fill();
    

  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (fabs(genParts->at(i).eta()) > 2.8) continue;
    double genIso = 0.;
    for (size_t j = 0; j < jGenJets.size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
      std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
      for (size_t k = 0; k < jconst.size(); k++) {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
        if (deltaR < 0.01) continue;
        if (abs(genParts->at(i).pdgId()) == 13 && deltaR > 0.4) continue;
        if (abs(genParts->at(i).pdgId()) == 11 && deltaR > 0.3) continue;
        genIso = genIso + jconst[k]->pt();
      }
    }
    genIso = genIso / genParts->at(i).pt();
    if (genIso > 0.15) continue;
    if (abs(genParts->at(i).pdgId()) == 13) {
    }
    if (abs(genParts->at(i).pdgId()) == 11) {
    }
    else if (abs(genParts->at(i).pdgId()) == 5) {
      //B-Jets
	  nPatGenB++;
      //double btagDisc = -1.;
      if (useDeepCSV_){
          //btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
              //jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
        }
      else{
           //btagDisc = jets->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        }
	  if (nPatGenB == 1){
      if (useDeepCSV_ )//&& btagDisc > 0.131) //|| (!useDeepCSV_ && btagDisc > 0.8484))  //0.6324
        {
        patGenB1.SetPtEtaPhiM(genJets->at(i).pt(),genJets->at(i).eta(), genJets->at(i).phi(), genJets->at(i).mass());
        }
      }
	  else if (nGenB == 2){
        patGenB2.SetPtEtaPhiM(genJets->at(i).pt(),genJets->at(i).eta(), genJets->at(i).phi(), genJets->at(i).mass());
        patGenHiggs2 = patGenB1+patGenB2;
        patGenBJetDble_Higgs_mass = patGenHiggs2.M();
        h_patGenbjetHiggsMass_->Fill(patGenBJetDble_Higgs_mass);

        }
      }

  }

  h_goodPatBJets_Higgs_n_->Fill(nPatGenB);

  
  // Jets
  size_t nGoodJets = 0;
  size_t nbGoodJets = 0;
  size_t nGoodLightJets = 0;
  size_t nbGoodLightJets = 0;
  size_t nGoodBtaggedJets = 0;
  size_t nbGoodBtaggedJets = 0;
  for (size_t i =0; i < jets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < elecs->size(); j++) {
      if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    for (size_t j = 0; j < muons->size(); j++) {
      if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;

    double btagDisc = -1.;
    if (useDeepCSV_)
        btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    else
        btagDisc = jets->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    
    h_allJets_pt_->Fill(jets->at(i).pt());
    h_allJets_phi_->Fill(jets->at(i).phi());
    h_allJets_eta_->Fill(jets->at(i).eta());
    h_allJets_csv_->Fill(btagDisc); 
    h_allJets_id_->Fill(0.);
    pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    retLoose.set(false);
    if (jetIDLoose_(jets->at(i), retLoose)) h_allJets_id_->Fill(1.);
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    if (jetIDTight_(jets->at(i), retTight)) h_allJets_id_->Fill(2.);

    if (jets->at(i).pt() < 30.) continue;
    if (fabs(jets->at(i).eta()) > 4.7) continue;
    if (!jetIDLoose_(jets->at(i), retLoose)) continue;
    h_goodJets_pt_->Fill(jets->at(i).pt());
    h_goodJets_phi_->Fill(jets->at(i).phi());
    h_goodJets_eta_->Fill(jets->at(i).eta());
    h_goodJets_csv_->Fill(btagDisc); 
    ++nGoodJets;
    //size_t nJetMother = jets->at(i).numberOfMothers();
    if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodJets;
    if ((useDeepCSV_ && btagDisc > 0.6324)
            || (!useDeepCSV_ && btagDisc > 0.8484)){  
      h_goodBJets_pt_->Fill(jets->at(i).pt());
      h_goodBJets_phi_->Fill(jets->at(i).phi());
      h_goodBJets_eta_->Fill(jets->at(i).eta());
      h_goodBJets_csv_->Fill(btagDisc);
      ++nGoodBtaggedJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodBtaggedJets;
    } else {
      h_goodLJets_pt_->Fill(jets->at(i).pt());
      h_goodLJets_phi_->Fill(jets->at(i).phi());
      h_goodLJets_eta_->Fill(jets->at(i).eta());
      h_goodLJets_csv_->Fill(btagDisc); 
      ++nGoodLightJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodLightJets;
    }
//Added Sam code
   // if (abs(jets->at(i).genParton()->pdgId()) == 5) {
   //   //B-Jets
   //   nrecoJetGenB++;
   //   //double btagDisc = -1.;
   //   if (useDeepCSV_){
   //       //btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
   //           //jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
   //     }
   //   else{
   //        //btagDisc = jets->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
   //     }
   //   if (nrecoJetGenB == 1){
   //   if (useDeepCSV_ )//&& btagDisc > 0.131) //|| (!useDeepCSV_ && btagDisc > 0.8484))  //0.6324
   //     {
   //     recoJetGenB1.SetPtEtaPhiM(genJets->at(i).pt(),genJets->at(i).eta(), genJets->at(i).phi(), genJets->at(i).mass());
   //     }
   //   }
   //   else if (nrecoJetGenB == 2){
   //     recoJetGenB2.SetPtEtaPhiM(genJets->at(i).pt(),genJets->at(i).eta(), genJets->at(i).phi(), genJets->at(i).mass());
   //     recoJetGenHiggs2 = recoJetGenB1+recoJetGenB2;
   //     recoJetGenBJetDble_Higgs_mass = recoJetGenHiggs2.M();
   //     h_recobjetHiggsMass_->Fill(recoJetGenBJetDble_Higgs_mass);

   //     }
   //   }


  }
  h_goodrecoJetBJets_Higgs_n_->Fill(nrecoJetGenB);
  h_goodLJets_n_->Fill(nGoodLightJets);
  h_goodLJets_nb_->Fill(nbGoodLightJets);
  h_goodBJets_n_->Fill(nGoodBtaggedJets);
  h_goodBJets_nb_->Fill(nbGoodBtaggedJets);
  h_goodJets_n_->Fill(nGoodJets);
  h_goodJets_nb_->Fill(nbGoodJets);
  h_allJets_n_->Fill(jets->size());
  
  // MET
  if (mets->size() > 0) {
    h_goodMET_pt_->Fill(mets->at(0).pt());
    h_goodMET_phi_->Fill(mets->at(0).phi());
  }

  // Photons
  int photonCounter = 0;

  float Rpt1 = 0, Rpt2 = 0;
  float Rpt1_raw = 0, Rpt2_raw = 0;

  recoPhoton1.SetPtEtaPhiM(0, 0, 0, 0);
  recoPhoton2.SetPtEtaPhiM(0, 0, 0, 0);
  recoPhoton1_raw.SetPtEtaPhiM(0, 0, 0, 0);
  recoPhoton2_raw.SetPtEtaPhiM(0, 0, 0, 0);
  
  for (int i=0; i<int(photons.size()); i++)
  {

    pat::Photon currentPhoton = photons.at(i);

    float photonEt       = currentPhoton.et();
    float photonPt       = currentPhoton.pt();
    float superClusterEt = (currentPhoton.superCluster()->energy())/(cosh(currentPhoton.superCluster()->position().eta()));
    float superClusterEta = currentPhoton.superCluster()->position().eta();
    float superClusterPhi = currentPhoton.superCluster()->position().phi();
    float phoIsoNeuHad = currentPhoton.neutralHadronIso();
    float phoIsoCharHad = currentPhoton.chargedHadronIso();
    float photonIso = currentPhoton.photonIso();
    float puppiPhoIsoNeuHad = currentPhoton.puppiNeutralHadronIso();
    float puppiPhoIsoCharHad = currentPhoton.puppiChargedHadronIso();
    float puppiPhotonIso = currentPhoton.puppiPhotonIso();




    // Only store photon candidates (SuperClusters) that pass some simple cuts
    bool passCuts = (              photonEt > minPhotonEt_ &&
       //             (      fabs(currentPhoton.eta()) > minPhotonAbsEta_ ) &&
				   fabs(currentPhoton.eta()) < maxPhotonAbsEta_);

      //         (          currentPhoton.r9() > minPhotonR9_     ) &&
      //                  ( currentPhoton.hadronicOverEm() < maxPhotonHoverE_ ) &&
      //               ( currentPhoton.sigmaIetaIeta() < maxIetaIeta_ ) &&
      //                  ( phoIsoCharHad < 1.295) &&
      //                 ( phoIsoNeuHad < 10.910+0.0148*photonPt+0.000017*TMath::Power(photonPt,2)) &&
      //                   ( photonIso <  3.630+0.0047*photonPt ) ;
                       

    if ( passCuts )
    {
      ///////////////////////////////////////////////////////
      //                fill histograms                    //
      ///////////////////////////////////////////////////////
      // PhotonID Variables
      /*
      h_isoEcalRecHit_->Fill(currentPhoton.ecalRecHitSumEtConeDR04());
      h_isoHcalRecHit_->Fill(currentPhoton.hcalTowerSumEtConeDR04());
      h_trk_pt_solid_ ->Fill(currentPhoton.trkSumPtSolidConeDR04());
      h_trk_pt_hollow_->Fill(currentPhoton.trkSumPtHollowConeDR04());
      h_ntrk_solid_->   Fill(currentPhoton.nTrkSolidConeDR04());
      h_ntrk_hollow_->  Fill(currentPhoton.nTrkHollowConeDR04());
      h_ebgap_->        Fill(currentPhoton.isEBGap());
      h_eeGap_->        Fill(currentPhoton.isEEGap());
      h_ebeeGap_->      Fill(currentPhoton.isEBEEGap());
      */

      //      tp_photonIso_nVtx_->Fill(photonIso,nVtx,1);

      // It passed photon cuts, mark it
      h_nPassingPho_->Fill(1.0);
      //      if(photons.size()>1){

      recoPhoton.SetPtEtaPhiM(photonPt,currentPhoton.eta(),currentPhoton.phi(),0.0);  
      bool matched = false;

      if (genPhoton1.DeltaR(recoPhoton) < 0.2) {
	matched = true;
	recoPhoton1.SetPtEtaPhiM(photonPt,currentPhoton.eta(),currentPhoton.phi(),0.0);  
	recoPhoton1_raw.SetPtEtaPhiM(superClusterEt,superClusterEta,superClusterPhi,0.0);  
	Rpt1 =  recoPhoton1.Pt()/genPhoton1.Pt();
	Rpt1_raw =  recoPhoton1_raw.Pt()/genPhoton1.Pt();
      }
      else if (genPhoton2.DeltaR(recoPhoton) < 0.2) {
	matched = true;
	recoPhoton2.SetPtEtaPhiM(photonPt,currentPhoton.eta(),currentPhoton.phi(),0.0);  
	recoPhoton2_raw.SetPtEtaPhiM(superClusterEt,superClusterEta,superClusterPhi,0.0);  
	Rpt2 =  recoPhoton2.Pt()/genPhoton2.Pt();
	Rpt2_raw =  recoPhoton2_raw.Pt()/genPhoton2.Pt();
      }

      if (matched){

	h_hadoverem_-> Fill(currentPhoton.hadronicOverEm());
	h_photonIetaIeta_->Fill(currentPhoton.sigmaIetaIeta());
	h_phoIsoNeuHad_->Fill(phoIsoNeuHad);
	h_phoIsoCharHad_->Fill(phoIsoCharHad);
	h_photonIso_->Fill(photonIso);
	h_puppiPhoIsoNeuHad_->Fill(puppiPhoIsoNeuHad);
	h_puppiPhoIsoCharHad_->Fill(puppiPhoIsoCharHad);
	h_puppiPhotonIso_->Fill(puppiPhotonIso);

	h_photonPt_->  Fill(currentPhoton.pt());
	h_photonEta_->  Fill(currentPhoton.eta());
	h_r9_->           Fill(currentPhoton.r9());


      }




      /*
            if((i+1)==int(photons.size())){
            //std::cout<<"The photon 1  Pt is: "<<firstPt<<" 2 photon "<<secondPt<<std::endl;
            recoHiggs = recoPhoton1 + recoPhoton2;
            h_recoPhotonHiggsMass->Fill(recoHiggs.M());
            //std::cout<<"RecoHiggs mass: "<<recoHiggs.M()<<std::endl;
            firstPt=0.0;
            secondPt=0.0;

            //Matching Delta R for efficiency and pt histograms. 
            for(int j=0; j<int(genPho1.size());j++){
            //std::cout<<"Reco DelR: "<<recoPhoton1.DeltaR(recoPhoton2)<<" Gen DelR: "<<genPho1[j].DeltaR(genPho2[j]);
            // Gen DelR: 2.96679Reco DelR: 2.86304 Gen DelR: 1.57869Reco DelR: 2.86304
            if(recoPhoton1.DeltaR(recoPhoton2) > genPho1[j].DeltaR(genPho2[j])-0.001 &&
                  recoPhoton1.DeltaR(recoPhoton2)< genPho1[j].DeltaR(genPho2[j])+0.001  ){
                h_genDeltaR->Fill(genPho1[j].DeltaR(genPho2[j]));
                h_recoDeltaR->Fill(recoPhoton1.DeltaR(recoPhoton2));
                h_matchPhotonPt->Fill(recoPhoton1.Pt()+recoPhoton2.Pt());
                h_matchPhotonPtGen->Fill(genPho1[j].Pt()+genPho2[j].Pt());
                tp_RecoGenEfficiencyPt->Fill(genPho1[j].Pt()+genPho2[j].Pt(),(recoPhoton1.Pt()+recoPhoton2.Pt())/(genPho1[j].Pt()+genPho2[j].Pt()));
            }


            //recoPhoton1 recoPhoton2;
            //genPhoton1 genPhoton2;
            }

            }

      */
    }
  }

  
  if (recoPhoton1.Pt() > 0 && recoPhoton2.Pt() > 0 && genPhoton1.Pt() > 0 && genPhoton2.Pt() > 0){
    /*
    cout << "recoPhoton1.Pt() = " << recoPhoton1.Pt() << " recoPhoton2.Pt() = " << recoPhoton2.Pt() << endl;
    cout << "genPhoton1.Pt() = " << genPhoton1.Pt() << " genPhoton2.Pt() = " << genPhoton2.Pt() << endl;

    cout << "recoPhoton1.Etat() = " << recoPhoton1.Eta() << " recoPhoton2.Eta() = " << recoPhoton2.Eta() << endl;
    cout << "genPhoton1.Eta() = " << genPhoton1.Eta() << " genPhoton2.Eta() = " << genPhoton2.Eta() << endl;

    cout << "recoPhoton1.Phit() = " << recoPhoton1.Phi() << " recoPhoton2.Phi() = " << recoPhoton2.Phi() << endl;
    cout << "genPhoton1.Phi() = " << genPhoton1.Phi() << " genPhoton2.Phi() = " << genPhoton2.Phi() << endl;


    cout << "genPhoton1.DeltaR(recoPhoton1) = " << genPhoton1.DeltaR(recoPhoton1) << " genPhoton1.DeltaR(recoPhoton2) = " << genPhoton1.DeltaR(recoPhoton2) << endl;
    cout << "genPhoton2.DeltaR(recoPhoton1) = " << genPhoton2.DeltaR(recoPhoton1) << " genPhoton2.DeltaR(recoPhoton2) = " << genPhoton2.DeltaR(recoPhoton2) << endl;
    */
    /*
      if (genPhoton1.DeltaR(recoPhoton1) < 0.2) Rpt1 =  recoPhoton1.Pt()/genPhoton1.Pt();
      else if (genPhoton1.DeltaR(recoPhoton2) < 0.2) Rpt1 =  recoPhoton2.Pt()/genPhoton1.Pt();

      if (genPhoton2.DeltaR(recoPhoton1) < 0.2) Rpt2 =  recoPhoton1.Pt()/genPhoton2.Pt();
      else if (genPhoton2.DeltaR(recoPhoton2) < 0.2) Rpt2 =  recoPhoton2.Pt()/genPhoton2.Pt();
    */
    //    cout << "Rpt1 = " << Rpt1 << " Rpt2 = " << Rpt2 << endl;      

      if (Rpt1 > 0.01 && Rpt2 > 0.01){
	recoHiggs = recoPhoton1 + recoPhoton2;
	recoHiggs_raw = recoPhoton1_raw + recoPhoton2_raw;

	h_PhotonPtr->Fill(genPhoton1.Pt(), Rpt1);
	h_PhotonPtr->Fill(genPhoton2.Pt(), Rpt2);
	h_PhotonPtr_raw->Fill(genPhoton1.Pt(), Rpt1_raw);
	h_PhotonPtr_raw->Fill(genPhoton2.Pt(), Rpt2_raw);

	double pT1 = recoPhoton1.Pt(), pT2 = recoPhoton2.Pt(), mgg = recoHiggs.M(), mHH = genHH.M();
	if (pT1 < pT2) pT2 = recoPhoton1.Pt(), pT1 = recoPhoton2.Pt();

	if (pT1 > 30 && pT2 > 20 && pT1 > mgg/3 && pT2 > mgg/4) {
	  h_recoPhotonHiggsMass->Fill(recoHiggs.M());
	  if(mHH>350) h_recoPhotonHiggsMass_HM->Fill(recoHiggs.M());
	  else if (mHH > 250 && mHH < 350)h_recoPhotonHiggsMass_LM->Fill(recoHiggs.M());

	  h_recoPhotonHiggsMass_raw->Fill(recoHiggs_raw.M());
	  if(mHH>350) h_recoPhotonHiggsMass_HM_raw->Fill(recoHiggs_raw.M());
	  else if (mHH > 250 && mHH < 350) h_recoPhotonHiggsMass_LM_raw->Fill(recoHiggs_raw.M());


	}
      // Photon Variables



      }
      
    }
    /* 
      ///////////////////////////////////////////////////////
      //                fill TTree (optional)              //
      ///////////////////////////////////////////////////////
      if ( createPhotonTTree_ ) {
        recPhoton.isolationEcalRecHit    = currentPhoton.ecalRecHitSumEtConeDR04();
        recPhoton.isolationHcalRecHit    = currentPhoton.hcalTowerSumEtConeDR04();
        recPhoton.isolationSolidTrkCone  = currentPhoton.trkSumPtSolidConeDR04();
        recPhoton.isolationHollowTrkCone = currentPhoton.trkSumPtHollowConeDR04();
        recPhoton.nTrkSolidCone          = currentPhoton.nTrkSolidConeDR04();
        recPhoton.nTrkHollowCone         = currentPhoton.nTrkHollowConeDR04();
        recPhoton.isEBGap                = currentPhoton.isEBGap();
        recPhoton.isEEGap                = currentPhoton.isEEGap();
        recPhoton.isEBEEGap              = currentPhoton.isEBEEGap();
        recPhoton.r9                     = currentPhoton.r9();
        recPhoton.pt                     = currentPhoton.pt();
        recPhoton.et                     = currentPhoton.et();
        recPhoton.eta                    = currentPhoton.eta();
        recPhoton.phi                    = currentPhoton.phi();
        recPhoton.hadronicOverEm         = currentPhoton.hadronicOverEm();
        recPhoton.ecalIso                = currentPhoton.ecalIso();
        recPhoton.hcalIso                = currentPhoton.hcalIso();
        recPhoton.trackIso               = currentPhoton.trackIso();

        // Fill the tree (this records all the recPhoton.* since
        // tree_PhotonAll_ was set to point at that.
        tree_PhotonAll_->Fill();
      }

      // Record whether it was near any module gap.
      // Very convoluted at the moment.
      bool inAnyGap = currentPhoton.isEBEEGap() || (currentPhoton.isEB()&&currentPhoton.isEBGap()) || (currentPhoton.isEE()&&currentPhoton.isEEGap());
      if (inAnyGap) {
        h_photonInAnyGap_->Fill(1.0);
      } else {
        h_photonInAnyGap_->Fill(0.0);
      }

      photonCounter++;
    }
    else
    {
      // This didn't pass photon cuts, mark it
      h_nPassingPho_->Fill(0.0);
    }

} // End Loop over photons
    */


  //Setting the total number of photons
  h_nPho_->Fill(photonCounter);



  //Setting the TProfile parameters
  //  tp_photonIso_nVtx_->GetXaxis()->SetTitle("# of Isolated Photons");
  //  tp_photonIso_nVtx_->GetYaxis()->SetTitle("# of Verticies");

//  tp_RecoGenEfficiencyPt->GetXaxis()->SetTitle("Pt_{Gen} [GeV]");
//  tp_RecoGenEfficiencyPt->GetYaxis()->SetTitle("Pt_{Reco}/Pt_{Gen} match to #Delta R");

}


// ------------ method check that an e passes loose id ----------------------------------
  bool BasicPatDistrib::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
  if (patEl.hcalOverEcal() > 6.741) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 73.76) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
  bool
BasicPatDistrib::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
  if (patEl.hcalOverEcal() > 7.371) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 22.6) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
  bool
BasicPatDistrib::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
  if (patEl.hcalOverEcal() > 4.492) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 18.26) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method to improve ME0 muon ID ----------------
  bool 
BasicPatDistrib::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = fabs(chamber->x - segment->x);
          deltaY   = fabs(chamber->y - segment->y);
          pullX    = fabs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = fabs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = fabs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

// ------------ method called once each job just before starting event loop  ------------
  void 
BasicPatDistrib::beginJob()
{
  // go to *OUR* rootfile
  //rootFile_->cd();

  // Book Histograms
  

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
BasicPatDistrib::endJob() 
{

  // go to *OUR* root file and store histograms
  ////rootFile_->cd();

  //// PhotonID Histograms
  //h_isoEcalRecHit_->Write();
  //h_isoHcalRecHit_->Write();
  //h_trk_pt_solid_-> Write();
  //h_trk_pt_hollow_->Write();
  //h_ntrk_solid_->   Write();
  //h_ntrk_hollow_->  Write();
  //h_ebgap_->     Write();
  //h_eeGap_->     Write();
  //h_ebeeGap_->   Write();
  //h_r9_->        Write();

  //// Photon Histograms
  //h_photonEt_->  Write();
  //h_photonEta_-> Write();
  //h_photonPhi_-> Write();
  //h_hadoverem_-> Write();

  //// Photon's SuperCluster Histograms
  //h_photonScEt_->      Write();
  //h_photonScEta_->     Write();
  //h_photonScPhi_->     Write();
  //h_photonScEtaWidth_->Write();

  //// Composite or Other Histograms
  //h_photonInAnyGap_->Write();
  //h_nPassingPho_->   Write();
  //h_nPho_->          Write();

  //// Write the root file (really writes the TTree)
  //rootFile_->Write();
  //rootFile_->Close();


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BasicPatDistrib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicPatDistrib);
