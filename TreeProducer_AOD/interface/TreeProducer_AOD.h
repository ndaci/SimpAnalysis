// -*- C++ -*-
//
// Package:    TreeProducer_AOD
// Class:      TreeProducer_AOD
// 
/**\class TreeProducer_AOD TreeProducer_AOD.cc SimpAnalysis/TreeProducer_AOD/src/TreeProducer_AOD.cc

 Description: EDAnalyzer produce flat trees from AOD for SimpAnalysis

*/

// C++ lib
#include <memory>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

// CMSSW standard lib
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// others
using namespace std;
int verbose=1;
const UInt_t nJ=8; // #jets
const UInt_t nGJ=4; // #genjets
const UInt_t nV=10; // #vertices
const UInt_t nT=10; // #tracks
const UInt_t nP=4; // #photons

//
// class declaration
//

class TreeProducer_AOD : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {
 public:
  explicit TreeProducer_AOD(const edm::ParameterSet&);
  ~TreeProducer_AOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool ptSorter(const reco::Track & i, const reco::Track & j);

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void Init();
	
  // ----------member data ---------------------------
  edm::InputTag _trigResultsTag;
//   edm::InputTag _trigResultsTag1;
//   edm::InputTag _trigResultsTag2;
//   edm::InputTag _prescalesTag;
//   edm::InputTag _prescalesTag2;
  edm::InputTag _METfilterTag;
  edm::InputTag _pfjetCollectionTag;
  edm::InputTag _genjetCollectionTag;
  edm::InputTag _vertexCollectionTag;
  edm::InputTag _METCollectionTag;
  edm::InputTag _photonCollectionTag;
  edm::InputTag _trackCollectionTag;
	edm::InputTag _phoLooseIdMapTag;
	edm::InputTag _phoMediumIdMapTag;
	edm::InputTag _phoTightIdMapTag;
  edm::InputTag _srcPFRhoTag;
  edm::InputTag _tagBS;
  edm::InputTag _tagCONV;
  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
//   edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken2;
//   edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken;
//   edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken2;
  edm::EDGetTokenT<edm::TriggerResults> _METfilterToken;
  edm::EDGetTokenT<vector<reco::PFJet> > _pfjetCollectionToken;
  edm::EDGetTokenT<vector<reco::GenJet> > _genjetCollectionToken;
  edm::EDGetTokenT<vector<reco::Vertex> > _vertexCollectionToken;
  edm::EDGetTokenT<vector<reco::PFMET> > _METCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Photon> > _photonCollectionToken;
  edm::EDGetTokenT<vector<reco::Track> > _trackCollectionToken;
	edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
	edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_;
	edm::EDGetTokenT<double> _pfRhoToken;
  
  bool _isData;
  
  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescale_;
  
  edm::EDGetTokenT<reco::ConversionCollection> tokenCOV;
  edm::EDGetTokenT<reco::BeamSpot> tokenBS;

  std::vector<std::string> triggerNames_;
  std::vector<unsigned int> triggerIndex_;
  
//   GlobalPoint vertexPosition;
  
  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nJet, _nJet_stored, _nGenJet, _nGenJet_stored, _nTrack, _nTrack_stored;

  // Vertices
  int _vtx_N, _vtx_N_stored, _vtx_ndof[nV], _vtx_nTracks[nV];
  double _vtx_x[nV], _vtx_y[nV], _vtx_z[nV];
  double _vtx_normalizedChi2[nV], _vtx_d0[nV];
	
	//Tracks
	int _track_fromPV[nT], _track_Nhits[nT], _track_NpixHits[nT], _track_purity[nT], _track_ndof[nT];
	double _track_eta[nT], _track_pt[nT], _track_phi[nT], _track_ptError[nT], _track_dxy[nT], _track_d0[nT], _track_dzError[nT], _track_dz[nT], _track_normalizedChi2[nT];	

  // Jets
  int _jet_mult_ch[nJ], _jet_mult_mu[nJ], _jet_mult_ne[nJ]; // multiplicities
  double _jet_vx[nJ], _jet_vy[nJ], _jet_vz[nJ];//vertex position
  double _jet_area[nJ];
  double _jet_eta[nJ], _jet_phi[nJ], _jet_pt[nJ], _jet_e[nJ], _jet_m[nJ];
  double _jet_efrac_ne_Had[nJ], _jet_efrac_ne_EM[nJ]; // neutral energy fractions
  double _jet_efrac_ch_Had[nJ], _jet_efrac_ch_EM[nJ], _jet_efrac_ch_Mu[nJ]; // charged energy fractions
  double _jet_efrac_photon[nJ];
  double _jet_unc[nJ], _jet_ptCor_up[nJ], _jet_ptCor_down[nJ]; //JEC uncertainties
  
  // GenJets
  double _genjet_vx[nGJ], _genjet_vy[nGJ], _genjet_vz[nGJ];//vertex position
  double _genjet_area[nGJ];
  double _genjet_eta[nGJ], _genjet_phi[nGJ], _genjet_pt[nGJ], _genjet_e[nGJ], _genjet_m[nGJ];
  double _genjet_efrac_ch[nGJ];// charged energy fractions
  
  //Photons
  int _nPhoton_stored;
  double _photon_pt[nP], _photon_eta[nP], _photon_phi[nP];
	int _passLooseId[nP], _passMediumId[nP], _passTightId[nP];
  int pc_matched[nP];
  double convtracks_pt[nP];

  // MET
  double _MET, _MET_phi;
  
  //Trigger
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  int _dijet_170_0p1, _dijet_220_0p3, _dijet_330_0p5, _dijet_430, _dijet_170, _singlejet_170_0p1;
  //prescales
  double  _pswgt_dijet_170,  _pswgt_singlejet_170_0p1;
  
  //MET filters
  std::vector<std::string>   filterPathsVector;
  std::map<std::string, int> filterPathsMap;
  int _HBHENoiseFlag, _HBHENoiseIsoFlag, _ECALFlag, _vertexFlag, _eeFlag, _beamhaloFlag;
  
//PFRho
  double _pfrho;
};

namespace reco {
	template<typename T> 
	class RecoPtSorter{
	public:
		bool operator ()(const T & i, const T & j) const {
			return (i->pt() > j->pt());
		}

	};
}

//
// constants, enums and typedefs
//
namespace reco {
  typedef std::vector<Track> TrackCollection; 
  typedef edm::Ref<TrackCollection> TrackRef; 
}
//
// static data member definitions
//
