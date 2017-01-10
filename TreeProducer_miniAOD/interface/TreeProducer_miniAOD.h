// -*- C++ -*-
//
// Package:    TreeProducer_miniAOD
// Class:      TreeProducer_miniAOD
// 
/**\class TreeProducer_miniAOD TreeProducer_miniAOD.cc SimpAnalysis/TreeProducer_miniAOD/src/TreeProducer_miniAOD.cc

 Description: EDAnalyzer produce flat trees from miniAOD for SimpAnalysis

*/

// C++ lib
#include <memory>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

// CMSSW standard lib
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// others
using namespace std;
int verbose=1;
const UInt_t nJ=8; // #jets
const UInt_t nGJ=4; // #genjets
const UInt_t nV=3; // #vertices
const UInt_t nT=10; // #tracks
const UInt_t nP=4; // #photons

//
// class declaration
//

class TreeProducer_miniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {
 public:
  explicit TreeProducer_miniAOD(const edm::ParameterSet&);
  ~TreeProducer_miniAOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool ptSorter(const pat::PackedCandidate & i, const pat::PackedCandidate & j);

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
  edm::InputTag _prescalesTag;
  edm::InputTag _prescalesTag2;
  edm::InputTag _METfilterTag;
  edm::InputTag _pfjetCollectionTag;
  edm::InputTag _genjetCollectionTag;
  edm::InputTag _vertexCollectionTag;
  edm::InputTag _METCollectionTag;
  edm::InputTag _photonCollectionTag;
  edm::InputTag _packedPFCollectionTag;
  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
//   edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken2;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken2;
  edm::EDGetTokenT<edm::TriggerResults> _METfilterToken;
  edm::EDGetTokenT<vector<pat::Jet> > _pfjetCollectionToken;
  edm::EDGetTokenT<vector<reco::GenJet> > _genjetCollectionToken;
  edm::EDGetTokenT<vector<reco::Vertex> > _vertexCollectionToken;
  edm::EDGetTokenT<vector<pat::MET> > _METCollectionToken;
  edm::EDGetTokenT<vector<pat::Photon> > _photonCollectionToken;
  edm::EDGetTokenT<vector<pat::PackedCandidate> > _packedPFCollectionToken;
//   GlobalPoint vertexPosition;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nJet, _nJet_stored, _nGenJet, _nGenJet_stored, _nTrack, _nTrack_stored;

  // Vertices
  int _vtx_N, _vtx_N_stored, _vtx_ndof[nV], _vtx_nTracks[nV];
  double _vtx_x[nV], _vtx_y[nV], _vtx_z[nV];
  double _vtx_normalizedChi2[nV], _vtx_d0[nV];
	
	//Photons
	
	
	//Tracks
	int _track_fromPV[nT], _track_Nhits[nT], _track_NpixHits[nT], _track_purity[nT], _track_ndof[nT];
	double _track_eta[nT], _track_pt[nT], _track_phi[nT], _track_ptError[nT], _track_dz[nT], _track_dzError[nT], _track_normalizedChi2[nT];	

  // Jets
  int _jet_mult_ch[nJ], _jet_mult_mu[nJ], _jet_mult_ne[nJ]; // multiplicities
  double _jet_vx[nJ], _jet_vy[nJ], _jet_vz[nJ];//vertex position
  double _jet_area[nJ];
  double _jet_eta[nJ], _jet_phi[nJ], _jet_pt[nJ], _jet_e[nJ], _jet_m[nJ];
  double _jet_efrac_ne_Had[nJ], _jet_efrac_ne_EM[nJ]; // neutral energy fractions
  double _jet_efrac_ch_Had[nJ], _jet_efrac_ch_EM[nJ], _jet_efrac_ch_Mu[nJ]; // charged energy fractions
  
  // GenJets
  double _genjet_vx[nGJ], _genjet_vy[nGJ], _genjet_vz[nGJ];//vertex position
  double _genjet_area[nGJ];
  double _genjet_eta[nGJ], _genjet_phi[nGJ], _genjet_pt[nGJ], _genjet_e[nGJ], _genjet_m[nGJ];
  double _genjet_efrac_ch[nGJ];// charged energy fractions
  
  //Photons
  int _nPhoton_stored;
  double _photon_pt[nP], _photon_eta[nP], _photon_phi[nP];

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

};

namespace pat {
	template<typename T> 
	class PatPtSorter{
	public:
		bool operator ()(const T & i, const T & j) const {
			return (i->pt() > j->pt());
		}

	};
}

//
// constants, enums and typedefs
//
namespace pat {
  typedef std::vector<PackedCandidate> PackedCollection; 
  typedef edm::Ref<PackedCollection> PackedCandidateRef; 
}
//
// static data member definitions
//
