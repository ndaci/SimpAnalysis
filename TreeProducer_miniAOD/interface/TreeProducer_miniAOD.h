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
// #include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// others
using namespace std;
int verbose=1;
const UInt_t nJ=4;
const UInt_t nV=3;

//
// class declaration
//

class TreeProducer_miniAOD : public edm::EDAnalyzer {
 public:
  explicit TreeProducer_miniAOD(const edm::ParameterSet&);
  ~TreeProducer_miniAOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


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
  edm::InputTag _trigResultsTag2;
  edm::InputTag _prescalesTag;
  edm::InputTag _prescalesTag2;
  edm::InputTag _METfilterTag;
  edm::InputTag _pfjetCollectionTag;
  edm::InputTag _vertexCollectionTag;
  edm::InputTag _METCollectionTag;
  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken2;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> _triggerPrescalesToken2;
  edm::EDGetTokenT<edm::TriggerResults> _METfilterToken;
  edm::EDGetTokenT<vector<pat::Jet> > _pfjetCollectionToken;
  edm::EDGetTokenT<vector<reco::Vertex> > _vertexCollectionToken;
  edm::EDGetTokenT<vector<pat::MET> > _METCollectionToken;
//   GlobalPoint vertexPosition;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nJet;

  // Vertices
  int _vtx_N, _vtx_N_stored;
  double _vtx_x[nV], _vtx_y[nV], _vtx_z[nV];
  double _vtx_normalizedChi2[nV], _vtx_ndof[nV], _vtx_nTracks[nV], _vtx_d0[nV];

  // Jets
  int _jet_mult_ch[nJ], _jet_mult_mu[nJ], _jet_mult_ne[nJ]; // multiplicities
  double _jet_vx[nJ], _jet_vy[nJ], _jet_vz[nJ];//vertex position
  double _jet_area[nJ];
  double _jet_eta[nJ], _jet_phi[nJ], _jet_pt[nJ], _jet_e[nJ], _jet_m[nJ];
  double _jet_efrac_ne_Had[nJ], _jet_efrac_ne_EM[nJ]; // neutral energy fractions
  double _jet_efrac_ch_Had[nJ], _jet_efrac_ch_EM[nJ], _jet_efrac_ch_Mu[nJ]; // charged energy fractions

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

//
// constants, enums and typedefs
//

//
// static data member definitions
//


