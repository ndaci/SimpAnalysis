// -*- C++ -*-
//
// Package:    TreeProducer
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc SimpAnalysis/TreeProducer/src/TreeProducer.cc

 Description: EDAnalyzer produce flat trees from EDM for SimpAnalysis

*/
// Original Author:  ndaci
//         Created:  Fri Jan  9 09:54:43 CET 2015
// $V0$
//

// C++ lib
#include <memory>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

// CMSSW standard lib
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// others
using namespace std;
int verbose=1;
const UInt_t nJ=3;

//
// class declaration
//

class TreeProducer : public edm::EDAnalyzer {
 public:
  explicit TreeProducer(const edm::ParameterSet&);
  ~TreeProducer();

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
  string   _hltProcessName;
  edm::InputTag _trigResultsLabel;
  edm::InputTag _pfjetCollection;

  TTree* _tree;

  // For branches
  int _nEvent, _nRun, _nLumi, _nJet;

  // Jets
  int _jet_mult_ch[nJ], _jet_mult_mu[nJ], _jet_mult_ne[nJ]; // multiplicities
  double _jet_eta[nJ], _jet_phi[nJ], _jet_pt[nJ], _jet_e[nJ], _jet_m[nJ];
  double _jet_efrac_ne_Had[nJ], _jet_efrac_ne_EM[nJ]; // neutral energy fractions
  double _jet_efrac_ch_Had[nJ], _jet_efrac_ch_EM[nJ], _jet_efrac_ch_Mu[nJ]; // charged energy fractions

  // Vertices

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


