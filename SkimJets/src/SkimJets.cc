// -*- C++ -*-
//
// Package:    SkimJets
// Class:      SkimJets
// 
/**\class SkimJets SkimJets.cc SimpAnalysis/SkimJets/src/SkimJets.cc

   Description: Skim input files for SIMP analysis

*/
//
// Original Author:  ndaci
//         Created:  Wed Jan  7 10:25:09 CET 2015
// $V0$

// C++ lib
#include <memory>

// CMSSW core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW specific
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//
// class declaration
//

class SkimJets : public edm::EDFilter {
public:
  explicit SkimJets(const edm::ParameterSet&);
  ~SkimJets();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  edm::InputTag jetCollection_;
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
SkimJets::SkimJets(const edm::ParameterSet& iConfig)
{
  jetCollection_ = iConfig.getParameter<edm::InputTag>("jetCollection");
}


SkimJets::~SkimJets()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SkimJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::PFJetCollection> H_pfjets;
  iEvent.getByLabel(jetCollection_ , H_pfjets);

  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SkimJets::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimJets::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
SkimJets::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SkimJets::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SkimJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SkimJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SkimJets);
