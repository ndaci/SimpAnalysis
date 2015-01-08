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
#include <algorithm> 
#include <vector>

// CMSSW core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW specific
#include "DataFormats/JetReco/interface/PFJetCollection.h"
// #include "DataFormats/JetReco/interface/PFJet.h"

#include "CommonTools/Utils/interface/PtComparator.h"

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
  edm::Handle<reco::PFJetCollection> H_pfjets;
  iEvent.getByLabel(jetCollection_ , H_pfjets);
  
  double jetpt1 = 0;
  double jetpt2 = 0;
  double jeteta1 = 0;
  double jeteta2 = 0;
  
  double ptcut = 250;
  double etacut = 2.5;
  
  int index = 0;
  
  for (reco::PFJetCollection::const_iterator jet = H_pfjets->begin(); jet != H_pfjets->end(); ++jet){
    if (index == 0){
      jetpt1 = jet->pt();
      jeteta1 = fabs(jet->eta());
    }
    if (index == 1){
      jetpt2 = jet->pt();
      jeteta2 = fabs(jet->eta());
    }
    index++;
  }
  
  if (jetpt1 > ptcut && jetpt2 > ptcut && jeteta1 < etacut && jeteta2 < etacut) return true;
  else return false;

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
