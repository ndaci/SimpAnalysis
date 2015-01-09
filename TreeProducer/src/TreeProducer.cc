#include "SimpAnalysis/TreeProducer/interface/TreeProducer.h"

//
// constructors and destructor
//
TreeProducer::TreeProducer(const edm::ParameterSet& pset):
  _hltProcessName(pset.getParameter<string>("hltProcessName")),
  _trigResultsLabel("TriggerResults", "", _hltProcessName),
  _pfjetCollection(pset.getParameter<edm::InputTag>("pfjetCollection"))
{

  // Initialize when class is created
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("SimpAnalysis","tree");
  
}


TreeProducer::~TreeProducer()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // HANDLES
  // Get collections
  edm::Handle<edm::TriggerResults> H_trig;
  iEvent.getByLabel(_trigResultsLabel, H_trig);

  edm::Handle<reco::PFJetCollection> H_pfjets;
  iEvent.getByLabel(_pfjetCollection , H_pfjets);

  // Check validity
  if(!H_trig.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trigResultsLabel << " ... skip entry !" << endl;
    return;
  }

  if(!H_pfjets.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _pfjetCollection << " ... skip entry !" << endl;
    return;
  }

  // STORE JET INFORMATION //
  //double jpt,jeta,jphi,je;
  //jpt=jeta=jphi=je=0;
  //UInt_t iJ=0;

  // Loop over PFJets where theJet is a pointer to a PFJet
  // loop only over 3 highest-pt jets

  /*
  // Kinematics
  jpt  = theJet->Pt();
  jeta = theJet->Eta();
  jphi = theJet->Phi();
  je   = theJet->E();
  jet_TLV[iJ].SetPtEtaPhiE(jpt,jeta,jphi,je);

  // Energy fractions
  jet_efrac_ne_Had[iJ] = theJet->neutralHadronEnergyFraction();
  jet_efrac_ne_EM[ iJ] = theJet->neutralEmEnergyFraction();
  jet_efrac_ch_Had[iJ] = theJet->chargedHadronEnergyFraction();
  jet_efrac_ch_EM[ iJ] = theJet->chargedEmEnergyFraction();
  jet_efrac_ch_Mu[ iJ] = theJet->chargedMuEnergyFraction();

  // Multiplicities
  jet_mult_ch[iJ] = theJet->chargedMultiplicity();
  jet_mult_mu[iJ] = theJet->muonMultiplicity();
  jet_mult_ne[iJ] = theJet->neutralMultiplicity();

  iJ++ ;
  if(iJ>=3) break;
  */

}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TreeProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TreeProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
