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
  
  // Declare tree's branches
  //_tree->Branch("",&,"");
  //
  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Jets
  _tree->Branch("nJet",&_nJet,"nJet/I");
  //
  _tree->Branch("jet_eta",&_jet_eta,"jet_eta/D");
  _tree->Branch("jet_phi",&_jet_phi,"jet_phi/D");
  _tree->Branch("jet_pt",&_jet_pt,"jet_pt/D");
  _tree->Branch("jet_e",&_jet_e,"jet_e/D");
  _tree->Branch("jet_m",&_jet_m,"jet_m/D");
  //
  _tree->Branch("jet_mult_ch",&_jet_mult_ch,"jet_mult_ch[nJet]/I");
  _tree->Branch("jet_mult_mu",&_jet_mult_mu,"jet_mult_mu[nJet]/I");
  _tree->Branch("jet_mult_ne",&_jet_mult_ne,"jet_mult_ne[nJet]/I");
  //
  _tree->Branch("jet_efrac_ne_Had", &_jet_efrac_ne_Had, "jet_efrac_ne_Had[nJet]/D");
  _tree->Branch("jet_efrac_ne_EM",  &_jet_efrac_ne_EM,  "jet_efrac_ne_EM[nJet]/D" );
  _tree->Branch("jet_efrac_ch_Had", &_jet_efrac_ch_Had, "jet_efrac_ch_Had[nJet]/D");
  _tree->Branch("jet_efrac_ch_EM",  &_jet_efrac_ch_EM,  "jet_efrac_ch_EM[nJet]/D" );
  _tree->Branch("jet_efrac_ch_Mu",  &_jet_efrac_ch_Mu,  "jet_efrac_ch_Mu[nJet]/D" );
  //

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

  // Initialize branches
  Init();

  // HANDLES //
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
  // Loop over PFJets where theJet is a pointer to a PFJet
  // loop only over 3 highest-pt jets
  //
  UInt_t iJ=0;
  //
  for (reco::PFJetCollection::const_iterator theJet = H_pfjets->begin(); theJet != H_pfjets->end(); ++theJet){

    // Kinematics
    _jet_pt[iJ]  = theJet->pt();
    _jet_eta[iJ] = theJet->eta();
    _jet_phi[iJ] = theJet->phi();
    _jet_e[iJ]   = theJet->energy();
    _jet_m[iJ]   = theJet->mass();

    // Energy fractions
    _jet_efrac_ne_Had[iJ] = theJet->neutralHadronEnergyFraction();
    _jet_efrac_ne_EM[ iJ] = theJet->neutralEmEnergyFraction();
    _jet_efrac_ch_Had[iJ] = theJet->chargedHadronEnergyFraction();
    _jet_efrac_ch_EM[ iJ] = theJet->chargedEmEnergyFraction();
    _jet_efrac_ch_Mu[ iJ] = theJet->chargedMuEnergyFraction();

    // Multiplicities
    _jet_mult_ch[iJ] = theJet->chargedMultiplicity();
    _jet_mult_mu[iJ] = theJet->muonMultiplicity();
    _jet_mult_ne[iJ] = theJet->neutralMultiplicity();

    iJ++ ;
    if(iJ>=nJ) break;
  }

  _nJet = nJ;

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

void
TreeProducer::Init()
{

  _nEvent = _nRun = _nLumi = 0;

  for(UInt_t i=0 ; i<nJ ; i++) {
    _jet_eta[i] = 0;
    _jet_phi[i] = 0;
    _jet_pt[i]  = 0;
    _jet_e[i]   = 0;
    _jet_m[i]   = 0;
    _jet_mult_ch[i] = 0;
    _jet_mult_mu[i] = 0;
    _jet_mult_ne[i] = 0;
    _jet_efrac_ne_Had[i] = 0;
    _jet_efrac_ne_EM[i] = 0;
    _jet_efrac_ch_Had[i] = 0; 
    _jet_efrac_ch_EM[i] = 0; 
    _jet_efrac_ch_Mu[i] = 0;
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
