#include "SimpAnalysis/TreeProducer_miniAOD/interface/TreeProducer_miniAOD.h"

//
// constructors and destructor
//
TreeProducer_miniAOD::TreeProducer_miniAOD(const edm::ParameterSet& pset):
  _trigResultsTag(pset.getParameter<edm::InputTag>("triggerResults")),
  _prescalesTag(pset.getParameter<edm::InputTag>("prescales")),
  _METfilterTag(pset.getParameter<edm::InputTag>("METfilter")),
  _pfjetCollectionTag(pset.getParameter<edm::InputTag>("pfjetCollection")),
  _vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  _METCollectionTag(pset.getParameter<edm::InputTag>("METCollection")),
  
  _trigResultsToken(consumes<edm::TriggerResults>(_trigResultsTag)),
  _triggerPrescalesToken(consumes<pat::PackedTriggerPrescales>(_prescalesTag)),
  _METfilterToken(consumes<edm::TriggerResults>(_METfilterTag)),
  _pfjetCollectionToken(consumes<vector<pat::Jet> >(_pfjetCollectionTag)),
  _vertexCollectionToken(consumes<vector<reco::Vertex> >(_vertexCollectionTag)),
  _METCollectionToken(consumes<vector<pat::MET> >(_METCollectionTag))
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
  // Vertices
  _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  _tree->Branch("vtx_N_stored",&_vtx_N_stored,"vtx_N_stored/I");
  _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[vtx_N_stored]/D");
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[vtx_N_stored]/D");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[vtx_N_stored]/D");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[vtx_N_stored]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[vtx_N_stored]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[vtx_N_stored]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[vtx_N_stored]/D");
  //
  // Jets
  _tree->Branch("nJet",&_nJet,"nJet/I");
  //
  _tree->Branch("jet_eta",&_jet_eta,"jet_eta[nJet]/D");
  _tree->Branch("jet_phi",&_jet_phi,"jet_phi[nJet]/D");
  _tree->Branch("jet_pt",&_jet_pt,"jet_pt[nJet]/D");
  _tree->Branch("jet_e",&_jet_e,"jet_e[nJet]/D");
  _tree->Branch("jet_m",&_jet_m,"jet_m[nJet]/D");
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
  
  //MET
  _tree->Branch("MET", &_MET, "MET/D");
  _tree->Branch("MET_phi", &_MET_phi, "MET_phi/D");
  
  //Trigger
  _tree->Branch("HLT_DiCentralPFJet170_CFMax0p1", &_dijet_170_0p1);
  _tree->Branch("HLT_DiCentralPFJet220_CFMax0p3", &_dijet_220_0p3);
  _tree->Branch("HLT_DiCentralPFJet330_CFMax0p5", &_dijet_330_0p5);
  _tree->Branch("HLT_DiCentralPFJet430", &_dijet_430);
  _tree->Branch("HLT_DiCentralPFJet170", &_dijet_170);
  _tree->Branch("HLT_SingleCentralPFJet170_CFMax0p1", &_singlejet_170_0p1);
  //prescales
  _tree->Branch("pswgt_dijet_170", &_pswgt_dijet_170, "pswgt_dijet_170/D");
  _tree->Branch("pswgt_singlejet_170_0p1", &_pswgt_singlejet_170_0p1, "pswgt_singlejet_170_0p1/D");
  
  
  //MET filters
  _tree->Branch("Flag_HBHENoiseFilter", &_HBHENoiseFlag);
  _tree->Branch("Flag_HBHENoiseIsoFilter", &_HBHENoiseIsoFlag);
  _tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &_ECALFlag);
  _tree->Branch("Flag_goodVertices", &_vertexFlag);
  _tree->Branch("Flag_eeBadScFilter", &_eeFlag);
  _tree->Branch("Flag_globalTightHalo2016Filter", &_beamhaloFlag);
}


TreeProducer_miniAOD::~TreeProducer_miniAOD()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer_miniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize branches
  Init();

  // HANDLES //
  // Get collections
  edm::Handle<edm::TriggerResults> H_METfilter;
  iEvent.getByToken(_METfilterToken, H_METfilter);
  
  edm::Handle<edm::TriggerResults> H_trig;
  iEvent.getByToken(_trigResultsToken, H_trig);
  
  edm::Handle<pat::PackedTriggerPrescales> H_prescale;
  iEvent.getByToken(_triggerPrescalesToken, H_prescale);

  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(_vertexCollectionToken, H_vert);

//   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
//   iEvent.getByType("offlineBeamSpot", recoBeamSpotHandle);
//   const reco::BeamSpot bs = *recoBeamSpotHandle;

  edm::Handle<vector<pat::Jet> > H_pfjets;
  iEvent.getByToken(_pfjetCollectionToken , H_pfjets);
  
  edm::Handle<vector<pat::MET> > H_MET;
  iEvent.getByToken(_METCollectionToken , H_MET);

  // Check validity
  if(!H_trig.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trigResultsTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_pfjets.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _pfjetCollectionTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!H_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _vertexCollectionTag << " ... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();

  // VERTICES //
  UInt_t vtx_counter=0;
  _vtx_N = H_vert->size();
  _vtx_N_stored = nV;
	
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(H_vert.product()) );

//   if(_vtx_N > 0) {
//     GlobalPoint local_vertexPosition(sortedVertices.front().position().x(),
// 				     sortedVertices.front().position().y(),
// 				     sortedVertices.front().position().z());
//     vertexPosition = local_vertexPosition;
//   }
//   else {
//     GlobalPoint local_vertexPosition(bs.position().x(),
// 				     bs.position().y(),
// 				     bs.position().z());
//     vertexPosition = local_vertexPosition;
//   }

  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
		
    _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
    _vtx_ndof[vtx_counter] = PV->ndof();
    _vtx_nTracks[vtx_counter] = PV->tracksSize();
    _vtx_d0[vtx_counter] = PV->position().Rho();
    _vtx_x[vtx_counter] = PV->x();
    _vtx_y[vtx_counter] = PV->y();
    _vtx_z[vtx_counter] = PV->z();
		
    vtx_counter++;
    if(vtx_counter >= nV) break;
  } // for loop on primary vertices


  // STORE JET INFORMATION //
  // Loop over PFJets where theJet is a pointer to a PFJet
  // loop only over 3 highest-pt jets
  //
  UInt_t iJ=0;
  //
  for (vector<pat::Jet>::const_iterator theJet = H_pfjets->begin(); theJet != H_pfjets->end(); ++theJet){

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
  
  //MET//
  const pat::MET &met = H_MET->front();
  _MET = met.pt();
  _MET_phi = met.phi();
  
  //Trigger//
  if (triggerPathsMap[triggerPathsVector[0]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[0]])) _dijet_170_0p1 = 1;
  if (triggerPathsMap[triggerPathsVector[1]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[1]])) _dijet_220_0p3 = 1;
  if (triggerPathsMap[triggerPathsVector[2]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[2]])) _dijet_330_0p5 = 1;
  if (triggerPathsMap[triggerPathsVector[3]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[3]])) _dijet_430 = 1;
  if (triggerPathsMap[triggerPathsVector[4]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[4]])) _dijet_170 = 1;
  if (triggerPathsMap[triggerPathsVector[5]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[5]])) _singlejet_170_0p1 = 1;
  //prescales
  const edm::TriggerNames &trignames = iEvent.triggerNames(*H_trig);
    for (size_t i = 0; i < H_trig->size(); i++) {
        if (trignames.triggerName(i).find("HLT_DiCentralPFJet170_v") != string::npos) _pswgt_dijet_170 = H_prescale->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_SingleCentralPFJet170_CFMax0p1_v") != string::npos) _pswgt_singlejet_170_0p1 = H_prescale->getPrescaleForIndex(i);
}

  //MET filters
  if (filterPathsMap[filterPathsVector[0]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[0]])) _HBHENoiseFlag = 1;
  if (filterPathsMap[filterPathsVector[1]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[1]])) _HBHENoiseIsoFlag = 1;
  if (filterPathsMap[filterPathsVector[2]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[2]])) _ECALFlag = 1;
  if (filterPathsMap[filterPathsVector[3]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[3]])) _vertexFlag = 1;
  if (filterPathsMap[filterPathsVector[4]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[4]])) _eeFlag = 1;
  if (filterPathsMap[filterPathsVector[5]] != -1 && H_METfilter->accept(filterPathsMap[filterPathsVector[5]])) _beamhaloFlag = 1;

  _tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer_miniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer_miniAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TreeProducer_miniAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  triggerPathsVector.push_back("HLT_DiCentralPFJet170_CFMax0p1_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet220_CFMax0p3_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet330_CFMax0p5_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet430_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet170_v");
  triggerPathsVector.push_back("HLT_SingleCentralPFJet170_CFMax0p1_v");
  
  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, _trigResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }
  
  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
        triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  } 
  
  filterPathsVector.push_back("Flag_HBHENoiseFilter");
  filterPathsVector.push_back("Flag_eeBadScFilter");
  filterPathsVector.push_back("Flag_HBHENoiseIsoFilter");
  filterPathsVector.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filterPathsVector.push_back("Flag_goodVertices");
  filterPathsVector.push_back("Flag_globalTightHalo2016Filter");
  
  HLTConfigProvider fltrConfig;
  fltrConfig.init(iRun, iSetup, _METfilterTag.process(), changedConfig);
  
  for (size_t i = 0; i < filterPathsVector.size(); i++) {
    filterPathsMap[filterPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < filterPathsVector.size(); i++){
    TPRegexp pattern(filterPathsVector[i]);
    for(size_t j = 0; j < fltrConfig.triggerNames().size(); j++){
      std::string pathName = fltrConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	filterPathsMap[filterPathsVector[i]] = j;
      }
    }
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
TreeProducer_miniAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TreeProducer_miniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TreeProducer_miniAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer_miniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TreeProducer_miniAOD::Init()
{

  _nEvent = _nRun = _nLumi = 0;

  //Trigger
  _dijet_170_0p1 = 0;      
  _dijet_220_0p3 = 0;      
  _dijet_330_0p5 = 0;       
  _dijet_430 = 0;          
  _dijet_170 = 0;           
  _singlejet_170_0p1 = 0;
  
  //prescales
  _pswgt_dijet_170 = 1;
  _pswgt_singlejet_170_0p1 = 1;

  //MET filters
  _HBHENoiseFlag = 0;
  _HBHENoiseIsoFlag = 0;
  _ECALFlag = 0; 
  _vertexFlag = 0;
  _eeFlag = 0;
  _beamhaloFlag = 0;
  
  //MET
  _MET = 0;
  _MET_phi = 0;
  
  // Vertices
  _vtx_N = 0; 
  _vtx_N_stored = 0;
  for(UInt_t iv=0;iv<nV;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }

  //Jets
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
DEFINE_FWK_MODULE(TreeProducer_miniAOD);
