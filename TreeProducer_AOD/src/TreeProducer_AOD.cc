#include "SimpAnalysis/TreeProducer_AOD/interface/TreeProducer_AOD.h"

//
// constructors and destructor
//
TreeProducer_AOD::TreeProducer_AOD(const edm::ParameterSet& pset):
  _trigResultsTag(pset.getParameter<edm::InputTag>("triggerResults")),
//   _trigResultsTag2(pset.getParameter<edm::InputTag>("triggerResults2")),
//   _prescalesTag(pset.getParameter<edm::InputTag>("prescales")),
//   _prescalesTag2(pset.getParameter<edm::InputTag>("prescales2")),
  _METfilterTag(pset.getParameter<edm::InputTag>("METfilter")),
  _pfjetCollectionTag(pset.getParameter<edm::InputTag>("pfjetCollection")),
  _genjetCollectionTag(pset.getParameter<edm::InputTag>("genjetCollection")),
  _vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  _METCollectionTag(pset.getParameter<edm::InputTag>("METCollection")),
  _photonCollectionTag(pset.getParameter<edm::InputTag>("photonCollection")),
  _trackCollectionTag(pset.getParameter<edm::InputTag>("trackCollection")),
  _phoLooseIdMapTag(pset.getParameter<edm::InputTag>("phoLooseIdMap")),
  _phoMediumIdMapTag(pset.getParameter<edm::InputTag>("phoMediumIdMap")),
  _phoTightIdMapTag(pset.getParameter<edm::InputTag>("phoTightIdMap")),
  
  _trigResultsToken(consumes<edm::TriggerResults>(_trigResultsTag)),
//   _trigResultsToken2(consumes<edm::TriggerResults>(_trigResultsTag2)),
//   _triggerPrescalesToken(consumes<pat::PackedTriggerPrescales>(_prescalesTag)),
//   _triggerPrescalesToken2(consumes<pat::PackedTriggerPrescales>(_prescalesTag2)),
  _METfilterToken(consumes<edm::TriggerResults>(_METfilterTag)),
  _pfjetCollectionToken(consumes<vector<reco::PFJet> >(_pfjetCollectionTag)),
  _genjetCollectionToken(consumes<vector<reco::GenJet> >(_genjetCollectionTag)),
  _vertexCollectionToken(consumes<vector<reco::Vertex> >(_vertexCollectionTag)),
  _METCollectionToken(consumes<vector<reco::PFMET> >(_METCollectionTag)),
  _photonCollectionToken(consumes<edm::View<reco::Photon> >(_photonCollectionTag)),
  _trackCollectionToken(consumes<vector<reco::Track> >(_trackCollectionTag)),
  phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(_phoLooseIdMapTag)),
  phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(_phoMediumIdMapTag)),
	phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(_phoTightIdMapTag)),
  
  _isData(pset.getUntrackedParameter<bool>("isData"))  
{

}


TreeProducer_AOD::~TreeProducer_AOD()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize branches
  Init();

  // HANDLES //
  // Get collections
  edm::Handle<edm::TriggerResults> H_METfilter;
	iEvent.getByToken(_METfilterToken, H_METfilter);
  
  edm::Handle<edm::TriggerResults> H_trig;//, H_trig1, H_trig2;
// 	try {iEvent.getByToken(_trigResultsToken, H_trig1);} catch (...) {;}
// 	try {iEvent.getByToken(_trigResultsToken2, H_trig2);} catch (...) {;}
// 	if(H_trig1.isValid()){
// 		H_trig = H_trig1;
// 		_trigResultsTag = _trigResultsTag1;
// 	}
// 	else if (H_trig2.isValid()){
// 		H_trig = H_trig2;
// 		_trigResultsTag = _trigResultsTag2;
// 	}
  iEvent.getByToken(_trigResultsToken, H_trig);
	
//   edm::Handle<pat::PackedTriggerPrescales> H_prescale, H_prescale1, H_prescale2;
//   try {iEvent.getByToken(_triggerPrescalesToken, H_prescale1);} catch (...) {;}
//   try {iEvent.getByToken(_triggerPrescalesToken2, H_prescale2);} catch (...) {;}
//   if(H_prescale1.isValid()) H_prescale = H_prescale1;
//   else if(H_prescale2.isValid()) H_prescale = H_prescale2;

  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(_vertexCollectionToken, H_vert);

  edm::Handle<vector<reco::PFJet> > H_pfjets;
  iEvent.getByToken(_pfjetCollectionToken , H_pfjets);
	
  edm::Handle<vector<reco::GenJet> > H_genjets;
  iEvent.getByToken(_genjetCollectionToken , H_genjets);
  
  edm::Handle<vector<reco::PFMET> > H_MET;
  iEvent.getByToken(_METCollectionToken , H_MET);
  
  edm::Handle<edm::View<reco::Photon> > H_photon;
  iEvent.getByToken(_photonCollectionToken , H_photon);
	
  edm::Handle<vector<reco::Track> > H_track;
  iEvent.getByToken(_trackCollectionToken , H_track);
	
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(phoTightIdMapToken_ ,tight_id_decisions);
  
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

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
  
  if(!H_photon.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _photonCollectionTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!H_track.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trackCollectionTag << " ... skip entry !" << endl;
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

  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
    _vtx_ndof[vtx_counter] = PV->ndof();
    _vtx_nTracks[vtx_counter] = PV->tracksSize();
//     _vtx_nTracks[vtx_counter] = PV->nTracks();
    _vtx_d0[vtx_counter] = PV->position().Rho();
    _vtx_x[vtx_counter] = PV->x();
    _vtx_y[vtx_counter] = PV->y();
    _vtx_z[vtx_counter] = PV->z();
		
    vtx_counter++;
    if(vtx_counter >= nV) break;
  } // for loop on primary vertices

  // TRACKS //
	vector<reco::TrackRef> trackRef;
	for (vector<reco::Track>::const_iterator theTrack = H_track->begin(); theTrack != H_track->end(); ++theTrack){
		reco::TrackRef ref(H_track, theTrack - H_track->begin());
			trackRef.push_back(ref);
	}
	
	reco::RecoPtSorter<reco::TrackRef> trackSorter;
	std::sort( trackRef.begin(), trackRef.end(), trackSorter);
	
  UInt_t iT=0;  
	for (size_t i = 0; i < trackRef.size(); i++) {
		_track_purity[iT] = trackRef[i]->highPurity;
		_track_Nhits[iT] = trackRef[i]->numberOfValidHits();
		_track_NpixHits[iT] = trackRef[i]->hitPattern().numberOfValidPixelHits();
		if (trackRef[i]->vz() == 0) _track_fromPV[iT] = 1;
		else _track_fromPV[iT] = 0;
		_track_pt[iT] = trackRef[i]->pt();
		_track_eta[iT] = trackRef[i]->eta();
		_track_phi[iT] = trackRef[i]->phi();
		_track_normalizedChi2[iT] = trackRef[i]->normalizedChi2();
		_track_ndof[iT] = trackRef[i]->ndof();
		_track_ptError[iT] = trackRef[i]->ptError();
		_track_dzError[iT] = trackRef[i]->dzError();
		_track_dz[iT] = trackRef[i]->dz();
		_track_dxy[iT] = trackRef[i]->dxy();
		_track_d0[iT] = trackRef[i]->d0();
		iT++ ;
    if(iT>=nT) break;
	}
		
  _nTrack = iT;
  _nTrack_stored = nT;

  // STORE JET INFORMATION //
  // Loop over PFJets where theJet is a pointer to a PFJet
  // loop only over 8 highest-pt jets
  // only count jets with pt > 20 GeV
  //
  UInt_t iJ=0;
  //
  for (vector<reco::PFJet>::const_iterator theJet = H_pfjets->begin(); theJet != H_pfjets->end(); ++theJet){
		if (theJet->pt() > 20){
			//Area
			_jet_area[iJ] = theJet->jetArea();  
				
			// Vertex
			_jet_vx[iJ] = theJet->vx();
			_jet_vy[iJ] = theJet->vy();
			_jet_vz[iJ] = theJet->vz();

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
      
      jecUnc->setJetEta(_jet_eta[iJ]);
      jecUnc->setJetPt(_jet_pt[iJ]); // here you must use the CORRECTED jet pt
      _jet_unc[iJ] = jecUnc->getUncertainty(true);
      _jet_ptCor_up[iJ] = (theJet->pt())*(1+_jet_unc[iJ]) ; // shift = +1(up), or -1(down)
      _jet_ptCor_down[iJ] = (theJet->pt())*(1-_jet_unc[iJ]) ;
			
			iJ++ ;
		}
    if(iJ>=nJ) break;
  }
  _nJet = iJ;
  _nJet_stored = nJ;
	
  UInt_t iGJ=0;
  //
	if(!_isData){
    for (vector<reco::GenJet>::const_iterator theGenJet = H_genjets->begin(); theGenJet != H_genjets->end(); ++theGenJet){
      if (theGenJet->pt() > 20){
        //Area
        _genjet_area[iJ] = theGenJet->jetArea();  
          
        // Vertex
        _genjet_vx[iGJ] = theGenJet->vx();
        _genjet_vy[iGJ] = theGenJet->vy();
        _genjet_vz[iGJ] = theGenJet->vz();

        // Kinematics
        _genjet_pt[iGJ]  = theGenJet->pt();
        _genjet_eta[iGJ] = theGenJet->eta();
        _genjet_phi[iGJ] = theGenJet->phi();
        _genjet_e[iGJ]   = theGenJet->energy();
        _genjet_m[iGJ]   = theGenJet->mass();

        // Energy fractions
        for (size_t i = 0; i < theGenJet->numberOfDaughters(); i++){
          const reco::Candidate * constituent = theGenJet->daughter(i);
          if (constituent->charge() != 0) _genjet_efrac_ch[iGJ] += constituent->energy()/theGenJet->energy();
        }
        
        iGJ++ ;
      }
      if(iGJ>=nGJ) break;
    }
    _nGenJet = iGJ;
    _nGenJet_stored = nGJ;
  }
  
  //MET//
  const reco::PFMET &met = H_MET->front();
  _MET = met.pt();
  _MET_phi = met.phi();
	
	// PHOTONS//
	UInt_t iP=0;
// 	for (vector<reco::Photon>::const_iterator thePhoton = H_photon->begin(); thePhoton != H_photon->end(); ++thePhoton){
  for (size_t i = 0; i < H_photon->size(); ++i){
		const auto thePhoton = H_photon->ptrAt(i);
		_photon_pt[iP] = thePhoton->pt();
		_photon_eta[iP] = thePhoton->eta();
		_photon_phi[iP] = thePhoton->phi();
		bool isPassLoose  = (*loose_id_decisions)[thePhoton];
    bool isPassMedium = (*medium_id_decisions)[thePhoton];
    bool isPassTight  = (*tight_id_decisions)[thePhoton];
    _passLooseId[iP] = ( (int)isPassLoose );
    _passMediumId[iP] = ( (int)isPassMedium);
		_passTightId[iP] = ( (int)isPassTight );
		iP++;
		if (iP>=nP) break;
	}
  _nPhoton_stored = nP;
  
  //TRIGGER//
  if (triggerPathsMap[triggerPathsVector[0]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[0]])) _dijet_170_0p1 = 1;
  if (triggerPathsMap[triggerPathsVector[1]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[1]])) _dijet_220_0p3 = 1;
  if (triggerPathsMap[triggerPathsVector[2]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[2]])) _dijet_330_0p5 = 1;
  if (triggerPathsMap[triggerPathsVector[3]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[3]])) _dijet_430 = 1;
  if (triggerPathsMap[triggerPathsVector[4]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[4]])) _dijet_170 = 1;
  if (triggerPathsMap[triggerPathsVector[5]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[5]])) _singlejet_170_0p1 = 1;
  //prescales
//   const edm::TriggerNames &trignames = iEvent.triggerNames(*H_trig);
// 	for (size_t i = 0; i < H_trig->size(); i++) {
// 			if (trignames.triggerName(i).find("HLT_DiCentralPFJet170_v") != string::npos) _pswgt_dijet_170 = H_prescale->getPrescaleForIndex(i);
// 			if (trignames.triggerName(i).find("HLT_SingleCentralPFJet170_CFMax0p1_v") != string::npos) _pswgt_singlejet_170_0p1 = H_prescale->getPrescaleForIndex(i);
// 	}

  //MET filters//
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
TreeProducer_AOD::beginJob()
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
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[vtx_N_stored]/I");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[vtx_N_stored]/I");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[vtx_N_stored]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[vtx_N_stored]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[vtx_N_stored]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[vtx_N_stored]/D");
  //
	// Tracks
	_tree->Branch("nTrack_stored",&_nTrack_stored,"nTrack_stored/I");
	_tree->Branch("nTrack",&_nTrack,"nTrack/I");
	
	_tree->Branch("track_pt",&_track_pt,"track_pt[nTrack_stored]/D");
	_tree->Branch("track_eta",&_track_eta,"track_eta[nTrack_stored]/D");
	_tree->Branch("track_phi",&_track_phi,"track_phi[nTrack_stored]/D");
	
	_tree->Branch("track_normalizedChi2",&_track_normalizedChi2,"track_normalizedChi2[nTrack_stored]/D");
	_tree->Branch("track_ndof",&_track_ndof,"track_ndof[nTrack_stored]/I");
	_tree->Branch("track_ptError",&_track_ptError,"track_ptError[nTrack_stored]/D");
	_tree->Branch("track_dzError",&_track_dzError,"track_dzError[nTrack_stored]/D");
	_tree->Branch("track_dz",&_track_dz,"track_dz[nTrack_stored]/D");
	_tree->Branch("track_fromPV",&_track_fromPV,"track_fromPV[nTrack_stored]/I");
	_tree->Branch("track_purity",&_track_purity,"track_purity[nTrack_stored]/I");
	_tree->Branch("track_nhits",&_track_Nhits,"track_nhits[nTrack_stored]/I");
	_tree->Branch("track_nPixHits",&_track_NpixHits,"track_nPixHits[nTrack_stored]/I");
	_tree->Branch("track_d0",&_track_d0,"track_d0[nTrack_stored]/D");
	_tree->Branch("track_dxy",&_track_dxy,"track_dxy[nTrack_stored]/D");
	//
  // Jets
  _tree->Branch("nJet_stored",&_nJet_stored,"nJet_stored/I");
  _tree->Branch("nJet",&_nJet,"nJet/I");
  //
  _tree->Branch("jetArea",&_jet_area,"jetArea[nJet_stored]/D");
  //
  _tree->Branch("jet_vx",&_jet_vx,"jet_vx[nJet_stored]/D");
  _tree->Branch("jet_vy",&_jet_vy,"jet_vy[nJet_stored]/D");
  _tree->Branch("jet_vz",&_jet_vz,"jet_vz[nJet_stored]/D");
  //
  _tree->Branch("jet_eta",&_jet_eta,"jet_eta[nJet_stored]/D");
  _tree->Branch("jet_phi",&_jet_phi,"jet_phi[nJet_stored]/D");
  _tree->Branch("jet_pt",&_jet_pt,"jet_pt[nJet_stored]/D");
  _tree->Branch("jet_e",&_jet_e,"jet_e[nJet_stored]/D");
  _tree->Branch("jet_m",&_jet_m,"jet_m[nJet_stored]/D");
  //
  _tree->Branch("jet_mult_ch",&_jet_mult_ch,"jet_mult_ch[nJet_stored]/I");
  _tree->Branch("jet_mult_mu",&_jet_mult_mu,"jet_mult_mu[nJet_stored]/I");
  _tree->Branch("jet_mult_ne",&_jet_mult_ne,"jet_mult_ne[nJet_stored]/I");
  //
  _tree->Branch("jet_efrac_ne_Had", &_jet_efrac_ne_Had, "jet_efrac_ne_Had[nJet_stored]/D");
  _tree->Branch("jet_efrac_ne_EM",  &_jet_efrac_ne_EM,  "jet_efrac_ne_EM[nJet_stored]/D" );
  _tree->Branch("jet_efrac_ch_Had", &_jet_efrac_ch_Had, "jet_efrac_ch_Had[nJet_stored]/D");
  _tree->Branch("jet_efrac_ch_EM",  &_jet_efrac_ch_EM,  "jet_efrac_ch_EM[nJet_stored]/D" );
  _tree->Branch("jet_efrac_ch_Mu",  &_jet_efrac_ch_Mu,  "jet_efrac_ch_Mu[nJet_stored]/D" );
  _tree->Branch("jet_unc",&_jet_unc,"jet_unc[nJet_stored]/D");
  _tree->Branch("_jet_ptCor_up",&_jet_ptCor_up,"_jet_ptCor_up[nJet_stored]/D");
  _tree->Branch("_jet_ptCor_down",&_jet_ptCor_down,"_jet_ptCor_down[nJet_stored]/D");
	//
  // GenJets
  _tree->Branch("nGenJet_stored",&_nGenJet_stored,"nGenJet_stored/I");
  _tree->Branch("nGenJet",&_nGenJet,"nGenJet/I");
  //
  _tree->Branch("genjetArea",&_genjet_area,"genjetArea[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_vx",&_genjet_vx,"genjet_vx[nGenJet_stored]/D");
  _tree->Branch("genjet_vy",&_genjet_vy,"genjet_vy[nGenJet_stored]/D");
  _tree->Branch("genjet_vz",&_genjet_vz,"genjet_vz[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_eta",&_genjet_eta,"genjet_eta[nGenJet_stored]/D");
  _tree->Branch("genjet_phi",&_genjet_phi,"genjet_phi[nGenJet_stored]/D");
  _tree->Branch("genjet_pt",&_genjet_pt,"genjet_pt[nGenJet_stored]/D");
  _tree->Branch("genjet_e",&_genjet_e,"genjet_e[nGenJet_stored]/D");
  _tree->Branch("genjet_m",&_genjet_m,"genjet_m[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_efrac_ch", &_genjet_efrac_ch, "genjet_efrac_ch[nGenJet_stored]/D");
	
	//Photons
  _tree->Branch("nPhoton_stored",&_nPhoton_stored,"nPhoton_stored/I");
  _tree->Branch("photon_eta",&_photon_eta,"photon_eta[nPhoton_stored]/D");
  _tree->Branch("photon_phi",&_photon_phi,"photon_phi[nPhoton_stored]/D");
  _tree->Branch("photon_pt",&_photon_pt,"photon_pt[nPhoton_stored]/D");
  _tree->Branch("photon_passLooseId",&_passLooseId, "photon_passLooseId[nPhoton_stored]/I");
  _tree->Branch("photon_passMediumId",&_passMediumId, "photon_passLooseId[nPhoton_stored]/I");
  _tree->Branch("photon_passTightId",&_passTightId, "photon_passLooseId[nPhoton_stored]/I");
  
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

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer_AOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TreeProducer_AOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
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
  bool changedConfig2 = false;
  fltrConfig.init(iRun, iSetup, _METfilterTag.process(), changedConfig2);
  
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
TreeProducer_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TreeProducer_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TreeProducer_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TreeProducer_AOD::Init()
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
  
  //Tracks
  _nTrack = 0;
	_nTrack_stored = 0;
  for(UInt_t it=0;it<nT;it++) {
		_track_eta[it] = 0;
		_track_fromPV[it] = 0;
		_track_ndof[it] = 0;
		_track_Nhits[it] = 0;
		_track_normalizedChi2[it] = 0;
		_track_NpixHits[it] = 0;
		_track_phi[it] = 0;
		_track_pt[it] = 0;
		_track_ptError[it] = 0;
		_track_dzError[it] = 0;
		_track_dz[it] = 0;
		_track_d0[it] = 0;
		_track_dxy[it] = 0;
		_track_purity[it] = 0;
	}  

  //Jets
  _nJet = 0;
	_nJet_stored = 0;
  for(UInt_t i=0 ; i<nJ ; i++) {
    _jet_vx[i]  = 0;
    _jet_vy[i]  = 0;
    _jet_vz[i]  = 0;
    _jet_area[i]= 0;
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
    _jet_unc[i]  = 0;
    _jet_ptCor_up[i]  = 0;
    _jet_ptCor_down[i]  = 0;
	}

  //GenJets
  _nGenJet = 0;
	_nGenJet_stored = 0;
  for(UInt_t i=0 ; i<nGJ ; i++) {
    _genjet_vx[i]  = 0;
    _genjet_vy[i]  = 0;
    _genjet_vz[i]  = 0;
    _genjet_area[i]= 0;
    _genjet_eta[i] = 0;
    _genjet_phi[i] = 0;
    _genjet_pt[i]  = 0;
    _genjet_e[i]   = 0;
    _genjet_m[i]   = 0;
    _genjet_efrac_ch[i] = 0;
	}
		
	//Photons
	for(UInt_t i=0 ; i<nP ; i++) {
		_photon_eta[i] = 0;
		_photon_phi[i] = 0;
		_photon_pt[i] = 0;
		_passLooseId[i] = 0;
		_passMediumId[i] = 0;
		_passTightId[i] = 0;
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
