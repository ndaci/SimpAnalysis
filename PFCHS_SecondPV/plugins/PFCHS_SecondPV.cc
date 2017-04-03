//-------------------------------------------------
//
//   Class: PFCHS_SecondPV
//
//   PFCHS_SecondPV
//
//
//   Author :
//   G. Flouris               VUB        Feb. 2016
//--------------------------------------------------

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ConsumesCollector.h>
#include <FWCore/Framework/interface/one/EDProducer.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
//vector<reco::Vertex>                  "offlinePrimaryVertices"    ""                "RECO"    


using namespace std;

class PFCHS_SecondPV: public edm::one::EDProducer<edm::one::SharedResources> {
public:
  PFCHS_SecondPV(const edm::ParameterSet & pset);
  ~PFCHS_SecondPV() {}
  void produce(edm::Event & e, const edm::EventSetup& c);
private:
  edm::EDGetToken m_vertices;

};


PFCHS_SecondPV::PFCHS_SecondPV(const edm::ParameterSet & pset) {
m_vertices      = consumes<vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("PV_Source"));

produces<vector<reco::Vertex> >("SecondOfflinePrimaryVerices");
}


void PFCHS_SecondPV::produce(edm::Event& e, const edm::EventSetup& c) {


  edm::Handle<vector<reco::Vertex> > H_vertices;
  e.getByToken(m_vertices, H_vertices);

  std::auto_ptr<vector<reco::Vertex> > out_vertices(new vector<reco::Vertex>);



    for (vector<reco::Vertex>::const_iterator thevtx = H_vertices->begin();
    thevtx != H_vertices->end(); ++thevtx){ 
      out_vertices->push_back(*thevtx);
    }

  if (out_vertices->size() > 1) swap((*out_vertices)[0], (*out_vertices)[1]);
  e.put(out_vertices,"SecondOfflinePrimaryVerices");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCHS_SecondPV);