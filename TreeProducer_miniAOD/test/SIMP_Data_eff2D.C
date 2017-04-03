#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TChain.h>
#include <TH2.h>

// #include "lists/list_JetHT_2016B.h"
// #include "lists/list_JetHT_2016C.h"
// #include "lists/list_JetHT_2016D.h"
// #include "lists/list_JetHT_2016E.h"
// #include "lists/list_JetHT_2016F.h"
// #include "lists/list_JetHT_2016G.h"
/*
#include "lists/list_JetHT_rereco_2016B.h"
#include "lists/list_JetHT_rereco_2016C.h"
#include "lists/list_JetHT_rereco_2016D.h"
#include "lists/list_JetHT_rereco_2016E.h"
#include "lists/list_JetHT_rereco_2016F.h"
#include "lists/list_JetHT_rereco_2016G.h"
#include "lists/list_JetHT_rereco_2016Hv2.h"
#include "lists/list_JetHT_rereco_2016Hv3.h"*/

#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016B.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016C.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016D.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016E.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016F.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016G.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016H2.h"
#include "../../TreeProducer_AOD/test/lists/list_JetHT_2016H3.h"

void SIMP_Data_eff2D(){
	
  TChain * chain = new TChain("treeCorr/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
	list_JetHT_2016H2(chain);
	list_JetHT_2016H3(chain);
  
// 	list_JetHT_rereco_2016B(chain);
// 	list_JetHT_rereco_2016C(chain);
// 	list_JetHT_rereco_2016D(chain);
// 	list_JetHT_rereco_2016E(chain);
// 	list_JetHT_rereco_2016F(chain);
// 	list_JetHT_rereco_2016G(chain);
// 	list_JetHT_rereco_2016Hv2(chain);
// 	list_JetHT_rereco_2016Hv3(chain);
	std::cout<<"TChain ready"<<std::endl;
  
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8];
	double track_pt[10], track_ptError[10];
	int nPixHits[10];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	int dijet_170;
	double pswgt_dijet_170;
	
	double chf_cuts[12] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_Data_Rereco_newAOD_lowstat.root", "RECREATE");
	TH2D* total = new TH2D("total", "total", 4, eta_bins, 9, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	TH2D* eff[12]    = {0,0,0,0,0,0,0,0,0,0,0,0};
	
	std::cout<<"CHF cuts: ";
	for(int j = 0; j < 12; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<" "<<cut<<" ";
		std::string title_passed = "passed_"+cut;
		std::string title_eff = "eff_"+cut;
		passed[j] = new TH2D(title_passed.c_str(), title_passed.c_str(), 4, eta_bins, 9, pt_bins);
		passed[j]->GetYaxis()->SetTitle("p_{T}");
		passed[j]->GetXaxis()->SetTitle("#eta");
		passed[j]->Sumw2();
		eff[j] = new TH2D(title_eff.c_str(), title_eff.c_str(), 4, eta_bins, 9, pt_bins);
		eff[j]->GetYaxis()->SetTitle("p_{T}");
		eff[j]->GetXaxis()->SetTitle("#eta");
		eff[j]->Sumw2();
	}
	std::cout<<std::endl;
		
	chain->SetBranchAddress("jet_pt", &jet_pt);
	chain->SetBranchAddress("jet_eta", &jet_eta);
	chain->SetBranchAddress("jet_phi", &jet_phi);
	chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
	chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
	chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
	chain->SetBranchAddress("track_nPixHits", &nPixHits);
	chain->SetBranchAddress("photon_pt", &photon_pt);
	chain->SetBranchAddress("photon_eta", &photon_eta);
	chain->SetBranchAddress("photon_phi", &photon_phi);
	chain->SetBranchAddress("photon_passLooseId",&passLooseId);
	chain->SetBranchAddress("photon_passMediumId",&passMediumId);
	chain->SetBranchAddress("photon_passTightId",&passTightId);
// 		chain->SetBranchAddress("track_pt", &track_pt);
// 		chain->SetBranchAddress("track_ptError", &track_ptError);
	
// 	Int_t Nentries = chain->GetEntries(); 
	Int_t Nentries = 5000000; 
	
	for(Int_t entry = 0; entry < Nentries; ++entry){
		if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
		chain->GetEntry(entry);
		
		double deltajet_phi = jet_phi[0] - jet_phi[1];
		if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
		if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
		deltajet_phi = fabs(deltajet_phi);
		
		double deltaphi_jet1photon = jet_phi[0] - photon_phi[0];
		if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
		if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
		double deltaphi_jet2photon = jet_phi[1] - photon_phi[0];
		if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
		if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
		
		double deltaeta_jet1photon = jet_eta[0] - photon_eta[0];
		double deltaeta_jet2photon = jet_eta[1] - photon_eta[0];
		
		double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
		double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
		
		for (int i = 0; i < 8; i++){
			CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
		} 
		
		output->cd();
		
		if (dijet_170 == 1 && jet_pt[0] > 250.0 && jet_pt[1] > 250.0 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && /*nPixHits > 0 &&*/ (passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
			if(CHEF_jet[0] > 0.5){
				total->Fill(fabs(jet_eta[1]), jet_pt[1], pswgt_dijet_170);
				for(int j = 0; j < 12; j++){
					if (CHEF_jet[1]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[1]), jet_pt[1], pswgt_dijet_170);
				}
			}
			if(CHEF_jet[1] > 0.5){
				total->Fill(fabs(jet_eta[0]), jet_pt[0], pswgt_dijet_170);
				for(int j = 0; j < 12; j++){
					if (CHEF_jet[0]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[0]), jet_pt[0], pswgt_dijet_170);
				}
			}
		}    
	}
	for(int j = 0; j < 12; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}