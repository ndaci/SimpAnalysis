#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TH2.h>

#include "SANtuple.h"
// #include "json_parser.cpp"

#include "lists/list_JetHT_2016B.h"
#include "lists/list_JetHT_2016C.h"
#include "lists/list_JetHT_2016D.h"
#include "lists/list_JetHT_2016E.h"
#include "lists/list_JetHT_2016F.h"
#include "lists/list_JetHT_2016G.h"
#include "lists/list_JetHT_2016H2.h"
#include "lists/list_JetHT_2016H3.h"

void SIMP_Data_eff2D(){
	  
  TChain* chain = new TChain("treeCorr/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
	list_JetHT_2016H2(chain);
	list_JetHT_2016H3(chain);
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
	list_JetHT_2016B(SPVchain);
	list_JetHT_2016C(SPVchain);
	list_JetHT_2016D(SPVchain);
	list_JetHT_2016E(SPVchain);
	list_JetHT_2016F(SPVchain);
	list_JetHT_2016G(SPVchain);
	list_JetHT_2016H2(SPVchain);
	list_JetHT_2016H3(SPVchain);
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
	
	double chf_cuts[23] = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[17] = {-2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
  
//   pt::ptree root;
//   pt::read_json("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt", root);
  
  TFile *output = new TFile("eff2D_Data_Rereco_AOD_withNvertexCut_Njets2_finerEtaBinning.root", "RECREATE");
	TH2D* total = new TH2D("total", "total", 16, eta_bins, 9, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	TH2D* eff[23]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	std::cout<<"CHF cuts: ";
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<" "<<cut<<" ";
		std::string title_passed = "passed_"+cut;
		std::string title_eff = "eff_"+cut;
		passed[j] = new TH2D(title_passed.c_str(), title_passed.c_str(), 16, eta_bins, 9, pt_bins);
		passed[j]->GetYaxis()->SetTitle("p_{T}");
		passed[j]->GetXaxis()->SetTitle("#eta");
		passed[j]->Sumw2();
		eff[j] = new TH2D(title_eff.c_str(), title_eff.c_str(), 16, eta_bins, 9, pt_bins);
		eff[j]->GetYaxis()->SetTitle("p_{T}");
		eff[j]->GetXaxis()->SetTitle("#eta");
	}
	std::cout<<std::endl;
		
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t entries = corrjets.fChain->GetEntries(); 
//   Int_t entries = 3000000; 
  std::cout<<"Processing "<<entries<<std::endl;
  
  for(Int_t entry = 0; entry < entries; ++entry){
    if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
    SPVjets.GetEntry(entry);
    corrjets.GetEntry(entry);
    
//     bool isInJSON = Check_Run_Lumi(corrjets.nRun, corrjets.nLumi, root);
//     if (!isInJSON) continue;
    
    double njets = 0;
    for (int k = 0; k < 8; ++k){
      if (corrjets.jet_pt[k] > 30) njets++;
    }
    
    double deltajet_phi = corrjets.jet_phi[0] - corrjets.jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
    deltajet_phi = fabs(deltajet_phi);
      
    double photonpt = 0;
    photon_nr = 0;
    for (int i = 0; i < 4; ++i){
      if(corrjets.photon_pt[i]>photonpt){
        photon_nr = i;
        photonpt = corrjets.photon_pt[i];
      }
    }
    
    double deltaphi_jet1photon = corrjets.jet_phi[0] - corrjets.photon_phi[photon_nr];
    if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
    if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
    double deltaphi_jet2photon = corrjets.jet_phi[1] - corrjets.photon_phi[photon_nr];
    if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
    if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
    
    double deltaeta_jet1photon = corrjets.jet_eta[0] - corrjets.photon_eta[photon_nr];
    double deltaeta_jet2photon = corrjets.jet_eta[1] - corrjets.photon_eta[photon_nr];
    
    double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
    double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
    
    for (int i = 0; i < 8; i++){
      CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
    } 
    
    output->cd();
    
    if (corrjets.HLT_DiCentralPFJet170 == 1 && corrjets.jet_pt[0] > 250.0 && corrjets.jet_pt[1] > 250.0 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && ( corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.vtx_N >= 2 && njets == 2){
      if(CHEF_corrjet[0] > 0.5){
//         std::cout<<corrjets.pswgt_dijet_170<<std::endl;
        total->Fill(corrjets.jet_eta[1], corrjets.jet_pt[1], corrjets.pswgt_dijet_170);
        for(int j = 0; j < 23; j++){
          if (CHEF_corrjet[1]<chf_cuts[j]) passed[j]->Fill(corrjets.jet_eta[1], corrjets.jet_pt[1], corrjets.pswgt_dijet_170);
        }
      }
      if(CHEF_corrjet[1] > 0.5){
        total->Fill(corrjets.jet_eta[0], corrjets.jet_pt[0], corrjets.pswgt_dijet_170);
        for(int j = 0; j < 23; j++){
          if (CHEF_corrjet[0]<chf_cuts[j]) passed[j]->Fill(corrjets.jet_eta[0], corrjets.jet_pt[0], corrjets.pswgt_dijet_170);
        }
      }
    }    
	}
	for(int j = 0; j < 23; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}