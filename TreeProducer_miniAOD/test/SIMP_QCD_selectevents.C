#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

// #include "list_QCD_300To500.h"
// #include "list_QCD_500To700.h"
// #include "list_QCD_700To1000.h"
// #include "list_QCD_1000To1500.h"
// #include "list_QCD_1500To2000.h"
// #include "list_QCD_2000ToInf.h"

#include "list_QCD_300To500_PUMoriond17.h"
#include "list_QCD_500To700_PUMoriond17.h"
#include "list_QCD_700To1000_PUMoriond17.h"
#include "list_QCD_1000To1500_PUMoriond17.h"
#include "list_QCD_1500To2000_PUMoriond17.h"
#include "list_QCD_2000ToInf_PUMoriond17.h"

void SIMP_QCD_selectevents(){
  
  TChain* chain0 = new TChain("tree/SimpAnalysis");
	list_QCD_300To500(chain0);
  TChain* chain1 = new TChain("tree/SimpAnalysis");
	list_QCD_500To700(chain1);
  TChain* chain2 = new TChain("tree/SimpAnalysis");
	list_QCD_700To1000(chain2);
  TChain* chain3 = new TChain("tree/SimpAnalysis");
	list_QCD_1000To1500(chain3);
  TChain* chain4 = new TChain("tree/SimpAnalysis");
	list_QCD_1500To2000(chain4);
  TChain* chain5 = new TChain("tree/SimpAnalysis");
	list_QCD_2000ToInf(chain5);
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
	
	std::cout<<"TChains ready"<<std::endl;
	
  int run, LS, event;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8];
	double CHEF_jet[8];
	double track_pt[10], track_ptError[10];
	
// 	std::cout<<std::endl;
	
// 	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
// 	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
// 	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};
// // 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
// 	double lumi = 20;
  
	for (int l = 2; l < 3; l++){
		std::cout<<"QCD bin "<<l<<"..."<<std::flush;
		TChain* chain = chains[l];
		chain->SetBranchAddress("nRun", &run);
		chain->SetBranchAddress("nLumi", &LS);
		chain->SetBranchAddress("nEvent", &event);
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		chain->SetBranchAddress("track_pt", &track_pt);
		chain->SetBranchAddress("track_ptError", &track_ptError);
		
		Int_t Nentries = chain->GetEntries(); 
		std::cout<<"processing "<<Nentries<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			chain->GetEntry(entry);
			
			double deltajet_phi = jet_phi[0] - jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			
			for (int i = 0; i < 8; i++){
				CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
			} 
			
			if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2){
// 				if(track_ptError[0]/track_pt[0] > 0.5) std::cout<<run<<":"<<LS<<":"<<event<<std::endl;
				if((CHEF_jet[0]<0.2 && CHEF_jet[1]<0.2)&&(CHEF_jet[1]<0.01 || CHEF_jet[0]<0.01)) std::cout<<run<<":"<<LS<<":"<<event<<std::endl;
// 				if((CHEF_jet[0]>0.5 && CHEF_jet[1]<0.01)||(CHEF_jet[1]>0.5 && CHEF_jet[0]<0.01)) std::cout<<run<<":"<<LS<<":"<<event<<",";
// 				if(CHEF_jet[0]>0.5 && CHEF_jet[1]>0.5) std::cout<<run<<":"<<LS<<":"<<event<<", "<<std::flush;
			}    
		}
		std::cout<<std::endl;
	}
}