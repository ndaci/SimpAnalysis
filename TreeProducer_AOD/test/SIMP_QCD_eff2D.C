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

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"

void SIMP_QCD_eff2D(){
  
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
//   TChain* chain5 = new TChain("tree/SimpAnalysis");
// 	list_QCD_2000ToInf(chain5);
	TChain* chains[5] = {chain0, chain1, chain2, chain3, chain4};
	std::cout<<"TChains ready"<<std::endl;
  
// 	TRandom3 r;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8];
	double track_pt[10], track_ptError[10], track_dzError[10], track_dz[10];
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_QCD_dzCut1_PUMoriond17_AOD.root", "RECREATE");
	TH2D* total = new TH2D("total", "total", 4, eta_bins, 9, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[11] = {0,0,0,0,0,0,0,0,0,0,0};
	TH2D* eff[11]    = {0,0,0,0,0,0,0,0,0,0,0};
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
	
	std::cout<<"CHF cuts: ";
	for(int j = 0; j < 11; j++){
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
	}
	std::cout<<std::endl;
		
	for (int l = 0; l < 5; l++){
		TChain* chain = chains[l];		
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		chain->SetBranchAddress("track_pt", &track_pt);
		chain->SetBranchAddress("track_ptError", &track_ptError);
		chain->SetBranchAddress("track_dzError", &track_dzError);
		chain->SetBranchAddress("track_dz", &track_dz);
		
		Int_t Nentries = chain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"QCD bin "<<l<<": Processing "<<Nentries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
			chain->GetEntry(entry);
			
// 			Double_t random1 = r.Uniform();
			
			double deltajet_phi = jet_phi[0] - jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			
			for (int i = 0; i < 8; i++){
				CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (jet_pt[0] > 250.0 && jet_pt[1] > 250.0 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && track_dzError[0]/track_dz[0] < 1){
// 				if (track_ptError[0]/track_pt[0] < 0.5){
					if(CHEF_jet[0] > 0.5){
						total->Fill(fabs(jet_eta[1]), jet_pt[1], weight);
						for(int j = 0; j < 11; j++){
							if (CHEF_jet[1]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[1]), jet_pt[1], weight);
						}
					}
					if(CHEF_jet[1] > 0.5){
						total->Fill(fabs(jet_eta[0]), jet_pt[0], weight);
						for(int j = 0; j < 11; j++){
							if (CHEF_jet[0]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[0]), jet_pt[0], weight);
						}
					}
// 				} else std::cout<<track_pt[0]<<" "<<track_ptError[0]<<" "<<CHEF_jet[0]<<" "<<CHEF_jet[1]<<std::endl;
			}    
		}
	}
	for(int j = 0; j < 11; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}