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

// #include "lists/list_QCD_300To500.h"
// #include "lists/list_QCD_500To700.h"
// #include "lists/list_QCD_700To1000.h"
// #include "lists/list_QCD_1000To1500.h"
// #include "lists/list_QCD_1500To2000.h"
// #include "lists/list_QCD_2000ToInf.h"

// #include "lists/list_QCD_300To500_PUMoriond17.h"
// #include "lists/list_QCD_500To700_PUMoriond17.h"
// #include "lists/list_QCD_700To1000_PUMoriond17.h"
// #include "lists/list_QCD_1000To1500_PUMoriond17.h"
// #include "lists/list_QCD_1500To2000_PUMoriond17.h"
// #include "lists/list_QCD_2000ToInf_PUMoriond17.h"

#include "lists/list_QCD_300To500_PUMoriond17_photoninfo.h"
#include "lists/list_QCD_500To700_PUMoriond17_photoninfo.h"
#include "lists/list_QCD_700To1000_PUMoriond17_photoninfo.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_photoninfo.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_photoninfo.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_photoninfo.h"

void SIMP_QCD_NF_eff2D(){
  
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
  
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ne_Had[8], jet_efrac_ne_EM[8], NEF_jet[8];
	double track_pt[10], track_ptError[10];
	int nPixHits[10];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	
	double nf_cuts[11] = {0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_QCD_nPixHits_photonVeto_NHEF_PUMoriond17.root", "RECREATE");
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
	
	std::cout<<"NF cuts: ";
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = nf_cuts[j];
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
		
	for (int l = 0; l < 6; l++){
		TChain* chain = chains[l];		
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ne_Had", &jet_efrac_ne_Had);
		chain->SetBranchAddress("jet_efrac_ne_EM", &jet_efrac_ne_EM);
		chain->SetBranchAddress("track_pt", &track_pt);
		chain->SetBranchAddress("track_ptError", &track_ptError);
		chain->SetBranchAddress("track_nPixHits", &nPixHits);
    chain->SetBranchAddress("photon_pt", &photon_pt);
    chain->SetBranchAddress("photon_eta", &photon_eta);
    chain->SetBranchAddress("photon_phi", &photon_phi);
		chain->SetBranchAddress("photon_passLooseId",&passLooseId);
		chain->SetBranchAddress("photon_passMediumId",&passMediumId);
		chain->SetBranchAddress("photon_passTightId",&passTightId);
		
		Int_t Nentries = chain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"QCD bin "<<l<<": Processing "<<Nentries<<" entries with weight "<<weight<<std::endl;
		
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
				NEF_jet[i] = jet_efrac_ne_Had[i];
			} 
			
			output->cd();
			
			if (jet_pt[0] > 250.0 && jet_pt[1] > 250.0 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && nPixHits > 0 && ( passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)) ){
// 				if (track_ptError[0]/track_pt[0] < 0.5){
					if(NEF_jet[0] < 0.5){
						total->Fill(fabs(jet_eta[1]), jet_pt[1], weight);
						for(int j = 0; j < 11; j++){
							if (NEF_jet[1]>nf_cuts[j]) passed[j]->Fill(fabs(jet_eta[1]), jet_pt[1], weight);
						}
					}
					if(NEF_jet[1] < 0.5){
						total->Fill(fabs(jet_eta[0]), jet_pt[0], weight);
						for(int j = 0; j < 11; j++){
							if (NEF_jet[0]>nf_cuts[j]) passed[j]->Fill(fabs(jet_eta[0]), jet_pt[0], weight);
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