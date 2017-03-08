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

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD.h"

void SIMP_QCD_eff2D(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_300To500(chain0);
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_500To700(chain1);
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_700To1000(chain2);
  TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain3);
  TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain4);
  TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain5);
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_300To500(SPVchain0);
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_500To700(SPVchain1);
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_700To1000(SPVchain2);
  TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain3);
  TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain4);
  TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain5);
	TChain* SPVchains[6] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5};
	std::cout<<"TChains ready"<<std::endl;
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_QCD_nopixcut_PUMoriond17_AOD.root", "RECREATE");
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
		
	for (int l = 0; l < 6; l++){
		TChain* chain = chains[l];
		TChain* SPVchain = SPVchains[l];		
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t Nentries = corrjets.fChain->GetEntries(); 
// 		Int_t Nentries = 1000; 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"QCD bin "<<l<<": Processing "<<Nentries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
			SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
			
// 			Double_t random1 = r.Uniform();
			
			double deltajet_phi = corrjets.jet_phi[0] - corrjets.jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			deltajet_phi = fabs(deltajet_phi);
			
			double deltaphi_jet1photon = corrjets.jet_phi[0] - corrjets.photon_phi[0];
			if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
			if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
			double deltaphi_jet2photon = corrjets.jet_phi[1] - corrjets.photon_phi[0];
			if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
			if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
			
			double deltaeta_jet1photon = corrjets.jet_eta[0] - corrjets.photon_eta[0];
			double deltaeta_jet2photon = corrjets.jet_eta[1] - corrjets.photon_eta[0];
			
			double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
			double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
			
			for (int i = 0; i < 8; i++){
				CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
				CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (corrjets.jet_pt[0] > 250.0 && corrjets.jet_pt[1] > 250.0 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2.0 /*&& npixhits[0] > 0*/ && ( corrjets.photon_passLooseId[0] == 0 || (corrjets.photon_passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
        if(CHEF_corrjet[0] > 0.5){
          total->Fill(fabs(corrjets.jet_eta[1]), corrjets.jet_pt[1], weight);
          for(int j = 0; j < 11; j++){
            if (/*CHEF_SPVjet[1]<chf_cuts[j] &&*/ CHEF_corrjet[1]<chf_cuts[j]) passed[j]->Fill(fabs(corrjets.jet_eta[1]), corrjets.jet_pt[1], weight);
          }
        }
        if(CHEF_corrjet[1] > 0.5){
          total->Fill(fabs(corrjets.jet_eta[0]), corrjets.jet_pt[0], weight);
          for(int j = 0; j < 11; j++){
            if (/*CHEF_SPVjet[0]<chf_cuts[j] && */CHEF_corrjet[0]<chf_cuts[j]) passed[j]->Fill(fabs(corrjets.jet_eta[0]), corrjets.jet_pt[0], weight);
          }
        }
			}    
		}
	}
	for(int j = 0; j < 11; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}