#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "SANtuple.h"

void SIMP_cutflow(){
	
	double mass[7] = {1, 10, 100, 200, 400, 700, 1000};
  
  TChain * chain0 = new TChain("treeCorr/SimpAnalysis");
  chain0->Add("SIMPs_PUMoriond17_AOD_M1.root");
  TChain * chain1 = new TChain("treeCorr/SimpAnalysis");
  chain1->Add("SIMPs_PUMoriond17_AOD_M10_1.root");
  chain1->Add("SIMPs_PUMoriond17_AOD_M10_2.root");
  chain1->Add("SIMPs_PUMoriond17_AOD_M10_3.root");
  TChain * chain2 = new TChain("treeCorr/SimpAnalysis");
  chain2->Add("SIMPs_PUMoriond17_AOD_M100.root");
  TChain * chain3 = new TChain("treeCorr/SimpAnalysis");
  chain3->Add("SIMPs_PUMoriond17_AOD_M200.root");
  TChain * chain4 = new TChain("treeCorr/SimpAnalysis");
  chain4->Add("SIMPs_PUMoriond17_AOD_M400.root");
  TChain * chain5 = new TChain("treeCorr/SimpAnalysis");
  chain5->Add("SIMPs_PUMoriond17_AOD_M700.root");
  TChain * chain6 = new TChain("treeCorr/SimpAnalysis");
  chain6->Add("SIMPs_PUMoriond17_AOD_M1000.root");
  
  TChain * SPVchain0 = new TChain("treeSPV/SimpAnalysis");
  SPVchain0->Add("SIMPs_PUMoriond17_AOD_M1.root");
  TChain * SPVchain1 = new TChain("treeSPV/SimpAnalysis");
  SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_1.root");
  SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_2.root");
  SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_3.root");
  TChain * SPVchain2 = new TChain("treeSPV/SimpAnalysis");
  SPVchain2->Add("SIMPs_PUMoriond17_AOD_M100.root");
  TChain * SPVchain3 = new TChain("treeSPV/SimpAnalysis");
  SPVchain3->Add("SIMPs_PUMoriond17_AOD_M200.root");
  TChain * SPVchain4 = new TChain("treeSPV/SimpAnalysis");
  SPVchain4->Add("SIMPs_PUMoriond17_AOD_M400.root");
  TChain * SPVchain5 = new TChain("treeSPV/SimpAnalysis");
  SPVchain5->Add("SIMPs_PUMoriond17_AOD_M700.root");
  TChain * SPVchain6 = new TChain("treeSPV/SimpAnalysis");
  SPVchain6->Add("SIMPs_PUMoriond17_AOD_M1000.root");
	
	TChain* chains[7] = {chain0, chain1, chain2, chain3, chain4, chain5, chain6};
	TChain* SPVchains[7] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5, SPVchain6};
//
	TFile* output = new TFile("SIMP_cutflow.root", "RECREATE");
	
  double CHEF_SPVjet[8], CHEF_corrjet[8];

	double passed_deltaphi[7] = {0,0,0,0,0,0,0};
	double passed_eta[7] = {0,0,0,0,0,0,0};
	double passed_photonveto[7] = {0,0,0,0,0,0,0};
	double passed_pt[7] = {0,0,0,0,0,0,0};
	double passed[7][12] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
	double passed_SPV[7][12] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
	double chf_cuts[12] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	
	double lumi = 33.095*1000;
	double QCD_xsec[7] = {4.461, 4.402, 2.553, 0.7903, 0.07434, 0.004846, 0.0005714}; 

	TH1D* pt_eff[7] = {0,0,0,0,0,0,0};
	TH1D* eta_eff[7] = {0,0,0,0,0,0,0};
	TH1D* deltaphi_eff[7] = {0,0,0,0,0,0,0};
	TH1D* SPV_eff[7] = {0,0,0,0,0,0,0};
	TH1D* photonVeto_eff[7] = {0,0,0,0,0,0,0};
	TH1D* ChF_eff[7] = {0,0,0,0,0,0,0};
			
	for (int l = 0; l < 7; l++){
		std::ostringstream strs;
		double dbl = mass[l];
		strs << dbl;
		std::string m = strs.str();
		
		pt_eff[l] = new TH1D(("M"+m+"_pt_eff").c_str(), "passed pt cuts", 100, 200, 300);
		eta_eff[l] = new TH1D(("M"+m+"_eta_eff").c_str(), "passed eta cuts", 100, 0, 5);
		deltaphi_eff[l] = new TH1D(("M"+m+"_deltaphi_eff").c_str(), "passed #Delta#phi cut", 100, 0, 5);
		SPV_eff[l] = new TH1D(("M"+m+"_SPV_eff").c_str(), "passed SPV ChF cuts", 100, 0, 0.5);
		photonVeto_eff[l] = new TH1D(("M"+m+"_photonVeto_eff").c_str(), "passed photon veto", 10, 0, 1);
		ChF_eff[l] = new TH1D(("M"+m+"_ChF_eff").c_str(), "passed ChF cuts", 100, 0, 0.5);
	}
  
	for (int l = 0; l < 7; l++){
		TChain* chain = chains[l];		
		TChain* SPVchain = SPVchains[l];
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t Nentries = corrjets.fChain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
		for(Int_t entry = 0; entry < Nentries; ++entry){
      SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
						
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
			
			if (fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0){
				passed_eta[l] += weight;
				eta_eff[l]->Fill(2.0, weight);
				if(corrjets.jet_pt[0] > 250 && corrjets.jet_pt[1] > 250){
					passed_pt[l] += weight;
					pt_eff[l]->Fill(250, weight);
					if(deltajet_phi > 2){
						passed_deltaphi[l] += weight;
						deltaphi_eff[l]->Fill(2, weight);
            if(corrjets.photon_passLooseId[0] == 0 || (corrjets.photon_passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)){
              passed_photonveto[l] += weight;
              photonVeto_eff[l]->Fill(0.1, weight);
              for(int j = 0; j < 12; j++){
                if (CHEF_corrjet[0]<chf_cuts[j] && CHEF_corrjet[1]<chf_cuts[j]){
                  passed[l][j] += weight;
                  ChF_eff[l]->Fill(chf_cuts[j], weight);
                  if (CHEF_SPVjet[0]<chf_cuts[j] && CHEF_SPVjet[1]<chf_cuts[j]){
                    passed_SPV[l][j] += weight;
                    SPV_eff[l]->Fill(chf_cuts[j], weight);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	std::cout<<"SIMPs & $m_{chi}=1$ GeV & $m_{chi}=10$ GeV & $m_{chi}=100$ GeV & $m_{chi}=200$ GeV & $m_{chi}=400$ GeV & $m_{chi}=700$ GeV & $m_{chi} = 1000$ GeV"<<std::endl; 
	std::cout<<"$|\\eta_{j1, j2}|<2.0$";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_eta[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>250$ GeV";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_pt[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_deltaphi[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_photonveto[l];
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 12; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2} <$ "<<cut;	
		for(int l = 0; l < 7; l++) std::cout<<" & "<<passed[l][j];
		std::cout<<" \\\\"<<std::endl;
		std::cout<<"ChF$_{j1, j2} <$ (both collections)"<<cut;	
		for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_SPV[l][j];
		std::cout<<" \\\\"<<std::endl;
	}	
	output->Write();
	output->Close();	
}