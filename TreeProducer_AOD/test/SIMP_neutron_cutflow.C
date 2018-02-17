#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "SANtuple.h"

void SIMP_neutron_cutflow(){
	
	double mass[8] = {1, 10, 100, 200, 400, 700, 1000, 2000};
  
  TChain * chain0 = new TChain("treeCorr/SimpAnalysis");
  chain0->Add("SIMPs_Neutrons_M1_AOD_1.root");
  chain0->Add("SIMPs_Neutrons_M1_AOD_2.root");
  chain0->Add("SIMPs_Neutrons_M1_AOD_3.root");
  TChain * chain1 = new TChain("treeCorr/SimpAnalysis");
  chain1->Add("SIMPs_Neutrons_M10_AOD_1.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_2.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_3.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_4.root");
  TChain * chain2 = new TChain("treeCorr/SimpAnalysis");
  chain2->Add("SIMPs_Neutrons_M100_AOD_1.root");
  chain2->Add("SIMPs_Neutrons_M100_AOD_2.root");
  chain2->Add("SIMPs_Neutrons_M100_AOD_3.root");
  TChain * chain3 = new TChain("treeCorr/SimpAnalysis");
  chain3->Add("SIMPs_Neutrons_M200_AOD.root");
  TChain * chain4 = new TChain("treeCorr/SimpAnalysis");
  chain4->Add("SIMPs_Neutrons_M400_AOD.root");
  TChain * chain5 = new TChain("treeCorr/SimpAnalysis");
  chain5->Add("SIMPs_Neutrons_M700_AOD.root");
  TChain * chain6 = new TChain("treeCorr/SimpAnalysis");
  chain6->Add("SIMPs_Neutrons_M1000_AOD.root");
  TChain * chain7 = new TChain("treeCorr/SimpAnalysis");
  chain7->Add("SIMPs_Neutrons_M2000_AOD.root");
  
  TChain * SPVchain0 = new TChain("treeSPV/SimpAnalysis");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_1.root");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_2.root");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_3.root");
  TChain * SPVchain1 = new TChain("treeSPV/SimpAnalysis");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_1.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_2.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_3.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_4.root");
  TChain * SPVchain2 = new TChain("treeSPV/SimpAnalysis");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_1.root");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_2.root");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_3.root");
  TChain * SPVchain3 = new TChain("treeSPV/SimpAnalysis");
  SPVchain3->Add("SIMPs_Neutrons_M200_AOD.root");
  TChain * SPVchain4 = new TChain("treeSPV/SimpAnalysis");
  SPVchain4->Add("SIMPs_Neutrons_M400_AOD.root");
  TChain * SPVchain5 = new TChain("treeSPV/SimpAnalysis");
  SPVchain5->Add("SIMPs_Neutrons_M700_AOD_new.root");
  TChain * SPVchain6 = new TChain("treeSPV/SimpAnalysis");
  SPVchain6->Add("SIMPs_Neutrons_M1000_AOD.root");
  TChain * SPVchain7 = new TChain("treeSPV/SimpAnalysis");
  SPVchain7->Add("SIMPs_Neutrons_M2000_AOD.root");
	
	TChain* chains[8] = {chain0, chain1, chain2, chain3, chain4, chain5, chain6, chain7};
	TChain* SPVchains[8] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5, SPVchain6, SPVchain7};
//
	TFile* output = new TFile("SIMP_Neutrons_cutflow_M2000.root", "RECREATE");
	
  double CHEF_SPVjet[8], CHEF_corrjet[8], EMF[8];
  int photon_nr;

	double passed_deltaphi[8] = {0,0,0,0,0,0,0,0};
	double total[8] = {0,0,0,0,0,0,0,0};
	double passed_eta[8] = {0,0,0,0,0,0,0,0};
	double passed_photonveto[8] = {0,0,0,0,0,0,0,0};
	double passed_conversionveto[8] = {0,0,0,0,0,0,0,0};
	double passed_njets[8] = {0,0,0,0,0,0,0,0};
	double passed_nvtx[8] = {0,0,0,0,0,0,0,0};
	double passed_filters[8] = {0,0,0,0,0,0,0,0};
	double passed_pt[8] = {0,0,0,0,0,0,0,0};
	double passed[8][12] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
	double passed_SPV[8][12] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
	double chf_cuts[12] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	
// 	double lumi = 33.095*1000;
	double lumi = 16.146*1000;
	double QCD_xsec[8] = {4.461, 4.402, 2.553, 0.7903, 0.07434, 0.004846, 0.0005714, 0.000001819}; 

	TH1D* pt_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* eta_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* deltaphi_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* SPV_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* photonVeto_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* conversionVeto_eff[8] = {0,0,0,0,0,0,0,0};
  TH1D* njets_eff[8] = {0,0,0,0,0,0,0,0};
  TH1D* nvtx_eff[8] = {0,0,0,0,0,0,0,0};
  TH1D* filters_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* ChF_eff[8] = {0,0,0,0,0,0,0,0};
			
	for (int l = 0; l < 8; l++){
		std::ostringstream strs;
		double dbl = mass[l];
		strs << dbl;
		std::string m = strs.str();
		
		pt_eff[l] = new TH1D(("M"+m+"_pt_eff").c_str(), "passed pt cuts", 100, 200, 300);
		eta_eff[l] = new TH1D(("M"+m+"_eta_eff").c_str(), "passed eta cuts", 100, 0, 5);
		deltaphi_eff[l] = new TH1D(("M"+m+"_deltaphi_eff").c_str(), "passed #Delta#phi cut", 100, 0, 5);
		SPV_eff[l] = new TH1D(("M"+m+"_SPV_eff").c_str(), "passed SPV ChF cuts", 100, 0, 0.5);
		photonVeto_eff[l] = new TH1D(("M"+m+"_photonVeto_eff").c_str(), "passed photon veto", 10, 0, 1);
		conversionVeto_eff[l] = new TH1D(("M"+m+"_conversionVeto_eff").c_str(), "passed conversion veto", 10, 0, 1);
    njets_eff[l] = new TH1D(("M"+m+"njets_eff").c_str(), "passed njet cut", 5, 0, 4);
    filters_eff[l] = new TH1D(("M"+m+"filters_eff").c_str(), "passed filters", 5, 0, 4);
    nvtx_eff[l] = new TH1D(("M"+m+"nvtx_eff").c_str(), "passed nvtx cut", 5, 0, 4);  
		ChF_eff[l] = new TH1D(("M"+m+"_ChF_eff").c_str(), "passed ChF cuts", 100, 0, 0.5);
	}
  
	for (int l = 0; l < 8; l++){
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
				CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
				CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
        EMF[i] = corrjets.jet_efrac_ch_EM[i] + corrjets.jet_efrac_ne_EM[i];
      }
      
      bool pass_conv_1 = true;
      bool pass_conv_2 = true;
      if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
      if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
			
			output->cd();
			
//       if (corrjets.vtx_N > 25){
        total[l] += weight;
          
        if(corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550){
          passed_pt[l] += weight;
          pt_eff[l]->Fill(550, weight);
            
          if (fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0){
            passed_eta[l] += weight;
            eta_eff[l]->Fill(2.0, weight);
                
            if (njets ==2) {
              passed_njets[l] += weight;
              njets_eff[l]->Fill(2.0, weight);
        
              if((corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9 ){
                passed_photonveto[l] += weight;
                photonVeto_eff[l]->Fill(0.1, weight);
              
                if(pass_conv_1 && pass_conv_2){
                  passed_conversionveto[l] += weight;
                  conversionVeto_eff[l]->Fill(0.3, weight);
                          
                  if (corrjets.vtx_N >= 2) {
                    passed_nvtx[l] += weight;
                    nvtx_eff[l]->Fill(2.0, weight);
        
                    if(deltajet_phi > 2){
                      passed_deltaphi[l] += weight;
                      deltaphi_eff[l]->Fill(2, weight);
                      
                      if (corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){
                        passed_filters[l] += weight;
                        filters_eff[l]->Fill(1, weight);
                      
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
          }
        }
//       }
		}
	}
	std::cout<<"SIMPs & $m_{\\chi}=1$ GeV & $m_{\\chi}=10$ GeV & $m_{\\chi}=100$ GeV & $m_{\\chi}=200$ GeV & $m_{\\chi}=400$ GeV & $m_{\\chi}=700$ GeV & $m_{\\chi} = 1000$ GeV \\\\"<<std::endl; 
	std::cout<<"total";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<total[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>550$ GeV";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_pt[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$|\\eta_{j1, j2}|<2.0$";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_eta[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"njets";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_njets[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_photonveto[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"conversion veto";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_conversionveto[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"nvtx";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_nvtx[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_deltaphi[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"filters";
	for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_filters[l];
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 12; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2} <$ "<<cut;	
		for(int l = 0; l < 8; l++) std::cout<<" & "<<passed[l][j];
		std::cout<<" \\\\"<<std::endl;
		std::cout<<"ChF$_{j1, j2} <$ (both collections)"<<cut;	
		for(int l = 0; l < 8; l++) std::cout<<" & "<<passed_SPV[l][j];
		std::cout<<" \\\\"<<std::endl;
	}	
	output->Write();
	output->Close();	
}
