#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "SANtuple.h"

#include "lists/list_QCD_1000To1500_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_filters.h"
// #include "lists/list_QCD_1000To1500_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1500To2000_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_newtrigger_2.h"

void SIMP_QCD_cutflow(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain0);
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain1);
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain2);
	TChain* chains[3] = {chain0, chain1, chain2};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain0);
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain1);
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain2);
	TChain* SPVchains[3] = {SPVchain0, SPVchain1, SPVchain2};
	std::cout<<"TChains ready"<<std::endl; 
  
	TFile* output = new TFile("QCD_cutflow_withHBHENoiseIsoFilter.root", "RECREATE");
	
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
	
  double passed_trigger = 0;
  double passed_filter = 0;
	double passed_deltaphi = 0;
	double passed_eta = 0;
	double passed_photonveto = 0;
	double passed_conversionveto = 0;
	double passed_NEMF = 0;
	double passed_pt = 0;
  double passed_njets = 0;
  double passed_nvtx = 0;
	double passed[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_SPV[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double chf_cuts[24] = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.0};
	
  TH1D* trigger_eff = new TH1D("trigger_eff", "passed control trigger", 3, 0, 2);
  TH1D* pt_eff = new TH1D("pt_eff", "passed pt cuts", 100, 200, 300);
  TH1D* eta_eff = new TH1D("eta_eff", "passed eta cuts", 100, 0, 5);
  TH1D* deltaphi_eff = new TH1D("deltaphi_eff", "passed #Delta#phi cut", 100, 0, 5);
  TH1D* photonVeto_eff = new TH1D("photonVeto_eff", "passed photon veto", 10, 0, 1);
  TH1D* conversionVeto_eff = new TH1D("conversionVeto_eff", "passed conversion veto", 10, 0, 1);
  TH1D* njets_eff = new TH1D("njets_eff", "passed njet cut", 5, 0, 4);
  TH1D* NEMF_eff = new TH1D("NEMF_eff", "passed NEMF cut", 5, 0, 4);
  TH1D* nvtx_eff = new TH1D("nvtx_eff", "passed nvtx cut", 5, 0, 4);
  TH1D* filter_eff = new TH1D("filter_eff", "passed all filters", 5, 0, 4);   
  TH1D* ChF_eff_1leg = new TH1D("ChF_eff_1leg", "passed ChF cuts", 100, 0, 0.5);
  TH1D* ChF_eff_2leg = new TH1D("ChF_eff_2leg", "passed ChF cuts", 100, 0, 0.5);
  
	TFile* efficiencies = new TFile("eff2D_QCD_filters.root", "READ");
	TH2D* eff_histos[23];
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
	
	double QCD_xsec[3] = {/*346400, 32010, 6842,*/ 1203, 120.1, 25.40}; 
	double lumi = 16.146*1000;
  
	for (int l = 0; l < 3; l++){
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
      if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
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
      }
        
      bool pass_conv_1 = true;
      bool pass_conv_2 = true;
      if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
      if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
      
      output->cd();
      
      if (corrjets.HLT_PFJet450 == 1){
        passed_trigger += weight;
        trigger_eff->Fill(1.0, weight);
          
        if(corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550){
          passed_pt += weight;
          pt_eff->Fill(550, weight);
      
          if (fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0){
            passed_eta += weight;
            eta_eff->Fill(2.0, weight);
                
            if (njets ==2) {                  
              passed_njets += weight;
              njets_eff->Fill(2.0, weight);
          
              if(corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)){              
                passed_photonveto += weight;
                photonVeto_eff->Fill(0.1, weight);
                
                if(pass_conv_1 && pass_conv_2){
                  passed_conversionveto += weight;
                  conversionVeto_eff->Fill(0.3, weight); 
                  
                  if (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9){
                    passed_NEMF += weight;
                    NEMF_eff->Fill(0.9, weight);
                  
                    if (corrjets.vtx_N >= 2) {                    
                      passed_nvtx += weight;
                      nvtx_eff->Fill(2.0, weight);
          
                      if(deltajet_phi > 2){
                        passed_deltaphi += weight;
                        deltaphi_eff->Fill(2, weight);
                        
                        if (corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){              
                          passed_filter += weight;
                          filter_eff->Fill(1.0, weight);  
                      
                          for(int j = 0; j < 23; j++){
                            double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(corrjets.jet_eta[0]), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
                            double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(corrjets.jet_eta[1]), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
    //                         passed[j] += weight*eff1*eff2;
                            if (CHEF_corrjet[1] < chf_cuts[j]) ChF_eff_1leg->Fill(chf_cuts[j], weight*eff1/2);
                            if (CHEF_corrjet[0] < chf_cuts[j]) ChF_eff_1leg->Fill(chf_cuts[j], weight*eff2/2);
                            ChF_eff_2leg->Fill(chf_cuts[j], weight*eff1*eff2);
                            if (CHEF_corrjet[0] < chf_cuts[j] && CHEF_corrjet[1] < chf_cuts[j] && CHEF_SPVjet[0] < chf_cuts[j] && CHEF_SPVjet[1] < chf_cuts[j]) passed[j] += weight;
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
    }
	}
	std::cout<<"trigger";
	 std::cout<<" & "<<passed_trigger;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>550$ GeV";
	 std::cout<<" & "<<passed_pt;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$|\\eta_{j1, j2}|<2.0$";
	 std::cout<<" & "<<passed_eta;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"njets";
	 std::cout<<" & "<<passed_njets;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto";
	 std::cout<<" & "<<passed_photonveto;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"conversion veto";
	 std::cout<<" & "<<passed_conversionveto;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"NEMF cut";
	 std::cout<<" & "<<passed_NEMF;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"nvtx";
	 std::cout<<" & "<<passed_nvtx;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$";
	 std::cout<<" & "<<passed_deltaphi;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"filters";
	 std::cout<<" & "<<passed_filter;
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2} <$ "<<cut;	
		 std::cout<<" & "<<passed[j];
		std::cout<<" \\\\"<<std::endl;
	}	
	output->Write();
	output->Close();	
}