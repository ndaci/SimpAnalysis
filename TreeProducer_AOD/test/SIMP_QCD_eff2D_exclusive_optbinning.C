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

#include "lists/list_QCD_1000To1500_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_filters.h"
// #include "lists/list_QCD_1000To1500_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1500To2000_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1000To1500_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_1500To2000_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_2000ToInf_ext_PUMoriond17_AOD_newtrigger.h"



void SIMP_QCD_eff2D_exclusive_optbinning(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain0);
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain1);
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain2);
//   TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1000To1500_ext(chain3);
//   TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1500To2000_ext(chain4);
//   TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_2000ToInf_ext(chain5);
	TChain* chains[3] = {chain0, chain1, chain2};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain0);
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain1);
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain2);
//   TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1000To1500_ext(SPVchain3);
//   TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1500To2000_ext(SPVchain4);
//   TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_2000ToInf_ext(SPVchain5);
	TChain* SPVchains[3] = {SPVchain0, SPVchain1, SPVchain2};
	std::cout<<"TChains ready"<<std::endl;
	  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
	
	double chf_cuts[25] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0};
	double pt_bins[7] = {550, 600, 700, 800, 950, 1200, 10000};
	double eta_bins[9] = {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
  
  TString name = "eff2D_QCD_exclusive_filters.root";
  TFile *output = new TFile(name, "RECREATE");
  std::cout<<"output file: "<<name<<std::endl;
	TH2D* total = new TH2D("total", "total", 8, eta_bins, 6, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	TH2D* eff[24]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //Spring16
	double QCD_xsec[6] = {/*346400, 32010, 6842,*/ 1203, 120.1, 25.40, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
	
	std::cout<<"CHF cuts: ";
	for(int j = 0; j < 24; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		double dbl2 = chf_cuts[j+1];
		strs << dbl << "To" << dbl2;
		std::string cut = strs.str();
		std::cout<<" "<<cut<<" ";
		std::string title_passed = "passed_"+cut;
		std::string title_eff = "eff_"+cut;
		passed[j] = new TH2D(title_passed.c_str(), title_passed.c_str(), 8, eta_bins, 6, pt_bins);
		passed[j]->GetYaxis()->SetTitle("p_{T}");
		passed[j]->GetXaxis()->SetTitle("#eta");
		passed[j]->Sumw2();
		eff[j] = new TH2D(title_eff.c_str(), title_eff.c_str(), 8, eta_bins, 6, pt_bins);
		eff[j]->GetYaxis()->SetTitle("p_{T}");
		eff[j]->GetXaxis()->SetTitle("#eta");
	}
	std::cout<<std::endl;
		
	for (int l = 0; l < 3; l++){
		TChain* chain = chains[l];
		TChain* SPVchain = SPVchains[l];		
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t entries = corrjets.fChain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/entries;
		std::cout<<"QCD bin "<<l<<": Processing "<<entries<<" entries with weight "<<weight<<std::endl;

    for(Int_t entry = 0; entry < entries; ++entry){
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
        CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
      } 
        
      bool pass_conv_1 = true;
      bool pass_conv_2 = true;
      if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
      if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
      
      output->cd();
      
      if (corrjets.HLT_PFJet450 == 1 && corrjets.jet_pt[0] > 550.0 && corrjets.jet_pt[1] > 550.0 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9 && corrjets.vtx_N >= 2 && njets == 2 && pass_conv_1 && pass_conv_2 && corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){      
        if(CHEF_corrjet[0] > 0.25){
          total->Fill(fabs(corrjets.jet_eta[1]), corrjets.jet_pt[1], weight);
          for(int j = 0; j < 24; j++){
            if (CHEF_corrjet[1]>=chf_cuts[j] && CHEF_corrjet[1]<chf_cuts[j+1]) passed[j]->Fill(fabs(corrjets.jet_eta[1]), corrjets.jet_pt[1], weight);
          }
        }
        if(CHEF_corrjet[1] > 0.25){
          total->Fill(fabs(corrjets.jet_eta[0]), corrjets.jet_pt[0], weight);
          for(int j = 0; j < 24; j++){
            if (CHEF_corrjet[0]>=chf_cuts[j] && CHEF_corrjet[0]<chf_cuts[j+1]) passed[j]->Fill(fabs(corrjets.jet_eta[0]), corrjets.jet_pt[0], weight);
          }
        }
      } 
    }
  }
	for(int j = 0; j < 24; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}