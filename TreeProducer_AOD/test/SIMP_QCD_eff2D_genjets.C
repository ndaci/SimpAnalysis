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

void SIMP_QCD_eff2D_genjets(){
  
    TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_300To500(chain0);
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_0.root");
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_1.root");
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_500To700(chain1);
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_0.root");
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_1.root");
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_700To1000(chain2);
  chain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17_ext.root");
  chain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1000To1500(chain3);
  chain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17_ext.root");
  chain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1500To2000(chain4);
  chain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17_ext.root");
  chain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_2000ToInf(chain5);
  chain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
  chain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17_ext.root");
//   chain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   chain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_300To500(SPVchain0);
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_0.root");
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_1.root");
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_500To700(SPVchain1);
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_0.root");
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_1.root");
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_700To1000(SPVchain2);
  SPVchain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17_ext.root");
  SPVchain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1000To1500(SPVchain3);
  SPVchain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17_ext.root");
  SPVchain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1500To2000(SPVchain4);
  SPVchain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17_ext.root");
  SPVchain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_2000ToInf(SPVchain5);
  SPVchain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
  SPVchain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17_ext.root");
//   SPVchain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   SPVchain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* SPVchains[6] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5};
	std::cout<<"TChains ready"<<std::endl;
  
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_QCD_GenJets_PUMoriond17_AOD_skim.root", "RECREATE");
	TH2D* total = new TH2D("total", "total", 4, eta_bins, 9, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[11] = {0,0,0,0,0,0,0,0,0,0,0};
	TH2D* eff[11]    = {0,0,0,0,0,0,0,0,0,0,0};
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
	double Nentries[6] = {36948578, 44029226, 29808140, 10360193, 7868538, 4047360};
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
		
		Int_t entries = corrjets.fChain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries[l];
		std::cout<<"QCD bin "<<l<<": Processing "<<entries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
			SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
			
// 			Double_t random1 = r.Uniform();
			
			double deltajet_phi = corrjets.genjet_phi[0] - corrjets.genjet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			deltajet_phi = fabs(deltajet_phi);
			
			output->cd();
			
			if (corrjets.genjet_pt[0] > 250.0 && corrjets.genjet_pt[1] > 250.0 && fabs(corrjets.genjet_eta[0]) < 2.0 && fabs(corrjets.genjet_eta[1]) < 2.0 && deltajet_phi > 2.0){
// 				if (track_ptError[0]/track_pt[0] < 0.5){
// 					if(CHEF_jet[0] > 0.5){
					if(corrjets.genjet_efrac_ch[0] > 0.5){
						total->Fill(fabs(corrjets.genjet_eta[1]), corrjets.genjet_pt[1], weight);
						for(int j = 0; j < 11; j++){
// 							if (CHEF_jet[1]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[1]), jet_pt[1], weight);
							if (corrjets.genjet_efrac_ch[1]<chf_cuts[j]) passed[j]->Fill(fabs(corrjets.genjet_eta[1]), corrjets.genjet_pt[1], weight);
						}
					}
// 					if(CHEF_jet[1] > 0.5){
					if(corrjets.genjet_efrac_ch[1] > 0.5){
						total->Fill(fabs(corrjets.genjet_eta[0]), corrjets.genjet_pt[0], weight);
						for(int j = 0; j < 11; j++){
// 							if (CHEF_jet[0]<chf_cuts[j]) passed[j]->Fill(fabs(jet_eta[0]), jet_pt[0], weight);
							if (corrjets.genjet_efrac_ch[0]<chf_cuts[j]) passed[j]->Fill(fabs(corrjets.genjet_eta[0]), corrjets.genjet_pt[0], weight);
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