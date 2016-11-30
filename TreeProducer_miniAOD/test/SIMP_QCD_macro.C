#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "list_QCD_300To500.h"
#include "list_QCD_500To700.h"
#include "list_QCD_700To1000.h"
#include "list_QCD_1000To1500.h"
#include "list_QCD_1500To2000.h"
#include "list_QCD_2000ToInf.h"

void SIMP_QCD_macro(){
  
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
  
  TFile *output = new TFile("plots_QCD_comparetodata.root", "RECREATE");
  
	int dijet_170, nJet;
	double pswgt_dijet_170;
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
// 	E[4], M[4];
  
	TH1F *njets = new TH1F("njets", "Number of jets", 11, -0.5, 10.5);
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 100, 0, 2000);
	TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 100, 0, 2000);
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *jetA_pt = new TH1F("jetA_pt", "ChF-leading jet pt", 100, 0, 2000);
	TH1F *jetB_pt = new TH1F("jetB_pt", "ChF-subleading jet pt", 100, 0, 2000);
	TH1F *jetA_eta = new TH1F("jetA_eta", "ChF-leading jet eta", 100, -3.14, 3.14);
	TH1F *jetB_eta = new TH1F("jetB_eta", "ChF-subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *chf = new TH1F("chf", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsCHF = new TH2D("CHFvsCHF", "Leading jet ChF vs. Subleading jet ChF", 100, 0, 1, 100, 0, 1);
	CHFvsCHF->GetXaxis()->SetTitle("jet2 ChF");
	CHFvsCHF->GetYaxis()->SetTitle("jet1 ChF");
	TH1D *CHF = new TH1D("CHF", "Jet ChF (both jets)", 100, 0, 1);
	
	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
  
	for (int j = 0; j < 6; j++){
		TChain* chain = chains[j];		
		chain->SetBranchAddress("nJet",&nJet);
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		
		Int_t Nentries = chain->GetEntries(); 
		std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
		
		double weight = QCD_xsec[j]*lumi/Nentries;
		std::cout<<"weight "<<j<<": "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			chain->GetEntry(entry);
			
			double deltajet_phi = jet_phi[0] - jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			
			for (int i = 0; i < 4; i++){
				CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2 && (CHEF_jet[0] > 0.5||CHEF_jet[1] > 0.5)){
				njets->Fill(nJet, weight);
				HT->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3], weight);
				HT_nowgt->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3]);
				jet1_pt->Fill(jet_pt[0], weight);
				jet2_pt->Fill(jet_pt[1], weight);
				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], weight);
				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], weight);
				jet1_eta->Fill(jet_eta[0], weight);
				jet2_eta->Fill(jet_eta[1], weight);
				jet1_phi->Fill(jet_phi[0], weight);
				jet2_phi->Fill(jet_phi[1], weight);
				deltaphi->Fill(deltajet_phi, weight);
				jet1_chf->Fill(CHEF_jet[0], weight);
				if (CHEF_jet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_jet[0], weight);
				jet2_chf->Fill(CHEF_jet[1], weight);
				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2], weight);
				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3], weight);
				CHFvsPT->Fill(jet_pt[0], CHEF_jet[0], weight);
				CHFvsPT->Fill(jet_pt[1], CHEF_jet[1], weight);
				CHF->Fill(CHEF_jet[0], weight);
				CHF->Fill(CHEF_jet[1], weight);
				CHFvsCHF->Fill(CHEF_jet[1], CHEF_jet[0], weight);
				if(CHEF_jet[0] > CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[0], weight);
					jetB_pt->Fill(jet_pt[1], weight);
					jetA_eta->Fill(jet_eta[0], weight);
					jetB_eta->Fill(jet_eta[1], weight);
					if (CHEF_jet[0] > 0.5) chf->Fill(CHEF_jet[1], weight);
				}
				else if(CHEF_jet[0] < CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[1], weight);
					jetB_pt->Fill(jet_pt[0], weight);
					jetA_eta->Fill(jet_eta[1], weight);
					jetB_eta->Fill(jet_eta[0], weight);
					if (CHEF_jet[1] > 0.1) chf->Fill(CHEF_jet[0], weight);
				}
			}    
		}
	}
  output->Write();
  output->Close();
}