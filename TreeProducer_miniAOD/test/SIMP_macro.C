#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "list_JetHT_2016B.h"
#include "list_JetHT_2016C.h"
#include "list_JetHT_2016D.h"
#include "list_JetHT_2016E.h"
#include "list_JetHT_2016F.h"
#include "list_JetHT_2016G.h"

void SIMP_macro(){
  
  TChain * chain = new TChain("tree/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
	std::cout<<"TChain ready"<<std::endl;
  
	int dijet_170, vtx_N;
	double pswgt_dijet_170;
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
  
  chain->SetBranchAddress("vtx_N", &vtx_N);
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_phi", &jet_phi);
  chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
  chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
  chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
  
  TFile *output = new TFile("control_region.root", "RECREATE");
	
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 100, 0, 2000);
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *nvtx_withps = new TH1F("nvtx_withps", "Number of vertices (weighted for prescale)", 51, -0.5, 50.5);
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet1_pt_withps = new TH1F("jet1_pt_withps", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt_withps = new TH1F("jet2_pt_withps", "Subleading jet pt", 100, 0, 2000);
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
	TH1F *jet1_chf_withps = new TH1F("jet1_chf_withps", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf_withps = new TH1F("jet2_chf_withps", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *chf = new TH1F("chf", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsCHF = new TH2D("CHFvsCHF", "Leading jet ChF vs. Subleading jet ChF", 100, 0, 1, 100, 0, 1);
	CHFvsCHF->GetXaxis()->SetTitle("jet2 ChF");
	CHFvsCHF->GetYaxis()->SetTitle("jet1 ChF");
	TH1D *CHF = new TH1D("CHF", "Jet ChF (both jets)", 100, 0, 1);
  
  Int_t Nentries = chain->GetEntries(); 
	std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
  for(Int_t entry = 0; entry < Nentries; ++entry){
    chain->GetEntry(entry);
    
    double deltajet_phi = jet_phi[0] - jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
    
		for (int i = 0; i < 4; i++){
			CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
		}
		
		output->cd();
	
		if (dijet_170 == 1 && jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2){
			if (CHEF_jet[0] > 0.5||CHEF_jet[1] > 0.5){
				HT->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3]);
				nvtx->Fill(vtx_N);
				nvtx_withps->Fill(vtx_N, pswgt_dijet_170);
				jet1_pt->Fill(jet_pt[0]);
				jet2_pt->Fill(jet_pt[1]);
				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2]);
				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3]);
				jet1_pt_withps->Fill(jet_pt[0], pswgt_dijet_170);
				jet2_pt_withps->Fill(jet_pt[1], pswgt_dijet_170);
				jet1_eta->Fill(jet_eta[0]);
				jet2_eta->Fill(jet_eta[1]);
				jet1_phi->Fill(jet_phi[0]);
				jet2_phi->Fill(jet_phi[1]);
				deltaphi->Fill(deltajet_phi);
				jet1_chf->Fill(CHEF_jet[0]);
				if (CHEF_jet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_jet[0]);
				jet2_chf->Fill(CHEF_jet[1]);
				jet1_chf_withps->Fill(CHEF_jet[0], pswgt_dijet_170);
				jet2_chf_withps->Fill(CHEF_jet[1], pswgt_dijet_170);
				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2]);
				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3]);
				CHFvsPT->Fill(jet_pt[0], CHEF_jet[0]);
				CHFvsPT->Fill(jet_pt[1], CHEF_jet[1]);
				CHF->Fill(CHEF_jet[0]);
				CHF->Fill(CHEF_jet[1]);
				CHFvsCHF->Fill(CHEF_jet[1], CHEF_jet[0]);
				if(CHEF_jet[0] > CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[0]);
					jetB_pt->Fill(jet_pt[1]);
					jetA_eta->Fill(jet_eta[0]);
					jetB_eta->Fill(jet_eta[1]);
					if (CHEF_jet[0] > 0.5) chf->Fill(CHEF_jet[1]);
				}
				else if(CHEF_jet[0] < CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[1]);
					jetB_pt->Fill(jet_pt[0]);
					jetA_eta->Fill(jet_eta[1]);
					jetB_eta->Fill(jet_eta[0]);
					if (CHEF_jet[1] > 0.5) chf->Fill(CHEF_jet[0]);
				}
			}
		}    
	}
  output->Write();
  output->Close();
}