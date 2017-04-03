#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "SANtuple.h"

#include "lists/list_JetHT_2016B.h"
#include "lists/list_JetHT_2016C.h"
#include "lists/list_JetHT_2016D.h"
#include "lists/list_JetHT_2016E.h"
#include "lists/list_JetHT_2016F.h"
#include "lists/list_JetHT_2016G.h"

void SIMP_Data_cutflow(){
	  
  TChain* chain = new TChain("treeCorr/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
	list_JetHT_2016B(SPVchain);
	list_JetHT_2016C(SPVchain);
	list_JetHT_2016D(SPVchain);
	list_JetHT_2016E(SPVchain);
	list_JetHT_2016F(SPVchain);
	list_JetHT_2016G(SPVchain);
//
	TFile* output = new TFile("Data_cutflow.root", "RECREATE");
	
  double CHEF_SPVjet[8], CHEF_corrjet[8];
	
	double passed_deltaphi = 0;
	double passed_eta = 0;
	double passed_photonveto = 0;
	double passed_pt = 0;
	double passed[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_SPV[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double chf_cuts[12] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	
  TH1D* pt_eff = new TH1D("pt_eff", "passed pt cuts", 100, 200, 300);
  TH1D* eta_eff = new TH1D("eta_eff", "passed eta cuts", 100, 0, 5);
  TH1D* deltaphi_eff = new TH1D("deltaphi_eff", "passed #Delta#phi cut", 100, 0, 5);
  TH1D* photonVeto_eff = new TH1D("photonVeto_eff", "passed photon veto", 10, 0, 1);
  TH1D* ChF_eff = new TH1D("ChF_eff", "passed ChF cuts", 100, 0, 0.5);
  
	TFile* efficiencies = new TFile("eff2D_Data_Rereco_AOD.root", "READ");
	TH2D* eff_histos[12];
	for(int j = 0; j < 12; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
  
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t Nentries = corrjets.fChain->GetEntries(); 
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
      passed_eta += corrjets.pswgt_dijet_170;
      eta_eff->Fill(2.0, corrjets.pswgt_dijet_170);
      if(corrjets.jet_pt[0] > 250 && corrjets.jet_pt[1] > 250){
        passed_pt += corrjets.pswgt_dijet_170;
        pt_eff->Fill(250, corrjets.pswgt_dijet_170);
        if(deltajet_phi > 2){
          passed_deltaphi += corrjets.pswgt_dijet_170;
          deltaphi_eff->Fill(2, corrjets.pswgt_dijet_170);
          if(corrjets.photon_passLooseId[0] == 0 || (corrjets.photon_passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)){
            passed_photonveto += corrjets.pswgt_dijet_170;
            photonVeto_eff->Fill(0.1, corrjets.pswgt_dijet_170);
            for(int j = 0; j < 12; j++){
              double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
              double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
              passed[j] += corrjets.pswgt_dijet_170*eff1*eff2;
            }
          }
        }
      }
    }
	}
	std::cout<<"$|\\eta_{j1, j2}|<2.0$";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_eta;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>250$ GeV";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_pt;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_deltaphi;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto";
	for(int l = 0; l < 7; l++) std::cout<<" & "<<passed_photonveto;
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2} <$ "<<cut;	
		for(int l = 0; l < 7; l++) std::cout<<" & "<<passed[j];
		std::cout<<" \\\\"<<std::endl;
	}	
	output->Write();
	output->Close();	
}