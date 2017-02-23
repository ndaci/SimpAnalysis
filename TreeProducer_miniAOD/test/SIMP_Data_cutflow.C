#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "lists/list_JetHT_rereco_2016B.h"
#include "lists/list_JetHT_rereco_2016C.h"
#include "lists/list_JetHT_rereco_2016D.h"
#include "lists/list_JetHT_rereco_2016E.h"
#include "lists/list_JetHT_rereco_2016F.h"
#include "lists/list_JetHT_rereco_2016G.h"

void SIMP_Data_cutflow(){
	
  TChain * chain = new TChain("tree/SimpAnalysis");
	list_JetHT_rereco_2016B(chain);
	list_JetHT_rereco_2016C(chain);
	list_JetHT_rereco_2016D(chain);
	list_JetHT_rereco_2016E(chain);
	list_JetHT_rereco_2016F(chain);
	list_JetHT_rereco_2016G(chain);
  
	TFile* output = new TFile("Data_Rereco_cutflow.root", "RECREATE");
	
	int nJet, dijet_170, dijet_170_0p1, dijet_220_0p3, dijet_330_0p5, dijet_430;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8], EMF_jet[8];
	double track_pt[10], track_ptError[10], track_dzError[10], track_dz[10];
	int track_fromPV[10], nhits[10], nPixHits[10], ndof[3];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
  double pswgt_dijet_170;
	
	double passed_deltaphi = 0;
	double passed_eta = 0;
	double passed_npix = 0;
	double passed_photonveto = 0;
	double passed_pt = 0;
	double passed[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
  
  TH1D* pt_eff = new TH1D("pt_eff", "passed pt cuts", 100, 200, 300);
  TH1D* eta_eff = new TH1D("eta_eff", "passed eta cuts", 100, 0, 5);
  TH1D* deltaphi_eff = new TH1D("deltaphi_eff", "passed #Delta#phi cut", 100, 0, 5);
  TH1D* nPix_eff = new TH1D("nPix_eff", "passed nPixHits cut", 10, 0, 10);
  TH1D* photonVeto_eff = new TH1D("photonVeto_eff", "passed photon veto", 10, 0, 1);
  TH1D* ChF_eff = new TH1D("ChF_eff", "passed ChF cuts", 100, 0, 0.5);
  
	TFile* efficiencies = new TFile("eff2D_photonVeto_Data_Rereco.root", "READ");
	TH2D* eff_histos[11];
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}

  chain->SetBranchAddress("nJet",&nJet);
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_phi", &jet_phi);
  chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
  chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
  chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("track_ptError", &track_ptError);
  chain->SetBranchAddress("track_nhits", &nhits);
  chain->SetBranchAddress("track_nPixHits", &nPixHits);
  chain->SetBranchAddress("track_pt", &track_pt);
  chain->SetBranchAddress("track_fromPV", &track_fromPV);
  chain->SetBranchAddress("track_dzError", &track_dzError);
  chain->SetBranchAddress("track_dz", &track_dz);
  chain->SetBranchAddress("photon_pt", &photon_pt);
  chain->SetBranchAddress("photon_eta", &photon_eta);
  chain->SetBranchAddress("photon_phi", &photon_phi);
  chain->SetBranchAddress("photon_passLooseId",&passLooseId);
  chain->SetBranchAddress("photon_passMediumId",&passMediumId);
  chain->SetBranchAddress("photon_passTightId",&passTightId);
  chain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &dijet_170_0p1);
  chain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &dijet_220_0p3);
  chain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &dijet_330_0p5);
  chain->SetBranchAddress("HLT_DiCentralPFJet430", &dijet_430);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
  
  Int_t Nentries = chain->GetEntries(); 
  std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
  for(Int_t entry = 0; entry < Nentries; ++entry){
    chain->GetEntry(entry);
    double weight = pswgt_dijet_170;
    
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
      CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
    } 
    
    output->cd();
    
    if (dijet_170 == 1 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0){
      passed_eta += weight;
      eta_eff->Fill(2.0, weight);
      if(jet_pt[0] > 250 && jet_pt[1] > 250){
        passed_pt += weight;
        pt_eff->Fill(250, weight);
        if(deltajet_phi > 2){
          passed_deltaphi += weight;
          deltaphi_eff->Fill(2, weight);
          if(nPixHits[0] > 0){
            passed_npix += weight;
            nPix_eff->Fill(0.0, weight);
            if(passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)){
              passed_photonveto += weight;
              photonVeto_eff->Fill(0.1, weight);
              for(int j = 0; j < 11; j++){
                double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
                double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
								passed[j] += weight*eff1*eff2;
								ChF_eff->Fill(chf_cuts[j], weight*eff1*eff2);
              }
            }
          }
        }
      }
    }
  }
	std::cout<<"$|\\eta_{j1, j2}|<2.0$ & "<<passed_eta;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>250$ GeV & "<<passed_pt;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$ & "<<passed_deltaphi;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$# pixel hits > 0$ & "<<passed_npix;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto & "<<passed_photonveto;
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2}$ < "<<cut;	
		std::cout<<" & "<<passed[j];
		std::cout<<" \\\\"<<std::endl;
	}
	
	output->Write();
	output->Close();	
}