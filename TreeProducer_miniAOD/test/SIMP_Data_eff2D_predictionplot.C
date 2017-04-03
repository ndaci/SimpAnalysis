#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TChain.h>
#include <TH2.h>

// #include "lists/list_JetHT_2016B.h"
// #include "lists/list_JetHT_2016C.h"
// #include "lists/list_JetHT_2016D.h"
// #include "lists/list_JetHT_2016E.h"
// #include "lists/list_JetHT_2016F.h"
// #include "lists/list_JetHT_2016G.h"

#include "lists/list_JetHT_rereco_2016B.h"
#include "lists/list_JetHT_rereco_2016C.h"
#include "lists/list_JetHT_rereco_2016D.h"
#include "lists/list_JetHT_rereco_2016E.h"
#include "lists/list_JetHT_rereco_2016F.h"
#include "lists/list_JetHT_rereco_2016G.h"
#include "lists/list_JetHT_rereco_2016Hv2.h"
#include "lists/list_JetHT_rereco_2016Hv3.h"

void SIMP_Data_eff2D_predictionplot(){
	
  TChain * chain = new TChain("tree/SimpAnalysis");
// 	list_JetHT_2016B(chain);
// 	list_JetHT_2016C(chain);
// 	list_JetHT_2016D(chain);
// 	list_JetHT_2016E(chain);
// 	list_JetHT_2016F(chain);
// 	list_JetHT_2016G(chain);
	list_JetHT_rereco_2016B(chain);
	list_JetHT_rereco_2016C(chain);
	list_JetHT_rereco_2016D(chain);
	list_JetHT_rereco_2016E(chain);
	list_JetHT_rereco_2016F(chain);
	list_JetHT_rereco_2016G(chain);
	list_JetHT_rereco_2016Hv2(chain);
	list_JetHT_rereco_2016Hv3(chain);
	std::cout<<"TChain ready"<<std::endl;
  
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8];
	double track_pt[10], track_ptError[10];
	int nPixHits[10];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	int dijet_170;
	double pswgt_dijet_170;
	
	double chf_cuts[201] = {0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 0.36, 0.365, 0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 0.4, 0.405, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435, 0.44, 0.445, 0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49, 0.495, 0.5, 0.505, 0.51, 0.515, 0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.6, 0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.635, 0.64, 0.645, 0.65, 0.655, 0.66, 0.665, 0.67, 0.675, 0.68, 0.685, 0.69, 0.695, 0.7, 0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.0};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  
  TFile *output = new TFile("eff2D_Data_Rereco_exclusiveChF.root", "RECREATE");
	TH2D* total = new TH2D("total", "total", 4, eta_bins, 9, pt_bins);
	total->GetYaxis()->SetTitle("p_{T}");
	total->GetXaxis()->SetTitle("#eta");
	total->Sumw2();
	TH2D* passed[200];
	TH2D* eff[200];
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
	
	std::cout<<"CHF cuts: ";
	for(int j = 0; j < 200; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		double dbl2 = chf_cuts[j+1];
		strs << dbl << "To" << dbl2;
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
		eff[j]->Sumw2();
	}
	std::cout<<std::endl;
		
	chain->SetBranchAddress("jet_pt", &jet_pt);
	chain->SetBranchAddress("jet_eta", &jet_eta);
	chain->SetBranchAddress("jet_phi", &jet_phi);
	chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
	chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
	chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
	chain->SetBranchAddress("track_nPixHits", &nPixHits);
	chain->SetBranchAddress("photon_pt", &photon_pt);
	chain->SetBranchAddress("photon_eta", &photon_eta);
	chain->SetBranchAddress("photon_phi", &photon_phi);
	chain->SetBranchAddress("photon_passLooseId",&passLooseId);
	chain->SetBranchAddress("photon_passMediumId",&passMediumId);
	chain->SetBranchAddress("photon_passTightId",&passTightId);
// 		chain->SetBranchAddress("track_pt", &track_pt);
// 		chain->SetBranchAddress("track_ptError", &track_ptError);
	
	Int_t Nentries = chain->GetEntries(); 
	
	for(Int_t entry = 0; entry < Nentries; ++entry){
		if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
		chain->GetEntry(entry);
		
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
		
		if (dijet_170 == 1 && jet_pt[0] > 250.0 && jet_pt[1] > 250.0 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && nPixHits > 0 && (passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
			if(CHEF_jet[0] > 0.5){
				total->Fill(fabs(jet_eta[1]), jet_pt[1], pswgt_dijet_170);
				for(int j = 0; j < 200; j++){
					if (chf_cuts[j] < CHEF_jet[1] && CHEF_jet[1]<chf_cuts[j+1]) passed[j]->Fill(fabs(jet_eta[1]), jet_pt[1], pswgt_dijet_170);
				}
			}
			if(CHEF_jet[1] > 0.5){
				total->Fill(fabs(jet_eta[0]), jet_pt[0], pswgt_dijet_170);
				for(int j = 0; j < 200; j++){
					if (chf_cuts[j] < CHEF_jet[0] && CHEF_jet[0]<chf_cuts[j+1]) passed[j]->Fill(fabs(jet_eta[0]), jet_pt[0], pswgt_dijet_170);
				}
			}
		}    
	}
	for(int j = 0; j < 200; j++){
		eff[j]->Divide(passed[j], total, 1, 1, "b");
	}
  output->Write();
  output->Close();
}