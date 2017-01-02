#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "list_QCD_300To500_PUMoriond17.h"
#include "list_QCD_500To700_PUMoriond17.h"
// #include "list_QCD_700To1000_PUMoriond17.h"
#include "list_QCD_1000To1500_PUMoriond17.h"
#include "list_QCD_1500To2000_PUMoriond17.h"
// #include "list_QCD_2000ToInf_PUMoriond17.h"

void SIMP_QCD_test(){
  
  TChain* chain0 = new TChain("tree/SimpAnalysis");
	list_QCD_300To500(chain0);
  TChain* chain1 = new TChain("tree/SimpAnalysis");
	list_QCD_500To700(chain1);
  TChain* chain2 = new TChain("tree/SimpAnalysis");
// 	list_QCD_700To1000(chain2);
  TChain* chain3 = new TChain("tree/SimpAnalysis");
	list_QCD_1000To1500(chain3);
  TChain* chain4 = new TChain("tree/SimpAnalysis");
	list_QCD_1500To2000(chain4);
  TChain* chain5 = new TChain("tree/SimpAnalysis");
// 	list_QCD_2000ToInf(chain5);
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
	std::cout<<"TChains ready"<<std::endl;
	
  int dijet_170;
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
// 	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
// 	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	
	double passed_Min0p5[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double passed_Max0p2[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_Min0p5[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_Max0p2[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double zero[11] = {0,0,0,0,0,0,0,0,0,0,0};
	
  TFile *output = new TFile("QCD_plot_lowpt_HT-500To700_PUMoriond17.root", "RECREATE");
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};//Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
  
	for (int l = 1; l < 2; l++){
		TChain* chain = chains[l];
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
		
		Int_t Nentries = chain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"QCD bin "<<l<<": Processing "<<Nentries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			chain->GetEntry(entry);
			
			double deltajet_phi = jet_phi[0] - jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			
			for (int i = 0; i < 4; i++){
				CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2){
				for(int j = 0; j < 11; j++){
					if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] > 0.5 && CHEF_jet[1] < chf_cuts[j]){
						passed_Min0p5[j] += weight;
						err_Min0p5[j] +=weight*weight;
					}
					if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] > 0.5 && CHEF_jet[0] < chf_cuts[j]){
						passed_Min0p5[j] += weight;
						err_Min0p5[j] +=weight*weight;
					}
					if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] < 0.2 && CHEF_jet[1] < chf_cuts[j]){
						passed_Max0p2[j] += weight;
						err_Max0p2[j] +=weight*weight;
					}
					if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] < 0.2 && CHEF_jet[0] < chf_cuts[j]){
						passed_Max0p2[j] += weight;
						err_Max0p2[j] +=weight*weight;
					}
				}
			}    
		}
	}
	for(int j = 0; j < 11; j++){
		err_Min0p5[j] = TMath::Sqrt(err_Min0p5[j]);
		err_Max0p2[j] = TMath::Sqrt(err_Max0p2[j]);
	}
	TCanvas *c1 = new TCanvas("Closure test", "Closure test");
	c1->cd();
	TGraphErrors *MC = new TGraphErrors(11, chf_cuts, passed_Min0p5, zero, err_Min0p5);
	MC->SetTitle("Closure test");
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->SetMarkerStyle(20);
	MC->Draw("AP");
	TGraphErrors *oneleg = new TGraphErrors(11, chf_cuts, passed_Max0p2, zero, err_Max0p2);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P");
	output->Append(c1);
  output->Write();
  output->Close();
}