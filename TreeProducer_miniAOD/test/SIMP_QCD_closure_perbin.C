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

void SIMP_QCD_closure_perbin(double pt_min, double pt_max, TString outputname){
// 	
// 	double pt_min = 275;
// 	double pt_max = 300;
  
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
  
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
// 	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
// 	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	
	double passed_eff[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff1[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff2[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double passed_effboth[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double passed_MCtruth[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_MCtruth[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_eff1[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_eff2[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_eff[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_effboth[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double zero[11] = {0,0,0,0,0,0,0,0,0,0,0};
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_QCD_random.root", "READ");
	TH2D* eff_histos[11];
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
	std::cout<<"done"<<std::endl;
	
  TFile *output = new TFile(outputname, "RECREATE");
	
	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
  
	for (int l = 0; l < 6; l++){
		TChain* chain = chains[l];
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		
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
				if(jet_pt[0] > pt_min && jet_pt[1] > pt_min && jet_pt[0] < pt_max && jet_pt[1] < pt_max){
					for(int j = 0; j < 11; j++){
						if (CHEF_jet[0]<chf_cuts[j] && CHEF_jet[1]<chf_cuts[j]){
							passed_MCtruth[j]+= weight;
							err_MCtruth[j] += pow(weight, 2);
						}
						double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
						double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
						if (CHEF_jet[1]<chf_cuts[j]) passed_eff1[j]+=weight*eff1;
						if (CHEF_jet[0]<chf_cuts[j]) passed_eff2[j]+=weight*eff2;
						passed_effboth[j]+=eff1*eff2*weight;	
					
// 					err_eff1[j] += ;
// 					err_eff2[j] += ;
// 					err_effboth[j] += ;
					}
				}
			}    
		}
	}
	for(int j = 0; j < 11; j++){
		err_MCtruth[j] = TMath::Sqrt(err_MCtruth[j]);
		err_eff[j] = (err_eff1[j]+err_eff2[j])/2.0 ;
		passed_eff[j] = (passed_eff1[j]+passed_eff2[j])/2;
	}
	TCanvas *c1 = new TCanvas("Closure test", "Closure test");
	c1->cd();
	TGraphErrors *MC = new TGraphErrors(11, chf_cuts, passed_MCtruth, zero, err_MCtruth);
	MC->SetTitle("Closure test");
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->SetMarkerStyle(20);
	MC->Draw("AP");
	TGraphErrors *oneleg = new TGraphErrors(11, chf_cuts, passed_eff, zero, err_eff);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P");
	TGraphErrors *twoleg = new TGraphErrors(11, chf_cuts, passed_effboth, zero, err_effboth);
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
	twoleg->Draw("P");
// 	c1->Print("closure.png");
	output->Append(c1);
  output->Write();
  output->Close();
}