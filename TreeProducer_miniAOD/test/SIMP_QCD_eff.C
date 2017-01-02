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

void SIMP_QCD_eff(){
  
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
  
  TFile *output = new TFile("eff_QCD.root", "RECREATE");
  
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
// 	E[4], M[4];
	
	double total_pt[9] = {0,0,0,0,0,0,0,0,0};
	double total_eta[4] = {0,0,0,0};
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	double passed_pt[9][11];
	double passed_eta[4][11];
	double eff_pt[9][11], red_pt[9][11], eff_eta[4][11], red_eta[4][11];
	double err_eff_pt[9][11], err_eff_eta[9][11], err_red_pt[4][11], err_red_eta[4][11];
	double zero[11] = {0,0,0,0,0,0,0,0,0,0,0};
	
	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 20;
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 9; i++){
			passed_pt[i][j] = 0;
		}
	}
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 4; i++){
			passed_eta[i][j] = 0;
		}
	}
  
	for (int j = 0; j < 6; j++){
		TChain* chain = chains[j];		
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
	//   chain->SetBranchAddress("jet_e", &E);
	//   chain->SetBranchAddress("jet_m", &M);
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
			
			if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2){
				
				if(CHEF_jet[0] > CHEF_jet[1] && CHEF_jet[0] > 0.5){
					for(int j = 0; j < 11; j++){
						for (int i = 0; i < 9; i++){
							if (jet_pt[1]>pt_bins[i] && jet_pt[1]<pt_bins[i+1]){
								total_pt[i] += weight;
								if (CHEF_jet[1]<chf_cuts[j]) passed_pt[i][j] += weight;
							}
						}
						for (int i = 0; i < 4; i++){
							if (jet_eta[1]>eta_bins[i] && jet_eta[1]<eta_bins[i+1]){
								total_eta[i] += weight;
								if (CHEF_jet[1]<chf_cuts[j]) passed_eta[i][j] += weight;
							}
						}
					}
				}
				else if(CHEF_jet[0] < CHEF_jet[1] && CHEF_jet[1] > 0.5){
					for(int j = 0; j < 11; j++){
						for (int i = 0; i < 9; i++){
							if (jet_pt[0]>pt_bins[i] && jet_pt[0]<pt_bins[i+1]){
								total_pt[i] += weight;
								if (CHEF_jet[0]<chf_cuts[j]) passed_pt[i][j] += weight;
							}
						}
						for (int i = 0; i < 4; i++){
							if (jet_eta[0]>eta_bins[i] && jet_eta[0]<eta_bins[i+1]){
								total_eta[i] += weight;
								if (CHEF_jet[0]<chf_cuts[j]) passed_eta[i][j] += weight;
							}
						}
					}
				}
			}    
		}
	}
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 9; i++){
			double tot = (double)total_pt[i];
			std::cout<<"pt "<<i<<": "<<tot<<std::endl;
			double passed = (double)passed_pt[i][j];
			eff_pt[i][j] = passed/tot;
			err_eff_pt[i][j] = pow(tot, -2)*TMath::Sqrt(pow(tot - passed, 2)*passed + (pow(passed, 2)*(tot - passed)));
			if (passed_pt[i][j] != 0){
				red_pt[i][j] = tot/(double)passed;
				err_red_pt[i][j] = pow(passed, -2)*TMath::Sqrt(pow(tot - passed, 2)*passed + (pow(passed, 2)*(tot - passed)));
			}
			else{
				red_pt[i][j] = 0;
				err_red_pt[i][j] = 0;
			}
// 			std::cout<<passed_pt[i][j]<<"  "<<eff_pt[i][j]<<std::endl;
// 			std::cout<<total<<"/"<<passed_pt[i][j]<<"  "<<red_pt[i][j]<<std::endl;
		}
	}
		for (int j = 0; j < 11; j++){
		for (int i = 0; i < 4; i++){
			double tot = (double)total_eta[i];
			std::cout<<"eta "<<i<<": "<<tot<<std::endl;
			double passed = (double)passed_eta[i][j];
			eff_eta[i][j] = passed/tot;
			err_eff_eta[i][j] = pow(tot, -2)*TMath::Sqrt(pow(tot - passed, 2)*passed + (pow(passed, 2)*(tot - passed)));
			if (passed_eta[i][j] != 0){
				red_eta[i][j] = tot/passed;
				err_red_eta[i][j] = pow(passed, -2)*TMath::Sqrt(pow(tot - passed, 2)*passed + (pow(passed, 2)*(tot - passed)));
			}
			else{
				red_eta[i][j] = 0;
				err_red_eta[i][j] = 0;
				
			}
		}
	}	
	
	TCanvas *c1 = new TCanvas("efficiency_pt", "efficiency_pt");
	c1->cd();
	TLegend *leg1 = new TLegend();
	TGraphErrors *firstptbin = new TGraphErrors(11, chf_cuts, eff_pt[0], zero, err_eff_pt[0]);
	firstptbin->SetTitle("CHF cut efficiency");
	firstptbin->GetXaxis()->SetTitle("CHF cut");
	firstptbin->GetYaxis()->SetTitle("efficiency");
	firstptbin->SetMarkerStyle(20);
	firstptbin->Draw("AP");
	leg1->AddEntry(firstptbin, "250 < pt < 400", "lep");
	for (int i = 1; i < 9; i++){
		TGraphErrors *otherptbin = new TGraphErrors(11, chf_cuts, eff_pt[i], zero, err_eff_pt[i]);
		otherptbin->SetMarkerStyle(20);
		otherptbin->SetMarkerColor(i+1);
		otherptbin->Draw("P");
		leg1->AddEntry(otherptbin, " < pt < ", "lep");
	}
	c1->cd();
	leg1->Draw();
	c1->Print("QCD_ptbinned_eff.png");
	
	TCanvas *c2 = new TCanvas("rejection_pt", "rejection_pt");
	c2->cd();
	TGraphErrors *firstptbin2 = new TGraphErrors(11, chf_cuts, red_pt[0], zero, err_red_pt[0]);
	firstptbin2->SetTitle("CHF cut rejection factor");
	firstptbin2->GetXaxis()->SetTitle("CHF cut");
	firstptbin2->GetYaxis()->SetTitle("rejection factor");
	firstptbin2->SetMarkerStyle(20);
	firstptbin2->Draw("AP");
	for (int i = 1; i < 9; i++){
		TGraphErrors *otherptbin2 = new TGraphErrors(11, chf_cuts, red_pt[i], zero, err_red_pt[i]);
		otherptbin2->SetMarkerStyle(20);
		otherptbin2->SetMarkerColor(i+1);
		otherptbin2->Draw("P");
	}
	c2->Print("QCD_ptbinned_red.png");
	
	TCanvas *c3 = new TCanvas("efficiency_eta", "efficiency_eta");
	c3->cd();
	TGraphErrors *firstetabin = new TGraphErrors(11, chf_cuts, eff_eta[0], zero, err_eff_eta[0]);
	firstetabin->SetTitle("CHF cut efficiency");
	firstetabin->GetXaxis()->SetTitle("CHF cut");
	firstetabin->GetYaxis()->SetTitle("efficiency");
	firstetabin->SetMarkerStyle(20);
	firstetabin->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *otheretabin = new TGraphErrors(11, chf_cuts, eff_eta[i], zero, err_eff_eta[i]);
		otheretabin->SetMarkerStyle(20);
		otheretabin->SetMarkerColor(i+1);
		otheretabin->Draw("P");
	}
	c3->Print("QCD_etabinned_eff.png");
	
	TCanvas *c4 = new TCanvas("rejection_eta", "rejection_eta");
	c4->cd();
	TGraphErrors *firstetabin2 = new TGraphErrors(11, chf_cuts, red_eta[0], zero, err_red_eta[0]);
	firstetabin2->SetTitle("CHF cut rejection factor");
	firstetabin2->GetXaxis()->SetTitle("CHF cut");
	firstetabin2->GetYaxis()->SetTitle("rejection factor");
	firstetabin2->SetMarkerStyle(20);
	firstetabin2->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *otheretabin2 = new TGraphErrors(11, chf_cuts, red_eta[i], zero, err_red_eta[i]);
		otheretabin2->SetMarkerStyle(20);
		otheretabin2->SetMarkerColor(i+1);
		otheretabin2->Draw("P");
	}
	c4->Print("QCD_etabinned_red.png");
	
	output->Append(c1);
	output->Append(c2);
	output->Append(c3);
	output->Append(c4);
	
// 	for (int j = 0; j < 11; j++){
// 		TH2D* histo = new TH2D("histo", "histo", 2, eta_bins, 9, pt_bins);
// 		for (int i = 0; i < 9; i++){
// 			histo->Fill(eff_eta[0][i]);
// 		}
// 		output->Append(histo);
// 	}
	
  output->Write();
  output->Close();
}