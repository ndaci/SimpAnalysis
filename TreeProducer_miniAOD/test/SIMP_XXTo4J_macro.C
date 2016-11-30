#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

void SIMP_XXTo4J_macro(){
  
  TChain * chain = new TChain("tree/SimpAnalysis");
	//ctau=1m M=3000GeV
	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M3000_new/161114_110710/0000/simptree_1.root");
  TFile *output = new TFile("XXTo4J_ctau1000_M3000_selection.root", "RECREATE");
	
	//ctau=2m M=50GeV
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M50/161110_132239/0000/simptree_1.root");
// 	TFile *output = new TFile("XXTo4J_ctau2000_M50_reallynocut.root", "RECREATE");	
	
	//ctau=2m M=3000GeV
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M3000/161110_132227/0000/simptree_1.root");
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M3000/161110_132227/0000/simptree_2.root");
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M3000/161110_132227/0000/simptree_3.root");
//   TFile *output = new TFile("XXTo4J_ctau2000_M3000_nocut.root", "RECREATE");
  
  
	int nJet, dijet_170, dijet_170_0p1, dijet_220_0p3, dijet_330_0p5, dijet_430;
	double pswgt_dijet_170;
  double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4];
	double CHEF_jet[4];
// 	E[4], M[4];
  
  chain->SetBranchAddress("nJet",&nJet);
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_phi", &jet_phi);
//   chain->SetBranchAddress("jet_e", &E);
//   chain->SetBranchAddress("jet_m", &M);
  chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
  chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
  chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &dijet_170_0p1);
  chain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &dijet_220_0p3);
  chain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &dijet_330_0p5);
  chain->SetBranchAddress("HLT_DiCentralPFJet430", &dijet_430);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
	
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
	TH1F *jet1_pt_withps = new TH1F("jet1_pt_withps", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt_withps = new TH1F("jet2_pt_withps", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH1F *jet1_chf_withps = new TH1F("jet1_chf_withps", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf_withps = new TH1F("jet2_chf_withps", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *trigger_dijet_170_0p1 = new TH1F("trigger_dijet_170_0p1", "Dijet_170_0p1", 2, -0.5, 1.5);
	TH1F *trigger_dijet_220_0p3 = new TH1F("trigger_dijet_220_0p3", "Dijet_220_0p3", 2, -0.5, 1.5);
	TH1F *trigger_dijet_330_0p5 = new TH1F("trigger_dijet_330_0p5", "Dijet_330_0p5", 2, -0.5, 1.5);
	TH1F *trigger_dijet_430 = new TH1F("trigger_dijet_430", "Dijet_430", 2, -0.5, 1.5);
	TH1F *njets = new TH1F("njets", "Number of jets", 2, -0.5, 10.5);
	
	int total = 0;
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[7] = {0, 400, 450, 500, 600, 800, 10000};
	double eta_bins[3] = {0, 1.5, 2.5};
	int passed_pt[6][11];
	int passed_eta[2][11];
	double eff_pt[6][11], red_pt[6][11], eff_eta[2][11], red_eta[2][11];
	double err_eff_pt[6][11], err_eff_eta[6][11], err_red_pt[2][11], err_red_eta[2][11];
	double zero[11] = {0,0,0,0,0,0,0,0,0,0,0};
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 6; i++){
			passed_pt[i][j] = 0;
		}
	}
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 2; i++){
			passed_eta[i][j] = 0;
		}
	}
  
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
		
		if (jet_phi[0] == 0 || jet_phi[1] == 0) continue;
		if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2){
			jet1_pt->Fill(jet_pt[0]);
			jet2_pt->Fill(jet_pt[1]);
			jet3_pt->Fill(jet_pt[2]);
			jet4_pt->Fill(jet_pt[3]);
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
			jet3_chf->Fill(CHEF_jet[2]);
			jet4_chf->Fill(CHEF_jet[3]);
			jet1_chf_withps->Fill(CHEF_jet[0], pswgt_dijet_170);
			jet2_chf_withps->Fill(CHEF_jet[1], pswgt_dijet_170);
			trigger_dijet_170_0p1->Fill(dijet_170_0p1);
			trigger_dijet_220_0p3->Fill(dijet_220_0p3);
			trigger_dijet_330_0p5->Fill(dijet_330_0p5);
			trigger_dijet_430->Fill(dijet_430);
			
			if(CHEF_jet[0] > CHEF_jet[1] && CHEF_jet[0] > 0.5){
				total++;
				for(int j = 0; j < 11; j++){
					for (int i = 0; i < 6; i++){
						if (jet_pt[1]>pt_bins[i] && jet_pt[1]<pt_bins[i+1] && CHEF_jet[1]<chf_cuts[j]) passed_pt[i][j]++;
					}
					for (int i = 0; i < 2; i++){
						if (jet_eta[1]>eta_bins[i] && jet_eta[1]<eta_bins[i+1] && CHEF_jet[1]<chf_cuts[j]) passed_eta[i][j]++;
					}
				}
			}
			else if(CHEF_jet[0] < CHEF_jet[1] && CHEF_jet[1] > 0.5){
				total++;
				for(int j = 0; j < 11; j++){
					for (int i = 0; i < 6; i++){
						if (jet_pt[0]>pt_bins[i] && jet_pt[0]<pt_bins[i+1] && CHEF_jet[0]<chf_cuts[j]) passed_pt[i][j]++;
					}
					for (int i = 0; i < 2; i++){
						if (jet_eta[0]>eta_bins[i] && jet_eta[0]<eta_bins[i+1] && CHEF_jet[0]<chf_cuts[j]) passed_eta[i][j]++;
					}
				}
			}
		}    
	}
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 6; i++){
			double tot = (double)total;
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
		for (int i = 0; i < 2; i++){
			double tot = (double)total;
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
	for (int i = 1; i < 6; i++){
		TGraphErrors *otherptbin = new TGraphErrors(11, chf_cuts, eff_pt[i], zero, err_eff_pt[i]);
		otherptbin->SetMarkerStyle(20);
		otherptbin->SetMarkerColor(i+1);
		otherptbin->Draw("P");
		leg1->AddEntry(otherptbin, " < pt < ", "lep");
	}
	leg1->Draw();
	c1->Print("ptbinned_eff.png");
	
	TCanvas *c2 = new TCanvas("rejection_pt", "rejection_pt");
	c2->cd();
	TGraphErrors *firstptbin2 = new TGraphErrors(11, chf_cuts, red_pt[0], zero, err_red_pt[0]);
	firstptbin2->SetTitle("CHF cut rejection factor");
	firstptbin2->GetXaxis()->SetTitle("CHF cut");
	firstptbin2->GetYaxis()->SetTitle("rejection factor");
	firstptbin2->SetMarkerStyle(20);
	firstptbin2->Draw("AP");
	for (int i = 1; i < 6; i++){
		TGraphErrors *otherptbin2 = new TGraphErrors(11, chf_cuts, red_pt[i], zero, err_red_pt[i]);
		otherptbin2->SetMarkerStyle(20);
		otherptbin2->SetMarkerColor(i+1);
		otherptbin2->Draw("P");
	}
	c2->Print("ptbinned_red.png");
	
	TCanvas *c3 = new TCanvas("efficiency_eta", "efficiency_eta");
	c3->cd();
	TGraphErrors *firstetabin = new TGraphErrors(11, chf_cuts, eff_eta[0], zero, err_eff_eta[0]);
	firstetabin->SetTitle("CHF cut efficiency");
	firstetabin->GetXaxis()->SetTitle("CHF cut");
	firstetabin->GetYaxis()->SetTitle("efficiency");
	firstetabin->SetMarkerStyle(20);
	firstetabin->Draw("AP");
	TGraphErrors *otheretabin = new TGraphErrors(11, chf_cuts, eff_eta[1], zero, err_eff_eta[1]);
	otheretabin->SetMarkerStyle(20);
	otheretabin->SetMarkerColor(2);
	otheretabin->Draw("P");
	c3->Print("etabinned_eff.png");
	
	TCanvas *c4 = new TCanvas("rejection_eta", "rejection_eta");
	c4->cd();
	TGraphErrors *firstetabin2 = new TGraphErrors(11, chf_cuts, red_eta[0], zero, err_red_eta[0]);
	firstetabin2->SetTitle("CHF cut rejection factor");
	firstetabin2->GetXaxis()->SetTitle("CHF cut");
	firstetabin2->GetYaxis()->SetTitle("rejection factor");
	firstetabin2->SetMarkerStyle(20);
	firstetabin2->Draw("AP");
	TGraphErrors *otheretabin2 = new TGraphErrors(11, chf_cuts, red_eta[1], zero, err_red_eta[1]);
	otheretabin2->SetMarkerStyle(20);
	otheretabin2->SetMarkerColor(2);
	otheretabin2->Draw("P");
	c4->Print("etabinned_red.png");
	
	output->Append(c1);
	output->Append(c2);
	output->Append(c3);
	output->Append(c4);
	
  output->Write();
  output->Close();
}