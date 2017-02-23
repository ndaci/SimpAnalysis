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
#include <TGraphErrors.h>
#include <TCanvas.h>

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

void SIMP_Data_closure(){
  
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
	std::cout<<"TChains ready"<<std::endl;
  
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8];
	double track_pt[10], track_ptError[10];
	int npixhits[10];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	int dijet_170;
	double pswgt_dijet_170;
	
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	
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
	
	double ratio[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio[11] = {0,0,0,0,0,0,0,0,0,0,0};
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_photonVeto_Data_Rereco_runH.root", "READ");
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
	
  TFile *output = new TFile("closure_test_photonVeto_Data_Rereco_noMCTruth_ratio_runH_test.root", "RECREATE");
  
	chain->SetBranchAddress("jet_pt", &jet_pt);
	chain->SetBranchAddress("jet_eta", &jet_eta);
	chain->SetBranchAddress("jet_phi", &jet_phi);
	chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
	chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
	chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
	chain->SetBranchAddress("track_nPixHits", &npixhits);
	chain->SetBranchAddress("photon_pt", &photon_pt);
	chain->SetBranchAddress("photon_eta", &photon_eta);
	chain->SetBranchAddress("photon_phi", &photon_phi);
	chain->SetBranchAddress("photon_passLooseId",&passLooseId);
	chain->SetBranchAddress("photon_passMediumId",&passMediumId);
	chain->SetBranchAddress("photon_passTightId",&passTightId);
// 		chain->SetBranchAddress("track_pt", &track_pt);
// 		chain->SetBranchAddress("track_ptError", &track_ptError);
	
	Int_t Nentries = chain->GetEntries(); 
// 		Int_t Nentries = 1000000;
	
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
		
		if (dijet_170 == 1 && jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2 && npixhits[0] > 0 && (passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
			for(int j = 0; j < 11; j++){
				if (CHEF_jet[0]<chf_cuts[j] && CHEF_jet[1]<chf_cuts[j]){
					passed_MCtruth[j]+= pswgt_dijet_170;
					err_MCtruth[j] += pow(pswgt_dijet_170, 2);
				}
				double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
				double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
				double erreff1 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
				double erreff2 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
// 					std::cout<<chf_cuts[j]<<std::endl;
// 					std::cout<<jet_eta[0]<<"  "<<jet_pt[0]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[0])<<std::endl;
// 					std::cout<<jet_eta[1]<<"  "<<jet_pt[1]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[1])<<std::endl;
// 					std::cout<<"eff1 = "<<eff1<<std::endl;
// 					std::cout<<"eff2 = "<<eff2<<std::endl;
				if (CHEF_jet[1]<chf_cuts[j]) passed_eff1[j]+=eff1*pswgt_dijet_170;
				if (CHEF_jet[0]<chf_cuts[j]) passed_eff2[j]+=eff2*pswgt_dijet_170;
				passed_effboth[j]+=eff1*eff2*pswgt_dijet_170;
// 				if (CHEF_jet[1]<chf_cuts[j]) passed_eff1[j]+=pswgt_dijet_170*eff1;
// 				if (CHEF_jet[0]<chf_cuts[j]) passed_eff2[j]+=pswgt_dijet_170*eff2;
// 				passed_effboth[j]+=eff1*eff2*pswgt_dijet_170;
				
				err_eff1[j] += pswgt_dijet_170*pswgt_dijet_170*erreff1*erreff1;
				err_eff2[j] += pswgt_dijet_170*pswgt_dijet_170*erreff2*erreff2;
				err_effboth[j] += (eff2*eff2*pswgt_dijet_170*pswgt_dijet_170*erreff1*erreff1) + (eff1*eff1*pswgt_dijet_170*pswgt_dijet_170*erreff2*erreff2);
			}
		}
	}
	for(int j = 0; j < 11; j++){
		err_MCtruth[j] = TMath::Sqrt(err_MCtruth[j]);
		err_effboth[j] = TMath::Sqrt(err_effboth[j]);
		err_eff1[j] = TMath::Sqrt(err_eff1[j]);
		err_eff2[j] = TMath::Sqrt(err_eff2[j]);
		err_eff[j] = (err_eff1[j]+err_eff2[j])/2.0 ;
		passed_eff[j] = (passed_eff1[j]+passed_eff2[j])/2;
		ratio[j] = passed_eff[j]/passed_effboth[j];
		err_ratio[j] = ratio[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j]));
		std::cout<<chf_cuts[j]<<" "<<passed_eff[j]<<" "<<passed_effboth[j]<<std::endl;
	}
	TCanvas *c1 = new TCanvas("Closure test", "Closure test");
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
// 	TGraphErrors *MC = new TGraphErrors(11, chf_cuts, passed_MCtruth, zero, err_MCtruth);
// 	MC->SetTitle("Closure test");
// 	MC->GetXaxis()->SetTitle("CHF cut");
// 	MC->GetYaxis()->SetTitle("# events");
// 	MC->SetMarkerStyle(20);
// 	MC->Draw("AP");
	
	TGraphErrors *oneleg = new TGraphErrors(11, chf_cuts, passed_eff, zero, err_eff);
	oneleg->SetTitle("Closure test");
	oneleg->GetXaxis()->SetTitle("CHF cut");
	oneleg->GetYaxis()->SetTitle("# events");
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
// 	oneleg->Draw("P");
	oneleg->Draw("AP");
// 	TGraphErrors *oneleg1 = new TGraphErrors(11, chf_cuts, passed_eff1, zero, zero);
// 	oneleg1->SetMarkerStyle(20);
// 	oneleg1->SetMarkerColor(2);
// 	oneleg1->Draw("P");
// 	TGraphErrors *oneleg2 = new TGraphErrors(11, chf_cuts, passed_eff2, zero, zero);
// 	oneleg2->SetMarkerStyle(20);
// 	oneleg2->SetMarkerColor(3);
// 	oneleg2->Draw("P");
	
	TGraphErrors *twoleg = new TGraphErrors(11, chf_cuts, passed_effboth, zero, err_effboth);
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
	twoleg->Draw("P");
	
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
	
	TGraphErrors *ratio_plot = new TGraphErrors(11, chf_cuts, ratio, zero, err_ratio);
	ratio_plot->SetMarkerStyle(20);
	ratio_plot->SetMarkerColor(2);
	ratio_plot->Draw("AP");
  ratio_plot->SetTitle("");
  ratio_plot->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio_plot->GetYaxis()->SetTitle("#frac{1-leg}{2-leg}");
  
  TLine *line = new TLine(0,1,0.55,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
	line->Draw();
	
	c1->cd();
	output->Append(c1);
  output->Write();
  output->Close();
}