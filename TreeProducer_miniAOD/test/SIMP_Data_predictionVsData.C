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

void SIMP_Data_predictionVsData(){
  
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
	
	double chf_cuts[21] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
  
  TH1F * prediction = new TH1F("prediction", "prediction", 20, chf_cuts);
  TH1F * data = new TH1F("data", "data", 20, chf_cuts);
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_Data_Rereco_exclusiveChF.root", "READ");
	TH2D* eff_histos[20];
	for(int j = 0; j < 20; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		double dbl2 = chf_cuts[j+1];
		strs << dbl << "To" << dbl2;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
	std::cout<<"done"<<std::endl;
	
  TFile *output = new TFile("predictionVsData_test.root", "RECREATE");
  
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
      
        if (CHEF_jet[0]>0.25) data->Fill(CHEF_jet[0], pswgt_dijet_170);
        else data->Fill(CHEF_jet[0], 0);
        if (CHEF_jet[1]>0.25) data->Fill(CHEF_jet[1], pswgt_dijet_170);
        else data->Fill(CHEF_jet[1], 0);
        
			for(int j = 0; j < 20; j++){        
				double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
				double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
				double erreff1 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[0]));
				double erreff2 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(jet_pt[1]));
        
        prediction->Fill(chf_cuts[j]+0.025, pswgt_dijet_170*eff1);
        prediction->Fill(chf_cuts[j]+0.025, pswgt_dijet_170*eff2);
			}
		}
	}
	TCanvas *c1 = new TCanvas("predictionVsData", "predictionVsData");
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
	
  data->SetTitle("");	
  data->SetStats(0);
	data->GetXaxis()->SetTitle("ChF");
	data->GetXaxis()->SetLabelSize(0.1);
	data->GetXaxis()->SetTitleSize(0.15);
	data->GetYaxis()->SetTitle("# jets");
	data->GetYaxis()->SetLabelSize(0.05);
	data->GetYaxis()->SetTitleSize(0.05);
  data->SetMarkerStyle(20);
  data->SetMarkerColor(1);
	data->DrawCopy("EP");
  prediction->SetLineColor(kTeal-5);
	prediction->Draw("C hist same");
	
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  
	data->GetYaxis()->SetLabelSize(0.1);
	data->GetYaxis()->SetTitleSize(0.1);
	data->GetYaxis()->SetTitleOffset(0.3);
	data->GetYaxis()->SetTitle("#frac{data}{prediction}");
  data->GetYaxis()->SetRangeUser(0.82, 1.18);
  data->Divide(prediction);
  data->Draw("EP");
	
// 	TGraphErrors *ratio_plot = new TGraphErrors(11, chf_cuts, ratio, zero, err_ratio);
// 	ratio_plot->SetMarkerStyle(20);
// 	ratio_plot->SetMarkerColor(2);
// 	ratio_plot->Draw("AP");
//   ratio_plot->SetTitle("");
//   ratio_plot->GetYaxis()->SetRangeUser(0.8, 1.2);
//   ratio_plot->GetYaxis()->SetTitle("#frac{1-leg}{2-leg}");
  
  TLine *line = new TLine(0,1,1,1);
  line->SetLineColor(kTeal-5);
  line->SetLineWidth(2);
	line->Draw();
	
	c1->cd();
	output->Append(c1);
  output->Write();
  output->Close();
}