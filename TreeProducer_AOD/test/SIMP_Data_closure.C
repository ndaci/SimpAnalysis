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
#include <TLine.h>

#include "SANtuple.h"
// #include "json_parser.cpp"

#include "lists/list_JetHT_2016B.h"
#include "lists/list_JetHT_2016C.h"
#include "lists/list_JetHT_2016D.h"
#include "lists/list_JetHT_2016E.h"
#include "lists/list_JetHT_2016F.h"
#include "lists/list_JetHT_2016G.h"
#include "lists/list_JetHT_2016H2.h"
#include "lists/list_JetHT_2016H3.h"

void SIMP_Data_closure(){
	  
  TChain* chain = new TChain("treeCorr/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
	list_JetHT_2016H2(chain);
	list_JetHT_2016H3(chain);
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
	list_JetHT_2016B(SPVchain);
	list_JetHT_2016C(SPVchain);
	list_JetHT_2016D(SPVchain);
	list_JetHT_2016E(SPVchain);
	list_JetHT_2016F(SPVchain);
	list_JetHT_2016G(SPVchain);
	list_JetHT_2016H2(SPVchain);
	list_JetHT_2016H3(SPVchain);
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
	
	double chf_cuts[23] = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
  double njet_cuts[4] = {2, 3, 4, 5};
// 	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
// 	double eta_bins[17] = {-2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};

//   pt::ptree root;
//   pt::read_json("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt", root);
  
	double passed_eff[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff1[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff2[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_effboth[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_MCtruth[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_MCtruth[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff1[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff2[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_effboth[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double zero[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	double ratio_1leg[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double ratio_2leg[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_1leg[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_Data_Rereco_AOD_withNvertexCut_Njets2_finerEtaBinning.root", "READ");
// 	TFile* efficiencies = new TFile("eff3D_Data_ChFMin0p5_AOD.root", "READ");
// 	TFile* efficiencies = new TFile("../../TreeProducer_miniAOD/test/eff2D_Data_Rereco_full2016.root", "READ");
	TH2D* eff_histos[4][23];
  for (int i = 0; i < 4; i++){
    for(int j = 0; j < 23; j++){
      std::ostringstream strs;
      double dbl = chf_cuts[j];
      strs << dbl;
      std::ostringstream strs2;
      double dbl2 = njet_cuts[i];
      strs2 << dbl2;
      std::string cut = strs.str();
      std::string cut2 = strs2.str();
//       std::string title_eff = "eff_njets_"+cut2+"_"+cut;
      std::string title_eff = "eff_"+cut;
      eff_histos[i][j] = (TH2D*) efficiencies->Get(title_eff.c_str());
    }
  }
	std::cout<<"done"<<std::endl;
	
  TFile *output = new TFile("closure_test_Data_nvertexcut_njets2_eff2DFinerEtaBinning.root", "RECREATE");
	  
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t Nentries = corrjets.fChain->GetEntries(); 
//   Int_t Nentries = 3000000; 
  
  for(Int_t entry = 0; entry < Nentries; ++entry){
    if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
    SPVjets.GetEntry(entry);
    corrjets.GetEntry(entry);
    
//     bool isInJSON = Check_Run_Lumi(corrjets.nRun, corrjets.nLumi, root);
//     if (!isInJSON) continue;
    
    double njets = 0;
    for (int k = 0; k < 8; ++k){
      if (corrjets.jet_pt[k] > 30) njets++;
    }
    
    double deltajet_phi =  corrjets.jet_phi[0] -  corrjets.jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
    deltajet_phi = fabs(deltajet_phi);
      
    double photonpt = 0;
    photon_nr = 0;
    for (int i = 0; i < 4; ++i){
      if(corrjets.photon_pt[i]>photonpt){
        photon_nr = i;
        photonpt = corrjets.photon_pt[i];
      }
    }
    
    double deltaphi_jet1photon =  corrjets.jet_phi[0] -  corrjets.photon_phi[photon_nr];
    if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
    if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
    double deltaphi_jet2photon =  corrjets.jet_phi[1] -  corrjets.photon_phi[photon_nr];
    if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
    if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
    
    double deltaeta_jet1photon =  corrjets.jet_eta[0] -  corrjets.photon_eta[photon_nr];
    double deltaeta_jet2photon =  corrjets.jet_eta[1] -  corrjets.photon_eta[photon_nr];
    
    double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
    double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
    
    for (int i = 0; i < 8; i++){
      CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
      CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
    } 
    
    output->cd();
    
    if (corrjets.HLT_DiCentralPFJet170 == 1 && corrjets.jet_pt[0] > 250 && corrjets.jet_pt[1] > 250 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.vtx_N >= 2 /*&& fabs(corrjets.jet_eta[0])>1.5 && fabs(corrjets.jet_eta[1])>1.5 && fabs(corrjets.jet_eta[0])<2.0 && fabs(corrjets.jet_eta[1])<2.0*/ && njets == 2){
      for(int j = 0; j < 23; j++){
        if (CHEF_SPVjet[0]<chf_cuts[j] && CHEF_SPVjet[1]<chf_cuts[j] && CHEF_corrjet[0]<chf_cuts[j] && CHEF_corrjet[1]<chf_cuts[j]){
          if (chf_cuts[j] > 0.25) passed_MCtruth[j]+= corrjets.pswgt_dijet_170;
          if (chf_cuts[j] > 0.25) err_MCtruth[j] += pow(corrjets.pswgt_dijet_170, 2);
        }
        int i = 0;
//         if (corrjets.nJet == 2) i = 0;
//         else if (corrjets.nJet == 3) i = 1;
//         else if (corrjets.nJet == 4) i = 2;
//         else if (corrjets.nJet >= 5) i = 3;
//         else std::cout<<"wrong nb of jets!!"<<std::endl;
        double eff1 = eff_histos[i][j]->GetBinContent(eff_histos[i][j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[i][j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
        double eff2 = eff_histos[i][j]->GetBinContent(eff_histos[i][j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[i][j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
        double erreff1 = eff_histos[i][j]->GetBinError(eff_histos[i][j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[i][j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
        double erreff2 = eff_histos[i][j]->GetBinError(eff_histos[i][j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[i][j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
// 					std::cout<<chf_cuts[j]<<std::endl;
// 					std::cout<<jet_eta[0]<<"  "<<jet_pt[0]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[0])<<std::endl;
// 					std::cout<<jet_eta[1]<<"  "<<jet_pt[1]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[1])<<std::endl;
// 					std::cout<<"eff1 = "<<eff1<<std::endl;
// 					std::cout<<"eff2 = "<<eff2<<std::endl;
        if (CHEF_corrjet[1]<chf_cuts[j]) passed_eff1[j]+=corrjets.pswgt_dijet_170*eff1;
        if (CHEF_corrjet[0]<chf_cuts[j]) passed_eff2[j]+=corrjets.pswgt_dijet_170*eff2;
        passed_effboth[j]+=eff1*eff2*corrjets.pswgt_dijet_170;
        
        err_eff1[j] += corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff1*erreff1;
        err_eff2[j] += corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff2*erreff2;
        err_effboth[j] += (eff2*eff2*corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff1*erreff1) + (eff1*eff1*corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff2*erreff2);
      }
    }
	}
	for(int j = 0; j < 23; j++){
		err_MCtruth[j] = TMath::Sqrt(err_MCtruth[j]);
		err_effboth[j] = TMath::Sqrt(err_effboth[j]);
		err_eff1[j] = TMath::Sqrt(err_eff1[j]);
		err_eff2[j] = TMath::Sqrt(err_eff2[j]);
		err_eff[j] = (err_eff1[j]+err_eff2[j])/2.0 ;
		passed_eff[j] = (passed_eff1[j]+passed_eff2[j])/2;
		ratio_1leg[j] = passed_MCtruth[j]/passed_eff[j];
		ratio_2leg[j] = passed_MCtruth[j]/passed_effboth[j];
		if (ratio_1leg[j] != 0) err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		if (ratio_2leg[j] != 0) err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
	}
	TCanvas *c1 = new TCanvas("Closure test", "Closure test",700,500);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
	TGraphErrors *MC = new TGraphErrors(23, chf_cuts, passed_MCtruth, zero, err_MCtruth);
	MC->SetTitle("");
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->GetYaxis()->SetTitleSize(0.06);
	MC->GetYaxis()->SetLabelSize(0.06);
	MC->GetYaxis()->SetTitleOffset(0.8);
	MC->SetMarkerStyle(20);
	MC->Draw("AP");
	
	TGraphErrors *oneleg = new TGraphErrors(23, chf_cuts, passed_eff, zero, err_eff);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P");
// 	TGraphErrors *oneleg1 = new TGraphErrors(23, chf_cuts, passed_eff1, zero, zero);
// 	oneleg1->SetMarkerStyle(20);
// 	oneleg1->SetMarkerColor(2);
// 	oneleg1->Draw("P");
// 	TGraphErrors *oneleg2 = new TGraphErrors(23, chf_cuts, passed_eff2, zero, zero);
// 	oneleg2->SetMarkerStyle(20);
// 	oneleg2->SetMarkerColor(3);
// 	oneleg2->Draw("P");
	
	TGraphErrors *twoleg = new TGraphErrors(23, chf_cuts, passed_effboth, zero, err_effboth);
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
	twoleg->Draw("P");
	
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->SetLogy();
  pad2->Draw();
  pad2->cd();
	
	TGraphErrors *ratio_oneleg = new TGraphErrors(23, chf_cuts, ratio_1leg, zero, err_ratio_1leg);
	ratio_oneleg->SetMarkerStyle(20);
	ratio_oneleg->SetMarkerColor(2);
	ratio_oneleg->Draw("AP");
  ratio_oneleg->SetTitle("");
  ratio_oneleg->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio_oneleg->GetYaxis()->SetTitle("#frac{data}{1/2-leg}");
	ratio_oneleg->GetYaxis()->SetTitleSize(0.11);
	ratio_oneleg->GetYaxis()->SetLabelSize(0.1);
	ratio_oneleg->GetYaxis()->SetTitleOffset(0.4);
  ratio_oneleg->GetXaxis()->SetTitle("ChF cut");
	ratio_oneleg->GetXaxis()->SetTitleSize(0.12);
	ratio_oneleg->GetXaxis()->SetLabelSize(0.1);
  
	TGraphErrors *ratio_twoleg = new TGraphErrors(23, chf_cuts, ratio_2leg, zero, err_ratio_2leg);
	ratio_twoleg->SetMarkerStyle(20);
	ratio_twoleg->SetMarkerColor(3);
	ratio_twoleg->Draw("P");
  
  TLine *line = new TLine(0,1,1,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
	line->Draw();
	
  c1->cd();
	output->Append(c1);
  output->Write();
  output->Close();
}
