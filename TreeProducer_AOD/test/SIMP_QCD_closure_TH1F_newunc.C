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

#include "lists/list_QCD_1000To1500_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_filters.h"
// #include "lists/list_QCD_1000To1500_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1500To2000_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1000To1500_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_1500To2000_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_2000ToInf_ext_PUMoriond17_AOD_newtrigger.h"

void SIMP_QCD_closure_TH1F_newunc(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain0);
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain1);
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain2);
//   TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1000To1500_ext(chain3);
//   TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1500To2000_ext(chain4);
//   TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_2000ToInf_ext(chain5);
	TChain* chains[3] = {chain0, chain1, chain2};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain0);
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain1);
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain2);
//   TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1000To1500_ext(SPVchain3);
//   TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1500To2000_ext(SPVchain4);
//   TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_2000ToInf_ext(SPVchain5);
	TChain* SPVchains[3] = {SPVchain0, SPVchain1, SPVchain2};
	std::cout<<"TChains ready"<<std::endl; 
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
  
  TH1D* histo_1leg = new TH1D("histo_1leg", "systematics from closure", 100, 0, 1.0);
  TH1D* histo_2leg = new TH1D("histo_2leg", "systematics from closure", 100, 0, 1.0);
	
	double chf_cuts[24] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9};
	double chf_cuts_tgraph[23] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9};
  double chf_bins[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.085, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925};
	
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
	double err_ratio_1leg_low[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg_low[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_1leg_high[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg_high[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  int N[3][23] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
  
  double weights[3];
  
  TH1D* MC_truth[3] = {new TH1D("MCtruth_1000To1500", "MCtruth_1000To1500", 23, chf_bins), new TH1D("MCtruth_1500To2000", "MCtruth_1500To2000", 23, chf_bins), new TH1D("MCtruth_2000ToInf", "MCtruth_2000ToInf", 23, chf_bins)};
  for (int i = 0; i < 3; i++){
    MC_truth[i]->Sumw2(kFALSE);
    MC_truth[i]->SetBinErrorOption(TH1::kPoisson);
  }
  
	std::cout<<"Getting the efficiency histos...";
// 	TFile* efficiencies = new TFile("eff2D_QCD_PUMoriond17_AOD_NvertexCut_Njets2_conversionCut.root", "READ");
	TFile* efficiencies = new TFile("eff2D_QCD_filters.root", "READ");
	TH2D* eff_histos[23];
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j+1];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
	std::cout<<"done"<<std::endl;
	
//   TFile *output = new TFile("closure_test_MCtruthSPVcut_PUMoriond17_AOD_conversionCut_lumiGH_correct.root", "RECREATE");
  TFile *output = new TFile("closure_test_QCD_filters_1ratioplot.root", "RECREATE");
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};//Spring16
  double QCD_xsec[6] = {/*346400, 32010, 6842,*/ 1203, 120.1, 25.40, 1203, 120.1, 25.40}; //PUMoriond17
// 	double lumi = 33.095*1000;
	double lumi = 16.146*1000;
  
	for (int l = 0; l < 3; l++){
		TChain* chain = chains[l];
		TChain* SPVchain = SPVchains[l];		
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
//     Int_t entries = 1000;
		Int_t entries = corrjets.fChain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/entries;
    weights[l] = QCD_xsec[l]*lumi/entries;
		std::cout<<"QCD bin "<<l<<": Processing "<<entries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
      SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
      
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
			
			for (int i = 0; i < 8; ++i){
				CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
				CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
			} 
        
      double deltaphi_jet1 = corrjets.jet_phi[0] - SPVjets.jet_phi[0];
      if(deltaphi_jet1 > TMath::Pi()) deltaphi_jet1 -= 2*TMath::Pi();
      if(deltaphi_jet1 < -TMath::Pi()) deltaphi_jet1 += 2*TMath::Pi();
      double deltaphi_jet2 = corrjets.jet_phi[1] - SPVjets.jet_phi[0];
      if(deltaphi_jet2 > TMath::Pi()) deltaphi_jet2 -= 2*TMath::Pi();
      if(deltaphi_jet2 < -TMath::Pi()) deltaphi_jet2 += 2*TMath::Pi();
      
      double deltaeta_jet1 = corrjets.jet_eta[0] - SPVjets.jet_eta[0];
      double deltaeta_jet2 = corrjets.jet_eta[1] - SPVjets.jet_eta[0];
      
      double dRjet1 = TMath::Sqrt(deltaphi_jet1*deltaphi_jet1 + deltaeta_jet1*deltaeta_jet1);
      double dRjet2 = TMath::Sqrt(deltaphi_jet2*deltaphi_jet2 + deltaeta_jet2*deltaeta_jet2); 
			
			output->cd();
      
			if (corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && ((corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9)) && corrjets.vtx_N >= 2 && njets == 2 && corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){
        bool pass_conv_1 = true;
        bool pass_conv_2 = true;
        if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
        if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
        if(pass_conv_1 && pass_conv_2){
        
          for(int j = 0; j < 23; j++){
            if (CHEF_SPVjet[0]<chf_cuts[j+1] && CHEF_SPVjet[1]<chf_cuts[j+1] && CHEF_corrjet[0]<chf_cuts[j+1] && CHEF_corrjet[1]<chf_cuts[j+1]){
              passed_MCtruth[j]+= weight;
              err_MCtruth[j] += pow(weight, 2);
              MC_truth[l]->Fill(chf_cuts[j+1]);
              N[l][j]++;
            }
            double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
            double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
            double erreff1 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
            double erreff2 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));

            if (CHEF_corrjet[1]<chf_cuts[j+1]) passed_eff1[j]+=weight*eff1;
            if (CHEF_corrjet[0]<chf_cuts[j+1]) passed_eff2[j]+=weight*eff2;
            passed_effboth[j]+=eff1*eff2*weight;
            
            err_eff1[j] += weight*weight*erreff1*erreff1;
            err_eff2[j] += weight*weight*erreff2*erreff2;
            err_effboth[j] += (eff2*eff2*weight*weight*erreff1*erreff1) + (eff1*eff1*weight*weight*erreff2*erreff2);
          }
        }
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
    std::cout<<chf_cuts_tgraph[j]<<": "<<passed_eff[j]<<" +- "<<err_eff[j]<<std::endl;
		ratio_1leg[j] = passed_MCtruth[j]/passed_eff[j];
		ratio_2leg[j] = passed_MCtruth[j]/passed_effboth[j];
	}
// 	for(int j = 0; j < 23; j++){
//     if(err_MCtruth[j] > fabs(passed_MCtruth[j]-passed_eff[j])) histo_1leg->Fill(chf_cuts[j+1], err_MCtruth[j]/passed_MCtruth[j]);
//     else histo_1leg->Fill(chf_cuts[j+1],fabs(passed_MCtruth[j]-passed_eff[j])/passed_eff[j]);
//     if(err_MCtruth[j] > fabs(passed_MCtruth[j]-passed_effboth[j])) histo_2leg->Fill(chf_cuts[j+1], err_MCtruth[j]/passed_MCtruth[j]);
//     else histo_2leg->Fill(chf_cuts[j+1], fabs(passed_MCtruth[j]-passed_effboth[j])/passed_effboth[j]);
//   }
  
	TCanvas *c1 = new TCanvas("Closure test", "Closure test",1400,1000);
//   TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
  
  TH1D *MC_TH1 = (TH1D*) MC_truth[0]->Clone();
  MC_TH1->Scale(weights[0]);
  MC_TH1->Add(MC_truth[1], weights[1]);
  MC_TH1->Add(MC_truth[2], weights[2]);
  
 	TGraphAsymmErrors *MC = new TGraphAsymmErrors(MC_TH1);

  std::cout<<"weight0="<<weights[0]<<" weight1="<<weights[1]<<" weight2="<<weights[2]<<std::endl;
	for(int j = 0; j < 23; j++){
    int N0 = N[0][j];
    int N1 = N[1][j];
    int N2 = N[2][j];
    double err0up = MC_truth[0]->GetBinErrorUp(MC_truth[0]->GetXaxis()->FindBin(chf_cuts[j+1]));
    double err1up = MC_truth[1]->GetBinErrorUp(MC_truth[1]->GetXaxis()->FindBin(chf_cuts[j+1]));
    double err2up = MC_truth[2]->GetBinErrorUp(MC_truth[2]->GetXaxis()->FindBin(chf_cuts[j+1]));
    double err0low = MC_truth[0]->GetBinErrorLow(MC_truth[0]->GetXaxis()->FindBin(chf_cuts[j+1]));
    double err1low = MC_truth[1]->GetBinErrorLow(MC_truth[1]->GetXaxis()->FindBin(chf_cuts[j+1]));
    double err2low = MC_truth[2]->GetBinErrorLow(MC_truth[2]->GetXaxis()->FindBin(chf_cuts[j+1]));
    std::cout<<"j="<<j<<" ChF cut="<<chf_cuts[j+1]<<std::endl;
    std::cout<<"N0="<<N0<<" N1="<<N1<<" N2="<<N2<<std::endl;
    std::cout<<"err0up="<<err0up<<" err1up="<<err1up<<" err2up="<<err2up<<std::endl;
    std::cout<<"err0low="<<err0low<<" err1low="<<err1low<<" err2low="<<err2low<<std::endl;
    if (N0 > 10){
      err0up = TMath::Sqrt(N0);
      err0low = TMath::Sqrt(N0);
    }
    if (N1 > 10){
      err1up = TMath::Sqrt(N1);
      err1low = TMath::Sqrt(N1);
    }
    if (N2 > 10){
      err2up = TMath::Sqrt(N2);
      err2low = TMath::Sqrt(N2);
    }
    std::cout<<"err0up="<<err0up<<" err1up="<<err1up<<" err2up="<<err2up<<std::endl;
    std::cout<<"err0low="<<err0low<<" err1low="<<err1low<<" err2low="<<err2low<<std::endl;
    MC->SetPointEXhigh(j, 0);
    MC->SetPointEXlow(j, 0);
    MC->SetPointEYhigh(j, TMath::Sqrt(pow(weights[0]*err0up,2) + pow(weights[1]*err1up,2) + pow(weights[2]*err2up,2)));
    MC->SetPointEYlow(j, TMath::Sqrt(pow(weights[0]*err0low,2) + pow(weights[1]*err1low,2) + pow(weights[2]*err2low,2)));
  }
	MC->SetTitle("");
	MC->SetMarkerStyle(20);
  MC->SetMarkerSize(1.2);
	MC->GetYaxis()->SetTitle("# events");
	MC->GetYaxis()->SetTitleSize(0.05);
	MC->GetYaxis()->SetLabelSize(0.05);
	MC->Draw("AP E0");
	
	TGraphErrors *oneleg = new TGraphErrors(23, chf_cuts_tgraph, passed_eff, zero, err_eff);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
  oneleg->SetMarkerSize(1.2);
	oneleg->Draw("P");
  
	TGraphErrors *twoleg = new TGraphErrors(23, chf_cuts_tgraph, passed_effboth, zero, err_effboth);
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
  twoleg->SetMarkerSize(1.2);
	twoleg->Draw("P");
  
  Double_t x[23], y[23];
	for(int j = 0; j < 23; j++){
    MC->GetPoint(j, x[j], y[j]);
    double MC_passed = y[j];
    double MC_err_high = MC->GetErrorYhigh(j);
    double MC_err_low = MC->GetErrorYlow(j);
    double MC_err = TMath::Max(MC_err_high, MC_err_low);
    if (chf_cuts[j+1] == 0.1) MC->SetPoint(j, 0.1, y[j]);
    if (MC_passed != 0){
//       if(MC_err > fabs(MC_passed - passed_eff[j])) histo_1leg->Fill(chf_cuts[j+1], MC_err/MC_passed);
//       else histo_1leg->Fill(chf_cuts[j+1],fabs(MC_passed-passed_eff[j])/passed_eff[j]);
//       if(MC_err > fabs(MC_passed-passed_effboth[j])) histo_2leg->Fill(chf_cuts[j+1], MC_err/MC_passed);
//       else histo_2leg->Fill(chf_cuts[j+1], fabs(MC_passed-passed_effboth[j])/passed_effboth[j]);
//       err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(MC_err*MC_err/MC_passed/MC_passed));
//       err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(MC_err*MC_err/MC_passed/MC_passed));
      
      double err_1leg = TMath::Max(MC_err/MC_passed,fabs(MC_passed-passed_eff[j])/passed_eff[j]);
      double err_2leg = TMath::Max(MC_err/MC_passed,fabs(MC_passed-passed_effboth[j])/passed_effboth[j]);
      histo_1leg->Fill(chf_cuts[j+1], err_1leg);
      histo_2leg->Fill(chf_cuts[j+1], err_2leg);
      
      err_ratio_1leg_high[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(MC_err_high*MC_err_high/MC_passed/MC_passed));
      err_ratio_1leg_low[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(MC_err_low*MC_err_low/MC_passed/MC_passed));
      err_ratio_2leg_low[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(MC_err_low*MC_err_low/MC_passed/MC_passed));
      err_ratio_2leg_high[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(MC_err_high*MC_err_high/MC_passed/MC_passed));
    }
  }
	
// 	c1->cd();
// 	TPad *pad2 = new TPad("pad2","pad2",0,0.23,1,0.4);
//   pad2->SetTopMargin(0);
//   pad2->SetBottomMargin(0);
//   pad2->Draw();
//   pad2->cd();
	
  TGraphAsymmErrors *ratio_oneleg = new TGraphAsymmErrors(23, chf_cuts_tgraph, ratio_1leg, zero, zero, err_ratio_1leg_low, err_ratio_1leg_high);
	ratio_oneleg->SetMarkerStyle(20);
	ratio_oneleg->SetMarkerColor(2);
  ratio_oneleg->SetTitle("");
  ratio_oneleg->GetYaxis()->SetRangeUser(0.81, 1.47);
  ratio_oneleg->GetYaxis()->SetTitle("#frac{MCtruth}{1/2-leg}");
  ratio_oneleg->GetYaxis()->CenterTitle(kTRUE);
  ratio_oneleg->GetYaxis()->SetTitleSize(0.16);
  ratio_oneleg->GetYaxis()->SetTitleOffset(0.3);
  ratio_oneleg->GetYaxis()->SetLabelSize(0.16);
  ratio_oneleg->SetMarkerSize(1.2);
// 	ratio_oneleg->Draw("AP");
	
	TGraphAsymmErrors *ratio_twoleg = new TGraphAsymmErrors(23, chf_cuts_tgraph, ratio_2leg, zero, zero, err_ratio_2leg_low, err_ratio_2leg_high);
	ratio_twoleg->SetMarkerStyle(20);
	ratio_twoleg->SetMarkerColor(3);
  ratio_twoleg->SetMarkerSize(1.2);
// 	ratio_twoleg->Draw("P");
  
  TLine *line = new TLine(0,1,0.99,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
// 	line->Draw();
  
	c1->cd();
// 	TPad *pad3 = new TPad("pad3","pad3",0,0,1,0.24);
	TPad *pad3 = new TPad("pad3","pad3",0,0,1,0.3);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.35);
  pad3->SetLogy(1);
  pad3->Draw();
  pad3->cd();
  
  
  ratio_oneleg->GetXaxis()->SetTitle("ChF cut");
  ratio_oneleg->GetXaxis()->SetTitleSize(0.12);
  ratio_oneleg->GetXaxis()->SetLabelSize(0.12);
	ratio_oneleg->Draw("AP");
	ratio_twoleg->Draw("P");
	line->Draw();
	
  c1->cd();
	output->Append(c1);
  histo_1leg->Write();
  histo_2leg->Write();
  oneleg->Write();
  twoleg->Write();
  output->Write();
  output->Close();
}
