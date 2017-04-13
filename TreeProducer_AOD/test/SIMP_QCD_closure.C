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

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD.h"

void SIMP_QCD_closure(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_300To500(chain0);
//   chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_500To700(chain1);
//   chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_700To1000(chain2);
//   chain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain3);
//   chain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain4);
//   chain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain5);
//   chain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_300To500(SPVchain0);
//   SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_500To700(SPVchain1);
//   SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_700To1000(SPVchain2);
//   SPVchain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain3);
//   SPVchain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain4);
//   SPVchain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain5);
//   SPVchain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
	TChain* SPVchains[6] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5};
	std::cout<<"TChains ready"<<std::endl;
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
  
  TH1D* histo_1leg = new TH1D("histo_1leg", "systematics from closure", 100, 0, 0.5);
  TH1D* histo_2leg = new TH1D("histo_2leg", "systematics from closure", 100, 0, 0.5);
	
	double chf_cuts[12] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.5};
	
	double passed_eff[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff1[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff2[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_effboth[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_MCtruth[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_MCtruth[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff1[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff2[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_effboth[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double zero[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	
	double ratio_1leg[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double ratio_2leg[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_1leg[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_QCD_PUMoriond17_AOD_NvertexCut_Njets2.root", "READ");
	TH2D* eff_histos[12];
	for(int j = 0; j < 12; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
	std::cout<<"done"<<std::endl;
	
  TFile *output = new TFile("closure_test_MCtruthSPVcut_PUMoriond17_AOD_pt350.root", "RECREATE");
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};//Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
  double Nentries[6] = {17035891, 18852895, 15629253, 4825904, 3970819, 1991645};
// 		Int_t Nentries[6] = {20000000, 2000000, 400000, 80000, 4000, 4000}; 
	double lumi = 33.095*1000;
  
	for (int l = 0; l < 6; l++){
		TChain* chain = chains[l];
		TChain* SPVchain = SPVchains[l];		
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
//     Int_t entries = 10000;
		Int_t entries = corrjets.fChain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/entries;
		std::cout<<"QCD bin "<<l<<": Processing "<<entries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
//       std::cout<<"processed "<<entry<<" events"<<std::endl;
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
      
			if (corrjets.jet_pt[0] > 350 && corrjets.jet_pt[1] > 350 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 /*&& corrjets.track_nPixHits[0] > 0*/ && (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.vtx_N >= 2 && njets == 2){
        
//         if (dRjet1 < 0.4 && CHEF_SPVjet[0] > 5*CHEF_corrjet[0] && CHEF_SPVjet[1] > 5*CHEF_corrjet[1]){
//           std::cout<<corrjets.nEvent<<" "<<corrjets.nLumi<<" "<<CHEF_corrjet[0]<<" "<<CHEF_SPVjet[0]<<" "<<CHEF_corrjet[1]<<" "<<CHEF_SPVjet[1]<<" "<<weight<<std::endl;
//           continue;			
//         }
//         if (dRjet2 < 0.4 && CHEF_SPVjet[0] > 5*CHEF_corrjet[1] && CHEF_SPVjet[1] > 5*CHEF_corrjet[0]){
//           std::cout<<corrjets.nEvent<<" "<<corrjets.nLumi<<" "<<CHEF_corrjet[0]<<" "<<CHEF_SPVjet[1]<<" "<<CHEF_corrjet[1]<<" "<<CHEF_SPVjet[0]<<" "<<weight<<std::endl;
//           continue;
//         }
      
        for(int j = 0; j < 12; j++){
          if (CHEF_SPVjet[0]<chf_cuts[j] && CHEF_SPVjet[1]<chf_cuts[j] && CHEF_corrjet[0]<chf_cuts[j] && CHEF_corrjet[1]<chf_cuts[j]){
            passed_MCtruth[j]+= weight;
            err_MCtruth[j] += pow(weight, 2);
//             if (j == 10) std::cout<<corrjets.nEvent<<" "<<corrjets.nLumi<<" "<<CHEF_corrjet[0]<<" "<<CHEF_SPVjet[0]<<" "<<CHEF_corrjet[1]<<" "<<CHEF_SPVjet[1]<<" "<<weight<<std::endl;
          }
          double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          double erreff1 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          double erreff2 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
// 					std::cout<<chf_cuts[j]<<std::endl;
// 					std::cout<<jet_eta[0]<<"  "<<jet_pt[0]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[0])<<std::endl;
// 					std::cout<<jet_eta[1]<<"  "<<jet_pt[1]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[1])<<std::endl;
// 					std::cout<<"eff1 = "<<eff1<<std::endl;
// 					std::cout<<"eff2 = "<<eff2<<std::endl;
          if (/*CHEF_SPVjet[1]<chf_cuts[j] && */CHEF_corrjet[1]<chf_cuts[j]) passed_eff1[j]+=weight*eff1;
          if (/*CHEF_SPVjet[0]<chf_cuts[j] && */CHEF_corrjet[0]<chf_cuts[j]) passed_eff2[j]+=weight*eff2;
          passed_effboth[j]+=eff1*eff2*weight;
          
          err_eff1[j] += weight*weight*erreff1*erreff1;
          err_eff2[j] += weight*weight*erreff2*erreff2;
          err_effboth[j] += (eff2*eff2*weight*weight*erreff1*erreff1) + (eff1*eff1*weight*weight*erreff2*erreff2);
        }
			}
		}
	}
	for(int j = 0; j < 12; j++){
		err_MCtruth[j] = TMath::Sqrt(err_MCtruth[j]);
		err_effboth[j] = TMath::Sqrt(err_effboth[j]);
		err_eff1[j] = TMath::Sqrt(err_eff1[j]);
		err_eff2[j] = TMath::Sqrt(err_eff2[j]);
		err_eff[j] = (err_eff1[j]+err_eff2[j])/2.0 ;
		passed_eff[j] = (passed_eff1[j]+passed_eff2[j])/2;
// 		ratio_1leg[j] = passed_eff[j]/passed_MCtruth[j];
// 		ratio_2leg[j] = passed_effboth[j]/passed_MCtruth[j];
// 		err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
// 		err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		ratio_1leg[j] = passed_MCtruth[j]/passed_eff[j];
		ratio_2leg[j] = passed_MCtruth[j]/passed_effboth[j];
		err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
	}
	for(int j = 0; j < 12; j++){
    histo_1leg->Fill(chf_cuts[j], ((passed_MCtruth[j]-passed_eff[j])/passed_MCtruth[j])+1);
    histo_2leg->Fill(chf_cuts[j], ((passed_MCtruth[j]-passed_effboth[j])/passed_MCtruth[j])+1);
  }
  
	TCanvas *c1 = new TCanvas("Closure test", "Closure test",1400,1000);
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
	TGraphErrors *MC = new TGraphErrors(12, chf_cuts, passed_MCtruth, zero, err_MCtruth);
	MC->SetTitle("");
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->SetMarkerStyle(20);
	MC->Draw("AP");
	
	TGraphErrors *oneleg = new TGraphErrors(12, chf_cuts, passed_eff, zero, err_eff);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P");
// 	TGraphErrors *oneleg1 = new TGraphErrors(12, chf_cuts, passed_eff1, zero, zero);
// 	oneleg1->SetMarkerStyle(20);
// 	oneleg1->SetMarkerColor(2);
// 	oneleg1->Draw("P");
// 	TGraphErrors *oneleg2 = new TGraphErrors(12, chf_cuts, passed_eff2, zero, zero);
// 	oneleg2->SetMarkerStyle(20);
// 	oneleg2->SetMarkerColor(3);
// 	oneleg2->Draw("P");
	
	TGraphErrors *twoleg = new TGraphErrors(12, chf_cuts, passed_effboth, zero, err_effboth);
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
	twoleg->Draw("P");
	
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.24,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->Draw();
  pad2->cd();
	
	TGraphErrors *ratio_oneleg = new TGraphErrors(12, chf_cuts, ratio_1leg, zero, err_ratio_1leg);
	ratio_oneleg->SetMarkerStyle(20);
	ratio_oneleg->SetMarkerColor(2);
  ratio_oneleg->SetTitle("");
  ratio_oneleg->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio_oneleg->GetYaxis()->SetTitle("#frac{MCtruth}{1/2-leg}");
  ratio_oneleg->GetYaxis()->SetTitleSize(0.06);
  ratio_oneleg->GetYaxis()->SetTitleOffset(0.4);
  ratio_oneleg->GetYaxis()->SetLabelSize(0.06);
	ratio_oneleg->Draw("AP");
	
	TGraphErrors *ratio_twoleg = new TGraphErrors(12, chf_cuts, ratio_2leg, zero, err_ratio_2leg);
	ratio_twoleg->SetMarkerStyle(20);
	ratio_twoleg->SetMarkerColor(3);
	ratio_twoleg->Draw("P");
  
  TLine *line = new TLine(0,1,0.55,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
	line->Draw();
  
	c1->cd();
	TPad *pad3 = new TPad("pad3","pad3",0,0,1,0.24);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.35);
  pad3->SetLogy(1);
  pad3->Draw();
  pad3->cd();
  
  ratio_oneleg->GetXaxis()->SetTitle("ChF cut");
  ratio_oneleg->GetXaxis()->SetTitleSize(0.06);
  ratio_oneleg->GetXaxis()->SetLabelSize(0.06);
	ratio_oneleg->Draw("AP");
	ratio_twoleg->Draw("P");
	line->Draw();
	
  c1->cd();
	output->Append(c1);
  histo_1leg->Write();
  histo_2leg->Write();
  output->Write();
  output->Close();
}
