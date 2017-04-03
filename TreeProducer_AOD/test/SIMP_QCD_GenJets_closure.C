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

#include "SANtuple.h"

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"

void SIMP_QCD_GenJets_closure(){
  
    TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_300To500(chain0);
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_0.root");
  chain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_1.root");
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_500To700(chain1);
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_0.root");
  chain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_1.root");
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_700To1000(chain2);
  chain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17_ext.root");
  chain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1000To1500(chain3);
  chain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17_ext.root");
  chain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_1500To2000(chain4);
  chain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17_ext.root");
  chain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
// 	list_QCD_2000ToInf(chain5);
  chain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
  chain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17_ext.root");
//   chain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   chain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_300To500(SPVchain0);
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17.root");
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_0.root");
  SPVchain0->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT300To500_PUMoriond17_ext_1.root");
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_500To700(SPVchain1);
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17.root");
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_0.root");
  SPVchain1->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT500To700_PUMoriond17_ext_1.root");
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_700To1000(SPVchain2);
  SPVchain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17_ext.root");
  SPVchain2->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT700To1000_PUMoriond17.root");
  TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1000To1500(SPVchain3);
  SPVchain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17_ext.root");
  SPVchain3->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1000To1500_PUMoriond17.root");
  TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_1500To2000(SPVchain4);
  SPVchain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17_ext.root");
  SPVchain4->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT1500To2000_PUMoriond17.root");
  TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
// 	list_QCD_2000ToInf(SPVchain5);
  SPVchain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
  SPVchain5->Add("/pnfs/iihe/cms/store/user/gflouris/SIMPS/QCD_Skimmed/SIMPs_QCD_HT2000ToInf_PUMoriond17_ext.root");
//   SPVchain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   SPVchain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* SPVchains[6] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5};
	std::cout<<"TChains ready"<<std::endl;
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
	
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
	
	double ratio_1leg[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double ratio_2leg[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_1leg[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg[11] = {0,0,0,0,0,0,0,0,0,0,0};
  
	std::cout<<"Getting the efficiency histos...";
	TFile* efficiencies = new TFile("eff2D_QCD_GenJets_PUMoriond17_AOD_skim.root", "READ");
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
	
  TFile *output = new TFile("closure_test_GenJets_PUMoriond17_AOD_skim.root", "RECREATE");
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91};//Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
	double Nentries[6] = {36948578, 44029226, 29808140, 10360193, 7868538, 4047360};
	double lumi = 33.095*1000;
  
	for (int l = 0; l < 6; l++){
		TChain* chain = chains[l];
		TChain* SPVchain = SPVchains[l];		
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t entries = chain->GetEntries(); 
// 		Int_t Nentries = 1000000;
		double weight = QCD_xsec[l]*lumi/Nentries[l];
		std::cout<<"QCD bin "<<l<<": Processing "<<entries<<" entries with weight "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"events"<<std::endl;
			chain->GetEntry(entry);
			
			double deltajet_phi = corrjets.genjet_phi[0] - corrjets.genjet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			deltajet_phi = fabs(deltajet_phi);
			
			output->cd();
			
			if (corrjets.genjet_pt[0] > 250 && corrjets.genjet_pt[1] > 250 && fabs(corrjets.genjet_eta[0]) < 2.0 && fabs(corrjets.genjet_eta[1]) < 2.0 && deltajet_phi > 2){
        for(int j = 0; j < 11; j++){
          if (corrjets.genjet_efrac_ch[0]<chf_cuts[j] && corrjets.genjet_efrac_ch[1]<chf_cuts[j]){
            passed_MCtruth[j]+= weight;
            err_MCtruth[j] += pow(weight, 2);
          }
          double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.genjet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.genjet_pt[0]));
          double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.genjet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.genjet_pt[1]));
          double erreff1 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          double erreff2 = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
// 					std::cout<<chf_cuts[j]<<std::endl;
// 					std::cout<<jet_eta[0]<<"  "<<jet_pt[0]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[0]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[0])<<std::endl;
// 					std::cout<<jet_eta[1]<<"  "<<jet_pt[1]<<" bin "<<eff_histos[j]->GetXaxis()->FindBin(fabs(jet_eta[1]))<<", "<<eff_histos[j]->GetYaxis()->FindBin(jet_pt[1])<<std::endl;
// 					std::cout<<"eff1 = "<<eff1<<std::endl;
// 					std::cout<<"eff2 = "<<eff2<<std::endl;
          if (corrjets.genjet_efrac_ch[1]<chf_cuts[j]) passed_eff1[j]+=weight*eff1;
          if (corrjets.genjet_efrac_ch[0]<chf_cuts[j]) passed_eff2[j]+=weight*eff2;
          passed_effboth[j]+=eff1*eff2*weight;
          
          err_eff1[j] += weight*weight*erreff1*erreff1;
          err_eff2[j] += weight*weight*erreff2*erreff2;
          err_effboth[j] += (eff2*eff2*weight*weight*erreff1*erreff1) + (eff1*eff1*weight*weight*erreff2*erreff2);
        }
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
// 		ratio_1leg[j] = passed_eff[j]/passed_MCtruth[j];
// 		ratio_2leg[j] = passed_effboth[j]/passed_MCtruth[j];
// 		err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
// 		err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		ratio_1leg[j] = passed_MCtruth[j]/passed_eff[j];
		ratio_2leg[j] = passed_MCtruth[j]/passed_effboth[j];
		err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
	}
	TCanvas *c1 = new TCanvas("Closure test", "Closure test",1400,1000);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
	TGraphErrors *MC = new TGraphErrors(11, chf_cuts, passed_MCtruth, zero, err_MCtruth);
	MC->SetTitle("GenJet Closure test");
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->SetMarkerStyle(20);
	MC->Draw("AP");
	
	TGraphErrors *oneleg = new TGraphErrors(11, chf_cuts, passed_eff, zero, err_eff);
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P");
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
  pad2->SetLogy();
  pad2->Draw();
  pad2->cd();
	
	TGraphErrors *ratio_oneleg = new TGraphErrors(11, chf_cuts, ratio_1leg, zero, err_ratio_1leg);
	ratio_oneleg->SetMarkerStyle(20);
	ratio_oneleg->SetMarkerColor(2);
	ratio_oneleg->Draw("AP");
  ratio_oneleg->SetTitle("");
  ratio_oneleg->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio_oneleg->GetYaxis()->SetTitle("#frac{MCtruth}{1/2-leg}");
	
	TGraphErrors *ratio_twoleg = new TGraphErrors(11, chf_cuts, ratio_2leg, zero, err_ratio_2leg);
	ratio_twoleg->SetMarkerStyle(20);
	ratio_twoleg->SetMarkerColor(3);
	ratio_twoleg->Draw("P");
  
  TLine *line = new TLine(0,1,0.55,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
	line->Draw();
	
  c1->cd();
	output->Append(c1);
  output->Write();
  output->Close();
}