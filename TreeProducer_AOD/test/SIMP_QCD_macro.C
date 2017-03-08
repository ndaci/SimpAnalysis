#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
// #include <TRandom3.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

#include "SANtuple.h"

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD.h"

void SIMP_QCD_macro(){
  
  TChain* chain0 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_300To500(chain0);
  TChain* chain1 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_500To700(chain1);
  TChain* chain2 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_700To1000(chain2);
  TChain* chain3 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1000To1500(chain3);
  TChain* chain4 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_1500To2000(chain4);
  TChain* chain5 = new TChain("treeCorr/SimpAnalysis");
	list_QCD_2000ToInf(chain5);
//   chain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   chain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* chains[6] = {chain0, chain1, chain2, chain3, chain4, chain5};
  
  TChain* SPVchain0 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_300To500(SPVchain0);
  TChain* SPVchain1 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_500To700(SPVchain1);
  TChain* SPVchain2 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_700To1000(SPVchain2);
  TChain* SPVchain3 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1000To1500(SPVchain3);
  TChain* SPVchain4 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_1500To2000(SPVchain4);
  TChain* SPVchain5 = new TChain("treeSPV/SimpAnalysis");
	list_QCD_2000ToInf(SPVchain5);
//   SPVchain5->Add("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root");
//   SPVchain5->Add("QCD_PUMoriond17_AOD_test.root");
	TChain* SPVchains[6] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5};
	std::cout<<"TChains ready"<<std::endl;
  
  TFile *output = new TFile("plots_QCD_lumi33p095_PUMoriond17_AOD.root", "RECREATE");
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int counter;
  
	TH1F *njets = new TH1F("njets", "Number of jets", 11, -0.5, 10.5);
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 100, 0, 2000);
	TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 100, 0, 2000);
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet1_pt_SPV = new TH1F("jet1_pt_SPV", "Leading jet pt SPV collection", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
// 	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
// 	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_eta_SPV = new TH1F("jet1_eta_SPV", "Leading jet eta SPV collection", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
// 	TH1F *jetA_pt = new TH1F("jetA_pt", "ChF-leading jet pt", 100, 0, 2000);
// 	TH1F *jetB_pt = new TH1F("jetB_pt", "ChF-subleading jet pt", 100, 0, 2000);
// 	TH1F *jetA_eta = new TH1F("jetA_eta", "ChF-leading jet eta", 100, -3.14, 3.14);
// 	TH1F *jetB_eta = new TH1F("jetB_eta", "ChF-subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet1_chf_SPV = new TH1F("jet1_chf_SPV", "Charged energy fraction (leading jet) SPV collection", 100, 0, 1);
	TH1F *jet2_chf_SPV = new TH1F("jet2_chf_SPV", "Charged energy fraction (subleading jet) SPV collection", 100, 0, 1);
	TH1F *ChfOverSPVChF = new TH1F("ChfOverSPVChF", "ChF/ChF(SPV)", 1000, 0, 10);
	TH1F *chf0p5 = new TH1F("chf_Min0p5", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
	TH1F *chf0p2 = new TH1F("chf_Max0p2", "Charged energy fraction (other jet ChF < 0.2)", 100, 0, 1);
// 	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
// 	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsPT_0p5 = new TH2D("CHFvsPT_0p5", "Jet ChF vs. pT (tag jet CHF > 0.5)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsPT_0p2 = new TH2D("CHFvsPT_0p2", "Jet ChF vs. pT (tag jet CHF < 0.2)", 100, 0, 5000, 100, 0, 1);
	TH2D *ChfVsSPVChf = new TH2D("ChfVsSPVChf", "jet ChF vs. SPV jet ChF", 100, 0, 1, 100, 0, 1);
	TH2D *CHFvsCHF = new TH2D("CHFvsCHF", "Leading jet ChF vs. Subleading jet ChF", 100, 0, 1, 100, 0, 1);
	CHFvsCHF->GetXaxis()->SetTitle("jet2 ChF");
	CHFvsCHF->GetYaxis()->SetTitle("jet1 ChF");
// 	TH1D *CHF = new TH1D("CHF", "Jet ChF (both jets)", 100, 0, 1);
// 	TH1D *ChFOtherJet = new TH1D("ChFOtherJet", "Jet ChF (other jet ChF > 0.5)", 100, 0, 1);
	TH2D *NVTXvsCHF = new TH2D("NVTXvsCHF", "Number of vertices vs. ChF (both jets)", 100, 0, 1, 51, -0.5, 50.5);
	TH2D *NVTXvsCHF1 = new TH2D("NVTXvsCHF1", "Number of vertices vs. ChF (leading jet)", 100, 0, 1, 51, -0.5, 50.5);
	TH2D *NVTXvsCHF2 = new TH2D("NVTXvsCHF2", "Number of vertices vs. ChF (subleading jet)", 100, 0, 1, 51, -0.5, 50.5);
	TH2D *Chi2vsCHF = new TH2D("Chi2vsCHF", "#chi^{2} vs. ChF (both jets)", 100, 0, 1, 100, 0, 100);
	TH2D *Chi2vsCHF1 = new TH2D("Chi2vsCHF1", "#chi^{2} vs. ChF (leading jet)", 100, 0, 1, 100, 0, 100);
	TH2D *Chi2vsCHF2 = new TH2D("Chi2vsCHF2", "#chi^{2} vs. ChF (subleading jet)", 100, 0, 1, 100, 0, 100);
	TH2D *ndofvsCHF = new TH2D("ndofvsCHF", "ndof vs. ChF (both jets)", 100, 0, 1, 501, -0.5, 500.5);
	TH2D *ndofvsCHF1 = new TH2D("ndofvsCHF1", "ndof vs. ChF (leading jet)", 100, 0, 1, 501, -0.5, 500.5);
	TH2D *ndofvsCHF2 = new TH2D("ndofvsCHF2", "ndof vs. ChF (subleading jet)", 100, 0, 1, 501, -0.5, 500.5);
	TH2D *Chi2vsNdof = new TH2D("Chi2vsNdof", "#chi^{2} vs. ndof", 101, -0.5, 100.5, 100, 0, 500);
	TH2D *track1_reldzErrorvsdz = new TH2D("track1_reldzErrorvsdz", "dzError/dz vs. dz (hardest track)", 100, 0, 100, 100, 0, 100);
	TH1F *track1_nhits = new TH1F("track1_nhits", "Number of hits for hardest track", 11, -0.5, 10.5);
	TH1F *track1_npixhits = new TH1F("track1_npixhits", "Number of pixel hits for hardest track", 11, -0.5, 10.5);
	TH1F *track1_pt = new TH1F("track1_pt", "pt of hardest track", 100, 0, 500);
	TH1F *track1_ptError = new TH1F("track1_ptError", "ptError of hardest track", 100, 0, 1500);
	TH1F *track1_dzError = new TH1F("track1_dzError", "dzError of hardest track", 100, 0, 100);
	TH1F *track1_dz = new TH1F("track1_dz", "dz of hardest track", 100, 0, 100);
	TH1F *track1_dzRelError = new TH1F("track1_dzRelError", "dzError/dz of hardest track", 100, 0, 100);
// 	TH1F *track1_fromPV = new TH1F("track1_fromPV", "Hardest track from PV?", 2, -0.5, 1.5);
// 	TH2D *Pt1vsCHF1 = new TH2D("Pt1vsCHF1", "p_{T} vs. ChF (leading jet)", 100, 0, 1, 100, 0, 2000);
// 	TH2D *eta1vsCHF1 = new TH2D("eta1vsCHF1", "#eta vs. ChF (leading jet)", 100, 0, 1, 100, -3.14, 3.14);
// 	TH2D *Pt2vsCHF2 = new TH2D("Pt2vsCHF2", "p_{T} vs. ChF (subleading jet)", 100, 0, 1, 100, 0, 2000);
// 	TH2D *eta2vsCHF2 = new TH2D("eta2vsCHF2", "#eta vs. ChF (subleading jet)", 100, 0, 1, 100, -3.14, 3.14);
// 	TH2D *Pt2vsCHF1 = new TH2D("Pt2vsCHF1", "p_{T} vs. ChF (leading jet)", 100, 0, 1, 100, 0, 2000);
// 	TH2D *eta2vsCHF1 = new TH2D("eta2vsCHF1", "#eta vs. ChF (leading jet)", 100, 0, 1, 100, -3.14, 3.14);
// 	TH2D *Pt1vsCHF2 = new TH2D("Pt1vsCHF2", "p_{T} vs. ChF (subleading jet)", 100, 0, 1, 100, 0, 2000);
// 	TH2D *eta1vsCHF2 = new TH2D("eta1vsCHF2", "#eta vs. ChF (subleading jet)", 100, 0, 1, 100, -3.14, 3.14);
	TH2D *CHFvsGenCHF1 = new TH2D("CHFvsGenCHF1", "ChF vs. ChF GenJet (leading jet)", 100, 0, 1, 100, 0, 1);
	TH2D *CHFvsGenCHF2 = new TH2D("CHFvsGenCHF2", "ChF vs. ChF GenJet (subleading jet)", 100, 0, 1, 100, 0, 1);
	TH1F *relChF1 = new TH1F("relChF1", "(ChF - ChF GenJet)/ChF GenJet (leading jet)", 1000, -10, 10);
	TH1F *relChF2 = new TH1F("relChF2", "(ChF - ChF GenJet)/ChF GenJet (subleading jet)", 1000, -10, 10);
	TH1F *genjet1_chf = new TH1F("genjet1_chf", "Charged energy fraction (leading genjet)", 100, 0, 1);
	TH1F *genjet2_chf = new TH1F("genjet2_chf", "Charged energy fraction (subleading genjet)", 100, 0, 1);
	TH2D *eta1vsnpixhits = new TH2D("eta1vsnpixhits", "eta vs. # pixel hits (leading jet)", 11, -0.5, 10.5, 100, -3.14, 3.14);
	TH2D *eta2vsnpixhits = new TH2D("eta2vsnpixhits", "eta vs. # pixel hits (subleading jet)", 11, -0.5, 10.5, 100, -3.14, 3.14);
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //in pb//Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
// 	double lumi = 0.062*1000; // in 1000 fb^-1 = pb^-1
  double lumi = 33.095*1000; //RunB-H
  double Nentries[6] = {54552852, 0, 0, 0, 11839357, 6039005};
  
	for (int j = 0; j < 6; j++){
		TChain* chain = chains[j];		
		TChain* SPVchain = SPVchains[j];
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t entries = corrjets.fChain->GetEntries();
// 		Int_t entries = 3000000;
		std::cout<<"Processing "<<entries<<"entries"<<std::endl;
		
// 		double weight = QCD_xsec[j]*lumi/Nentries[j];
		double weight = QCD_xsec[j]*lumi/entries;
		std::cout<<"weight "<<j<<": "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
      SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
			
// 			Double_t random1 = r.Uniform();
			
			double deltajet_phi = corrjets.jet_phi[0] - corrjets.jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			deltajet_phi = fabs(deltajet_phi);
			
			double deltaphi_jet1photon = corrjets.jet_phi[0] - corrjets.photon_phi[0];
			if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
			if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
			double deltaphi_jet2photon = corrjets.jet_phi[1] - corrjets.photon_phi[0];
			if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
			if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
			
			double deltaeta_jet1photon = corrjets.jet_eta[0] - corrjets.photon_eta[0];
			double deltaeta_jet2photon = corrjets.jet_eta[1] - corrjets.photon_eta[0];
			
			double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
			double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
			
			for (int i = 0; i < 8; i++){
				CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
				CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (corrjets.jet_pt[0] > 250 && corrjets.jet_pt[1] > 250 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && corrjets.track_nPixHits[0] > 0 && (corrjets.photon_passLooseId[0] == 0 || (corrjets.photon_passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
				njets->Fill(corrjets.nJet, weight);
				nvtx->Fill(corrjets.vtx_N, weight);
				HT->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3], weight);
				HT_nowgt->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]);
				
				jet1_pt->Fill(corrjets.jet_pt[0], weight);
        jet1_pt_SPV->Fill(SPVjets.jet_pt[0], weight);
				jet2_pt->Fill(corrjets.jet_pt[1], weight);
// 				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], weight);
// 				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], weight);
				
				jet1_eta->Fill(corrjets.jet_eta[0], weight);
				jet1_eta_SPV->Fill(SPVjets.jet_eta[0], weight);
				jet2_eta->Fill(corrjets.jet_eta[1], weight);
				jet1_phi->Fill(corrjets.jet_phi[0], weight);
				jet2_phi->Fill(corrjets.jet_phi[1], weight);
				deltaphi->Fill(deltajet_phi, weight);
				
				jet1_chf->Fill(CHEF_corrjet[0], weight);
				jet2_chf->Fill(CHEF_corrjet[1], weight);
				if (CHEF_corrjet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_corrjet[0], weight);
				jet1_chf_SPV->Fill(CHEF_SPVjet[0], weight);
				jet2_chf_SPV->Fill(CHEF_SPVjet[1], weight);
				ChfOverSPVChF->Fill(CHEF_corrjet[0]/CHEF_SPVjet[0], weight);
				ChfOverSPVChF->Fill(CHEF_corrjet[1]/CHEF_SPVjet[1], weight);
        
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
        if (dRjet1 < 0.4){
          if (CHEF_corrjet[1]<0.05) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[0], weight);
          if (CHEF_corrjet[0]<0.05)ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[1], weight);
        }else if(dRjet2 < 0.4){
          if (CHEF_corrjet[1]<0.05) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[1], weight);
          if (CHEF_corrjet[0]<0.05) ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[0], weight);          
        }else counter++;
// 				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2], weight);
// 				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3], weight);
// 				CHF->Fill(CHEF_jet[0], weight);
// 				CHF->Fill(CHEF_jet[1], weight);
				CHFvsCHF->Fill(CHEF_corrjet[1], CHEF_corrjet[0], weight);				
				CHFvsPT->Fill(corrjets.jet_pt[0], CHEF_corrjet[0], weight);
				CHFvsPT->Fill(corrjets.jet_pt[1], CHEF_corrjet[1], weight);
				
				genjet1_chf->Fill(corrjets.genjet_efrac_ch[0], weight);
				genjet2_chf->Fill(corrjets.genjet_efrac_ch[1], weight);
				CHFvsGenCHF1->Fill(corrjets.genjet_efrac_ch[0], CHEF_corrjet[0], weight);
				CHFvsGenCHF2->Fill(corrjets.genjet_efrac_ch[1], CHEF_corrjet[1], weight);
				relChF1->Fill((CHEF_corrjet[0]-corrjets.genjet_efrac_ch[0])/corrjets.genjet_efrac_ch[0], weight);
				relChF2->Fill((CHEF_corrjet[1]-corrjets.genjet_efrac_ch[1])/corrjets.genjet_efrac_ch[1], weight);
				
				eta1vsnpixhits->Fill(corrjets.track_nPixHits[0], corrjets.jet_eta[0], weight);
				eta2vsnpixhits->Fill(corrjets.track_nPixHits[0], corrjets.jet_eta[1], weight);
				
				NVTXvsCHF->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
				NVTXvsCHF->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);
				NVTXvsCHF1->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
				NVTXvsCHF2->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);
				
				Chi2vsCHF->Fill(CHEF_corrjet[0], corrjets.track_normalizedChi2[0], weight);
				Chi2vsCHF->Fill(CHEF_corrjet[1], corrjets.track_normalizedChi2[0], weight);
				Chi2vsCHF1->Fill(CHEF_corrjet[0], corrjets.track_normalizedChi2[0], weight);
				Chi2vsCHF2->Fill(CHEF_corrjet[1], corrjets.track_normalizedChi2[0], weight);
				
				ndofvsCHF->Fill(CHEF_corrjet[0], corrjets.track_ndof[0], weight);
				ndofvsCHF->Fill(CHEF_corrjet[1], corrjets.track_ndof[0], weight);
				ndofvsCHF1->Fill(CHEF_corrjet[0], corrjets.track_ndof[0], weight);
				ndofvsCHF2->Fill(CHEF_corrjet[1], corrjets.track_ndof[0], weight);
				
// 				Pt1vsCHF1->Fill(CHEF_jet[0], jet_pt[0], weight);
// 				Pt2vsCHF2->Fill(CHEF_jet[1], jet_pt[1], weight);
// 				eta1vsCHF1->Fill(CHEF_jet[0], jet_eta[0], weight);
// 				eta2vsCHF2->Fill(CHEF_jet[1], jet_eta[1], weight);				
// 				Pt2vsCHF1->Fill(CHEF_jet[0], jet_pt[1], weight);
// 				Pt1vsCHF2->Fill(CHEF_jet[1], jet_pt[0], weight);
// 				eta2vsCHF1->Fill(CHEF_jet[0], jet_eta[1], weight);
// 				eta1vsCHF2->Fill(CHEF_jet[1], jet_eta[0], weight);
				
				track1_nhits->Fill(corrjets.track_nhits[0], weight);
				track1_npixhits->Fill(corrjets.track_nPixHits[0], weight);
				track1_pt->Fill(corrjets.track_pt[0], weight);
				track1_ptError->Fill(corrjets.track_ptError[0], weight);
				track1_dzError->Fill(corrjets.track_dzError[0], weight);
				track1_dz->Fill(corrjets.track_dz[0], weight);
				track1_dzRelError->Fill(corrjets.track_dzError[0]/corrjets.track_dz[0], weight);
				track1_reldzErrorvsdz->Fill(corrjets.track_dz[0], corrjets.track_dzError[0]/corrjets.track_dz[0], weight);
// 				track1_fromPV->Fill(track_fromPV[0], weight);
				
// 				if (ndof[0] == 0) std::cout<<vtx_x[0]<<"\t"<<vtx_y[0]<<"\t"<<vtx_z[0]<<"\t"<<vtx_d0[0]<<std::endl;
				
				Chi2vsNdof->Fill(corrjets.track_ndof[0], corrjets.track_normalizedChi2[0], weight);
				
// 				if(CHEF_jet[0] > CHEF_jet[1]){
// 					jetA_pt->Fill(jet_pt[0], weight);
// 					jetB_pt->Fill(jet_pt[1], weight);
// 					jetA_eta->Fill(jet_eta[0], weight);
// 					jetB_eta->Fill(jet_eta[1], weight);
// 				}
// 				else if(CHEF_jet[0] < CHEF_jet[1]){
// 					jetA_pt->Fill(jet_pt[1], weight);
// 					jetB_pt->Fill(jet_pt[0], weight);
// 					jetA_eta->Fill(jet_eta[1], weight);
// 					jetB_eta->Fill(jet_eta[0], weight);
// 				}
// 				if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] > 0.5) chf0p5->Fill(CHEF_jet[1], weight);
// 				if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] > 0.5) chf0p5->Fill(CHEF_jet[0], weight);
// 				if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] < 0.2) chf0p2->Fill(CHEF_jet[1], weight);
// 				if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] < 0.2) chf0p2->Fill(CHEF_jet[0], weight);
// 				
// 				if (CHEF_jet[1] > 0.5) CHFvsPT_0p5->Fill(jet_pt[0], CHEF_jet[0], weight);
// 				if (CHEF_jet[0] > 0.5) CHFvsPT_0p5->Fill(jet_pt[1], CHEF_jet[1], weight);
// 				if (CHEF_jet[1] < 0.2) CHFvsPT_0p2->Fill(jet_pt[0], CHEF_jet[0], weight);
// 				if (CHEF_jet[0] < 0.2) CHFvsPT_0p2->Fill(jet_pt[1], CHEF_jet[1], weight);
				
// 				if(random1 <= 0.5 && CHEF_jet[0] > 0.5) ChFOtherJet->Fill(CHEF_jet[1], weight);
// 				else if (random1 > 0.5 && CHEF_jet[1] > 0.5) ChFOtherJet->Fill(CHEF_jet[0], weight);
			}    
		}
	}
  output->Write();
  output->Close();
  std::cout<<"counter: "<<counter<<std::endl;
}
