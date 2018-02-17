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

#include "lists/list_QCD_1000To1500_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD_filters.h"
#include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_filters.h"
// #include "lists/list_QCD_1000To1500_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1500To2000_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_2000ToInf_PUMoriond17_AOD_newtrigger_2.h"
// #include "lists/list_QCD_1000To1500_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_1500To2000_ext_PUMoriond17_AOD_newtrigger.h"
// #include "lists/list_QCD_2000ToInf_ext_PUMoriond17_AOD_newtrigger.h"

void SIMP_QCD_macro_newtrigger(){
  
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
  
  TFile *output = new TFile("plots_QCD_PUMoriond17_AOD_filters_noDeltaphi.root", "RECREATE");
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
  int counter;
  
	TH1F *nJets = new TH1F("nJets", "Number of jets", 11, -0.5, 10.5);
	TH1F *nJets_NEMFMin0p9 = new TH1F("nJets_NEMFMin0p9", "Number of jets (NEMF > 0.9)", 11, -0.5, 10.5);
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 150, 1000, 4000);
	TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 150, 1000, 4000);
  TH1F *METOverHT = new TH1F("METOverHT", "MET / HT(4 jets)", 100, 0, 1);	
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
	TH1F *deltaphi_ChFMax0p5 = new TH1F("DeltaPhi_ChFMax0p5", "DeltaPhi (ChF < 0.5)", 100, 0, 3.14);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet_chf_pt550To600 = new TH1F("jet_chf_pt550To600", "Charged energy fraction (pt 550-600 GeV)", 100, 0, 1);
	TH1F *jet_chf_pt600To700 = new TH1F("jet_chf_pt600To700", "Charged energy fraction (pt 600-700 GeV)", 100, 0, 1);
	TH1F *jet_chf_pt700To800 = new TH1F("jet_chf_pt700To800", "Charged energy fraction (pt 700-800 GeV)", 100, 0, 1);
	TH1F *jet_chf_pt800To950 = new TH1F("jet_chf_pt800To950", "Charged energy fraction (pt 800-950 GeV)", 100, 0, 1);
	TH1F *jet_chf_pt950To1200 = new TH1F("jet_chf_pt950To1200", "Charged energy fraction (pt 950-1200 GeV)", 100, 0, 1);
	TH1F *jet_chf_pt1200 = new TH1F("jet_chf_pt1200", "Charged energy fraction (pt>1200 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt550To600 = new TH1F("genjet_chf_pt550To600", "GenJet Charged energy fraction (pt 550-600 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt600To700 = new TH1F("jgenet_chf_pt600To700", "GenJet Charged energy fraction (pt 600-700 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt700To800 = new TH1F("genjet_chf_pt700To800", "GenJet Charged energy fraction (pt 700-800 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt800To950 = new TH1F("genjetchf_pt800To950", "GenJet Charged energy fraction (pt 800-950 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt950To1200 = new TH1F("genjet_chf_pt950To1200", "GenJet Charged energy fraction (pt 950-1200 GeV)", 100, 0, 1);
	TH1F *genjet_chf_pt1200 = new TH1F("genjet_chf_pt1200", "GenJet Charged energy fraction (pt>1200 GeV)", 100, 0, 1);
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
	double QCD_xsec[6] = {/*346400, 32010, 6842,*/ 1203, 120.1, 25.40, 1203, 120.1, 25.40}; //PUMoriond17
// 	double lumi = 0.062*1000; // in 1000 fb^-1 = pb^-1
//   double lumi = 33.095*1000; //RunB-H
  double lumi = 16.146*1000; //RunG-H
  
	for (int j = 0; j < 3; j++){
		TChain* chain = chains[j];		
		TChain* SPVchain = SPVchains[j];
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t entries = corrjets.fChain->GetEntries();
// 		Int_t entries = 3000000;
		std::cout<<"Processing "<<entries<<"entries"<<std::endl;
		
		double weight = QCD_xsec[j]*lumi/entries;
		std::cout<<"weight "<<j<<": "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < entries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
      SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
      
      double njets = 0;
      for (int k = 0; k < 8; ++k){
        if (corrjets.jet_pt[k] > 30) njets++;
      }
			
			double deltajet_phi = corrjets.jet_phi[0] - corrjets.jet_phi[1];
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
			
			double deltaphi_jet1photon = corrjets.jet_phi[0] - corrjets.photon_phi[photon_nr];
			if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
			if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
			double deltaphi_jet2photon = corrjets.jet_phi[1] - corrjets.photon_phi[photon_nr];
			if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
			if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
			
			double deltaeta_jet1photon = corrjets.jet_eta[0] - corrjets.photon_eta[photon_nr];
			double deltaeta_jet2photon = corrjets.jet_eta[1] - corrjets.photon_eta[photon_nr];
			
			double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
			double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
			
			for (int i = 0; i < 8; i++){
				CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
				CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
			}
			
      bool pass_conv_1 = true;
      bool pass_conv_2 = true;
      if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
      if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
			
			output->cd();
			
//       std::cout<<"pt1: "<< corrjets.jet_pt[0]<<" pt2: "<<corrjets.jet_pt[1]<<" eta1: "<<corrjets.jet_eta[0]<<" eta2: "<<corrjets.jet_eta[1]<<" deltaphi: "<<deltajet_phi<<std::endl;
      
			if (corrjets.HLT_PFJet450 == 1 && corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 /*&& deltajet_phi > 2 */&& (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9) && corrjets.vtx_N >= 2 && njets == 2 && pass_conv_1 && pass_conv_2 && corrjets.Flag_HBHENoiseFilter == 1 /*&& corrjets.Flag_HBHENoiseIsoFilter == 1 */&& corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){				
				jet1_eta->Fill(corrjets.jet_eta[0], weight);
				jet1_eta_SPV->Fill(SPVjets.jet_eta[0], weight);
				jet2_eta->Fill(corrjets.jet_eta[1], weight);
        
        nJets->Fill(njets, weight);
//         if (corrjets.jet_efrac_photon[0] >= 0.9 && corrjets.jet_efrac_photon[1] >= 0.9) nJets_NEMFMin0p9->Fill(njets, weight);
        nvtx->Fill(corrjets.vtx_N, weight);
        HT->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3], weight);
        HT_nowgt->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]);
        METOverHT->Fill(corrjets.MET/(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]), weight);
        
        jet1_pt->Fill(corrjets.jet_pt[0], weight);
        jet1_pt_SPV->Fill(SPVjets.jet_pt[0], weight);
        jet2_pt->Fill(corrjets.jet_pt[1], weight);
// 				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], weight);
// 				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], weight);
        jet1_phi->Fill(corrjets.jet_phi[0], weight);
        jet2_phi->Fill(corrjets.jet_phi[1], weight);
        deltaphi->Fill(deltajet_phi, weight);
        if (CHEF_corrjet[0] < 0.5 && CHEF_corrjet[1] < 0.5 && CHEF_SPVjet[0] < 0.5 && CHEF_SPVjet[1] < 0.5) deltaphi_ChFMax0p5->Fill(deltajet_phi, weight);
        
        jet1_chf->Fill(CHEF_corrjet[0], weight);
        if (corrjets.jet_pt[0] > 550 && corrjets.jet_pt[0] < 600) jet_chf_pt550To600->Fill(CHEF_corrjet[0], weight);
        else if (corrjets.jet_pt[0] > 600 && corrjets.jet_pt[0] < 700) jet_chf_pt600To700->Fill(CHEF_corrjet[0], weight);
        else if (corrjets.jet_pt[0] > 700 && corrjets.jet_pt[0] < 800) jet_chf_pt700To800->Fill(CHEF_corrjet[0], weight);
        else if (corrjets.jet_pt[0] > 800 && corrjets.jet_pt[0] < 950) jet_chf_pt800To950->Fill(CHEF_corrjet[0], weight);
        else if (corrjets.jet_pt[0] > 950 && corrjets.jet_pt[0] < 1200) jet_chf_pt950To1200->Fill(CHEF_corrjet[0], weight);
        else if (corrjets.jet_pt[0] > 1200) jet_chf_pt1200->Fill(CHEF_corrjet[0], weight);
        jet2_chf->Fill(CHEF_corrjet[1], weight);
        if (corrjets.jet_pt[1] > 550 && corrjets.jet_pt[1] < 600) jet_chf_pt550To600->Fill(CHEF_corrjet[1], weight);
        else if (corrjets.jet_pt[1] > 600 && corrjets.jet_pt[1] < 700) jet_chf_pt600To700->Fill(CHEF_corrjet[1], weight);
        else if (corrjets.jet_pt[1] > 700 && corrjets.jet_pt[1] < 800) jet_chf_pt700To800->Fill(CHEF_corrjet[1], weight);
        else if (corrjets.jet_pt[1] > 800 && corrjets.jet_pt[1] < 950) jet_chf_pt800To950->Fill(CHEF_corrjet[1], weight);
        else if (corrjets.jet_pt[1] > 950 && corrjets.jet_pt[1] < 1200) jet_chf_pt950To1200->Fill(CHEF_corrjet[1], weight);
        else if (corrjets.jet_pt[1] > 1200) jet_chf_pt1200->Fill(CHEF_corrjet[1], weight);
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
//           ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[0], weight);
          if (CHEF_corrjet[0]>CHEF_corrjet[1]) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[0], weight);
          if (CHEF_corrjet[1]>CHEF_corrjet[0]) ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[1], weight);
//           std::cout<<"corrjet1: "<<CHEF_corrjet[0]<<" SPVjet1: "<<CHEF_SPVjet[0]<<" corrjet2: "<<CHEF_corrjet[1]<<" SPVjet2: "<<CHEF_SPVjet[1]<<std::endl;
        }else if(dRjet2 < 0.4){
//           ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[1], weight);
          if (CHEF_corrjet[0]>CHEF_corrjet[1]) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[1], weight);
          if (CHEF_corrjet[1]>CHEF_corrjet[0]) ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[0], weight);   
//           std::cout<<"corrjet1: "<<CHEF_corrjet[0]<<" SPVjet1: "<<CHEF_SPVjet[0]<<" corrjet2: "<<CHEF_corrjet[1]<<" SPVjet2: "<<CHEF_SPVjet[1]<<std::endl;       
        }else counter++;
// 				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2], weight);
// 				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3], weight);
// 				CHF->Fill(CHEF_jet[0], weight);
// 				CHF->Fill(CHEF_jet[1], weight);
        CHFvsCHF->Fill(CHEF_corrjet[1], CHEF_corrjet[0], weight);				
        CHFvsPT->Fill(corrjets.jet_pt[0], CHEF_corrjet[0], weight);
        CHFvsPT->Fill(corrjets.jet_pt[1], CHEF_corrjet[1], weight);
        
        genjet1_chf->Fill(corrjets.genjet_efrac_ch[0], weight);
        
        if (corrjets.genjet_pt[0] > 550 && corrjets.genjet_pt[0] < 600) genjet_chf_pt550To600->Fill(corrjets.genjet_efrac_ch[0], weight);
        else if (corrjets.genjet_pt[0] > 600 && corrjets.genjet_pt[0] < 700) genjet_chf_pt600To700->Fill(corrjets.genjet_efrac_ch[0], weight);
        else if (corrjets.genjet_pt[0] > 700 && corrjets.genjet_pt[0] < 800) genjet_chf_pt700To800->Fill(corrjets.genjet_efrac_ch[0], weight);
        else if (corrjets.genjet_pt[0] > 800 && corrjets.genjet_pt[0] < 950) genjet_chf_pt800To950->Fill(corrjets.genjet_efrac_ch[0], weight);
        else if (corrjets.genjet_pt[0] > 950 && corrjets.genjet_pt[0] < 1200) genjet_chf_pt950To1200->Fill(corrjets.genjet_efrac_ch[0], weight);
        else if (corrjets.genjet_pt[0] > 1200) genjet_chf_pt1200->Fill(corrjets.genjet_efrac_ch[0], weight);
        
        genjet2_chf->Fill(corrjets.genjet_efrac_ch[1], weight);
        
        if (corrjets.genjet_pt[1] > 550 && corrjets.genjet_pt[1] < 600) genjet_chf_pt550To600->Fill(corrjets.genjet_efrac_ch[1], weight);
        else if (corrjets.genjet_pt[1] > 600 && corrjets.genjet_pt[1] < 700) genjet_chf_pt600To700->Fill(corrjets.genjet_efrac_ch[1], weight);
        else if (corrjets.genjet_pt[1] > 700 && corrjets.genjet_pt[1] < 800) genjet_chf_pt700To800->Fill(corrjets.genjet_efrac_ch[1], weight);
        else if (corrjets.genjet_pt[1] > 800 && corrjets.genjet_pt[1] < 950) genjet_chf_pt800To950->Fill(corrjets.genjet_efrac_ch[1], weight);
        else if (corrjets.genjet_pt[1] > 950 && corrjets.genjet_pt[1] < 1200) genjet_chf_pt950To1200->Fill(corrjets.genjet_efrac_ch[1], weight);
        else if (corrjets.genjet_pt[1] > 1200) genjet_chf_pt1200->Fill(corrjets.genjet_efrac_ch[1], weight);
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
