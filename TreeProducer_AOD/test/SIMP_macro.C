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

void SIMP_macro(int j){
	
	int mass[7] = {1, 10, 100, 200, 400, 700, 1000};
  TFile* outputfiles[7];
  
//   TChain * chain0 = new TChain("treeCorr/SimpAnalysis");
//   chain0->Add("SIMPs_PUMoriond17_AOD_M1_1_filters.root");
//   chain0->Add("SIMPs_PUMoriond17_AOD_M1_2_filters.root");
//   TChain * chain1 = new TChain("treeCorr/SimpAnalysis");
//   chain1->Add("SIMPs_PUMoriond17_AOD_M10_1_filters.root");
//   chain1->Add("SIMPs_PUMoriond17_AOD_M10_2_filters.root");
//   chain1->Add("SIMPs_PUMoriond17_AOD_M10_3_filters.root");
//   TChain * chain2 = new TChain("treeCorr/SimpAnalysis");
//   chain2->Add("SIMPs_PUMoriond17_AOD_M100_filters.root");
//   TChain * chain3 = new TChain("treeCorr/SimpAnalysis");
//   chain3->Add("SIMPs_PUMoriond17_AOD_M200_filters.root");
//   TChain * chain4 = new TChain("treeCorr/SimpAnalysis");
//   chain4->Add("SIMPs_PUMoriond17_AOD_M400_filters.root");
//   TChain * chain5 = new TChain("treeCorr/SimpAnalysis");
//   chain5->Add("SIMPs_PUMoriond17_AOD_M700_filters.root");
//   TChain * chain6 = new TChain("treeCorr/SimpAnalysis");
//   chain6->Add("SIMPs_PUMoriond17_AOD_M1000_filters.root");
//   
//   TChain * SPVchain0 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain0->Add("SIMPs_PUMoriond17_AOD_M1_1_filters.root");
//   SPVchain0->Add("SIMPs_PUMoriond17_AOD_M1_2_filters.root");
//   TChain * SPVchain1 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_1_filters.root");
//   SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_2_filters.root");
//   SPVchain1->Add("SIMPs_PUMoriond17_AOD_M10_3_filters.root");
//   TChain * SPVchain2 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain2->Add("SIMPs_PUMoriond17_AOD_M100_filters.root");
//   TChain * SPVchain3 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain3->Add("SIMPs_PUMoriond17_AOD_M200_filters.root");
//   TChain * SPVchain4 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain4->Add("SIMPs_PUMoriond17_AOD_M400_filters.root");
//   TChain * SPVchain5 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain5->Add("SIMPs_PUMoriond17_AOD_M700_filters.root");
//   TChain * SPVchain6 = new TChain("treeSPV/SimpAnalysis");
//   SPVchain6->Add("SIMPs_PUMoriond17_AOD_M1000_filters.root");
  
  TChain * chain0 = new TChain("treeCorr/SimpAnalysis");
  chain0->Add("SIMPs_Neutrons_M1_AOD_1.root");
  chain0->Add("SIMPs_Neutrons_M1_AOD_2.root");
  chain0->Add("SIMPs_Neutrons_M1_AOD_3.root");
  TChain * chain1 = new TChain("treeCorr/SimpAnalysis");
  chain1->Add("SIMPs_Neutrons_M10_AOD_1.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_2.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_3.root");
  chain1->Add("SIMPs_Neutrons_M10_AOD_4.root");
  TChain * chain2 = new TChain("treeCorr/SimpAnalysis");
  chain2->Add("SIMPs_Neutrons_M100_AOD_1.root");
  chain2->Add("SIMPs_Neutrons_M100_AOD_2.root");
  chain2->Add("SIMPs_Neutrons_M100_AOD_3.root");
  TChain * chain3 = new TChain("treeCorr/SimpAnalysis");
  chain3->Add("SIMPs_Neutrons_M200_AOD.root");
  TChain * chain4 = new TChain("treeCorr/SimpAnalysis");
  chain4->Add("SIMPs_Neutrons_M400_AOD.root");
  TChain * chain5 = new TChain("treeCorr/SimpAnalysis");
  chain5->Add("SIMPs_Neutrons_M700_AOD.root");
  TChain * chain6 = new TChain("treeCorr/SimpAnalysis");
  chain6->Add("SIMPs_Neutrons_M1000_AOD.root");
  
  TChain * SPVchain0 = new TChain("treeSPV/SimpAnalysis");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_1.root");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_2.root");
  SPVchain0->Add("SIMPs_Neutrons_M1_AOD_3.root");
  TChain * SPVchain1 = new TChain("treeSPV/SimpAnalysis");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_1.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_2.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_3.root");
  SPVchain1->Add("SIMPs_Neutrons_M10_AOD_4.root");
  TChain * SPVchain2 = new TChain("treeSPV/SimpAnalysis");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_1.root");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_2.root");
  SPVchain2->Add("SIMPs_Neutrons_M100_AOD_3.root");
  TChain * SPVchain3 = new TChain("treeSPV/SimpAnalysis");
  SPVchain3->Add("SIMPs_Neutrons_M200_AOD.root");
  TChain * SPVchain4 = new TChain("treeSPV/SimpAnalysis");
  SPVchain4->Add("SIMPs_Neutrons_M400_AOD.root");
  TChain * SPVchain5 = new TChain("treeSPV/SimpAnalysis");
  SPVchain5->Add("SIMPs_Neutrons_M700_AOD_new.root");
  TChain * SPVchain6 = new TChain("treeSPV/SimpAnalysis");
  SPVchain6->Add("SIMPs_Neutrons_M1000_AOD.root");
	
	TChain* chains[7] = {chain0, chain1, chain2, chain3, chain4, chain5, chain6};
	TChain* SPVchains[7] = {SPVchain0, SPVchain1, SPVchain2, SPVchain3, SPVchain4, SPVchain5, SPVchain6};
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];	
	
	double SIMP_xsec[7] = {4.461, 4.402, 2.553, 0.7903, 0.07434, 0.004846, 0.0005714}; //1/pb
//   double lumi = 33.095*1000; //RunB-H pb
  double lumi = 16.146*1000; //RunG-H pb
  
// 	for (int j = 0; j < 7; j++){
		std::ostringstream strs;
		double dbl = mass[j];
		strs << dbl;
		std::string m = strs.str();
    outputfiles[j] = new TFile(("plots_SIMP_M"+m+"_neutrons_noDeltaphi.root").c_str(), "RECREATE");
//   }
  
// 	for (int j = 0; j < 7; j++){

    TFile *output = outputfiles[j];
    TH1F *njets = new TH1F("njets", "Number of jets", 11, -0.5, 10.5);
    TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
    TH1F *vtx_ntracks = new TH1F("vtx_ntracks", "# tracks PV", 101, -0.5, 100.5);
    TH1F *HT = new TH1F("HT", "HT (4 jets)", 150, 1000, 4000);
    TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 150, 1000, 4000);
    TH1F *METOverHT = new TH1F("METOverHT", "MET / HT(4 jets)", 100, 0, 1);	
    TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
    TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
  // 	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
  // 	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
    TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
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
    TH1F *jet1_nhf = new TH1F("jet1_nhf", "Neutral hadron energy fraction (leading jet)", 100, 0, 1);
    TH1F *jet1_nemf = new TH1F("jet1_nemf", "Neutral em energy fraction (leading jet)", 100, 0, 1);
    TH1F *jet1_cemf = new TH1F("jet1_cemf", "Charged em energy fraction (leading jet)", 100, 0, 1);
    TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
    TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
    TH1F *chf0p5 = new TH1F("chf_Min0p5", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
    TH1F *chf0p2 = new TH1F("chf_Max0p2", "Charged energy fraction (other jet ChF < 0.2)", 100, 0, 1);
  // 	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
  // 	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
    TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
    TH2D *CHFvsPT_0p5 = new TH2D("CHFvsPT_0p5", "Jet ChF vs. pT (tag jet CHF > 0.5)", 100, 0, 5000, 100, 0, 1);
    TH2D *CHFvsPT_0p2 = new TH2D("CHFvsPT_0p2", "Jet ChF vs. pT (tag jet CHF < 0.2)", 100, 0, 5000, 100, 0, 1);
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
  
		TChain* chain = chains[j];		
		TChain* SPVchain = SPVchains[j];
    SANtuple corrjets;
    corrjets.Init(chain);
    SANtuple SPVjets;
    SPVjets.Init(SPVchain);
		
		Int_t Nentries = corrjets.fChain->GetEntries(); 
		std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
		
		double weight = SIMP_xsec[j]*lumi/Nentries;
		std::cout<<"weight "<<j<<": "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
      SPVjets.GetEntry(entry);
      corrjets.GetEntry(entry);
    
      double njet = 0;
      for (int k = 0; k < 8; ++k){
        if (corrjets.jet_pt[k] > 30) njet++;
      }
			
// 			Double_t random1 = r.Uniform();
			
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
			
			if (corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 /*&& deltajet_phi > 2*/ && (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9) && corrjets.vtx_N >= 2 && njet == 2 && pass_conv_1 && pass_conv_2 && corrjets.Flag_HBHENoiseFilter == 1 /*&& corrjets.Flag_HBHENoiseIsoFilter == 1*/ && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){
				njets->Fill(njet, weight);
				nvtx->Fill(corrjets.vtx_N, weight);
				vtx_ntracks->Fill(corrjets.vtx_nTracks[0], weight);
				HT->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3], weight);
				METOverHT->Fill(corrjets.MET/(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]), weight);
				HT_nowgt->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]);
				
				jet1_pt->Fill(corrjets.jet_pt[0], weight);
				jet2_pt->Fill(corrjets.jet_pt[1], weight);
// 				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], weight);
// 				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], weight);
				
				jet1_eta->Fill(corrjets.jet_eta[0], weight);
				jet2_eta->Fill(corrjets.jet_eta[1], weight);
				jet1_phi->Fill(corrjets.jet_phi[0], weight);
				jet2_phi->Fill(corrjets.jet_phi[1], weight);
				deltaphi->Fill(deltajet_phi, weight);
        if (CHEF_corrjet[0] < 0.5 && CHEF_corrjet[1] < 0.5 && CHEF_SPVjet[0] < 0.5 && CHEF_SPVjet[1] < 0.5) deltaphi_ChFMax0p5->Fill(deltajet_phi, weight);
				
				jet1_chf->Fill(CHEF_corrjet[0], weight);
				jet1_nhf->Fill(corrjets.jet_efrac_ne_Had[0], weight);
				jet1_nemf->Fill(corrjets.jet_efrac_ne_EM[0], weight);
				jet1_cemf->Fill(corrjets.jet_efrac_ch_EM[0], weight);
				if (CHEF_corrjet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_corrjet[0], weight);
				jet2_chf->Fill(CHEF_corrjet[1], weight);
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
				
// 				eta1vsnpixhits->Fill(nPixHits[0], jet_eta[0], weight);
// 				eta2vsnpixhits->Fill(nPixHits[0], jet_eta[1], weight);
				
// 				NVTXvsCHF->Fill(CHEF_jet[0], vtx_N, weight);
// 				NVTXvsCHF->Fill(CHEF_jet[1], vtx_N, weight);
// 				NVTXvsCHF1->Fill(CHEF_jet[0], vtx_N, weight);
// 				NVTXvsCHF2->Fill(CHEF_jet[1], vtx_N, weight);
// 				
// 				Chi2vsCHF->Fill(CHEF_jet[0], chi2[0], weight);
// 				Chi2vsCHF->Fill(CHEF_jet[1], chi2[0], weight);
// 				Chi2vsCHF1->Fill(CHEF_jet[0], chi2[0], weight);
// 				Chi2vsCHF2->Fill(CHEF_jet[1], chi2[0], weight);
// 				
// 				ndofvsCHF->Fill(CHEF_jet[0], ndof[0], weight);
// 				ndofvsCHF->Fill(CHEF_jet[1], ndof[0], weight);
// 				ndofvsCHF1->Fill(CHEF_jet[0], ndof[0], weight);
// 				ndofvsCHF2->Fill(CHEF_jet[1], ndof[0], weight);
				
// 				Pt1vsCHF1->Fill(CHEF_jet[0], jet_pt[0], weight);
// 				Pt2vsCHF2->Fill(CHEF_jet[1], jet_pt[1], weight);
// 				eta1vsCHF1->Fill(CHEF_jet[0], jet_eta[0], weight);
// 				eta2vsCHF2->Fill(CHEF_jet[1], jet_eta[1], weight);				
// 				Pt2vsCHF1->Fill(CHEF_jet[0], jet_pt[1], weight);
// 				Pt1vsCHF2->Fill(CHEF_jet[1], jet_pt[0], weight);
// 				eta2vsCHF1->Fill(CHEF_jet[0], jet_eta[1], weight);
// 				eta1vsCHF2->Fill(CHEF_jet[1], jet_eta[0], weight);
				
// 				track1_nhits->Fill(nhits[0], weight);
// 				track1_npixhits->Fill(nPixHits[0], weight);
// 				track1_pt->Fill(track_pt[0], weight);
// 				track1_ptError->Fill(track_ptError[0], weight);
// 				track1_dzError->Fill(track_dzError[0], weight);
// 				track1_dz->Fill(track_dz[0], weight);
// 				track1_dzRelError->Fill(track_dzError[0]/track_dz[0], weight);
// 				track1_reldzErrorvsdz->Fill(track_dz[0], track_dzError[0]/track_dz[0], weight);
// 				track1_fromPV->Fill(track_fromPV[0], weight);
				
// 				if (ndof[0] == 0) std::cout<<vtx_x[0]<<"\t"<<vtx_y[0]<<"\t"<<vtx_z[0]<<"\t"<<vtx_d0[0]<<std::endl;
				
// 				Chi2vsNdof->Fill(ndof[0], chi2[0], weight);
				
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
  output->Write();
  output->Close();
// 	}
}
