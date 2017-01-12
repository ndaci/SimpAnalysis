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

#include "lists/list_QCD_300To500_PUMoriond17_AOD.h"
#include "lists/list_QCD_500To700_PUMoriond17_AOD.h"
#include "lists/list_QCD_700To1000_PUMoriond17_AOD.h"
#include "lists/list_QCD_1000To1500_PUMoriond17_AOD.h"
#include "lists/list_QCD_1500To2000_PUMoriond17_AOD.h"

void SIMP_QCD_macro(){
  
  TChain* chain0 = new TChain("tree/SimpAnalysis");
	list_QCD_300To500(chain0);
  TChain* chain1 = new TChain("tree/SimpAnalysis");
	list_QCD_500To700(chain1);
  TChain* chain2 = new TChain("tree/SimpAnalysis");
	list_QCD_700To1000(chain2);
  TChain* chain3 = new TChain("tree/SimpAnalysis");
	list_QCD_1000To1500(chain3);
  TChain* chain4 = new TChain("tree/SimpAnalysis");
	list_QCD_1500To2000(chain4);
	TChain* chains[5] = {chain0, chain1, chain2, chain3, chain4};
	std::cout<<"TChains ready"<<std::endl;
  
  TFile *output = new TFile("QCD_plots_prescaledlumi0p062_PUMoriond17_AOD.root", "RECREATE");
  
// 	TRandom3 r;
	int nJet, vtx_N;
	double pswgt_dijet_170;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], chi2[3], vtx_x[3], vtx_y[3], vtx_z[3], vtx_d0[3];
	double CHEF_jet[8], genjet_efrac_ch[4];
	double track_pt[10], track_dzError[10], track_dz[10], track_ptError[10];
	int track_fromPV[10], nhits[10], nPixHits[10], ndof[3];
// 	E[4], M[4];
  
	TH1F *njets = new TH1F("njets", "Number of jets", 11, -0.5, 10.5);
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 100, 0, 2000);
	TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 100, 0, 2000);
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
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
// 	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
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
	
// 	double QCD_xsec[6] = {343500, 32050, 6791, 1214, 118.7, 24.91}; //in pb//Spring16
	double QCD_xsec[6] = {346400, 32010, 6842, 1203, 120.1, 25.40}; //PUMoriond17
// 	double QCD_events[5] = {16830696, 19199088, 15621634, 4980387, 3846616};
	double lumi = 0.062*1000; // in 1000 fb^-1 = pb^-1
  
	for (int j = 0; j < 5; j++){
		TChain* chain = chains[j];		
		chain->SetBranchAddress("vtx_N", &vtx_N);
		chain->SetBranchAddress("nJet",&nJet);
		chain->SetBranchAddress("jet_pt", &jet_pt);
		chain->SetBranchAddress("jet_eta", &jet_eta);
		chain->SetBranchAddress("jet_phi", &jet_phi);
		chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
		chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
		chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
		chain->SetBranchAddress("genjet_efrac_ch", &genjet_efrac_ch);
		chain->SetBranchAddress("vtx_normalizedChi2", &chi2);
// 		chain->SetBranchAddress("vtx_d0", &vtx_d0);
// 		chain->SetBranchAddress("vtx_x", &vtx_x);
// 		chain->SetBranchAddress("vtx_y", &vtx_y);
// 		chain->SetBranchAddress("vtx_z", &vtx_z);
		chain->SetBranchAddress("vtx_ndof", &ndof);
		chain->SetBranchAddress("track_ptError", &track_ptError);
		chain->SetBranchAddress("track_dzError", &track_dzError);
		chain->SetBranchAddress("track_dz", &track_dz);
		chain->SetBranchAddress("track_nhits", &nhits);
		chain->SetBranchAddress("track_nPixHits", &nPixHits);
		chain->SetBranchAddress("track_pt", &track_pt);
		chain->SetBranchAddress("track_fromPV", &track_fromPV);
		
		Int_t Nentries = chain->GetEntries(); 
		std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
		
		double weight = QCD_xsec[j]*lumi/Nentries;
		std::cout<<"weight "<<j<<": "<<weight<<std::endl;
		
		for(Int_t entry = 0; entry < Nentries; ++entry){
			if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
			chain->GetEntry(entry);
			
// 			Double_t random1 = r.Uniform();
			
			double deltajet_phi = jet_phi[0] - jet_phi[1];
			if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
			if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
			
			for (int i = 0; i < 8; i++){
				CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
			} 
			
			output->cd();
			
			if (jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2 /*&& nPixHits[0] > 0*/){
				njets->Fill(nJet, weight);
				nvtx->Fill(vtx_N, weight);
				HT->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3], weight);
				HT_nowgt->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3]);
				
				jet1_pt->Fill(jet_pt[0], weight);
				jet2_pt->Fill(jet_pt[1], weight);
// 				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], weight);
// 				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], weight);
				
				jet1_eta->Fill(jet_eta[0], weight);
				jet2_eta->Fill(jet_eta[1], weight);
				jet1_phi->Fill(jet_phi[0], weight);
				jet2_phi->Fill(jet_phi[1], weight);
				deltaphi->Fill(deltajet_phi, weight);
				
				jet1_chf->Fill(CHEF_jet[0], weight);
// 				if (CHEF_jet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_jet[0], weight);
				jet2_chf->Fill(CHEF_jet[1], weight);
// 				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2], weight);
// 				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3], weight);
// 				CHF->Fill(CHEF_jet[0], weight);
// 				CHF->Fill(CHEF_jet[1], weight);
				CHFvsCHF->Fill(CHEF_jet[1], CHEF_jet[0], weight);				
				CHFvsPT->Fill(jet_pt[0], CHEF_jet[0], weight);
				CHFvsPT->Fill(jet_pt[1], CHEF_jet[1], weight);
				
				genjet1_chf->Fill(genjet_efrac_ch[0], weight);
				genjet2_chf->Fill(genjet_efrac_ch[1], weight);
				CHFvsGenCHF1->Fill(genjet_efrac_ch[0], CHEF_jet[0], weight);
				CHFvsGenCHF2->Fill(genjet_efrac_ch[1], CHEF_jet[1], weight);
				relChF1->Fill((CHEF_jet[0]-genjet_efrac_ch[0])/genjet_efrac_ch[0], weight);
				relChF2->Fill((CHEF_jet[1]-genjet_efrac_ch[1])/genjet_efrac_ch[1], weight);
				
				eta1vsnpixhits->Fill(nPixHits[0], jet_eta[0], weight);
				eta2vsnpixhits->Fill(nPixHits[0], jet_eta[1], weight);
				
				NVTXvsCHF->Fill(CHEF_jet[0], vtx_N, weight);
				NVTXvsCHF->Fill(CHEF_jet[1], vtx_N, weight);
				NVTXvsCHF1->Fill(CHEF_jet[0], vtx_N, weight);
				NVTXvsCHF2->Fill(CHEF_jet[1], vtx_N, weight);
				
				Chi2vsCHF->Fill(CHEF_jet[0], chi2[0], weight);
				Chi2vsCHF->Fill(CHEF_jet[1], chi2[0], weight);
				Chi2vsCHF1->Fill(CHEF_jet[0], chi2[0], weight);
				Chi2vsCHF2->Fill(CHEF_jet[1], chi2[0], weight);
				
				ndofvsCHF->Fill(CHEF_jet[0], ndof[0], weight);
				ndofvsCHF->Fill(CHEF_jet[1], ndof[0], weight);
				ndofvsCHF1->Fill(CHEF_jet[0], ndof[0], weight);
				ndofvsCHF2->Fill(CHEF_jet[1], ndof[0], weight);
				
// 				Pt1vsCHF1->Fill(CHEF_jet[0], jet_pt[0], weight);
// 				Pt2vsCHF2->Fill(CHEF_jet[1], jet_pt[1], weight);
// 				eta1vsCHF1->Fill(CHEF_jet[0], jet_eta[0], weight);
// 				eta2vsCHF2->Fill(CHEF_jet[1], jet_eta[1], weight);				
// 				Pt2vsCHF1->Fill(CHEF_jet[0], jet_pt[1], weight);
// 				Pt1vsCHF2->Fill(CHEF_jet[1], jet_pt[0], weight);
// 				eta2vsCHF1->Fill(CHEF_jet[0], jet_eta[1], weight);
// 				eta1vsCHF2->Fill(CHEF_jet[1], jet_eta[0], weight);
				
				track1_nhits->Fill(nhits[0], weight);
				track1_npixhits->Fill(nPixHits[0], weight);
				track1_pt->Fill(track_pt[0], weight);
				track1_ptError->Fill(track_ptError[0], weight);
				track1_dzError->Fill(track_dzError[0], weight);
				track1_dz->Fill(track_dz[0], weight);
				track1_dzRelError->Fill(track_dzError[0]/track_dz[0], weight);
				track1_reldzErrorvsdz->Fill(track_dz[0], track_dzError[0]/track_dz[0], weight);
// 				track1_fromPV->Fill(track_fromPV[0], weight);
				
// 				if (ndof[0] == 0) std::cout<<vtx_x[0]<<"\t"<<vtx_y[0]<<"\t"<<vtx_z[0]<<"\t"<<vtx_d0[0]<<std::endl;
				
				Chi2vsNdof->Fill(ndof[0], chi2[0], weight);
				
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
				if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] > 0.5) chf0p5->Fill(CHEF_jet[1], weight);
				if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] > 0.5) chf0p5->Fill(CHEF_jet[0], weight);
				if (jet_pt[1] < 300 && fabs(jet_eta[1]) < 1.0 && CHEF_jet[0] < 0.2) chf0p2->Fill(CHEF_jet[1], weight);
				if (jet_pt[0] < 300 && fabs(jet_eta[0]) < 1.0 && CHEF_jet[1] < 0.2) chf0p2->Fill(CHEF_jet[0], weight);
				
				if (CHEF_jet[1] > 0.5) CHFvsPT_0p5->Fill(jet_pt[0], CHEF_jet[0], weight);
				if (CHEF_jet[0] > 0.5) CHFvsPT_0p5->Fill(jet_pt[1], CHEF_jet[1], weight);
				if (CHEF_jet[1] < 0.2) CHFvsPT_0p2->Fill(jet_pt[0], CHEF_jet[0], weight);
				if (CHEF_jet[0] < 0.2) CHFvsPT_0p2->Fill(jet_pt[1], CHEF_jet[1], weight);
				
// 				if(random1 <= 0.5 && CHEF_jet[0] > 0.5) ChFOtherJet->Fill(CHEF_jet[1], weight);
// 				else if (random1 > 0.5 && CHEF_jet[1] > 0.5) ChFOtherJet->Fill(CHEF_jet[0], weight);
			}    
		}
	}
  output->Write();
  output->Close();
}
