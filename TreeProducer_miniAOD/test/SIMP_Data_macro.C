#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

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

void SIMP_Data_macro(){
  
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
	std::cout<<"TChain ready"<<std::endl;

  TFile *output = new TFile("plots_CR_rereco_photonVeto_withprescale_runH.root", "RECREATE");
  
	TRandom3 *r = new TRandom3();
	int dijet_170, vtx_N;
	double pswgt_dijet_170, MET;
//   double jet_pt[4], jet_eta[4], jet_phi[4], jet_efrac_ch_Had[4], jet_efrac_ch_EM[4], jet_efrac_ch_Mu[4], CHEF_jet[4];
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ne_Had[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], jet_efrac_ne_EM[8], chi2[3], vtx_x[3], vtx_y[3], vtx_z[3], vtx_d0[3], CHEF_jet[8], EMF_jet[8];
	double track_pt[10], track_ptError[10], track_dzError[10], track_dz[10], genjet_efrac_ch[4];
	int track_fromPV[10], nhits[10], nPixHits[10], ndof[3];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	int selection = 0;
	int CR_region = 0;
  
	chain->SetBranchAddress("MET", &MET);
  chain->SetBranchAddress("vtx_N", &vtx_N);
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_phi", &jet_phi);
	chain->SetBranchAddress("genjet_efrac_ch", &genjet_efrac_ch);
  chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
  chain->SetBranchAddress("jet_efrac_ne_Had", &jet_efrac_ne_Had);
  chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
	chain->SetBranchAddress("jet_efrac_ne_EM", &jet_efrac_ne_EM);
  chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
	chain->SetBranchAddress("vtx_normalizedChi2", &chi2);
	chain->SetBranchAddress("vtx_d0", &vtx_d0);
	chain->SetBranchAddress("vtx_x", &vtx_x);
	chain->SetBranchAddress("vtx_y", &vtx_y);
	chain->SetBranchAddress("vtx_z", &vtx_z);
	chain->SetBranchAddress("vtx_ndof", &ndof);
	chain->SetBranchAddress("track_ptError", &track_ptError);
	chain->SetBranchAddress("track_nhits", &nhits);
	chain->SetBranchAddress("track_nPixHits", &nPixHits);
	chain->SetBranchAddress("track_pt", &track_pt);
	chain->SetBranchAddress("track_fromPV", &track_fromPV);
	chain->SetBranchAddress("track_dzError", &track_dzError);
	chain->SetBranchAddress("track_dz", &track_dz);
	chain->SetBranchAddress("photon_pt", &photon_pt);
	chain->SetBranchAddress("photon_eta", &photon_eta);
	chain->SetBranchAddress("photon_phi", &photon_phi);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
  chain->SetBranchAddress("pswgt_dijet_170", &pswgt_dijet_170);
  chain->SetBranchAddress("photon_passLooseId",&passLooseId);
  chain->SetBranchAddress("photon_passMediumId",&passMediumId);
  chain->SetBranchAddress("photon_passTightId",&passTightId);
	
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 100, 0, 2000);
	TH1F *METOverHT = new TH1F("METOverHT", "MET / HT(4 jets)", 100, 0, 1);	
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *nvtx_withps = new TH1F("nvtx_withps", "Number of vertices (weighted for prescale)", 51, -0.5, 50.5);
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet1_pt_withps = new TH1F("jet1_pt_withps", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt_withps = new TH1F("jet2_pt_withps", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *jetA_pt = new TH1F("jetA_pt", "ChF-leading jet pt", 100, 0, 2000);
	TH1F *jetB_pt = new TH1F("jetB_pt", "ChF-subleading jet pt", 100, 0, 2000);
	TH1F *jetA_eta = new TH1F("jetA_eta", "ChF-leading jet eta", 100, -3.14, 3.14);
	TH1F *jetB_eta = new TH1F("jetB_eta", "ChF-subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *jet1_emf = new TH1F("jet1_emf", "EM energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_emf = new TH1F("jet2_emf", "EM energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet1_chhf = new TH1F("jet1_chhf", "Charged hadronic energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chhf = new TH1F("jet2_chhf", "Charged hadronic energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet1_nhf = new TH1F("jet1_nhf", "Neutral hadronic energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_nhf = new TH1F("jet2_nhf", "Neutral hadronic energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet1_chf_withps = new TH1F("jet1_chf_withps", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf_withps = new TH1F("jet2_chf_withps", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *chf = new TH1F("chf", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsCHF = new TH2D("CHFvsCHF", "Leading jet ChF vs. Subleading jet ChF", 100, 0, 1, 100, 0, 1);
	CHFvsCHF->GetXaxis()->SetTitle("jet2 ChF");
	CHFvsCHF->GetYaxis()->SetTitle("jet1 ChF");
	TH1D *CHF = new TH1D("CHF", "Jet ChF (both jets)", 100, 0, 1);
	TH1D *ChFOtherJet = new TH1D("ChFOtherJet", "Jet ChF (other jet ChF > 0.5)", 100, 0, 1);
  TH1F *photon1_pt = new TH1F("photon1_pt", "Photon p_{T} (leading photon)", 100, 0, 500);
  TH1F *photon2_pt = new TH1F("photon2_pt", "Photon p_{T} (subleading photon)", 100, 0, 500);
  TH1F *photon1_eta = new TH1F("photon1_eta", "Photon #eta (leading photon)", 100, -3.14, 3.14);
  TH1F *photon2_eta = new TH1F("photon2_eta", "Photon #eta (subleading photon)", 100, -3.14, 3.14);
  TH1F *photon1_phi = new TH1F("photon1_phi", "Photon #phi (leading photon)", 100, -3.14, 3.14);
  TH1F *photon2_phi = new TH1F("photon2_phi", "Photon #phi (subleading photon)", 100, -3.14, 3.14);
  TH1F *dR_1 = new TH1F("dR_1", "dR_{#gamma, j1}", 100, -1, 4);
  TH1F *dR_2 = new TH1F("dR_2", "dR_{#gamma, j2}", 100, -1, 4);
  TH2D *dRvsdR = new TH2D("dRvsdR", "dR_{#gamma, j1} vs. dR_{#gamma, j2}", 100, -1, 4, 100, -1, 4);
	dRvsdR->GetXaxis()->SetTitle("dR_{#gamma, j2}");
	dRvsdR->GetYaxis()->SetTitle("dR_{#gamma, j1}");
  TH2D *dRvsChF1 = new TH2D("dRvsChF1", "dR_{#gamma, j1} vs. ChF_{j1}", 100, 0, 1, 100, -1, 4);
  TH2D *dRvsChF2 = new TH2D("dRvsChF2", "dR_{#gamma, j2} vs. Chf_{j2}", 100, 0, 1, 100, -1, 4);
  
	std::cout<<"Getting entries....";
  Int_t Nentries = chain->GetEntries(); 
	std::cout<<"processing "<<Nentries<<"entries"<<std::endl;
  for(Int_t entry = 0; entry < Nentries; ++entry){
		if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
    chain->GetEntry(entry);
			
		Double_t random1 = r->Uniform();
    
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
			EMF_jet[i] = jet_efrac_ch_EM[i]+jet_efrac_ne_EM[i];
		}
		
		output->cd();
	
		if (dijet_170 == 1 && jet_pt[0] > 250 && jet_pt[1] > 250 && fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0 && deltajet_phi > 2 && nPixHits[0] > 0 && (passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
			selection ++;
			if (CHEF_jet[0] > 0.5||CHEF_jet[1] > 0.5){
				CR_region++;
				HT->Fill(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3], pswgt_dijet_170);
				METOverHT->Fill(MET/(jet_pt[0]+jet_pt[1]+jet_pt[2]+jet_pt[3]), pswgt_dijet_170);
				nvtx->Fill(vtx_N, pswgt_dijet_170);
				
				jet1_pt->Fill(jet_pt[0], pswgt_dijet_170);
				jet2_pt->Fill(jet_pt[1], pswgt_dijet_170);
				if (jet_pt[2] > 20) jet3_pt->Fill(jet_pt[2], pswgt_dijet_170);
				if (jet_pt[3] > 20) jet4_pt->Fill(jet_pt[3], pswgt_dijet_170);
				
				jet1_eta->Fill(jet_eta[0], pswgt_dijet_170);
				jet2_eta->Fill(jet_eta[1], pswgt_dijet_170);
				jet1_phi->Fill(jet_phi[0], pswgt_dijet_170);
				jet2_phi->Fill(jet_phi[1], pswgt_dijet_170);
				deltaphi->Fill(deltajet_phi, pswgt_dijet_170);
				
				jet1_emf->Fill(EMF_jet[0], pswgt_dijet_170);
				jet2_emf->Fill(EMF_jet[1], pswgt_dijet_170);
				jet1_chf->Fill(CHEF_jet[0], pswgt_dijet_170);
				jet2_chf->Fill(CHEF_jet[1], pswgt_dijet_170);
				jet1_nhf->Fill(jet_efrac_ne_Had[0], pswgt_dijet_170);
				jet2_nhf->Fill(jet_efrac_ne_Had[1], pswgt_dijet_170);
				jet1_chhf->Fill(jet_efrac_ch_Had[0], pswgt_dijet_170);
				jet2_chhf->Fill(jet_efrac_ch_Had[1], pswgt_dijet_170);
				if (CHEF_jet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_jet[0], pswgt_dijet_170);
				if (jet_pt[2] > 20) jet3_chf->Fill(CHEF_jet[2], pswgt_dijet_170);
				if (jet_pt[3] > 20) jet4_chf->Fill(CHEF_jet[3], pswgt_dijet_170);
				CHFvsPT->Fill(jet_pt[0], CHEF_jet[0], pswgt_dijet_170);
				CHFvsPT->Fill(jet_pt[1], CHEF_jet[1], pswgt_dijet_170);
				CHF->Fill(CHEF_jet[0], pswgt_dijet_170);
				CHF->Fill(CHEF_jet[1], pswgt_dijet_170);
				CHFvsCHF->Fill(CHEF_jet[1], CHEF_jet[0], pswgt_dijet_170);
				if(CHEF_jet[0] > CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[0], pswgt_dijet_170);
					jetB_pt->Fill(jet_pt[1], pswgt_dijet_170);
					jetA_eta->Fill(jet_eta[0], pswgt_dijet_170);
					jetB_eta->Fill(jet_eta[1], pswgt_dijet_170);
					if (CHEF_jet[0] > 0.5) chf->Fill(CHEF_jet[1], pswgt_dijet_170);
				}
				else if(CHEF_jet[0] < CHEF_jet[1]){
					jetA_pt->Fill(jet_pt[1], pswgt_dijet_170);
					jetB_pt->Fill(jet_pt[0], pswgt_dijet_170);
					jetA_eta->Fill(jet_eta[1], pswgt_dijet_170);
					jetB_eta->Fill(jet_eta[0], pswgt_dijet_170);
					if (CHEF_jet[1] > 0.5) chf->Fill(CHEF_jet[0], pswgt_dijet_170);
				}
				if(random1 <= 0.5 && CHEF_jet[0] > 0.5) ChFOtherJet->Fill(CHEF_jet[1], pswgt_dijet_170);
				else if (random1 > 0.5 && CHEF_jet[1] > 0.5) ChFOtherJet->Fill(CHEF_jet[0], pswgt_dijet_170);
				
				
				if (passLooseId[0] == 1){
					dR_1->Fill(dR1, pswgt_dijet_170);
					dR_2->Fill(dR2, pswgt_dijet_170);
					dRvsChF1->Fill(CHEF_jet[0], dR1, pswgt_dijet_170);
					dRvsChF2->Fill(CHEF_jet[1], dR2, pswgt_dijet_170);
					dRvsdR->Fill(dR2, dR1, pswgt_dijet_170);
				}else{
					dR_1->Fill(-1, pswgt_dijet_170);
					dR_2->Fill(-1, pswgt_dijet_170);
					dRvsChF1->Fill(CHEF_jet[0], -1, pswgt_dijet_170);
					dRvsChF2->Fill(CHEF_jet[1], -1, pswgt_dijet_170);
					dRvsdR->Fill(-1, -1, pswgt_dijet_170);
				}
				
				if(photon_pt[0] > 10 && passLooseId[0] == 1) photon1_pt->Fill(photon_pt[0], pswgt_dijet_170);
				if(photon_pt[1] > 10 && passLooseId[1] == 1) photon2_pt->Fill(photon_pt[1], pswgt_dijet_170);
				if(photon_pt[0] > 10 && passLooseId[0] == 1) photon1_eta->Fill(photon_eta[0], pswgt_dijet_170);
				if(photon_pt[1] > 10 && passLooseId[1] == 1) photon2_eta->Fill(photon_eta[1], pswgt_dijet_170);
				if(photon_pt[0] > 10 && passLooseId[0] == 1) photon1_phi->Fill(photon_phi[0], pswgt_dijet_170);
				if(photon_pt[1] > 10 && passLooseId[1] == 1) photon2_phi->Fill(photon_phi[1], pswgt_dijet_170);
				
			}
		}
	}
	std::cout<<"selected: "<<selection<<std::endl;
	std::cout<<"CR region: "<<CR_region<<std::endl;
  output->Write();
  output->Close();
}