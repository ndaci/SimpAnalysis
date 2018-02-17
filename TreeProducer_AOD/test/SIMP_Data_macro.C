#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

#include "SANtuple.h"

// #include "lists/list_JetHT_2016B_new.h"
// #include "lists/list_JetHT_2016C_new.h"
// #include "lists/list_JetHT_2016D_new.h"
// #include "lists/list_JetHT_2016E_new.h"
// #include "lists/list_JetHT_2016F_new.h"
// #include "lists/list_JetHT_2016G_new.h"
// #include "lists/list_JetHT_2016H2_new.h"
// #include "lists/list_JetHT_2016H3_new.h"

// #include "lists/list_JetHT_2016G_newtrigger.h"
// #include "lists/list_JetHT_2016H2_newtrigger.h"
// #include "lists/list_JetHT_2016H3_newtrigger.h"

#include "lists/list_JetHT_2016G_newtrigger_2.h"
#include "lists/list_JetHT_2016H2_newtrigger_2.h"
#include "lists/list_JetHT_2016H3_newtrigger_2.h"

// #include "lists/list_GJets.h"

void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void SIMP_Data_macro(){
	  
  TChain* chain = new TChain("treeCorr/SimpAnalysis");
// 	list_JetHT_2016B(chain);
// 	list_JetHT_2016C(chain);
// 	list_JetHT_2016D(chain);
// 	list_JetHT_2016E(chain);
// 	list_JetHT_2016F(chain);
	list_JetHT_2016G(chain);
	list_JetHT_2016H2(chain);
	list_JetHT_2016H3(chain);
// 	list_GJets(chain);
//   chain->Add("/pnfs/iihe/cms/store/user/isdebruy/JetHT/SIMPs_JetHT_2016G_rereco_AOD_test/170802_172813/0000/Data_AOD_newtrigger_1.root");
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
//   list_JetHT_2016B(SPVchain);
//   list_JetHT_2016C(SPVchain);
//   list_JetHT_2016D(SPVchain);
//   list_JetHT_2016E(SPVchain);
//   list_JetHT_2016F(SPVchain);
	list_JetHT_2016G(SPVchain);
	list_JetHT_2016H2(SPVchain);
	list_JetHT_2016H3(SPVchain); 
// 	list_GJets(SPVchain); 
//   SPVchain->Add("/pnfs/iihe/cms/store/user/isdebruy/JetHT/SIMPs_JetHT_2016G_rereco_AOD_test/170802_172813/0000/Data_AOD_newtrigger_1.root");
  
  TFile *output = new TFile("plots_Data_RunGH_AOD_unblinded.root", "RECREATE");
//   TFile *output = new TFile("plots_GJets_nophotonveto.root", "RECREATE");
  
  double chf_cuts[25] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0};
  double chf_bins[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.06, 0.085, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.95};
  double CHEF_SPVjet[8], CHEF_corrjet[8], HoverE[8], EMF[8];
  int photon_nr;
  int counter = 0;
  
	TH1F *nJets = new TH1F("nJets", "Number of jets", 11, -0.5, 10.5);
	TH2F *nJetsVsNEMF1 = new TH2F("nJetsVsNEMF1", "Number of jets vs. NEMF1", 11, -0.5, 10.5, 100, 0, 1);
	TH2F *nJetsVsNEMF2 = new TH2F("nJetsVsNEMF2", "Number of jets vs. NEMF2", 11, -0.5, 10.5, 100, 0, 1);
	TH1F *nvtx = new TH1F("nvtx", "Number of vertices", 51, -0.5, 50.5);
	TH1F *HT = new TH1F("HT", "HT (4 jets)", 150, 1000, 4000);
	TH1F *HT_nowgt = new TH1F("HT_noweight", "HT (4 jets) without correct weights", 150, 1000, 4000);
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet1_pt_SPV = new TH1F("jet1_pt_SPV", "Leading jet pt SPV collection", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_eta_SPV = new TH1F("jet1_eta_SPV", "Leading jet eta SPV collection", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *deltaphi_ChFMax0p5 = new TH1F("DeltaPhi_ChFMax0p5", "DeltaPhi (ChF < 0.5)", 100, 0, 3.14);
	TH2F *deltaphiVsChF1 = new TH2F("DeltaPhiVsChF1", "DeltaPhiVsChF1", 100, 0, 3.14, 100, 0, 1);
	TH2F *deltaphiVsChF2 = new TH2F("DeltaPhiVsChF2", "DeltaPhiVsChF2", 100, 0, 3.14, 100, 0, 1);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet_chf_passed = new TH1F("jet_chf_passed", "Charged energy fraction (passing MET filters)", 100, 0, 1);  
	TH1F *jet_chf_failed = new TH1F("jet_chf_failed", "Charged energy fraction (failing MET filters)", 100, 0, 1);  
	TH1F *jet_chf_passed_HBHENoise = new TH1F("jet_chf_passed_HBHENoise", "Charged energy fraction (passing HBHENoise filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_HBHENoise = new TH1F("jet_chf_failed_HBHENoise", "Charged energy fraction (failing HBHENoise filters)", 100, 0, 1); 
	TH1F *jet_chf_passed_HBHENoiseIso = new TH1F("jet_chf_passed_HBHENoiseIso", "Charged energy fraction (passing HBHENoiseIso filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_HBHENoiseIso = new TH1F("jet_chf_failed_HBHENoiseIso", "Charged energy fraction (failing HBHENoiseIso filters)", 100, 0, 1); 
	TH1F *jet_chf_passed_EcalDeadCellTriggerPrimitiveFilter = new TH1F("jet_chf_passed_EcalDeadCellTriggerPrimitiveFilter", "Charged energy fraction (passing EcalDeadCellTriggerPrimitiveFilter filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_EcalDeadCellTriggerPrimitiveFilter = new TH1F("jet_chf_failed_EcalDeadCellTriggerPrimitiveFilter", "Charged energy fraction (failing EcalDeadCellTriggerPrimitiveFilter filters)", 100, 0, 1); 
	TH1F *jet_chf_passed_goodVertices = new TH1F("jet_chf_passed_goodVertices", "Charged energy fraction (passing goodVertices filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_goodVertices = new TH1F("jet_chf_failed_goodVertices", "Charged energy fraction (failing goodVertices filters)", 100, 0, 1); 
	TH1F *jet_chf_passed_eeBadScFilter = new TH1F("jet_chf_passed_eeBadScFilter", "Charged energy fraction (passing eeBadScFilter filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_eeBadScFilter = new TH1F("jet_chf_failed_eeBadScFilter", "Charged energy fraction (failing eeBadScFilter filters)", 100, 0, 1); 
	TH1F *jet_chf_passed_globalTightHalo2016Filter = new TH1F("jet_chf_passed_globalTightHalo2016Filter", "Charged energy fraction (passing globalTightHalo2016Filter filters)", 100, 0, 1);  
	TH1F *jet_chf_failed_globalTightHalo2016Filter = new TH1F("jet_chf_failed_globalTightHalo2016Filter", "Charged energy fraction (failing globalTightHalo2016Filter filters)", 100, 0, 1);  
  TH1F *jet1_HOverE = new TH1F("jet1_HOverE", "Leading jet H/E", 100, 0, 1);
  TH1F *jet2_HOverE = new TH1F("jet2_HOverE", "Subleading jet H/E", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet1_chf_SPV = new TH1F("jet1_chf_SPV", "Charged energy fraction (leading jet) SPV collection", 100, 0, 1);
	TH1F *jet2_chf_SPV = new TH1F("jet2_chf_SPV", "Charged energy fraction (subleading jet) SPV collection", 100, 0, 1);
	TH1F *ChfOverSPVChF = new TH1F("ChfOverSPVChF", "ChF/ChF(SPV)", 1000, 0, 10);
	TH1F *chf0p5 = new TH1F("chf_Min0p5", "Charged energy fraction (other jet ChF > 0.5)", 100, 0, 1);
	TH1F *chf0p2 = new TH1F("chf_Max0p2", "Charged energy fraction (other jet ChF < 0.2)", 100, 0, 1);
	TH2D *CHFvsPT = new TH2D("CHFvsPT", "Jet ChF vs. pT (both jets)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsPT_0p5 = new TH2D("CHFvsPT_0p5", "Jet ChF vs. pT (tag jet CHF > 0.5)", 100, 0, 5000, 100, 0, 1);
	TH2D *CHFvsPT_0p2 = new TH2D("CHFvsPT_0p2", "Jet ChF vs. pT (tag jet CHF < 0.2)", 100, 0, 5000, 100, 0, 1);
	TH2D *ChfVsSPVChf = new TH2D("ChfVsSPVChf", "jet ChF vs. SPV jet ChF", 100, 0, 1, 100, 0, 1);
	TH2D *CHFvsCHF = new TH2D("CHFvsCHF", "Leading jet ChF vs. Subleading jet ChF", 100, 0, 1, 100, 0, 1);
  TH2D* ChFvsNvtx = new TH2D("ChFvsNvtx", "jet ChF vs. #vtx", 100, 0, 1, 100, 0, 50);
  TH2D* ChFvsNvtx1 = new TH2D("ChFvsNvtx1", "leading jet ChF vs. #vtx", 100, 0, 1, 100, 0, 50);
  TH2D* ChFvsNvtx2 = new TH2D("ChFvsNvtx2", "subleading jet ChF vs. #vtx", 100, 0, 1, 100, 0, 50);
	CHFvsCHF->GetXaxis()->SetTitle("jet2 ChF");
	CHFvsCHF->GetYaxis()->SetTitle("jet1 ChF");
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
	TH2D *CHFvsGenCHF1 = new TH2D("CHFvsGenCHF1", "ChF vs. ChF GenJet (leading jet)", 100, 0, 1, 100, 0, 1);
	TH2D *CHFvsGenCHF2 = new TH2D("CHFvsGenCHF2", "ChF vs. ChF GenJet (subleading jet)", 100, 0, 1, 100, 0, 1);
	TH1F *relChF1 = new TH1F("relChF1", "(ChF - ChF GenJet)/ChF GenJet (leading jet)", 1000, -10, 10);
	TH1F *relChF2 = new TH1F("relChF2", "(ChF - ChF GenJet)/ChF GenJet (subleading jet)", 1000, -10, 10);
	TH1F *genjet1_chf = new TH1F("genjet1_chf", "Charged energy fraction (leading genjet)", 100, 0, 1);
	TH1F *genjet2_chf = new TH1F("genjet2_chf", "Charged energy fraction (subleading genjet)", 100, 0, 1);
	TH2D *eta1vsnpixhits = new TH2D("eta1vsnpixhits", "eta vs. # pixel hits (leading jet)", 11, -0.5, 10.5, 100, -3.14, 3.14);
	TH2D *eta2vsnpixhits = new TH2D("eta2vsnpixhits", "eta vs. # pixel hits (subleading jet)", 11, -0.5, 10.5, 100, -3.14, 3.14);
  TH1F *photon1_ID = new TH1F("photon1_ID", "Photon loose ID (leading photon)", 2, -0.5, 1.5);
  TH1F *photon1_pt = new TH1F("photon1_pt", "Photon p_{T} (leading photon)", 100, 0, 500);
  TH1F *photon1_eta = new TH1F("photon1_eta", "Photon #eta (leading photon)", 100, -3.14, 3.14);
  TH1F *photon1_phi = new TH1F("photon1_phi", "Photon #phi (leading photon)", 100, -3.14, 3.14);
  TH1F *photon1_dR1 = new TH1F("photon1_dR1", "Photon dR1 (leading photon)", 100, 0, 3.14);
  TH1F *photon1_dR2 = new TH1F("photon1_dR2", "Photon dR2 (leading photon)", 100, 0, 3.14);
  TH1F *photon2_ID = new TH1F("photon2_ID", "Photon loose ID (subleading photon)", 2, -0.5, 1.5);
  TH1F *photon2_pt = new TH1F("photon2_pt", "Photon p_{T} (subleading photon)", 100, 0, 500);
  TH1F *photon2_eta = new TH1F("photon2_eta", "Photon #eta (subleading photon)", 100, -3.14, 3.14);
  TH1F *photon2_phi = new TH1F("photon2_phi", "Photon #phi (subleading photon)", 100, -3.14, 3.14);
  TH1F *ntracks = new TH1F("ntracks", "total number of tracks", 201, -0.5, 2000.5);
  TH1F *nhits = new TH1F("nhits", "total nhits for first 10 tracks", 251, -0.5, 250.5);
  TH1F *nhitspertrack = new TH1F("nhitspertrack", "nhits per track", 51, -0.5, 50.5);
  
  TH2F *jet_photon_frac = new TH2F("jet_photon_frac", "jet photon fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14, 100, 0, 1);
  TH2F *jet1_photon_frac = new TH2F("jet1_photon_frac", "jet 1 photon fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14, 100, 0, 1);
  TH2F *jet2_photon_frac = new TH2F("jet2_photon_frac", "jet 2 photon fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14, 100, 0, 1);
  
  TH2F *photon1_ptVsID = new TH2F("photon1_ptVsID", "Leading photon pt vs. ID", 2, -0.5, 1.5, 100, 0, 1000);
  
  TProfile *photon_frac = new TProfile("photon_frac", "jet photon fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *photon_frac_0p3 = new TProfile("photon_frac_0p3", "jet photon fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
  
  TProfile *ch_had_frac = new TProfile("ch_had_frac", "jet charged hadronic fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *ch_had_frac_0p3 = new TProfile("ch_had_frac_0p3", "jet charged hadronic fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
  
  TProfile *ch_EM_frac = new TProfile("ch_EM_frac", "jet charged EM fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *ch_EM_frac_0p3 = new TProfile("ch_EM_frac_0p3", "jet charged EM fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
  
  TProfile *ch_Mu_frac = new TProfile("ch_Mu_frac", "jet charged muon fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *ch_Mu_frac_0p3 = new TProfile("ch_Mu_frac_0p3", "jet charged muon fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
  
  TProfile *ne_Had_frac = new TProfile("ne_Had_frac", "jet neutral hadronic fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *ne_Had_frac_0p3 = new TProfile("ne_Had_frac_0p3", "jet neutral hadronic fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
  
  TProfile *ne_EM_frac = new TProfile("ne_EM_frac", "jet neutral EM fraction vs. eta", 100, -3.14, 3.14);  
  TProfile *ne_EM_frac_0p3 = new TProfile("ne_EM_frac_0p3", "jet neutral EM fraction vs. eta (ChF < 0.3)", 100, -3.14, 3.14);
	
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t entries = corrjets.fChain->GetEntries();
// 		Int_t entries = 3000000;
  std::cout<<"Processing "<<entries<<"entries"<<std::endl;
//   double weight = 93.47*16146/entries;//GJets
  
  for(Int_t entry = 0; entry < entries; ++entry){
    if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
    SPVjets.GetEntry(entry);
    corrjets.GetEntry(entry);
    
    double njets = 0;
    for (int k = 0; k < 8; ++k){
      if (corrjets.jet_pt[k] > 30) njets++;
    }
  
//     double weight = corrjets.pswgt_dijet_170;
    double weight = 1;
    
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
    int arr[4] = {0,1,2,3};
    
    //bubbleSort
    for (int i = 0; i < 4; i++){     
       // Last i elements are already in place   
      for (int j = 0; j < 3-i; j++){
        if (corrjets.photon_pt[arr[j]] < corrjets.photon_pt[arr[j+1]]) swap(&arr[j], &arr[j+1]);
      }
    }
    
//     std::cout<<corrjets.photon_pt[0]<<" "<<corrjets.photon_pt[1]<<" "<<corrjets.photon_pt[2]<<" "<<corrjets.photon_pt[3]<<std::endl;
//     std::cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<std::endl;
    
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
      HoverE[i] = (corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ne_Had[i])/(corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ne_EM[i]);
      EMF[i] = corrjets.jet_efrac_ch_EM[i] + corrjets.jet_efrac_ne_EM[i];
    } 
      
    bool pass_conv_1 = true;
    bool pass_conv_2 = true;
    if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
    if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
    
    output->cd();
    
    if (corrjets.HLT_PFJet450 == 1 && corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1]> 550 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && (corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9) && corrjets.vtx_N >= 2 && njets == 2 && pass_conv_1 && pass_conv_2){
      if (corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){
      nJets->Fill(njets, weight);
      nJetsVsNEMF1->Fill(njets, corrjets.jet_efrac_photon[0], weight);
      nJetsVsNEMF2->Fill(njets, corrjets.jet_efrac_photon[1], weight);
      nvtx->Fill(corrjets.vtx_N, weight);
      HT->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3], weight);
      HT_nowgt->Fill(corrjets.jet_pt[0]+corrjets.jet_pt[1]+corrjets.jet_pt[2]+corrjets.jet_pt[3]);
      
      jet1_pt->Fill(corrjets.jet_pt[0], weight);
      jet1_pt_SPV->Fill(SPVjets.jet_pt[0], weight);
      jet2_pt->Fill(corrjets.jet_pt[1], weight);
      
      jet1_eta->Fill(corrjets.jet_eta[0], weight);
      jet1_eta_SPV->Fill(SPVjets.jet_eta[0], weight);
      jet2_eta->Fill(corrjets.jet_eta[1], weight);
      jet1_phi->Fill(corrjets.jet_phi[0], weight);
      jet2_phi->Fill(corrjets.jet_phi[1], weight);
      deltaphi->Fill(deltajet_phi, weight);
      if (deltajet_phi < 2.5) std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<corrjets.nEvent<<std::endl;
      if (CHEF_corrjet[0] < 0.5 && CHEF_corrjet[1] < 0.5 && CHEF_SPVjet[0] < 0.5 && CHEF_SPVjet[1] < 0.5) deltaphi_ChFMax0p5->Fill(deltajet_phi, weight);
      if (CHEF_corrjet[1] > 0.25) deltaphiVsChF1->Fill(deltajet_phi, CHEF_corrjet[0], weight);
      if (CHEF_corrjet[0] > 0.25) deltaphiVsChF2->Fill(deltajet_phi, CHEF_corrjet[1], weight);
      
      jet1_HOverE->Fill(HoverE[0], weight);
      jet2_HOverE->Fill(HoverE[1], weight);
      
      if (CHEF_corrjet[0] >= 0.0){
        jet1_chf->Fill(CHEF_corrjet[0], weight);
        jet1_chf_SPV->Fill(CHEF_SPVjet[0], weight);
        ChfOverSPVChF->Fill(CHEF_corrjet[0]/CHEF_SPVjet[0], weight);
        CHFvsPT->Fill(corrjets.jet_pt[0], CHEF_corrjet[0], weight);
        ChFvsNvtx1->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
      }
      
      if (CHEF_corrjet[1] >= 0.0){
        jet2_chf->Fill(CHEF_corrjet[1], weight);
        if (CHEF_corrjet[1] > 0.5){
          jet1_chf_jet2_0p5->Fill(CHEF_corrjet[0], weight);
//           if (CHEF_corrjet[0] < 0.01){
//             if (corrjets.nEvent < 0){  
//               UInt_t event = corrjets.nEvent;
//               std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"/*<<corrjets.nEvent<<":"*/<<event<<std::endl; 
//             }else std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<corrjets.nEvent<<std::endl;
//           }
        }
        jet2_chf_SPV->Fill(CHEF_SPVjet[1], weight);
        ChfOverSPVChF->Fill(CHEF_corrjet[1]/CHEF_SPVjet[1], weight);
        CHFvsPT->Fill(corrjets.jet_pt[1], CHEF_corrjet[1], weight);
        ChFvsNvtx2->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);
      }
      
      for(int j = 0; j < 24; j++){
        if (CHEF_SPVjet[0]<chf_cuts[j+1] && CHEF_SPVjet[1]<chf_cuts[j+1] && CHEF_corrjet[0]<chf_cuts[j+1] && CHEF_corrjet[1]<chf_cuts[j+1] && (CHEF_corrjet[1]>=chf_cuts[j] || CHEF_corrjet[0]>=chf_cuts[j]) && chf_cuts[j] >= 0.05){
          if (corrjets.Flag_HBHENoiseFilter == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_HBHENoise->Fill(chf_bins[j], weight);
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_HBHENoise->Fill(chf_bins[j], weight);            
          }
          if (corrjets.Flag_HBHENoiseIsoFilter == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_HBHENoiseIso->Fill(chf_bins[j], weight);            
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_HBHENoiseIso->Fill(chf_bins[j], weight);             
          }
          if (corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_EcalDeadCellTriggerPrimitiveFilter->Fill(chf_bins[j], weight);
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_EcalDeadCellTriggerPrimitiveFilter->Fill(chf_bins[j], weight);            
          }
          if (corrjets.Flag_goodVertices == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_goodVertices->Fill(chf_bins[j], weight);
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_goodVertices->Fill(chf_bins[j], weight);            
          }
          if (corrjets.Flag_eeBadScFilter == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_eeBadScFilter->Fill(chf_bins[j], weight);
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_eeBadScFilter->Fill(chf_bins[j], weight);            
          }
          if (corrjets.Flag_globalTightHalo2016Filter == 1){
            jet_chf_passed->Fill(chf_bins[j], weight);
            jet_chf_passed_globalTightHalo2016Filter->Fill(chf_bins[j], weight);
          }else{
            jet_chf_failed->Fill(chf_bins[j], weight);
            jet_chf_failed_globalTightHalo2016Filter->Fill(chf_bins[j], weight);            
          }
        }
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
      
      if (CHEF_corrjet[0] >= 0.2 && CHEF_corrjet[1] >= 0.2){
        if (dRjet1 < 0.4){
  //           ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[0], weight);
          if (CHEF_corrjet[0]>CHEF_corrjet[1]) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[0], weight);
          if (CHEF_corrjet[1]>CHEF_corrjet[0]) ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[1], weight);
        }else if(dRjet2 < 0.4){
  //           ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[1], weight);
          if (CHEF_corrjet[0]>CHEF_corrjet[1]) ChfVsSPVChf->Fill(CHEF_corrjet[0], CHEF_SPVjet[1], weight);
          if (CHEF_corrjet[1]>CHEF_corrjet[0]) ChfVsSPVChf->Fill(CHEF_corrjet[1], CHEF_SPVjet[0], weight);   
        }
        CHFvsCHF->Fill(CHEF_corrjet[1], CHEF_corrjet[0], weight);
        ChFvsNvtx->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
        ChFvsNvtx->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);				
      }
      
//       eta1vsnpixhits->Fill(corrjets.track_nPixHits[0], corrjets.jet_eta[0], weight);
//       eta2vsnpixhits->Fill(corrjets.track_nPixHits[0], corrjets.jet_eta[1], weight);
//       
//       NVTXvsCHF->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
//       NVTXvsCHF->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);
//       NVTXvsCHF1->Fill(CHEF_corrjet[0], corrjets.vtx_N, weight);
//       NVTXvsCHF2->Fill(CHEF_corrjet[1], corrjets.vtx_N, weight);
//       
//       Chi2vsCHF->Fill(CHEF_corrjet[0], corrjets.track_normalizedChi2[0], weight);
//       Chi2vsCHF->Fill(CHEF_corrjet[1], corrjets.track_normalizedChi2[0], weight);
//       Chi2vsCHF1->Fill(CHEF_corrjet[0], corrjets.track_normalizedChi2[0], weight);
//       Chi2vsCHF2->Fill(CHEF_corrjet[1], corrjets.track_normalizedChi2[0], weight);
//       
//       ndofvsCHF->Fill(CHEF_corrjet[0], corrjets.track_ndof[0], weight);
//       ndofvsCHF->Fill(CHEF_corrjet[1], corrjets.track_ndof[0], weight);
//       ndofvsCHF1->Fill(CHEF_corrjet[0], corrjets.track_ndof[0], weight);
//       ndofvsCHF2->Fill(CHEF_corrjet[1], corrjets.track_ndof[0], weight);
      
      track1_nhits->Fill(corrjets.track_nhits[0], weight);
      track1_npixhits->Fill(corrjets.track_nPixHits[0], weight);
      track1_pt->Fill(corrjets.track_pt[0], weight);
      track1_ptError->Fill(corrjets.track_ptError[0], weight);
      track1_dzError->Fill(corrjets.track_dzError[0], weight);
      track1_dz->Fill(corrjets.track_dz[0], weight);
      track1_dzRelError->Fill(corrjets.track_dzError[0]/corrjets.track_dz[0], weight);
      track1_reldzErrorvsdz->Fill(corrjets.track_dz[0], corrjets.track_dzError[0]/corrjets.track_dz[0], weight);
      
      ntracks->Fill(corrjets.nTrack, weight);
      int totnhits = 0;
      for (int i = 0; i < 10; ++i){
        nhitspertrack->Fill(corrjets.track_nhits[i], weight);
        totnhits += corrjets.track_nhits[i];
      }
      nhits->Fill(totnhits, weight);
      
      Chi2vsNdof->Fill(corrjets.track_ndof[0], corrjets.track_normalizedChi2[0], weight);
      
//       if(corrjets.photon_pt[photon_nr] > 10 && corrjets.photon_passLooseId[photon_nr] == 1) 
      photon1_pt->Fill(corrjets.photon_pt[photon_nr], weight);
//       if(corrjets.photon_pt[photon_nr] > 10 && corrjets.photon_passLooseId[photon_nr] == 1) 
      photon1_eta->Fill(corrjets.photon_eta[photon_nr], weight);
//       if(corrjets.photon_pt[photon_nr] > 10 && corrjets.photon_passLooseId[photon_nr] == 1) 
      photon1_phi->Fill(corrjets.photon_phi[photon_nr], weight);
      photon1_ID->Fill(corrjets.photon_passLooseId[photon_nr], weight);
      photon1_dR1->Fill(dR1, weight);
      photon1_dR2->Fill(dR2, weight);
      photon1_ptVsID->Fill(corrjets.photon_passLooseId[photon_nr], corrjets.photon_pt[photon_nr], weight);
      
      photon2_pt->Fill(corrjets.photon_pt[arr[1]], weight);
      photon2_eta->Fill(corrjets.photon_eta[arr[1]], weight);
      photon2_phi->Fill(corrjets.photon_phi[arr[1]], weight);
      photon2_ID->Fill(corrjets.photon_passLooseId[arr[1]], weight);
//       photon2_dR1->Fill(dR1, weight);
//       photon2_dR2->Fill(dR2, weight);
      
      photon_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_photon[0], weight);
      photon_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_photon[1], weight);
      
      ch_had_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_Had[0], weight);
      ch_had_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_Had[1], weight);
      
      ch_EM_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_EM[0], weight);
      ch_EM_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_EM[1], weight);
      
      ch_Mu_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_Mu[0], weight);
      ch_Mu_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_Mu[1], weight);
      
      ne_Had_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ne_Had[0], weight);
      ne_Had_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ne_Had[1], weight);
      
      ne_EM_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ne_EM[0], weight);
      ne_EM_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ne_EM[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) photon_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_photon[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) photon_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_photon[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) ch_had_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_Had[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) ch_had_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_Had[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) ch_EM_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_EM[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) ch_EM_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_EM[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) ch_Mu_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ch_Mu[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) ch_Mu_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ch_Mu[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) ne_Had_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ne_Had[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) ne_Had_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ne_Had[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) ne_EM_frac_0p3->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_ne_EM[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) ne_EM_frac_0p3->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_ne_EM[1], weight);
      
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) jet_photon_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_photon[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) jet_photon_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_photon[1], weight);
      if (CHEF_SPVjet[0]<0.3 && CHEF_corrjet[0]<0.3) jet1_photon_frac->Fill(corrjets.jet_eta[0], corrjets.jet_efrac_photon[0], weight);
      if (CHEF_SPVjet[1]<0.3 && CHEF_corrjet[1]<0.3) jet2_photon_frac->Fill(corrjets.jet_eta[1], corrjets.jet_efrac_photon[1], weight);
    }else counter++;
	}
  }
  output->Write();
  output->Close();
  std::cout<<"counter: "<<counter<<std::endl;
}