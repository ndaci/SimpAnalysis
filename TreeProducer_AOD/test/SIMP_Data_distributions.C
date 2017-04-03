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

void SIMP_Data_distributions(){
	  
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
	
  TFile *output = new TFile("Data_distributions.root", "RECREATE");
  
  TH1F *data_jet1_pt = new TH1F("data_jet1_pt", "Leading jet pt", 100, 0, 2000);
  TH1F *data_jet1_eta = new TH1F("data_jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
  TH1F *data_jet2_pt = new TH1F("data_jet2_pt", "Subleading jet pt", 100, 0, 2000);
 	TH1F *data_jet2_eta = new TH1F("data_jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *data_nvtx = new TH1F("data_nvtx", "Number of vertices", 51, -0.5, 50.5);
  
  TH1F *oneleg_jet1_pt = new TH1F("oneleg_jet1_pt", "Leading jet pt", 100, 0, 2000);
  TH1F *oneleg_jet1_eta = new TH1F("oneleg_jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
  TH1F *oneleg_jet2_pt = new TH1F("oneleg_jet2_pt", "Subleading jet pt", 100, 0, 2000);
 	TH1F *oneleg_jet2_eta = new TH1F("oneleg_jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *oneleg_nvtx = new TH1F("oneleg_nvtx", "Number of vertices", 51, -0.5, 50.5);
  
  TH1F *twoleg_jet1_pt = new TH1F("twoleg_jet1_pt", "Leading jet pt", 100, 0, 2000);
  TH1F *twoleg_jet1_eta = new TH1F("twoleg_jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
  TH1F *twoleg_jet2_pt = new TH1F("twoleg_jet2_pt", "Subleading jet pt", 100, 0, 2000);
 	TH1F *twoleg_jet2_eta = new TH1F("twoleg_jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *twoleg_nvtx = new TH1F("twoleg_nvtx", "Number of vertices", 51, -0.5, 50.5);
  
  data_jet1_pt->Sumw2();
  data_jet1_eta->Sumw2();
  data_jet2_pt->Sumw2();
  data_jet2_eta->Sumw2();
  data_nvtx->Sumw2();
  
  oneleg_jet1_pt->Sumw2();
  oneleg_jet1_eta->Sumw2();
  oneleg_jet2_pt->Sumw2();
  oneleg_jet2_eta->Sumw2();
  oneleg_nvtx->Sumw2();
  
  twoleg_jet1_pt->Sumw2();
  twoleg_jet1_eta->Sumw2();
  twoleg_jet2_pt->Sumw2();
  twoleg_jet2_eta->Sumw2();
  twoleg_nvtx->Sumw2();
  
	TFile* efficiencies = new TFile("eff2D_Data_Rereco_AOD.root", "READ");
	TH2D* eff_histo = (TH2D*) efficiencies->Get("eff_0.3");
	  
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
    
    double deltajet_phi = corrjets.jet_phi[0] - corrjets.jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
    deltajet_phi = fabs(deltajet_phi);
    
    double deltaphi_jet1photon =  corrjets.jet_phi[0] -  corrjets.photon_phi[0];
    if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
    if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
    double deltaphi_jet2photon =  corrjets.jet_phi[1] -  corrjets.photon_phi[0];
    if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
    if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
    
    double deltaeta_jet1photon =  corrjets.jet_eta[0] -  corrjets.photon_eta[0];
    double deltaeta_jet2photon =  corrjets.jet_eta[1] -  corrjets.photon_eta[0];
    
    double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
    double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
    
    for (int i = 0; i < 8; i++){
      CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
      CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
    } 
    
    output->cd();
    
    if (corrjets.HLT_DiCentralPFJet170 == 1 && corrjets.jet_pt[0] > 250 && corrjets.jet_pt[1] > 250 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && (corrjets.photon_passLooseId[0] == 0 || (corrjets.photon_passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1))){
      if (CHEF_SPVjet[0]<0.3 && CHEF_SPVjet[1]<0.3 && CHEF_corrjet[0]<0.3 && CHEF_corrjet[1]<0.3){
        data_jet1_pt->Fill(corrjets.jet_pt[0], corrjets.pswgt_dijet_170);
        data_jet1_eta->Fill(corrjets.jet_eta[0], corrjets.pswgt_dijet_170);
        data_jet2_pt->Fill(corrjets.jet_pt[1], corrjets.pswgt_dijet_170);
        data_jet2_eta->Fill(corrjets.jet_eta[1], corrjets.pswgt_dijet_170);
        data_nvtx->Fill(corrjets.vtx_N, corrjets.pswgt_dijet_170);
//         err_MCtruth[j] += pow(corrjets.pswgt_dijet_170, 2);
      }
      
      double eff1 = eff_histo->GetBinContent(eff_histo->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histo->GetYaxis()->FindBin(corrjets.jet_pt[0]));
      double eff2 = eff_histo->GetBinContent(eff_histo->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histo->GetYaxis()->FindBin(corrjets.jet_pt[1]));
      double erreff1 = eff_histo->GetBinError(eff_histo->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histo->GetYaxis()->FindBin(corrjets.jet_pt[0]));
      double erreff2 = eff_histo->GetBinError(eff_histo->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histo->GetYaxis()->FindBin(corrjets.jet_pt[1]));

      if (CHEF_corrjet[0]<0.3){
        oneleg_jet1_pt->Fill(corrjets.jet_pt[0], eff2*corrjets.pswgt_dijet_170);
        oneleg_jet1_eta->Fill(corrjets.jet_eta[0], eff2*corrjets.pswgt_dijet_170);
        oneleg_jet2_pt->Fill(corrjets.jet_pt[1], eff2*corrjets.pswgt_dijet_170);
        oneleg_jet2_eta->Fill(corrjets.jet_eta[1],eff2*corrjets.pswgt_dijet_170);
        oneleg_nvtx->Fill(corrjets.vtx_N, eff2*corrjets.pswgt_dijet_170);
      }
      if (CHEF_corrjet[1]<0.3){
        oneleg_jet1_pt->Fill(corrjets.jet_pt[0], eff1*corrjets.pswgt_dijet_170);
        oneleg_jet1_eta->Fill(corrjets.jet_eta[0], eff1*corrjets.pswgt_dijet_170);
        oneleg_jet2_pt->Fill(corrjets.jet_pt[1], eff1*corrjets.pswgt_dijet_170);
        oneleg_jet2_eta->Fill(corrjets.jet_eta[1],eff1*corrjets.pswgt_dijet_170);
        oneleg_nvtx->Fill(corrjets.vtx_N, eff1*corrjets.pswgt_dijet_170);
      }
      
      twoleg_jet1_pt->Fill(corrjets.jet_pt[0], eff1*eff2*corrjets.pswgt_dijet_170);
      twoleg_jet1_eta->Fill(corrjets.jet_eta[0], eff1*eff2*corrjets.pswgt_dijet_170);
      twoleg_jet2_pt->Fill(corrjets.jet_pt[1], eff1*eff2*corrjets.pswgt_dijet_170);
      twoleg_jet2_eta->Fill(corrjets.jet_eta[1], eff1*eff2*corrjets.pswgt_dijet_170);
      twoleg_nvtx->Fill(corrjets.vtx_N, eff1*eff2*corrjets.pswgt_dijet_170);
        
//       err_eff1[j] += corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff1*erreff1;
//       err_eff2[j] += corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff2*erreff2;
//       err_effboth[j] += (eff2*eff2*corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff1*erreff1) + (eff1*eff1*corrjets.pswgt_dijet_170*corrjets.pswgt_dijet_170*erreff2*erreff2);
    }
	}
  oneleg_jet1_pt->Scale(0.5);
  oneleg_jet1_eta->Scale(0.5);
  oneleg_jet2_pt->Scale(0.5);
  oneleg_jet2_eta->Scale(0.5);
  oneleg_nvtx->Scale(0.5);
  output->Write();
  output->Close();
}
