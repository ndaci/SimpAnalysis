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

// #include "runbunchmap_50bx.C"
#include "SANtuple.h"
// #include "checkLumi.cpp"
// #include "json_parser.cpp"

// #include "lists/list_JetHT_2016B_local.h"
// #include "lists/list_JetHT_2016C_local.h"
// #include "lists/list_JetHT_2016D_local.h"
// #include "lists/list_JetHT_2016E_local.h"
// #include "lists/list_JetHT_2016F_local.h"

// #include "lists/list_JetHT_2016B_new.h"
// #include "lists/list_JetHT_2016C_new.h"
// #include "lists/list_JetHT_2016D_new.h"
// #include "lists/list_JetHT_2016E_new.h"
// #include "lists/list_JetHT_2016F_new.h"
// #include "lists/list_JetHT_2016G_new.h"
// #include "lists/list_JetHT_2016H2_new.h"
// #include "lists/list_JetHT_2016H3_new.h"

#include "lists/list_JetHT_2016B_newtrigger_2.h"
#include "lists/list_JetHT_2016C_newtrigger_2.h"
#include "lists/list_JetHT_2016D_newtrigger_2.h"
#include "lists/list_JetHT_2016E_newtrigger_2.h"
#include "lists/list_JetHT_2016F_newtrigger_2.h"
#include "lists/list_JetHT_2016G_newtrigger_2.h"
#include "lists/list_JetHT_2016H2_newtrigger_2.h"
#include "lists/list_JetHT_2016H3_newtrigger_2.h"

// #include "lists/list_JetHT_2016B_allConversions.h"
// #include "lists/list_JetHT_2016C_allConversions.h"
// #include "lists/list_JetHT_2016D_allConversions.h"
// // #include "lists/list_JetHT_2016E_allConversions.h"
// #include "lists/list_JetHT_2016F_allConversions.h"
// #include "lists/list_JetHT_2016G_allConversions.h"
// #include "lists/list_JetHT_2016H2_allConversions.h"
// #include "lists/list_JetHT_2016H3_allConversions.h"

// #include "lists/list_JetHT_2016B_bx.h"
// #include "lists/list_JetHT_2016C_bx.h"
// #include "lists/list_JetHT_2016D_bx.h"
// #include "lists/list_JetHT_2016E_bx.h"
// #include "lists/list_JetHT_2016F_bx.h"
// #include "lists/list_JetHT_2016G_bx.h"
// #include "lists/list_JetHT_2016H2_bx.h"
// #include "lists/list_JetHT_2016H3_bx.h"

void SIMP_Data_closure_exclusive_optbinning(){
	  
  TChain* chain = new TChain("treeCorr/SimpAnalysis");
	list_JetHT_2016B(chain);
	list_JetHT_2016C(chain);
	list_JetHT_2016D(chain);
	list_JetHT_2016E(chain);
	list_JetHT_2016F(chain);
// 	list_JetHT_2016G(chain);
// 	list_JetHT_2016H2(chain);
// 	list_JetHT_2016H3(chain);
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
	list_JetHT_2016B(SPVchain);
	list_JetHT_2016C(SPVchain);
	list_JetHT_2016D(SPVchain);
	list_JetHT_2016E(SPVchain);
	list_JetHT_2016F(SPVchain);
// 	list_JetHT_2016G(SPVchain);
// 	list_JetHT_2016H2(SPVchain);
// 	list_JetHT_2016H3(SPVchain);
  
//   fillrunbunchmap();
  
  double CHEF_SPVjet[8], CHEF_corrjet[8];
  int photon_nr;
	
// 	double chf_cuts[23] = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	double chf_cuts[25] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0};
	double chf_bins[24] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.06, 0.085, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.95};
	double zero[24] = {0.005,0.005,0.005,0.005,0.005,0.01,0.015,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025, 0.025, 0.05};
// 	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
// 	double eta_bins[17] = {-2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};

//   pt::ptree root;
//   pt::read_json("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt", root);

  double eff1_low = 0;
  double eff2_low = 0;
  double erreff1_low = 0;
  double erreff2_low = 0;
  double eff1_high = 0;
  double eff2_high = 0;
  double erreff1_high = 0;
  double erreff2_high = 0;
  int N_1[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int N_2[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int N_3[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int N_4[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
	double passed_eff[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff1[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff2[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff3[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_eff4[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_effboth[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_MCtruth[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_MCtruth[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff1[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff2[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff3[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff4[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_eff[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_effboth[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	double ratio_1leg[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double ratio_2leg[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_1leg[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double err_ratio_2leg[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
// 	TFile* closure = new TFile("closure_test_MCtruthSPVcut_PUMoriond17_AOD_conversionCut_lumiGH_correct.root", "READ");
//   TH1D* histo_closure = (TH1D*) closure->Get("histo_2leg");
  
	std::cout<<"Getting the efficiency histos...";
// 	TFile* efficiencies = new TFile("eff2D_Data_RerecoGH_AOD_conversionVeto.root", "READ");
// 	TFile* efficiencies = new TFile("eff2D_Data_RerecoAF_AOD_withNvertexCut_Njets2_finerEtaBinning.root", "READ");
// 	TFile* efficiencies = new TFile("eff2D_Data_AOD_RunBToF_first50bx.root", "READ");
	TFile* efficiencies = new TFile("eff2D_Data_AOD_RunBF_filters.root", "READ");
// 	TFile* efficiencies = new TFile("eff2D_Data_AOD_RunGH_6ptbins8etabins.root", "READ");
	TH2D* eff_histos[23];
// 	double closure_unc[24];
  for(int j = 0; j < 23; j++){
    std::ostringstream strs;
    double dbl = chf_cuts[j+1];
		strs << dbl;
    std::string cut = strs.str();
    std::string title_eff = "eff_"+cut;
    eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
//     closure_unc[j] = histo_closure->GetBinContent(histo_closure->GetXaxis()->FindBin(chf_cuts[j]));
  }
	TFile* exclusive_efficiencies = new TFile("eff2D_Data_AOD_RunBF_exclusive_filters.root", "READ");
// 	TFile* exclusive_efficiencies = new TFile("eff2D_Data_AOD_RunGH_exclusive_6ptbins8etabins.root", "READ");
	TH2D* excl_eff_histos[24];
// 	double closure_unc[24];
  for(int j = 0; j < 24; j++){
    std::ostringstream strs;
    double dbl = chf_cuts[j];
		double dbl2 = chf_cuts[j+1];
		strs << dbl << "To" << dbl2;
    std::string cut = strs.str();
    std::string title_eff = "eff_"+cut;
    excl_eff_histos[j] = (TH2D*) exclusive_efficiencies->Get(title_eff.c_str());
//     closure_unc[j] = histo_closure->GetBinContent(histo_closure->GetXaxis()->FindBin(chf_cuts[j]));
  }
	std::cout<<"done"<<std::endl;
	
//   TFile *output = new TFile("closure_test_Data_RunGH_signaltrigger.root", "RECREATE");
  TFile *output = new TFile("closure_test_Data_RunBF_exclusive_unblinded.root", "RECREATE");
	  
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t Nentries = corrjets.fChain->GetEntries(); 
  std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
//   Nentries = 100000; 
  
  int counter = 0;
  
  for(Int_t entry = 0; entry < Nentries; ++entry){
    if(entry%1000000==0) std::cout<<"processed "<<entry/1000000<<"M events"<<std::endl;
    SPVjets.GetEntry(entry);
    corrjets.GetEntry(entry);
    
//     auto runit = runbunchmap.find(corrjets.nRun);
//     if(runit == runbunchmap.end()) continue;
//     
//     std::vector<int> &bxvec = runbunchmap[corrjets.nRun];
//     std::vector<int>::iterator it;
//     it = std::find( bxvec.begin(), bxvec.end(), corrjets.nBx);
//     if(it == bxvec.end()) continue;
//     
//     counter++;
    
//     if (corrjets.nRun > 278820) continue;
//     if (corrjets.nRun < 279479) continue;    
    
//     bool isInLowLumi = isLowLumi(corrjets.nRun);
//     if (!isInLowLumi) continue;
    
//     bool isInHighLumi = isHighLumi(corrjets.nRun);
//     if (!isInHighLumi) continue;
    
//     bool isInJSON = Check_Run_Lumi(corrjets.nRun, corrjets.nLumi, root);
//     if (!isInJSON) continue;
    
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
    
    for (int i = 0; i < 8; i++){
      CHEF_SPVjet[i] = SPVjets.jet_efrac_ch_Had[i]+SPVjets.jet_efrac_ch_EM[i]+SPVjets.jet_efrac_ch_Mu[i];
      CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
    } 
      
    bool pass_conv_1 = true;
    bool pass_conv_2 = true;
    if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
    if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
    
    output->cd();
    
    if (corrjets.HLT_PFJet450 == 1 && corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2 && ((corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9)) && corrjets.vtx_N >= 2 && njets == 2 && pass_conv_1 && pass_conv_2 && corrjets.Flag_HBHENoiseFilter == 1 && corrjets.Flag_HBHENoiseIsoFilter == 1 && corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){
      for(int j = 0; j < 24; j++){
//         if(corrjets.HLT_DiCentralPFJet170_CFMax0p1 == 1 || corrjets.HLT_DiCentralPFJet220_CFMax0p3 == 1 || corrjets.HLT_DiCentralPFJet330_CFMax0p5 == 1 || corrjets.HLT_DiCentralPFJet430 == 1){
          if (CHEF_SPVjet[0]<chf_cuts[j+1] && CHEF_SPVjet[1]<chf_cuts[j+1] && CHEF_corrjet[0]<chf_cuts[j+1] && CHEF_corrjet[1]<chf_cuts[j+1] && (CHEF_corrjet[1]>=chf_cuts[j] || CHEF_corrjet[0]>=chf_cuts[j])){
            /*if (chf_cuts[j] >= 0.05)*/ passed_MCtruth[j]+= /*corrjets.pswgt_dijet_170*/1;
            /*if (chf_cuts[j] >= 0.05)*/ err_MCtruth[j] += pow(/*corrjets.pswgt_dijet_170*/1, 2);
//             if (chf_cuts[j] >= 0.2) passed_MCtruth[j]++;
//             if (chf_cuts[j] >= 0.2) err_MCtruth[j] ++;
          }
//         }
//         if (corrjets.HLT_DiCentralPFJet170 == 1){
          if (j == 0) eff1_low = 0;
          else eff1_low = eff_histos[j-1]->GetBinContent(eff_histos[j-1]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j-1]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          if (j == 23) eff1_high = 1;
          else eff1_high = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          if (j == 0) erreff1_low = 0;
          else erreff1_low = eff_histos[j-1]->GetBinError(eff_histos[j-1]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j-1]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          if (j == 23) erreff1_high = 0;
          else erreff1_high = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          
          if (j == 0) eff2_low = 0;
          else eff2_low = eff_histos[j-1]->GetBinContent(eff_histos[j-1]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j-1]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          if (j == 23) eff2_high = 1;
          else eff2_high = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          if (j == 0) erreff2_low = 0;
          else erreff2_low = eff_histos[j-1]->GetBinError(eff_histos[j-1]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j-1]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          if (j == 23) erreff2_high = 0;
          else erreff2_high = eff_histos[j]->GetBinError(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          
          double excl_eff1 = excl_eff_histos[j]->GetBinContent(excl_eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), excl_eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          double excl_eff2 = excl_eff_histos[j]->GetBinContent(excl_eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), excl_eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
          double excl_erreff1 = excl_eff_histos[j]->GetBinError(excl_eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), excl_eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
          double excl_erreff2 = excl_eff_histos[j]->GetBinError(excl_eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), excl_eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));

          if (CHEF_corrjet[0]<chf_cuts[j+1] && CHEF_corrjet[0]>=chf_cuts[j]){
            passed_eff1[j]+=/*corrjets.pswgt_dijet_170*/1*eff2_high;
            N_1[j]++;
            err_eff1[j] += /*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*erreff2_high*erreff2_high;
          }
          if (CHEF_corrjet[0]<chf_cuts[j]){
            passed_eff2[j]+=/*corrjets.pswgt_dijet_170*/1*excl_eff2;
            N_2[j]++;
            err_eff2[j] += /*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*excl_erreff2*excl_erreff2;
          }
          if (CHEF_corrjet[1]<chf_cuts[j+1] && CHEF_corrjet[1]>=chf_cuts[j]){
            passed_eff3[j]+=/*corrjets.pswgt_dijet_170*/1*eff1_high;
            N_3[j]++;
            err_eff3[j] += /*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*erreff1_high*erreff1_high;
          }
          if (CHEF_corrjet[1]<chf_cuts[j]){
            passed_eff4[j]+=/*corrjets.pswgt_dijet_170*/1*excl_eff1;
            N_4[j]++;
            err_eff4[j] += /*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*excl_erreff1*excl_erreff1;
          }
          passed_effboth[j]+=((excl_eff1*eff2_high)+(eff1_low*excl_eff2))*/*corrjets.pswgt_dijet_170*/1;
          
          err_effboth[j] += (eff2_high*eff2_high*/*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*excl_erreff1*excl_erreff1) + (excl_eff1*excl_eff1*/*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*erreff2_high*erreff2_high) +  (excl_eff2*excl_eff2*/*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*erreff1_low*erreff1_low) + (eff1_low*eff1_low*/*corrjets.pswgt_dijet_170*/1*/*corrjets.pswgt_dijet_170*/1*excl_erreff2*excl_erreff2);
//         }
      }
    }
	}
	for(int j = 0; j < 24; j++){
//     std::cout<<chf_cuts[j]<<" "<<passed_effboth[j]<<" +- "<<TMath::Sqrt(err_effboth[j])<<" (stat) +- "<<closure_unc[j]<<" (syst)"<<std::endl;
		err_MCtruth[j] = TMath::Sqrt(err_MCtruth[j]);
//     err_effboth[j] += closure_unc[j]*closure_unc[j];
//     err_eff1[j] += closure_unc[j]*closure_unc[j];
//     err_eff2[j] += closure_unc[j]*closure_unc[j];
		err_effboth[j] = TMath::Sqrt(err_effboth[j]);
		err_eff1[j]+= pow(passed_eff1[j], 2)/N_1[j];
		err_eff2[j]+= pow(passed_eff2[j], 2)/N_2[j];
		err_eff3[j]+= pow(passed_eff3[j], 2)/N_3[j];
		err_eff4[j]+= pow(passed_eff4[j], 2)/N_4[j];
		err_eff[j] = TMath::Sqrt(err_eff1[j]+err_eff2[j]+err_eff3[j]+err_eff4[j])/2.0 ;
		passed_eff[j] = (passed_eff1[j]+passed_eff2[j]+passed_eff3[j]+passed_eff4[j])/2;
// 		ratio_1leg[j] = passed_MCtruth[j]/passed_eff[j];
		ratio_2leg[j] = passed_MCtruth[j]/passed_effboth[j];
// 		if (ratio_1leg[j] != 0) err_ratio_1leg[j] = ratio_1leg[j]*TMath::Sqrt((err_eff[j]*err_eff[j]/passed_eff[j]/passed_eff[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
		if (ratio_2leg[j] != 0) err_ratio_2leg[j] = ratio_2leg[j]*TMath::Sqrt((err_effboth[j]*err_effboth[j]/passed_effboth[j]/passed_effboth[j])+(err_MCtruth[j]*err_MCtruth[j]/passed_MCtruth[j]/passed_MCtruth[j]));
	
	}
	
	std::cout<<"# events in considered bx: "<<counter<<std::endl;
	
	TCanvas *c1 = new TCanvas("Closure test", "Closure test",700,500);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy(1);
  pad1->Draw();
  pad1->cd();
// 	TGraphErrors *MC = new TGraphErrors(24, chf_bins, passed_MCtruth, zero, err_MCtruth);  
	TH1D *MC = new TH1D("data", "data", 24, chf_cuts);
  for (int j = 0; j < 24; j++){
    MC->Fill(chf_bins[j], passed_MCtruth[j]/**0.025/zero[j]*/);
    MC->SetBinError(MC->GetXaxis()->FindBin(chf_bins[j]), err_MCtruth[j]/**0.025/zero[j]*/);
  }
	MC->SetTitle("");
  MC->SetStats(0);
	MC->GetXaxis()->SetTitle("CHF cut");
	MC->GetYaxis()->SetTitle("# events");
	MC->GetYaxis()->SetTitleSize(0.06);
	MC->GetYaxis()->SetLabelSize(0.06);
	MC->GetYaxis()->SetTitleOffset(0.8);
	MC->SetMarkerStyle(20);
	MC->Draw("P");
	
// 	TGraphErrors *oneleg = new TGraphErrors(24, chf_bins, passed_eff, zero, err_eff);
  TH1D *oneleg = new TH1D("1leg", "1-leg prediction", 24, chf_cuts);
  for (int j = 0; j < 24; j++){
    oneleg->Fill(chf_bins[j], passed_eff[j]/**0.025/zero[j]*/);
    oneleg->SetBinError(oneleg->GetXaxis()->FindBin(chf_bins[j]), err_eff[j]/**0.025/zero[j]*/);
  } 
	oneleg->SetMarkerStyle(20);
	oneleg->SetMarkerColor(2);
	oneleg->Draw("P SAME");
  
// 	TGraphErrors *twoleg = new TGraphErrors(24, chf_bins, passed_effboth, zero, err_effboth);
  TH1D *twoleg = new TH1D("2leg", "2-leg prediction", 24, chf_cuts);
  for (int j = 0; j < 24; j++){
    twoleg->Fill(chf_bins[j], passed_effboth[j]/**0.025/zero[j]*/);
    twoleg->SetBinError(twoleg->GetXaxis()->FindBin(chf_bins[j]), err_effboth[j]/**0.025/zero[j]*/);
  }
	twoleg->SetMarkerStyle(20);
	twoleg->SetMarkerColor(3);
	twoleg->Draw("P SAME");
	
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
//   pad2->SetLogy();
  pad2->Draw();
  pad2->cd();
	
// 	TGraphErrors *ratio_oneleg = new TGraphErrors(24, chf_cuts, ratio_1leg, zero, err_ratio_1leg);
// 	ratio_oneleg->SetMarkerStyle(20);
// 	ratio_oneleg->SetMarkerColor(2);
// 	ratio_oneleg->Draw("AP");
//   ratio_oneleg->SetTitle("");
//   ratio_oneleg->GetYaxis()->SetRangeUser(0.8, 1.2);
//   ratio_oneleg->GetYaxis()->SetTitle("#frac{data}{1/2-leg}");
//   ratio_oneleg->GetYaxis()->CenterTitle(kTRUE);
// 	ratio_oneleg->GetYaxis()->SetTitleSize(0.11);
// 	ratio_oneleg->GetYaxis()->SetLabelSize(0.1);
// 	ratio_oneleg->GetYaxis()->SetTitleOffset(0.4);
//   ratio_oneleg->GetXaxis()->SetTitle("ChF cut");
// 	ratio_oneleg->GetXaxis()->SetTitleSize(0.12);
// 	ratio_oneleg->GetXaxis()->SetLabelSize(0.1);
  
  TH1D *ratio_oneleg = new TH1D("1leg_ratio", "1leg_ratio", 24, chf_cuts);
  ratio_oneleg->Divide(MC,oneleg);
  ratio_oneleg->SetMarkerColor(2);
	ratio_oneleg->SetMarkerStyle(20);
  ratio_oneleg->SetStats(0);
  ratio_oneleg->SetTitle("");
  ratio_oneleg->GetYaxis()->SetRangeUser(0.55, 1.45);
  ratio_oneleg->GetYaxis()->SetTitle("#frac{data}{1/2-leg}");
  ratio_oneleg->GetYaxis()->CenterTitle(kTRUE);
	ratio_oneleg->GetYaxis()->SetTitleSize(0.11);
	ratio_oneleg->GetYaxis()->SetLabelSize(0.1);
	ratio_oneleg->GetYaxis()->SetTitleOffset(0.4);
  ratio_oneleg->GetXaxis()->SetTitle("ChF bin");
	ratio_oneleg->GetXaxis()->SetTitleSize(0.12);
	ratio_oneleg->GetXaxis()->SetLabelSize(0.1);
	ratio_oneleg->Draw("P");
  
  TH1D *ratio_twoleg = new TH1D("2leg_ratio", "2leg_ratio", 24, chf_cuts);
  ratio_twoleg->Divide(MC,twoleg);
	ratio_twoleg->SetMarkerStyle(20);
  ratio_twoleg->SetMarkerColor(3);
	ratio_twoleg->Draw("P SAME");
  
  
// 	TGraphErrors *ratio_twoleg = new TGraphErrors(24, chf_bins, ratio_2leg, zero, err_ratio_2leg);
// 	ratio_twoleg->SetMarkerStyle(20);
// 	ratio_twoleg->SetMarkerColor(3);
// // 	ratio_twoleg->Draw("P");
// 	ratio_twoleg->Draw("AP");
//   ratio_twoleg->SetTitle("");
//   ratio_twoleg->GetYaxis()->SetRangeUser(0.8, 1.2);
//   ratio_twoleg->GetYaxis()->SetTitle("#frac{data}{1/2-leg}");
//   ratio_twoleg->GetYaxis()->CenterTitle(kTRUE);
// 	ratio_twoleg->GetYaxis()->SetTitleSize(0.11);
// 	ratio_twoleg->GetYaxis()->SetLabelSize(0.1);
// 	ratio_twoleg->GetYaxis()->SetTitleOffset(0.4);
//   ratio_twoleg->GetXaxis()->SetTitle("ChF cut");
// 	ratio_twoleg->GetXaxis()->SetTitleSize(0.12);
// 	ratio_twoleg->GetXaxis()->SetLabelSize(0.1);
  
  TLine *line = new TLine(0,1,1,1);
  line->SetLineColor(kViolet+5);
  line->SetLineWidth(2);
	line->Draw();
	
  c1->cd();
	output->Append(c1);
  output->Write();
  output->Close();
}
