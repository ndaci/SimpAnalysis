#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

#include "SANtuple.h"

// #include "lists/list_JetHT_2016B_new.h"
// #include "lists/list_JetHT_2016C_new.h"
// #include "lists/list_JetHT_2016D_new.h"
// #include "lists/list_JetHT_2016E_new.h"
// #include "lists/list_JetHT_2016F_new.h"
// #include "lists/list_JetHT_2016G_new.h"
// #include "lists/list_JetHT_2016H2_new.h"
// #include "lists/list_JetHT_2016H3_new.h"

#include "lists/list_JetHT_2016G_newtrigger_2.h"
#include "lists/list_JetHT_2016H2_newtrigger_2.h"
#include "lists/list_JetHT_2016H3_newtrigger_2.h"

#include "lists/list_GJets.h"
#include "lists/list_WJets.h"

void SIMP_Data_cutflow(){
	  
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
// 	list_WJets(chain);
  
  TChain* SPVchain = new TChain("treeSPV/SimpAnalysis");
// 	list_JetHT_2016B(SPVchain);
// 	list_JetHT_2016C(SPVchain);
// 	list_JetHT_2016D(SPVchain);
// 	list_JetHT_2016E(SPVchain);
// 	list_JetHT_2016F(SPVchain);
	list_JetHT_2016G(SPVchain);
	list_JetHT_2016H2(SPVchain);
	list_JetHT_2016H3(SPVchain);
// 	list_GJets(SPVchain);
// 	list_WJets(SPVchain);
//
// 	TFile* output = new TFile("Data_cutflow_GH_controlTrigger_conversionVeto.root", "RECREATE");
	TFile* output = new TFile("Data_cutflow_GH_noHBHENoiseIsoFilter_unblinded.root", "RECREATE");
	
  double CHEF_SPVjet[8], CHEF_corrjet[8], EMF[8];
  int photon_nr;
  
  int HBHENoise = 0;
  int HBHENoiseIso = 0;
  int EcalDeadCell = 0;
  int goodVertices = 0;
  int eeBadScFilter = 0;
  int globalTightHalo2016 = 0;
	
  double passed_trigger = 0;
	double passed_deltaphi = 0;
	double passed_eta = 0;
	double passed_photonveto = 0;
	double passed_conversionveto = 0;
	double passed_pt = 0;
  double passed_njets = 0;
  double passed_nvtx = 0;
  double passed_filter = 0;
  double lowChF = 0;
  double highChF = 0;
	double passed[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_real[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double passed_SPV[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double chf_cuts[23] = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01};
	
  TH1D* trigger_eff = new TH1D("trigger_eff", "passed control trigger", 3, 0, 2);
  TH1D* pt_eff = new TH1D("pt_eff", "passed pt cuts", 100, 200, 300);
  TH1D* eta_eff = new TH1D("eta_eff", "passed eta cuts", 100, 0, 5);
  TH1D* deltaphi_eff = new TH1D("deltaphi_eff", "passed #Delta#phi cut", 100, 0, 5);
  TH1D* photonVeto_eff = new TH1D("photonVeto_eff", "passed photon veto", 10, 0, 1);
  TH1D* conversionVeto_eff = new TH1D("conversionVeto_eff", "passed conversion veto", 10, 0, 1);
  TH1D* njets_eff = new TH1D("njets_eff", "passed njet cut", 5, 0, 4);
  TH1D* nvtx_eff = new TH1D("nvtx_eff", "passed nvtx cut", 5, 0, 4);  
  TH1D* filter_eff = new TH1D("filter_eff", "passed all filters", 5, 0, 4);  
  TH1D* ChF_eff_1leg = new TH1D("ChF_eff_1leg", "passed ChF cuts (1-leg pred)", 100, 0, 0.5);
  TH1D* ChF_eff_2leg = new TH1D("ChF_eff_2leg", "passed ChF cuts (2-leg pred)", 100, 0, 0.5);
  TH1D* ChF_passed = new TH1D("ChF_passed", "passed ChF cuts", 100, 0, 0.5);
  
	TFile* efficiencies = new TFile("eff2D_Data_AOD_RunGH_filters.root", "READ");
	TH2D* eff_histos[23];
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::string title_eff = "eff_"+cut;
		eff_histos[j] = (TH2D*) efficiencies->Get(title_eff.c_str());
	}
  
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t Nentries = corrjets.fChain->GetEntries(); 
//   double weight = 93.47*16146/Nentries;//GJets
//   double weight = 95.14*16146/Nentries;//WJets
  double weight = 1;
  std::cout<<"Processing "<<Nentries<<"entries"<<" entries with weight "<<weight<<std::endl;
  
  for(Int_t entry = 0; entry < Nentries; ++entry){
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
      EMF[i] = corrjets.jet_efrac_ch_EM[i] + corrjets.jet_efrac_ne_EM[i];
    }
      
    bool pass_conv_1 = true;
    bool pass_conv_2 = true;
    if (corrjets.jet_efrac_photon[0] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR1 < 0.2) pass_conv_1 = false;
    if (corrjets.jet_efrac_photon[1] > 0.8 && corrjets.convtracks_pt[photon_nr]/corrjets.photon_pt[photon_nr] > 0.3 && dR2 < 0.2) pass_conv_2 = false;
    
    output->cd();
    
//     if (corrjets.HLT_DiCentralPFJet170 == 1){
    if (corrjets.HLT_PFJet450 == 1){
      passed_trigger += weight;
      trigger_eff->Fill(1.0, weight);
        
      if(corrjets.jet_pt[0] > 550 && corrjets.jet_pt[1] > 550){
        passed_pt += weight;
        pt_eff->Fill(550, weight);
    
        if (fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0){
          passed_eta += weight;
          eta_eff->Fill(2.0, weight);
              
          if (njets ==2) {                  
            passed_njets += weight;
            njets_eff->Fill(2.0, weight);
        
            if((corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && (corrjets.jet_efrac_photon[0] < 0.9 && corrjets.jet_efrac_photon[1] < 0.9)){              
              passed_photonveto += weight;
              photonVeto_eff->Fill(0.1, weight);
              
              if(pass_conv_1 && pass_conv_2){
                passed_conversionveto += weight;
                conversionVeto_eff->Fill(0.3, weight); 
                
                if (corrjets.vtx_N >= 2) {                    
                  passed_nvtx += weight;
                  nvtx_eff->Fill(2.0, weight);
        
                  if(deltajet_phi > 2){
                    passed_deltaphi += weight;
                    deltaphi_eff->Fill(2, weight);
                    
                    if (corrjets.Flag_HBHENoiseFilter == 1 && /*corrjets.Flag_HBHENoiseIsoFilter == 1 &&*/ corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && corrjets.Flag_goodVertices == 1 && corrjets.Flag_eeBadScFilter == 1 && corrjets.Flag_globalTightHalo2016Filter == 1){              
                      passed_filter += weight;
                      filter_eff->Fill(1.0, weight);            
                    
//                     if (corrjets.photon_passLooseId[photon_nr] == 0){
//                       if (corrjets.nEvent < 0){  
//                         UInt_t event = corrjets.nEvent;
//                         std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"/*<<corrjets.nEvent<<":"*/<<event<<std::endl; 
//                       }else std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<corrjets.nEvent<<std::endl;
//                     }
                    
                      for(int j = 0; j < 23; j++){
                        double eff1 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[0])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[0]));
                        double eff2 = eff_histos[j]->GetBinContent(eff_histos[j]->GetXaxis()->FindBin(fabs(corrjets.jet_eta[1])), eff_histos[j]->GetYaxis()->FindBin(corrjets.jet_pt[1]));
  //                       passed[j] += weight*eff1*eff2;
                        if (CHEF_corrjet[1] < chf_cuts[j]){
                          ChF_eff_1leg->Fill(chf_cuts[j], weight*eff1/2);
                          passed[j] += weight*eff1/2;
                        }
                        if (CHEF_corrjet[0] < chf_cuts[j]){
                          ChF_eff_1leg->Fill(chf_cuts[j], weight*eff2/2);
                          passed[j] += weight*eff2/2;
                        }
                        ChF_eff_2leg->Fill(chf_cuts[j], weight*eff1*eff2);
                        if (CHEF_corrjet[0] < chf_cuts[j] && CHEF_corrjet[1] < chf_cuts[j] && CHEF_SPVjet[0] < chf_cuts[j] && CHEF_SPVjet[1] < chf_cuts[j]){
                          passed_real[j] += weight; // only for Wjets and Gjets!!!
  //                        if (chf_cuts[j] == 0.1){
  //                          if (corrjets.nEvent < 0){  
  //                            UInt_t event = corrjets.nEvent;
  //                            std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"/*<<corrjets.nEvent<<":"*/<<event<<std::endl; 
  //                          }else std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<corrjets.nEvent<<std::endl;
  //                        }   
                        }
                      }
                    }
                  }
//                   else{
//                     if (CHEF_corrjet[0] < 0.05 && CHEF_corrjet[1] < 0.05){
//                       lowChF++;
// //                       if (corrjets.nEvent < 0){
// //                         UInt_t event = corrjets.nEvent;
// //                         std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<event<<std::endl;
// //                       }else std::cout<<corrjets.nRun<<":"<<corrjets.nLumi<<":"<<corrjets.nEvent<<std::endl;
//                       if (corrjets.Flag_HBHENoiseFilter == 0) HBHENoise++;
//                       if (corrjets.Flag_HBHENoiseIsoFilter == 0) HBHENoiseIso++;
//                       if (corrjets.Flag_EcalDeadCellTriggerPrimitiveFilter == 0) EcalDeadCell++;
//                       if (corrjets.Flag_goodVertices == 0) goodVertices++;
//                       if (corrjets.Flag_eeBadScFilter == 0) eeBadScFilter++;
//                       if (corrjets.Flag_globalTightHalo2016Filter == 0) globalTightHalo2016++;
//                     }
//                     else highChF++;
//                   }
                }
              }
            }
          }
        }
      }
    }
	}
	std::cout<<"trigger";
	 std::cout<<" & "<<passed_trigger;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>550$ GeV";
	 std::cout<<" & "<<passed_pt;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$|\\eta_{j1, j2}|<2.0$";
	 std::cout<<" & "<<passed_eta;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"njets";
	 std::cout<<" & "<<passed_njets;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto";
	 std::cout<<" & "<<passed_photonveto;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"conversion veto";
	 std::cout<<" & "<<passed_conversionveto;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"nvtx";
	 std::cout<<" & "<<passed_nvtx;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$";
	 std::cout<<" & "<<passed_deltaphi;
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"filters";
	 std::cout<<" & "<<passed_filter;
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 23; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"PRED ChF$_{j1, j2} <$ "<<cut;	
		 std::cout<<" & "<<passed[j];
		std::cout<<" \\\\"<<std::endl;
		std::cout<<"REAL ChF$_{j1, j2} <$ "<<cut;	
		 std::cout<<" & "<<passed_real[j];
		std::cout<<" \\\\"<<std::endl;
	}	
// 	std::cout<<"highChF: "<<highChF<<std::endl;	
// 	std::cout<<"lowChF: "<<lowChF<<std::endl;	
//   std::cout<<"HBHENoise"<<HBHENoise<<std::endl;
//   std::cout<<"HBHENoiseIso"<<HBHENoiseIso<<std::endl;
//   std::cout<<"EcalDeadCell"<<EcalDeadCell<<std::endl;
//   std::cout<<"goodVertices"<<goodVertices<<std::endl;
//   std::cout<<"eeBadScFilter"<<eeBadScFilter<<std::endl;
//   std::cout<<"globalTightHalo2016"<<globalTightHalo2016<<std::endl;
	output->Write();
	output->Close();	
}