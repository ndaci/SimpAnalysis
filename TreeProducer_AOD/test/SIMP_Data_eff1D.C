
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
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TH2.h>

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

void SIMP_Data_eff1D(){
	  
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
	int photon_nr;
  
// 	TRandom3 r;
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	double pt_bins[10] = {250, 275, 300, 350, 400, 450, 550, 700, 900, 10000};
	double eta_bins[5] = {0, 0.5, 1.0, 1.5, 2.0};
  double njets_bins[5] = {1.5, 2.5, 3.5, 4.5, 10.5};
	double total_pt[9] = {0,0,0,0,0,0,0,0,0};
	double total_eta[4] = {0,0,0,0};
	double total_njets[4] = {0,0,0,0};
	double passed_pt[9][11];
	double passed_eta[4][11];
	double passed_njets[4][11];
	double err_pt_passed[9][11];
	double err_eta_passed[4][11];
	double err_njets_passed[4][11];
	double err_pt_failed[9][11];
	double err_eta_failed[4][11];
	double err_njets_failed[4][11];
	double eff_pt[9][11], red_pt[9][11], eff_eta[4][11], red_eta[4][11], eff_njets[4][11], red_njets[4][11];
	double err_eff_pt[9][11], err_eff_eta[9][11], err_eff_njets[9][11], err_red_pt[4][11], err_red_eta[4][11], err_red_njets[4][11];
	double zero[11] = {0,0,0,0,0,0,0,0,0,0,0};
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 9; i++){
			passed_pt[i][j] = 0;
			err_pt_passed[i][j] = 0;
			err_pt_failed[i][j] = 0;
		}
	}
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 4; i++){
			passed_eta[i][j] = 0;
			err_eta_passed[i][j] = 0;
			err_eta_failed[i][j] = 0;
			passed_njets[i][j] = 0;
			err_njets_passed[i][j] = 0;
			err_njets_failed[i][j] = 0;
		}
	}
	
  TFile *output = new TFile("eff1D_Data_AOD.root", "RECREATE");
  
  SANtuple corrjets;
  corrjets.Init(chain);
  SANtuple SPVjets;
  SPVjets.Init(SPVchain);
  
  Int_t Nentries = corrjets.fChain->GetEntries(); 
  std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
  
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
      CHEF_corrjet[i] = corrjets.jet_efrac_ch_Had[i]+corrjets.jet_efrac_ch_EM[i]+corrjets.jet_efrac_ch_Mu[i];
    } 
    
    double weight = corrjets.pswgt_dijet_170;
    
    output->cd();
    
    if (corrjets.HLT_DiCentralPFJet170 == 1 && corrjets.jet_pt[0] > 250.0 && corrjets.jet_pt[1] > 250.0 && fabs(corrjets.jet_eta[0]) < 2.0 && fabs(corrjets.jet_eta[1]) < 2.0 && deltajet_phi > 2.0 && ( corrjets.photon_passLooseId[photon_nr] == 0 || (corrjets.photon_passLooseId[photon_nr] == 1 && dR1 > 0.1 && dR2 > 0.1)) && corrjets.vtx_N >= 2){ 
      for (int i = 0; i < 4; i++){
        if (njets>njets_bins[i] && njets<njets_bins[i+1]){
          total_njets[i] += weight;
          for(int j = 0; j < 11; j++){
            if (CHEF_corrjet[1]<chf_cuts[j]){
              passed_njets[i][j] += weight;
              err_njets_passed[i][j] += weight*weight;
            }else err_njets_failed[i][j] += weight*weight;
          }
        }
      }      
      if(CHEF_corrjet[0] > 0.5){
        for (int i = 0; i < 9; i++){
          if (corrjets.jet_pt[1]>pt_bins[i] && corrjets.jet_pt[1]<pt_bins[i+1]){
            total_pt[i] += weight;
            for(int j = 0; j < 11; j++){
              if (CHEF_corrjet[1]<chf_cuts[j]){
                passed_pt[i][j] += weight;
                err_pt_passed[i][j] += weight*weight;
              }else err_pt_failed[i][j] += weight*weight;
            }
          }
        }
        for (int i = 0; i < 4; i++){
          if (corrjets.jet_eta[1]>eta_bins[i] && corrjets.jet_eta[1]<eta_bins[i+1]){
            total_eta[i] += weight;
            for(int j = 0; j < 11; j++){
              if (CHEF_corrjet[1]<chf_cuts[j]){
                passed_eta[i][j] += weight;
                err_eta_passed[i][j] += weight*weight;
              }else err_eta_failed[i][j] += weight*weight;
            }
          }
        }
      }
      else if(CHEF_corrjet[1] > 0.5){
        for (int i = 0; i < 9; i++){
          if (corrjets.jet_pt[0]>pt_bins[i] && corrjets.jet_pt[0]<pt_bins[i+1]){
            total_pt[i] += weight;
            for(int j = 0; j < 11; j++){
              if (CHEF_corrjet[0]<chf_cuts[j]){
                passed_pt[i][j] += weight;
                err_pt_passed[i][j] += weight*weight;
              }else err_pt_failed[i][j] += weight*weight;
            }
          }
        }
        for (int i = 0; i < 4; i++){
          if (corrjets.jet_eta[0]>eta_bins[i] && corrjets.jet_eta[0]<eta_bins[i+1]){
            total_eta[i] += weight;
            for(int j = 0; j < 11; j++){
              if (CHEF_corrjet[0]<chf_cuts[j]){
                passed_eta[i][j] += weight;
                err_eta_passed[i][j] += weight*weight;
              }else err_eta_failed[i][j] += weight*weight;
            }
          }
        }
      }
    }    
  }
	
	for (int j = 0; j < 11; j++){
		for (int i = 0; i < 9; i++){
			double tot = (double)total_pt[i];
			double passed = (double)passed_pt[i][j];
			err_pt_passed[i][j] = TMath::Sqrt(err_pt_passed[i][j]);
			err_pt_failed[i][j] = TMath::Sqrt(err_pt_failed[i][j]);
			eff_pt[i][j] = passed/tot;
			err_eff_pt[i][j] = pow(tot, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_pt_passed[i][j], 2) + (pow(passed, 2)*pow(err_pt_failed[i][j], 2)));
			if (passed_pt[i][j] != 0){
				red_pt[i][j] = tot/(double)passed;
				err_red_pt[i][j] = pow(passed, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_pt_passed[i][j], 2) + (pow(passed, 2)*pow(err_pt_failed[i][j], 2)));
			}
			else{
				red_pt[i][j] = 0;
				err_red_pt[i][j] = 0;
			}
		}
	}
  for (int j = 0; j < 11; j++){
		for (int i = 0; i < 4; i++){
			double tot = (double)total_eta[i];
			double passed = (double)passed_eta[i][j];
			err_eta_passed[i][j] = TMath::Sqrt(err_eta_passed[i][j]);
			err_eta_failed[i][j] = TMath::Sqrt(err_eta_failed[i][j]);
			eff_eta[i][j] = passed/tot;
			err_eff_eta[i][j] = pow(tot, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_eta_passed[i][j], 2) + (pow(passed, 2)*pow(err_eta_failed[i][j], 2)));
			if (passed_eta[i][j] != 0){
				red_eta[i][j] = tot/passed;
				err_red_eta[i][j] = pow(passed, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_eta_passed[i][j], 2) + (pow(passed, 2)*pow(err_eta_failed[i][j], 2)));
			}
			else{
				red_eta[i][j] = 0;
				err_red_eta[i][j] = 0;
			}
		}
		for (int i = 0; i < 4; i++){
			double tot = (double)total_njets[i];
			double passed = (double)passed_njets[i][j];
			err_njets_passed[i][j] = TMath::Sqrt(err_njets_passed[i][j]);
			err_njets_failed[i][j] = TMath::Sqrt(err_njets_failed[i][j]);
			eff_njets[i][j] = passed/tot;
			err_eff_njets[i][j] = pow(tot, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_njets_passed[i][j], 2) + (pow(passed, 2)*pow(err_njets_failed[i][j], 2)));
			if (passed_njets[i][j] != 0){
				red_njets[i][j] = tot/passed;
				err_red_njets[i][j] = pow(passed, -2)*TMath::Sqrt(pow(tot - passed, 2)*pow(err_njets_passed[i][j], 2) + (pow(passed, 2)*pow(err_njets_failed[i][j], 2)));
			}
			else{
				red_njets[i][j] = 0;
				err_red_njets[i][j] = 0;
			}
		}
	}	
	
	TCanvas *c1 = new TCanvas("efficiency_pt", "efficiency_pt");
	c1->cd();
	TGraphErrors *firstptbin = new TGraphErrors(11, chf_cuts, eff_pt[0], zero, err_eff_pt[0]);
	firstptbin->SetTitle("");
	firstptbin->GetXaxis()->SetTitle("CHF cut");
	firstptbin->GetYaxis()->SetTitle("efficiency");
	firstptbin->SetMarkerStyle(20);
	firstptbin->Draw("AP");
	for (int i = 1; i < 9; i++){
		TGraphErrors *otherptbin = new TGraphErrors(11, chf_cuts, eff_pt[i], zero, err_eff_pt[i]);
		otherptbin->SetMarkerStyle(20);
		otherptbin->SetMarkerColor(i+1);
		otherptbin->Draw("P");
	}
	
	TCanvas *c2 = new TCanvas("rejection_pt", "rejection_pt");
	c2->cd();
	TGraphErrors *firstptbin2 = new TGraphErrors(11, chf_cuts, red_pt[0], zero, err_red_pt[0]);
	firstptbin2->SetTitle("");
	firstptbin2->GetXaxis()->SetTitle("CHF cut");
	firstptbin2->GetYaxis()->SetTitle("rejection factor");
	firstptbin2->SetMarkerStyle(20);
	firstptbin2->Draw("AP");
	for (int i = 1; i < 9; i++){
		TGraphErrors *otherptbin2 = new TGraphErrors(11, chf_cuts, red_pt[i], zero, err_red_pt[i]);
		otherptbin2->SetMarkerStyle(20);
		otherptbin2->SetMarkerColor(i+1);
		otherptbin2->Draw("P");
	}
	
	TCanvas *c3 = new TCanvas("efficiency_eta", "efficiency_eta");
	c3->cd();
	TGraphErrors *firstetabin = new TGraphErrors(11, chf_cuts, eff_eta[0], zero, err_eff_eta[0]);
	firstetabin->SetTitle("");
	firstetabin->GetXaxis()->SetTitle("CHF cut");
	firstetabin->GetYaxis()->SetTitle("efficiency");
	firstetabin->SetMarkerStyle(20);
	firstetabin->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *otheretabin = new TGraphErrors(11, chf_cuts, eff_eta[i], zero, err_eff_eta[i]);
		otheretabin->SetMarkerStyle(20);
		otheretabin->SetMarkerColor(i+1);
		otheretabin->Draw("P");
	}
	
	TCanvas *c4 = new TCanvas("rejection_eta", "rejection_eta");
	c4->cd();
	TGraphErrors *firstetabin2 = new TGraphErrors(11, chf_cuts, red_eta[0], zero, err_red_eta[0]);
	firstetabin2->SetTitle("");
	firstetabin2->GetXaxis()->SetTitle("CHF cut");
	firstetabin2->GetYaxis()->SetTitle("rejection factor");
	firstetabin2->SetMarkerStyle(20);
	firstetabin2->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *otheretabin2 = new TGraphErrors(11, chf_cuts, red_eta[i], zero, err_red_eta[i]);
		otheretabin2->SetMarkerStyle(20);
		otheretabin2->SetMarkerColor(i+1);
		otheretabin2->Draw("P");
	}
	
	TCanvas *c5 = new TCanvas("efficiency_njets", "efficiency_njets");
	c5->cd();
	TGraphErrors *firstnjetbin = new TGraphErrors(11, chf_cuts, eff_njets[0], zero, err_eff_njets[0]);
	firstnjetbin->SetTitle("");
	firstnjetbin->GetXaxis()->SetTitle("CHF cut");
	firstnjetbin->GetYaxis()->SetTitle("efficiency");
	firstnjetbin->SetMarkerStyle(20);
	firstnjetbin->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *othernjetbin = new TGraphErrors(11, chf_cuts, eff_njets[i], zero, err_eff_njets[i]);
		othernjetbin->SetMarkerStyle(20);
		othernjetbin->SetMarkerColor(i+1);
		othernjetbin->Draw("P");
	}
	
	TCanvas *c6 = new TCanvas("rejection_njets", "rejection_njets");
	c6->cd();
	TGraphErrors *firstnjetbin2 = new TGraphErrors(11, chf_cuts, red_njets[0], zero, err_red_njets[0]);
	firstnjetbin2->SetTitle("");
	firstnjetbin2->GetXaxis()->SetTitle("CHF cut");
	firstnjetbin2->GetYaxis()->SetTitle("rejection factor");
	firstnjetbin2->SetMarkerStyle(20);
	firstnjetbin2->Draw("AP");
	for (int i = 1; i < 4; i++){
		TGraphErrors *othernjetbin2 = new TGraphErrors(11, chf_cuts, red_njets[i], zero, err_red_njets[i]);
		othernjetbin2->SetMarkerStyle(20);
		othernjetbin2->SetMarkerColor(i+1);
		othernjetbin2->Draw("P");
	}
	
	output->Append(c1);
	output->Append(c2);
	output->Append(c3);
	output->Append(c4);
	output->Append(c5);
	output->Append(c6);
	
  output->Write();
  output->Close();
}