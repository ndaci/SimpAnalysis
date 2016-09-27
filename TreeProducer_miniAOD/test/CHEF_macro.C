#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>

void CHEF_macro(TString sampleName){
  
  TString sample = sampleName + ".root";
  
  TFile *input = new TFile(sample, "READ");
  TTree *tree = input->Get("tree/SimpAnalysis");
  
  double Pt[3], Eta[3], Phi[3], Ch_Had[3], Ch_EM[3], Ch_Mu[3], E[3], M[3];
  
  tree->SetBranchAddress("jet_pt", &Pt);
  tree->SetBranchAddress("jet_eta", &Eta);
  tree->SetBranchAddress("jet_phi", &Phi);
  tree->SetBranchAddress("jet_e", &E);
  tree->SetBranchAddress("jet_m", &M);
  tree->SetBranchAddress("jet_efrac_ch_Had", &Ch_Had);
  tree->SetBranchAddress("jet_efrac_ch_EM", &Ch_EM);
  tree->SetBranchAddress("jet_efrac_ch_Mu", &Ch_Mu);
  
  TString outputfile = sampleName+"_CHEF.root";
  std::cout<<"output: "<<outputfile<<std::endl;
  TFile *output = new TFile(outputfile, "RECREATE");
  
  double CHEF_jet1, CHEF_jet2, CHEF_jet3, Pt_jet1, Pt_jet2, Pt_jet3, Eta_jet1, Eta_jet2, Eta_jet3, DeltaPhi, Energy_jet1, Energy_jet2, Energy_jet3;
  double CHEF_sel_jet1, CHEF_sel_jet2, CHEF_sel_jet3, Pt_sel_jet1, Pt_sel_jet2, Pt_sel_jet3, Eta_sel_jet1, Eta_sel_jet2, Eta_sel_jet3, DeltaPhi_sel, Energy_sel_jet1, Energy_sel_jet2, Energy_sel_jet3;
  
  TTree *newtree = new TTree("CHEF_All", "CHEF_All");
  newtree->Branch("CHEF_jet1", &CHEF_jet1);
  newtree->Branch("CHEF_jet2", &CHEF_jet2);
  newtree->Branch("CHEF_jet3", &CHEF_jet3);
  newtree->Branch("Pt_jet1", &Pt_jet1);
  newtree->Branch("Pt_jet2", &Pt_jet2);
  newtree->Branch("Pt_jet3", &Pt_jet3);
  newtree->Branch("Eta_jet1", &Eta_jet1);
  newtree->Branch("Eta_jet2", &Eta_jet2);
  newtree->Branch("Eta_jet3", &Eta_jet3);
  newtree->Branch("DeltaPhi", &DeltaPhi);
  newtree->Branch("Energy_jet1", &Energy_jet1);
  newtree->Branch("Energy_jet2", &Energy_jet2);
  newtree->Branch("Energy_jet3", &Energy_jet3);
    
  TTree *seltree = new TTree("CHEF_selection", "CHEF_selection");
  seltree->Branch("CHEF_jet1_selection", &CHEF_sel_jet1);
  seltree->Branch("CHEF_jet2_selection", &CHEF_sel_jet2);
  seltree->Branch("CHEF_jet3_selection", &CHEF_sel_jet3);
  seltree->Branch("Pt_jet1_selection", &Pt_sel_jet1);
  seltree->Branch("Pt_jet2_selection", &Pt_sel_jet2);
  seltree->Branch("Pt_jet3_selection", &Pt_sel_jet3);
  seltree->Branch("Eta_jet1_selection", &Eta_sel_jet1);
  seltree->Branch("Eta_jet2_selection", &Eta_sel_jet2);
  seltree->Branch("Eta_jet3_selection", &Eta_sel_jet3);
  seltree->Branch("DeltaPhi_selection", &DeltaPhi_sel);
  seltree->Branch("Energy_jet1_selection", &Energy_sel_jet1);
  seltree->Branch("Energy_jet2_selection", &Energy_sel_jet2);
  seltree->Branch("Energy_jet3_selection", &Energy_sel_jet3);
  
  Int_t Nentries = tree->GetEntries(); 
  for(Int_t entry = 0; entry < Nentries; ++entry){
    tree->GetEntry(entry);
    
    double deltaPhi = Phi[0] - Phi[1];
    if(deltaPhi > TMath::Pi()) deltaPhi -= 2*TMath::Pi();
    if(deltaPhi < -TMath::Pi()) deltaPhi += 2*TMath::Pi();
    
    CHEF_jet1 = Ch_Had[0]+Ch_EM[0]+Ch_Mu[0];
    CHEF_jet2 = Ch_Had[1]+Ch_EM[1]+Ch_Mu[1];
    CHEF_jet3 = Ch_Had[2]+Ch_EM[2]+Ch_Mu[2];    
    
    Pt_jet1 = Pt[0];
    Pt_jet2 = Pt[1];
    Pt_jet3 = Pt[2];
    
    Eta_jet1 = Eta[0];
    Eta_jet2 = Eta[1];
    Eta_jet3 = Eta[2];
    
    DeltaPhi = deltaPhi;
    
    Energy_jet1 = E[0];
    Energy_jet2 = E[1];
    Energy_jet3 = E[2];
    
    newtree->Fill();
    
    if(fabs(Eta[0]) < 2.0 && fabs(Eta[1]) < 2.0 && Pt[0] > 350 && Pt[1] > 350 && fabs(deltaPhi) > 2){
      CHEF_sel_jet1 = Ch_Had[0]+Ch_EM[0]+Ch_Mu[0];
      CHEF_sel_jet2 = Ch_Had[1]+Ch_EM[1]+Ch_Mu[1];
      CHEF_sel_jet3 = Ch_Had[2]+Ch_EM[2]+Ch_Mu[2];      
    
      Pt_sel_jet1 = Pt[0];
      Pt_sel_jet2 = Pt[1];
      Pt_sel_jet3 = Pt[2];
      
      Eta_sel_jet1 = Eta[0];
      Eta_sel_jet2 = Eta[1];
      Eta_sel_jet3 = Eta[2];
    
      DeltaPhi_sel = deltaPhi;
      
      Energy_sel_jet1 = E[0];
      Energy_sel_jet2 = E[1];
      Energy_sel_jet3 = E[2];
    
      seltree->Fill();
    }
  }
  output->Write();
  output->Close();
}