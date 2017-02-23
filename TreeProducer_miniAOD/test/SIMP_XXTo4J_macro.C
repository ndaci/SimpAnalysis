#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

void SIMP_XXTo4J_macro(){
  
  TChain * chain = new TChain("tree/SimpAnalysis");
	
	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M300/170125_151018/0000/XXTo4J_1.root");

// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M500/170124_155759/0000/QCD_PUMoriond17_test_1.root");
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M500/170124_155759/0000/QCD_PUMoriond17_test_2.root");
// 
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M700/170124_155734/0000/QCD_PUMoriond17_test_1.root");
// 
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1000/170124_155709/0000/QCD_PUMoriond17_test_1.root");
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1000/170124_155709/0000/QCD_PUMoriond17_test_2.root");
// 
// 	chain->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1500/170124_155644/0000/QCD_PUMoriond17_test_1.root");
// 
// 	chain->Add("lists/XXTo4J_ctau-2m_M-3000.root");
	
// 	chain->Add("SIMPs_PUMoriond17_AOD_M1.root");
	
  TFile *output = new TFile("XXTo4J_ctau2000_M300.root", "RECREATE");
//   TFile *output = new TFile("XXTo4J_ctau2000_M500.root", "RECREATE");
//   TFile *output = new TFile("XXTo4J_ctau2000_M700.root", "RECREATE");
//   TFile *output = new TFile("XXTo4J_ctau2000_M1000.root", "RECREATE");
//   TFile *output = new TFile("XXTo4J_ctau2000_M1500.root", "RECREATE");
//   TFile *output = new TFile("XXTo4J_ctau2000_M3000.root", "RECREATE");
//   TFile *output = new TFile("SIMP_M1.root", "RECREATE");
//   
	int nJet, dijet_170, dijet_170_0p1, dijet_220_0p3, dijet_330_0p5, dijet_430;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8], EMF_jet[8];
	double track_pt[10], track_ptError[10], track_dzError[10], track_dz[10];
	int track_fromPV[10], nhits[10], nPixHits[10], ndof[3];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	int passed_deltaphi = 0, passed_eta = 0, passed_npix = 0, passed_photonveto = 0, passed_pt = 0;
	int passed[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
  
  chain->SetBranchAddress("nJet",&nJet);
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_phi", &jet_phi);
  chain->SetBranchAddress("jet_efrac_ch_Had", &jet_efrac_ch_Had);
  chain->SetBranchAddress("jet_efrac_ch_EM", &jet_efrac_ch_EM);
  chain->SetBranchAddress("jet_efrac_ch_Mu", &jet_efrac_ch_Mu);
	chain->SetBranchAddress("track_ptError", &track_ptError);
	chain->SetBranchAddress("track_nhits", &nhits);
	chain->SetBranchAddress("track_nPixHits", &nPixHits);
	chain->SetBranchAddress("track_pt", &track_pt);
	chain->SetBranchAddress("track_fromPV", &track_fromPV);
	chain->SetBranchAddress("track_dzError", &track_dzError);
	chain->SetBranchAddress("track_dz", &track_dz);
// 	chain->SetBranchAddress("photon_pt", &photon_pt);
// 	chain->SetBranchAddress("photon_eta", &photon_eta);
// 	chain->SetBranchAddress("photon_phi", &photon_phi);
//   chain->SetBranchAddress("photon_passLooseId",&passLooseId);
//   chain->SetBranchAddress("photon_passMediumId",&passMediumId);
//   chain->SetBranchAddress("photon_passTightId",&passTightId);
  chain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &dijet_170_0p1);
  chain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &dijet_220_0p3);
  chain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &dijet_330_0p5);
  chain->SetBranchAddress("HLT_DiCentralPFJet430", &dijet_430);
  chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
	
	TH1F *jet1_pt = new TH1F("jet1_pt", "Leading jet pt", 100, 0, 2000);
	TH1F *jet2_pt = new TH1F("jet2_pt", "Subleading jet pt", 100, 0, 2000);
	TH1F *jet3_pt = new TH1F("jet3_pt", "3rd jet pt", 100, 0, 2000);
	TH1F *jet4_pt = new TH1F("jet4_pt", "4th jet pt", 100, 0, 2000);
	TH1F *jet1_eta = new TH1F("jet1_eta", "Leading jet eta", 100, -3.14, 3.14);
	TH1F *jet2_eta = new TH1F("jet2_eta", "Subleading jet eta", 100, -3.14, 3.14);
	TH1F *jet1_phi = new TH1F("jet1_phi", "Leading jet phi", 100, -3.14, 3.14);
	TH1F *jet2_phi = new TH1F("jet2_phi", "Subleading jet phi", 100, -3.14, 3.14);
	TH1F *deltaphi = new TH1F("DeltaPhi", "DeltaPhi", 100, 0, 3.14);
	TH1F *jet1_chf = new TH1F("jet1_chf", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet1_chf_jet2_0p5 = new TH1F("jet1_chf_jet2_0p5", "Leading jet charged energy fraction (jet2_chf > 0.5)", 100, 0, 1);
	TH1F *jet2_chf = new TH1F("jet2_chf", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *jet3_chf = new TH1F("jet3_chf", "Charged energy fraction (3rd jet)", 100, 0, 1);
	TH1F *jet4_chf = new TH1F("jet4_chf", "Charged energy fraction (4th jet)", 100, 0, 1);
	TH1F *jet1_chf_withps = new TH1F("jet1_chf_withps", "Charged energy fraction (leading jet)", 100, 0, 1);
	TH1F *jet2_chf_withps = new TH1F("jet2_chf_withps", "Charged energy fraction (subleading jet)", 100, 0, 1);
	TH1F *trigger_dijet_170_0p1 = new TH1F("trigger_dijet_170_0p1", "Dijet_170_0p1", 2, -0.5, 1.5);
	TH1F *trigger_dijet_220_0p3 = new TH1F("trigger_dijet_220_0p3", "Dijet_220_0p3", 2, -0.5, 1.5);
	TH1F *trigger_dijet_330_0p5 = new TH1F("trigger_dijet_330_0p5", "Dijet_330_0p5", 2, -0.5, 1.5);
	TH1F *trigger_dijet_430 = new TH1F("trigger_dijet_430", "Dijet_430", 2, -0.5, 1.5);
	TH1F *njets = new TH1F("njets", "Number of jets", 2, -0.5, 10.5);
  
  Int_t Nentries = chain->GetEntries(); 
	std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
  for(Int_t entry = 0; entry < Nentries; ++entry){
    chain->GetEntry(entry);
    
    double deltajet_phi = jet_phi[0] - jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();
		deltajet_phi = fabs(deltajet_phi);
			
// 		double deltaphi_jet1photon = jet_phi[0] - photon_phi[0];
// 		if(deltaphi_jet1photon > TMath::Pi()) deltaphi_jet1photon -= 2*TMath::Pi();
// 		if(deltaphi_jet1photon < -TMath::Pi()) deltaphi_jet1photon += 2*TMath::Pi();
// 		double deltaphi_jet2photon = jet_phi[1] - photon_phi[0];
// 		if(deltaphi_jet2photon > TMath::Pi()) deltaphi_jet2photon -= 2*TMath::Pi();
// 		if(deltaphi_jet2photon < -TMath::Pi()) deltaphi_jet2photon += 2*TMath::Pi();
// 		
// 		double deltaeta_jet1photon = jet_eta[0] - photon_eta[0];
// 		double deltaeta_jet2photon = jet_eta[1] - photon_eta[0];
// 		
// 		double dR1 = TMath::Sqrt(deltaphi_jet1photon*deltaphi_jet1photon + deltaeta_jet1photon*deltaeta_jet1photon);
// 		double dR2 = TMath::Sqrt(deltaphi_jet2photon*deltaphi_jet2photon + deltaeta_jet2photon*deltaeta_jet2photon);
    
		for (int i = 0; i < 8; i++){
			CHEF_jet[i] = jet_efrac_ch_Had[i]+jet_efrac_ch_EM[i]+jet_efrac_ch_Mu[i];
		} 
		
		output->cd();
		
		if (fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0){
			passed_eta++;
			if(jet_pt[0] > 250 && jet_pt[1] > 250){
				passed_pt++;
				if(deltajet_phi > 2){
					passed_deltaphi++;
					if(nPixHits[0] > 0){
						passed_npix++;
// 						if(passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)){
							passed_photonveto++;
							jet1_pt->Fill(jet_pt[0]);
							jet2_pt->Fill(jet_pt[1]);
							jet3_pt->Fill(jet_pt[2]);
							jet4_pt->Fill(jet_pt[3]);
							jet1_eta->Fill(jet_eta[0]);
							jet2_eta->Fill(jet_eta[1]);
							jet1_phi->Fill(jet_phi[0]);
							jet2_phi->Fill(jet_phi[1]);
							deltaphi->Fill(deltajet_phi);
							jet1_chf->Fill(CHEF_jet[0]);
							if (CHEF_jet[1] > 0.5) jet1_chf_jet2_0p5->Fill(CHEF_jet[0]);
							jet2_chf->Fill(CHEF_jet[1]);
							jet3_chf->Fill(CHEF_jet[2]);
							jet4_chf->Fill(CHEF_jet[3]);
							trigger_dijet_170_0p1->Fill(dijet_170_0p1);
							trigger_dijet_220_0p3->Fill(dijet_220_0p3);
							trigger_dijet_330_0p5->Fill(dijet_330_0p5);
							trigger_dijet_430->Fill(dijet_430);
							for(int j = 0; j < 11; j++){
								if (CHEF_jet[0]<chf_cuts[j] && CHEF_jet[1]<chf_cuts[j]) passed[j]++;
							}
// 						}
					}
				}
			}
		}
	}
	std::cout<<"passed eta cut: "<<passed_eta<<std::endl;
	std::cout<<"passed pt cut: "<<passed_pt<<std::endl;
	std::cout<<"passed deltaphi cut: "<<passed_deltaphi<<std::endl;
	std::cout<<"passed nPixHits cut: "<<passed_npix<<std::endl;
	std::cout<<"passed photon veto: "<<passed_photonveto<<std::endl;
	for(int j = 0; j < 11; j++) std::cout<<passed[j]<<" ";
	
// 	TCanvas *c1 = new TCanvas("Signal yield", "Signal yield");
// 	c1->cd();
// 	TGraphErrors *XXTo4J = new TGraphErrors(11, chf_cuts, passed, zero, zero);
// 	XXTo4J->SetTitle("Signal yield");
// 	XXTo4J->GetXaxis()->SetTitle("CHF cut");
// 	XXTo4J->GetYaxis()->SetTitle("# events");
// 	XXTo4J->SetMarkerStyle(20);
// 	XXTo4J->Draw("AP");
	
  output->Write();
  output->Close();
}