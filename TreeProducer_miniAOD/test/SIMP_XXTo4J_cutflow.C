#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include <sstream> 
#include <stdio.h>
#include <TMath.h>

void SIMP_XXTo4J_cutflow(){
	
	double mass[8] = {50, 100, 300, 500, 700, 1000, 1500, 3000};
  
  TChain * chain0 = new TChain("tree/SimpAnalysis");
  TChain * chain1 = new TChain("tree/SimpAnalysis");
  TChain * chain2 = new TChain("tree/SimpAnalysis");
  TChain * chain3 = new TChain("tree/SimpAnalysis");
  TChain * chain4 = new TChain("tree/SimpAnalysis");
  TChain * chain5 = new TChain("tree/SimpAnalysis");
  TChain * chain6 = new TChain("tree/SimpAnalysis");
  TChain * chain7 = new TChain("tree/SimpAnalysis");
  
  //ctau = 2m
//   chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M50/170206_152723/0000/QCD_1.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_1.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_2.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_3.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_4.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_5.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_6.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_7.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M100/170130_091407/0000/QCD_8.root");
// 	chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M300/170125_151018/0000/XXTo4J_1.root");
// 	chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M500/170125_151048/0000/XXTo4J_1.root");
// 	chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M500/170125_151048/0000/XXTo4J_2.root");
// 	chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M700/170125_151108/0000/XXTo4J_1.root");
// 	chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1000/170125_151140/0000/XXTo4J_1.root");
// 	chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1000/170125_151140/0000/XXTo4J_2.root");
// 	chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau2000_M1500/170125_151158/0000/XXTo4J_1.root");
// 	chain7->Add("lists/XXTo4J_ctau-2m_M-3000.root");
  
  //ctau = 1m
//   chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M50/170206_141728/0000/QCD_1.root");
//   chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M50/170206_141728/0000/QCD_2.root");
//   chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M50/170206_141728/0000/QCD_3.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M100/170206_141810/0000/QCD_1.root");
//   chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M300/170206_141835/0000/QCD_1.root");
//   chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M500/170206_141931/0000/QCD_1.root");
//   chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M700/170206_141950/0000/QCD_1.root");
//   chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M700/170206_141950/0000/QCD_2.root");
//   chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M1000/170206_142028/0000/QCD_1.root");
//   chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M1500/170206_142056/0000/QCD_1.root");
//   chain7->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau1000_M3000/170206_142117/0000/QCD_1.root");
  
  //ctau = 30cm
//   chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M50/170206_144704/0000/QCD_1.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M100/170206_144247/0000/QCD_1.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M100/170206_144247/0000/QCD_2.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M100/170206_144247/0000/QCD_3.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M100/170206_144247/0000/QCD_4.root");
//   chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M100/170206_144247/0000/QCD_5.root");
//   chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M300/170206_144226/0000/QCD_1.root");
//   chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M300/170206_144226/0000/QCD_2.root");
//   chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M500/170206_144202/0000/QCD_1.root");
//   chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M700/170206_144141/0000/QCD_1.root");
//   chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M1000/170206_144118/0000/QCD_1.root");
//   chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M1000/170206_144118/0000/QCD_2.root");
//   chain6->Add(" /pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M1500/170206_144058/0000/QCD_1.root");
//   chain6->Add(" /pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-300mm_Tu/store/mc/RunIISpring16MiniAODv2/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/127F1E62-253B-E611-A51D-001E673D3629.rootneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M1500/170206_144058/0000/QCD_2.root");
//   chain7->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau300_M3000/170206_144039/0000/QCD_1.root");
  
  //ctau = 10cm
  chain0->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-50_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M50/170206_154235/0000/QCD_1.root");
  chain1->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-100_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M100/170206_154329/0000/QCD_1.root");
	chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M300/170206_154400/0000/QCD_1.root");
	chain2->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-300_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M300/170206_154400/0000/QCD_2.root");
	chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M500/170206_154420/0000/QCD_1.root");
	chain3->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M500/170206_154420/0000/QCD_2.root");
	chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M700/170206_154513/0000/QCD_1.root");
	chain4->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-700_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M700/170206_154513/0000/QCD_2.root");
	chain5->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1000_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M1000/170206_154539/0000/QCD_1.root");
	chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M1500/170206_154559/0000/QCD_1.root");
	chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M1500/170206_154559/0000/QCD_2.root");
	chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M1500/170206_154559/0000/QCD_3.root");
	chain6->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-1500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M1500/170206_154559/0000/QCD_4.root");
	chain7->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M3000/170206_154618/0000/QCD_1.root");
	chain7->Add("/pnfs/iihe/cms/store/user/isdebruy/XXTo4J_M-3000_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/SIMPs_XXTo4J_ctau100_M3000/170206_154618/0000/QCD_2.root");
	
	TChain* chains[8] = {chain0, chain1, chain2, chain3, chain4, chain5, chain6, chain7};
//
	TFile* output = new TFile("XXTo4J_ctau100_cutflow.root", "RECREATE");
	
	int nJet, dijet_170, dijet_170_0p1, dijet_220_0p3, dijet_330_0p5, dijet_430;
  double jet_pt[8], jet_eta[8], jet_phi[8], jet_efrac_ch_Had[8], jet_efrac_ch_EM[8], jet_efrac_ch_Mu[8], CHEF_jet[8], EMF_jet[8];
	double track_pt[10], track_ptError[10], track_dzError[10], track_dz[10];
	int track_fromPV[10], nhits[10], nPixHits[10], ndof[3];
	double photon_eta[4], photon_phi[4], photon_pt[4];
	int passLooseId[4], passMediumId[4], passTightId[4];
	
	double passed_deltaphi[8] = {0,0,0,0,0,0,0,0};
	double passed_eta[8] = {0,0,0,0,0,0,0,0};
	double passed_npix[8] = {0,0,0,0,0,0,0,0};
	double passed_photonveto[8] = {0,0,0,0,0,0,0,0};
	double passed_pt[8] = {0,0,0,0,0,0,0,0};
	double passed[8][11] = {{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0}};
	double chf_cuts[11] = {0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01};
	
	double lumi = 23.015;
  // in fb
// 	double QCD_xsec[8] = {0, 199.6, 3.156, 0.3376, 0.06179, 0.007885, 0.0004316, 0.0000002283};//ctau = 2m
// 	double QCD_xsec[8] = {0, 198.9, 3.182, 0.3373, 0.06228, 0.007855, 0.0004334, 0.0000002296};// ctau = 1m
// 	double QCD_xsec[8] = {802.3, 200.5, 3.219, 0.3384, 0.06193, 0.007848, 0.0004370, 0.0000002280}; //ctau = 30cm
	double QCD_xsec[8] = {805.6, 199.4, 3.183, 0.3348, 0.06225, 0.007819, 0.0004371, 0.0000002288}; //ctau = 10cm

	TH1D* pt_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* eta_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* deltaphi_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* nPix_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* photonVeto_eff[8] = {0,0,0,0,0,0,0,0};
	TH1D* ChF_eff[8] = {0,0,0,0,0,0,0,0};
			
	for (int l = 0; l < 8; l++){
		std::ostringstream strs;
		double dbl = mass[l];
		strs << dbl;
		std::string m = strs.str();
		
		pt_eff[l] = new TH1D(("M"+m+"_pt_eff").c_str(), "passed pt cuts", 100, 200, 300);
		eta_eff[l] = new TH1D(("M"+m+"_eta_eff").c_str(), "passed eta cuts", 100, 0, 5);
		deltaphi_eff[l] = new TH1D(("M"+m+"_deltaphi_eff").c_str(), "passed #Delta#phi cut", 100, 0, 5);
		nPix_eff[l] = new TH1D(("M"+m+"_nPix_eff").c_str(), "passed nPixHits cut", 10, 0, 10);
		photonVeto_eff[l] = new TH1D(("M"+m+"_photonVeto_eff").c_str(), "passed photon veto", 10, 0, 1);
		ChF_eff[l] = new TH1D(("M"+m+"_ChF_eff").c_str(), "passed ChF cuts", 100, 0, 0.5);
	}
  
	for (int l = 0; l < 8; l++){
		TChain* chain = chains[l];
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
		chain->SetBranchAddress("photon_pt", &photon_pt);
		chain->SetBranchAddress("photon_eta", &photon_eta);
		chain->SetBranchAddress("photon_phi", &photon_phi);
		chain->SetBranchAddress("photon_passLooseId",&passLooseId);
		chain->SetBranchAddress("photon_passMediumId",&passMediumId);
		chain->SetBranchAddress("photon_passTightId",&passTightId);
		chain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &dijet_170_0p1);
		chain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &dijet_220_0p3);
		chain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &dijet_330_0p5);
		chain->SetBranchAddress("HLT_DiCentralPFJet430", &dijet_430);
		chain->SetBranchAddress("HLT_DiCentralPFJet170", &dijet_170);
		
		Int_t Nentries = chain->GetEntries(); 
		double weight = QCD_xsec[l]*lumi/Nentries;
		std::cout<<"Processing "<<Nentries<<"entries"<<std::endl;
		for(Int_t entry = 0; entry < Nentries; ++entry){
			chain->GetEntry(entry);
			
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
			} 
			
			output->cd();
			
			if (fabs(jet_eta[0]) < 2.0 && fabs(jet_eta[1]) < 2.0){
				passed_eta[l] += weight;
				eta_eff[l]->Fill(2.0, weight);
				if(jet_pt[0] > 250 && jet_pt[1] > 250){
					passed_pt[l] += weight;
					pt_eff[l]->Fill(250, weight);
					if(deltajet_phi > 2){
						passed_deltaphi[l] += weight;
						deltaphi_eff[l]->Fill(2, weight);
						if(nPixHits[0] > 0){
							passed_npix[l] += weight;
							nPix_eff[l]->Fill(0.0, weight);
							if(passLooseId[0] == 0 || (passLooseId[0] == 1 && dR1 > 0.1 && dR2 > 0.1)){
								passed_photonveto[l] += weight;
								photonVeto_eff[l]->Fill(0.1, weight);
								for(int j = 0; j < 11; j++){
									if (CHEF_jet[0]<chf_cuts[j] && CHEF_jet[1]<chf_cuts[j]){
										passed[l][j] += weight;
										ChF_eff[l]->Fill(chf_cuts[j], weight);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	std::cout<<"XXTo4J & $M=300$ GeV & $M=500$ GeV & $M=700$ GeV & $M=1000$ GeV & $M=1500$ GeV & $M=3000$ GeV"<<std::endl; 
	std::cout<<"$|\\eta_{j1, j2}|<2.0$ & ";
	for(int l = 0; l < 6; l++) std::cout<<" & "<<passed_eta[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$p_T^{j1, j2}>250$ GeV & ";
	for(int l = 0; l < 6; l++) std::cout<<" & "<<passed_pt[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$\\Delta\\phi_{j1, j2}>2.0$ & ";
	for(int l = 0; l < 6; l++) std::cout<<" & "<<passed_deltaphi[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"$# pixel hits > 0$ & ";
	for(int l = 0; l < 6; l++) std::cout<<" & "<<passed_npix[l];
	std::cout<<" \\\\"<<std::endl;
	std::cout<<"photon veto & ";
	for(int l = 0; l < 6; l++) std::cout<<" & "<<passed_photonveto[l];
	std::cout<<" \\\\"<<std::endl;
	for(int j = 0; j < 11; j++){
		std::ostringstream strs;
		double dbl = chf_cuts[j];
		strs << dbl;
		std::string cut = strs.str();
		std::cout<<"ChF$_{j1, j2}$ < "<<cut;	
		for(int l = 0; l < 6; l++) std::cout<<" & "<<passed[l][j];
		std::cout<<" \\\\"<<std::endl;
	}
	
	output->Write();
	output->Close();	
}