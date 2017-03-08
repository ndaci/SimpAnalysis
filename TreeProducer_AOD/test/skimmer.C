#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TBranch.h"
#include "TString.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include <algorithm>
#include "TH1F.h"
#include "TMath.h"

#include <vector>

#include "SANtuple.h"

int skimmer(){

  TChain *oldtreeCorr = new TChain("treeCorr/SimpAnalysis");
  TChain *oldtreeSPV = new TChain("treeSPV/SimpAnalysis");

  //2000ToInf 67 ext 151 
  //1500To2000 170 ext 318
  //1000To1500 191 ext 365
  //700To1000 462 ext 878
  //500To700 571 ext 999,1345
  //300To500 494 ext 999,1040
  
  for(int i=1; i<=67; i++){
    stringstream a;
    a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT2000ToInf_PUMoriond17_AOD/170301_175146/0000/QCD_PUMoriond17_AOD_test_";
//     a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT1500To2000_PUMoriond17_AOD/170301_174832/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT1000To1500_PUMoriond17_AOD/170301_173348/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT700To1000_PUMoriond17_AOD/170301_173021/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT500To700_PUMoriond17_AOD/170301_163332/0000/QCD_PUMoriond17_AOD_test_";
// 		a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT300To500_PUMoriond17_AOD/170301_163036/0000/QCD_PUMoriond17_AOD_test_";

    a<<i; a<<".root";
    oldtreeCorr->Add((a.str()).c_str());
    oldtreeSPV->Add((a.str()).c_str());
  }
  
  
  for(int i=1; i<=151; i++){
    stringstream a;
    a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT2000ToInf_ext_PUMoriond17_AOD/170301_172407/0000/QCD_PUMoriond17_AOD_test_";
//     a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT1500To2000_ext_PUMoriond17_AOD/170301_172021/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT1000To1500_ext_PUMoriond17_AOD/170301_163332/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT700To1000_ext_PUMoriond17_AOD/170301_163036/0000/QCD_PUMoriond17_AOD_test_";
		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT500To700_ext_PUMoriond17_AOD/170301_162829/0000/QCD_PUMoriond17_AOD_test_";
// 		a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT300To500_ext_PUMoriond17_AOD/170301_162539/0000/QCD_PUMoriond17_AOD_test_";

    a<<i; a<<".root";
    oldtreeCorr->Add((a.str()).c_str());
    oldtreeSPV->Add((a.str()).c_str());
  }
  
//   for(int i=1000; i<=1040; i++){
//     stringstream a;
// 		//a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT500To700_ext_PUMoriond17_AOD/170301_162829/0001/QCD_PUMoriond17_AOD_test_";
// 		a<<"srm://maite.iihe.ac.be:8443//pnfs/iihe/cms/store/user/isdebruy/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SIMPs_QCD_HT300To500_ext_PUMoriond17_AOD/170301_162539/0001/QCD_PUMoriond17_AOD_test_";
// 
//     a<<i; a<<".root";
//     oldtreeCorr->Add((a.str()).c_str());
//     oldtreeSPV->Add((a.str()).c_str());
//   }

  SANtuple c_treeCorr;
  c_treeCorr.Init(oldtreeCorr);
  SANtuple c_treeSPV;
  c_treeSPV.Init(oldtreeSPV);

  Long64_t nentries = oldtreeCorr->GetEntries();

  oldtreeCorr->SetBranchStatus("*",1);
  oldtreeSPV->SetBranchStatus("*",1);

  const double eta_cut = 2.0, dphi_cut= 2.0, pt_cut = 250.;

   //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT2000ToInf_PUMoriond17.root","recreate");
//   TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT1500To2000_PUMoriond17.root","recreate");
//   TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT1000To1500_PUMoriond17.root","recreate");
//   TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT700To1000_PUMoriond17.root","recreate");
//   TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT500To700_PUMoriond17.root","recreate");
//   TFile *newfile = new TFile("ROOTFiles/SIMPs_QCD_HT300To500_PUMoriond17.root","recreate");
  TDirectory *dir1 = newfile->mkdir("treeCorr");
  dir1->cd();
  TTree *newtreeCorr = oldtreeCorr->CloneTree(0);
  TDirectory *dir2 = newfile->mkdir("treeSPV");
  dir2->cd();
  TTree *newtreeSPV = oldtreeSPV->CloneTree(0);

  for (Long64_t i=0;i<nentries; i++) {

    c_treeCorr.GetEntry(i);
    c_treeSPV.GetEntry(i);
    if(i%100000 == 0) cout<<i<<" / "<<nentries<<"\r"<<flush;
      ///Object Selection

    double deltajet_phi = c_treeCorr.jet_phi[0] - c_treeCorr.jet_phi[1];
    if(deltajet_phi > TMath::Pi()) deltajet_phi -= 2*TMath::Pi();
    if(deltajet_phi < -TMath::Pi()) deltajet_phi += 2*TMath::Pi();

    if(c_treeCorr.jet_pt[0]>pt_cut && c_treeCorr.jet_pt[1]>pt_cut && fabs(c_treeCorr.jet_eta[0])<eta_cut && fabs(c_treeCorr.jet_eta[1])<eta_cut && fabs(deltajet_phi)>dphi_cut){
      newtreeCorr->Fill();
      newtreeSPV->Fill();
    }
  }

  newtreeCorr->Print();
  newtreeSPV->Print();
  newfile->Write("", TObject::kOverwrite);
  newfile->Close();
   //delete oldfile;
   //delete newfile;

  return 0;
}
