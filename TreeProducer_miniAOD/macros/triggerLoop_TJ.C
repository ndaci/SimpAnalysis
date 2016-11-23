#include "myIncludes.h"
#include "myMath.h"

bool DEBUG=false;

Int_t quick(TString dir="v0_test", TString UsrCut="Loose",
	    TString era="2016B", UInt_t run1=274748, UInt_t run2=999999,
	    Int_t nRequest=-1)
{

  // Define input tree //////////////////////////////////////////////
  // Input dir
  TString dirIn="/user/ndaci/Data/XMET/TriggerTrees_07Oct2016/";

  // Define TChain
  cout << "- Define TChain for era " << era << " ";
  TChain* chain;
  chain = new TChain("tree/tree");
  //
  if(era.Contains("All")) chain->Add(dirIn+"/ntuple_SingleMuon_*.root");
  else                    chain->Add(dirIn+"/ntuple_SingleMuon_"+era+".root");
  cout << " Done!" << endl;

  // Branch variables
  uint32_t event, run, lumi;  
  uint32_t nvtx;
  uint8_t hltJet170CF, hltDiJet170, hltDiJet170CF, hltDiJet220CF, hltDiJet330CF, hltDiJet430;
  double  ps_DiJet170, ps_Jet170Chf0p1;

  std::vector<double> *combinejetpt      = new vector<double> ();
  std::vector<double> *combinejeteta     = new vector<double> ();
  std::vector<double> *combinejetphi     = new vector<double> ();
  std::vector<double> *combinejetCHfrac  = new vector<double> ();
  std::vector<double> *combinejetNHfrac  = new vector<double> ();
  std::vector<double> *combinejetEMfrac  = new vector<double> ();
  std::vector<double> *combinejetCEMfrac = new vector<double> ();

  // Branches
  chain->SetBranchAddress("event"                , &event                ); // , "event/i");
  chain->SetBranchAddress("run"                  , &run                  ); // , "run/i");
  chain->SetBranchAddress("lumi"                 , &lumi                 ); // , "lumi/i");
  chain->SetBranchAddress("nvtx"                 , &nvtx                 ); // , "nvtx/i");
  chain->SetBranchAddress("hltJet170CF"  , &hltJet170CF); // , "hltJet170CF/b");
  chain->SetBranchAddress("hltDiJet170"  , &hltDiJet170); // , "hltDiJet170/b");
  chain->SetBranchAddress("hltDiJet170CF", &hltDiJet170CF); // , "hltDiJet170CF/b");
  chain->SetBranchAddress("hltDiJet220CF", &hltDiJet220CF); // , "hltDiJet220CF/b");
  chain->SetBranchAddress("hltDiJet330CF", &hltDiJet330CF); // , "hltDiJet330CF/b");
  chain->SetBranchAddress("hltDiJet430"  , &hltDiJet430); // , "hltDiJet430/b");
  chain->SetBranchAddress("pswgt_DiJet170"          , &ps_DiJet170);
  chain->SetBranchAddress("pswgt_Jet170Chf0p1"      , &ps_Jet170Chf0p1);
  chain->SetBranchAddress("combinejetpt" ,     &combinejetpt);
  chain->SetBranchAddress("combinejeteta",     &combinejeteta);
  chain->SetBranchAddress("combinejetphi",     &combinejetphi);
  chain->SetBranchAddress("combinejetCHfrac",  &combinejetCHfrac);
  chain->SetBranchAddress("combinejetNHfrac",  &combinejetNHfrac);
  chain->SetBranchAddress("combinejetEMfrac",  &combinejetEMfrac);
  chain->SetBranchAddress("combinejetCEMfrac", &combinejetCEMfrac);

  ////////////////////////////////////////////////////////////////////

  // Define HLT paths
  const UInt_t nP=5;
  TString nameP[  nP]={"hltDiJet170","hltDiJet170CF","hltDiJet220CF","hltDiJet330CF","hltDiJet430"};
  
  // Define x-axis variables
  const UInt_t nV=5;  

  TString nameVar[nV] = {"leadingjetpt", "leadingjetChF","secondjetpt", "secondjetChF",
			 "highestCHF"};

  TString nameAxis[  nV] = {"Leading jet p_{T} [GeV]", "Leading jet ChF"
			    "Second jet p_{T} [GeV]" , "Second jet ChF",
			    "Max ChF"};

  const UInt_t nBinsPt=25;
  Float_t bins_pt[nBinsPt]  = { 50, 100, 120, 140, 160, 180, 
			       200, 220, 240, 260, 280, 300,
			       320, 340, 360, 380, 400, 440,
			       480, 520, 560, 600, 700, 800,
			       1000};

  const UInt_t nBinsChf=27;
  float bins_chf[nBinsChf] = {0.00, 0.02, 0.04, 0.06, 0.08, 
			      0.10, 0.12, 0.14, 0.16, 0.18,
			      0.20, 0.22, 0.24, 0.26, 0.28,
			      0.30, 0.35, 0.40, 0.45, 0.50,
			      0.55, 0.60, 0.65, 0.70, 0.80,
			      0.90, 1.00};

  Int_t nBins[nV]     = {nBinsPt, nBinsChf, nBinsPt, nBinsChf, nBinsChf};
  Float_t* v_xlow[nV] = {bins_pt, bins_chf, bins_pt, bins_chf, bins_chf};

  /* // Regular binning
  Int_t   nBins[  nV] = {100};
  Float_t xmin[   nV] = {0};
  Float_t xmax[   nV] = {1000};
  */

  // Output dir
  TString dirOut="/user/ndaci/Results/TracklessJets/Trigger/Efficiency/CMSSW_808/";
  
  // Declare output histograms
  cout << "- Declare output histograms" << endl;
  TH1F *hNum[nP][nV];
  TH1F *hDen[nP][nV];
  TEfficiency *teff[nP][nV];

  // Define histograms
  cout << "- Define output histograms" << endl;
  TString name_hNum[nP][nV];
  TString name_hDen[nP][nV];
  TString title_hNum[nP][nV];
  TString title_hDen[nP][nV];
  TString title_teff[nP][nV];
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    
    cout << "-- Era: " << era << " Variable: " << nameVar[iV] << endl;
    
    // Numerator histos: 1 per path #iP
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      
      cout << "--- Path: " << nameP[iP] << endl;
      
      // Denominator histo for requested era and variable #iV
      name_hDen[iP][iV]  = "hDen_"+nameVar[iV]+"_"+era+"_"+nameP[iP];
      title_hDen[iP][iV] = "Denominator: "+nameVar[iV]+" "+era+"_"+nameP[iP];
      //
      //hDen[iP][iV] = new TH1F(name_hDen[iP][iV], title_hDen[iP][iV], nBins[iV], xmin[iV], xmax[iV]);
      hDen[iP][iV] = new TH1F(name_hDen[iP][iV], title_hDen[iP][iV], nBins[iV]-1, v_xlow[iV]);
      //
      cout << "-- Denominator: " << name_hDen[iP][iV] << endl;
      
      name_hNum[ iP][iV]="hNum_"+nameVar[iV]+"_"+era+"_"+nameP[iP];
      title_hNum[iP][iV]="Numerator: "+nameVar[iV]+" "+era+" "+nameP[iP];
      title_teff[iP][iV]="Efficiency: "+nameP[iP]+" ("+era+")";
      //
      hNum[iP][iV] = new TH1F(name_hNum[iP][iV], title_hNum[iP][iV], nBins[iV]-1, v_xlow[iV]);
      //
      cout << "--- Numerator: " << name_hNum[iP][iV] << endl;
    }
  }

  // Define 2D Histograms
  const UInt_t nV2D = 3;
  TString nameVar2D[   nV2D] = {"jetpt", "jetChF", "extremaChF"};
  TString nameAxisX_2D[nV2D] = {"Leading jet p_{T} [GeV]", "Leading jet ChF", "Min jet ChF"};
  TString nameAxisY_2D[nV2D] = {"Second jet p_{T} [GeV]", "Second jet ChF", "Max jet ChF"};

  // Binning
  Int_t nBins2D[nV2D]      = {nBinsPt, nBinsChf, nBinsChf};
  Float_t* v_2D_xlow[nV2D] = {bins_pt, bins_chf, bins_chf};
  Float_t* v_2D_ylow[nV2D] = {bins_pt, bins_chf, bins_chf};

  // Declare output histograms
  cout << "- Declare output histograms" << endl;
  TH2F *hNum_2D[nP][nV2D];
  TH2F *hDen_2D[nP][nV2D];
  TH2F *hEff_2D[nP][nV2D];

  // Define histograms
  cout << "- Define output histograms" << endl;
  TString name_hNum_2D[ nP][nV2D];
  TString name_hDen_2D[ nP][nV2D];
  TString name_hEff_2D[nP][nV2D];
  TString title_hNum_2D[nP][nV2D];
  TString title_hDen_2D[nP][nV2D];
  TString title_hEff_2D[nP][nV2D];
  //
  for(UInt_t iV=0 ; iV<nV2D ; iV++) {
    
    cout << "-- Era: " << era << " Variable: " << nameVar2D[iV] << endl;
    
    // Numerator histos: 1 per path #iP
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      
      cout << "--- Path: " << nameP[iP] << endl;
      
      // Denominator histo for requested era and variable #iV
      name_hDen_2D[iP][iV]  = "hDen_2D_"+nameVar2D[iV]+"_"+era+"_"+nameP[iP];
      title_hDen_2D[iP][iV] = "Denominator: "+nameVar2D[iV]+" "+era+"_"+nameP[iP];
      //
      hDen_2D[iP][iV] = new TH2F(name_hDen_2D[iP][iV], title_hDen_2D[iP][iV], 
				 nBins2D[iV]-1, v_2D_xlow[iV], nBins2D[iV]-1, v_2D_ylow[iV]);
      //
      cout << "-- Denominator: " << name_hDen_2D[iP][iV] << endl;

      // Numerator histo      
      name_hNum_2D[ iP][iV]="hNum_2D_"+nameVar2D[iV]+"_"+era+"_"+nameP[iP];
      title_hNum_2D[iP][iV]="Numerator: "+nameVar2D[iV]+" "+era+" "+nameP[iP];
      //
      name_hEff_2D[iP][iV]="hEff_2D_"+nameVar2D[iV]+"_"+era+"_"+nameP[iP];
      title_hEff_2D[iP][iV]="Efficiency: "+nameP[iP]+" ("+era+")";
      //
      hNum_2D[iP][iV] = new TH2F(name_hNum_2D[iP][iV], title_hNum_2D[iP][iV], 
				 nBins2D[iV]-1, v_2D_xlow[iV], nBins2D[iV]-1, v_2D_ylow[iV]);

      hEff_2D[iP][iV] = new TH2F(name_hEff_2D[iP][iV], title_hEff_2D[iP][iV], 
				 nBins2D[iV]-1, v_2D_xlow[iV], nBins2D[iV]-1, v_2D_ylow[iV]);
      //
      hNum_2D[iP][iV]->SetXTitle( nameAxisX_2D[iV] );
      hNum_2D[iP][iV]->SetYTitle( nameAxisY_2D[iV] );
      //
      hEff_2D[iP][iV]->SetXTitle( nameAxisX_2D[iV] );
      hEff_2D[iP][iV]->SetYTitle( nameAxisY_2D[iV] );
      //
      cout << "--- Numerator: " << name_hNum_2D[iP][iV] << endl;
    }
  }

  //////// MAIN LOOP ///////
  // Prepare utility variables
  bool   print;
  double jetpt1, jetphi1, jetchf1;
  double jetpt2, jetphi2, jetchf2;
  double chf_max, chf_min;
  double pswgt;

  // Loop over chain entries //
  UInt_t nEntries = chain->GetEntries();
  UInt_t nProcess = nEntries;
  if( nRequest > 0 && (UInt_t)nRequest < nProcess ) nProcess = (UInt_t) nRequest;

  cout << "Start loop over: " << nEntries << " entries" << endl;
  //
  for(UInt_t iE=0 ; iE<nProcess ; iE++) {

    // Printout boolean
    if(iE%10000==0) print=true;
    else            print=false;

    // Printout current entry
    if(print) cout << "- Start processing entry # " << iE << " / " << nProcess << endl;

    // Initialize
    jetpt1 = jetphi1 = jetchf1 = 0;
    jetpt2 = jetphi2 = jetchf2 = 0;
    chf_max = chf_min = 0;
    event = run = lumi = nvtx = 0;
    hltJet170CF = hltDiJet170 = hltDiJet170CF = hltDiJet220CF = hltDiJet330CF = hltDiJet430 = 0;
    ps_DiJet170 = ps_Jet170Chf0p1 = pswgt = 1;
    //
    combinejetpt      ->clear();
    combinejeteta     ->clear();
    combinejetphi     ->clear();
    combinejetCHfrac  ->clear();
    combinejetNHfrac  ->clear();
    combinejetEMfrac  ->clear();
    combinejetCEMfrac ->clear();

    // Get entry iE    
    chain->GetEntry(iE);
    
    // Use utility variables
    if(combinejetpt->size()>0) {
      jetpt1  = (*combinejetpt )[0];
      jetphi1 = (*combinejetphi)[0];
      jetchf1 = (*combinejetCHfrac)[0]+(*combinejetCEMfrac)[0];
      chf_min = chf_max = jetchf1;
    }
    //
    if(combinejetpt->size()>1) {
      jetpt2  = (*combinejetpt )[1];
      jetphi2 = (*combinejetphi)[1];
      jetchf2 = (*combinejetCHfrac)[1]+(*combinejetCEMfrac)[1];
      chf_min = TMath::Min( jetchf1 , jetchf2 );
      chf_max = TMath::Max( jetchf1 , jetchf2 );
    }

    if(print) cout << "-- Check ChF: chf_min=" << chf_min 
		   << " chf_max=" << chf_max
		   << " (" << combinejetpt->size() << " jets)" << endl;
  
    // Denominator selection //
    // DeltaPhi >= 2 (back-to-back dijet)
    if( TMath::Abs( MyMath::DeltaPhi(jetphi1,jetphi2) ) < 2.0 ) continue;
    if( run<274748 ) continue; // triggers start to exist
    if( run<run1   ) continue; // input argument run1
    if( run>run2   ) continue; // input argument run2
    
    // Path/Variable-specific denominators
    bool cut_DenPath[nV][nP];
    //
    bool cut_Loose[nV][nP]={ {true, 
			      (jetchf1)<0.1 && (jetchf2)<0.1, 
			      (jetchf1)<0.3 && (jetchf2)<0.3, 
			      (jetchf1)<0.5 && (jetchf2)<0.5, 
			      true} ,
			     {jetpt1>170 && jetpt2>170, 
			      jetpt1>170 && jetpt2>170, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>330 && jetpt2>330, 
			      jetpt1>430 && jetpt2>430} ,
			     {true, 
			      (jetchf1)<0.1 && (jetchf2)<0.1, 
			      (jetchf1)<0.3 && (jetchf2)<0.3, 
			      (jetchf1)<0.5 && (jetchf2)<0.5, 
			      true} ,
			     {jetpt1>170 && jetpt2>170, 
			      jetpt1>170 && jetpt2>170, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>330 && jetpt2>330, 
			      jetpt1>430 && jetpt2>430} ,
			     {jetpt1>170 && jetpt2>170, 
			      jetpt1>170 && jetpt2>170, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>330 && jetpt2>330, 
			      jetpt1>430 && jetpt2>430} };
    //
    bool cut_Tight[nV][nP]={ {true, 
			      (jetchf1)<0.05 && (jetchf2)<0.05, 
			      (jetchf1)<0.2 && (jetchf2)<0.2, 
			      (jetchf1)<0.4 && (jetchf2)<0.4, 
			      true} ,
			     {jetpt1>220 && jetpt2>220, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>270 && jetpt2>270, 
			      jetpt1>380 && jetpt2>380, 
			      jetpt1>480 && jetpt2>480} ,
			     {true, 
			      (jetchf1)<0.05 && (jetchf2)<0.05, 
			      (jetchf1)<0.2 && (jetchf2)<0.2, 
			      (jetchf1)<0.4 && (jetchf2)<0.4, 
			      true} ,
			     {jetpt1>220 && jetpt2>220, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>270 && jetpt2>270, 
			      jetpt1>380 && jetpt2>380, 
			      jetpt1>480 && jetpt2>480} ,
			     {jetpt1>220 && jetpt2>220, 
			      jetpt1>220 && jetpt2>220, 
			      jetpt1>270 && jetpt2>270, 
			      jetpt1>380 && jetpt2>380, 
			      jetpt1>480 && jetpt2>480} };
    //
    for(UInt_t iP=0; iP<nP; iP++) {
      for(UInt_t iV=0; iV<nV; iV++) {
	if(UsrCut=="Tight") cut_DenPath[iV][iP] = cut_Tight[iV][iP];
	else                cut_DenPath[iV][iP] = cut_Loose[iV][iP];
      }
    }

    bool cut_NumPath[nP] = { hltDiJet170>0 , hltDiJet170CF>0 , hltDiJet220CF>0 , hltDiJet330CF>0 , hltDiJet430>0 };
    
    double var1D[nV] = {jetpt1, jetchf1, jetpt2, jetchf2, chf_max };

    // Fill histograms
    if(print && DEBUG) cout << "- Fill histograms" << endl;
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      if(print && DEBUG) cout << "-- Era: " << era << " Variable: " << nameVar[iV] << endl;

      for(UInt_t iP=0 ; iP<nP ; iP++) {

	if(iP==0) pswgt = ps_DiJet170;
	else pswgt = 1;

	if( cut_DenPath[iV][iP] ) {
	  hDen[iP][iV]->Fill( var1D[iV] , 1/pswgt );
	  if( cut_DenPath[iV][iP] && cut_NumPath[iP] ) {
	    hNum[iP][iV]->Fill( var1D[iV] );
	  }
	}

      } // end loop: iP

    } // end loop: iV

    // Fill 2D histograms
    for(UInt_t iP=0 ; iP<nP ; iP++) {

      if(iP==0) {
	// if(ps_DiJet170==0) pswgt=0;
	// else               pswgt = 1/ps_DiJet170;
	pswgt = ps_DiJet170;
      }
      else pswgt = 1;

      // x/y axes variables: Jet pT
      if( cut_DenPath[0][iP] && cut_DenPath[2][iP] ) {
	hDen_2D[iP][0]->Fill( jetpt1 , jetpt2 );
	if(cut_NumPath[iP]) {
	  hNum_2D[iP][0]->Fill( jetpt1 , jetpt2 , pswgt);
	}
      }
      
      // x/y axes variables: Jet ChF
      if( cut_DenPath[1][iP] && cut_DenPath[3][iP] ) {
	hDen_2D[iP][1]->Fill( jetchf1 , jetchf2 );
	if(cut_NumPath[iP]) {
	  hNum_2D[iP][1]->Fill( jetchf1 , jetchf2 , pswgt  );
	}
      }

      // x/y axes variables: min/max ChF
      if( cut_DenPath[1][iP] && cut_DenPath[3][iP] ) {
	hDen_2D[iP][2]->Fill( chf_min , chf_max );
	if(cut_NumPath[iP]) {
	  hNum_2D[iP][2]->Fill( chf_min , chf_max , pswgt );
	}
      }

    } // end loop: iP



  } // end loop: iE (chain entries)
  ////////// END LOOP //////////////////////

  // Create TCanvas
  TCanvas c("c","c",0,0,600,600); 
  //setTDRStyle();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(5);
  gPad->Update();

  // Create TEfficiencies
  for(UInt_t iP=0 ; iP<nP ; iP++) {

    for(UInt_t iV=0 ; iV<nV ; iV++) {	

      if( TEfficiency::CheckConsistency( *hNum[iP][iV] , *hDen[iP][iV] , "w" ) ) {

	teff[iP][iV] = new TEfficiency( *hNum[iP][iV] , *hDen[iP][iV] );

	TString theName = "teff_"+name_hNum[iP][iV];
	teff[iP][iV]->SetNameTitle( theName, title_teff[iP][iV]+";"+nameAxis[iV]+";Efficiency");

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gPad->Update();
	teff[iP][iV]->Draw("AP");

	gStyle->SetOptStat(0); 
	gStyle->SetOptFit(0); 
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(5);
	gPad->Update();

	c.Print(dirOut+"/"+dir+"/pdf/"+theName+".pdf","pdf");
	c.Print(dirOut+"/"+dir+"/png/"+theName+".png","png");
	//c.Print(dirOut+"/"+dir+"/"+theName+".C","C");
	//c.Print(dirOut+"/"+dir+"/"+theName+".root","root");
      }

    } //end loop:nV
  } // end loop:nP

  // Create 2D Efficiency plots
  // Create TEfficiencies
  for(UInt_t iP=0 ; iP<nP ; iP++) {

    for(UInt_t iV=0 ; iV<nV2D ; iV++) {	

      if(!hNum_2D[iP][iV] || !hDen_2D[iP][iV]) continue;

      hEff_2D[iP][iV] = (TH2F*) hNum_2D[iP][iV]->Clone(name_hEff_2D[iP][iV]);
      hEff_2D[iP][iV] -> Divide( hDen_2D[iP][iV] );

      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      gPad->Update();
      hEff_2D[iP][iV]->Draw("colz");

      c.Print(dirOut+"/"+dir+"/pdf/"+name_hEff_2D[iP][iV]+".pdf","pdf");
      c.Print(dirOut+"/"+dir+"/png/"+name_hEff_2D[iP][iV]+".png","png");
      //c.Print(dirOut+"/"+dir+"/"+theName+".C","C");
      //c.Print(dirOut+"/"+dir+"/"+theName+".root","root");

    } //end loop:nV
  } // end loop:nP

  return 0;
}
