#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include <sys/stat.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"

//Include BTagCalibrationStandalone class === Incomplete ===


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "csvEventVars.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;          

//*****************************************************************************

void csvReweight(bool isHF=1, int insample=1, int reportEvery=100000) {

  int maxNentries=-1;
  int Njobs=1;
  int jobN=1;
  double intLumi= 36460.;
  bool verbose = false;	
  string taggerName = "csv";
  int iSys = 0;

 //=================================================================================================================================
 //Load the input root file(s)
 //=================================================================================================================================
  
  std::cout << "===> Loading the input root file(s)..." << std::endl;
  
  double mySample_xSec_ = 1.;
  double mySample_nGen_ = 1.;
  std::string mySample_sampleName_ = "default";
  std::string mySample_inputDir_ = "/afs/cern.ch/work/s/spmondal/public/BTVCSVReweight/Samples/";
  
  //******MC*****
  if( insample==2500 ){
    mySample_xSec_ = 87.31;
    mySample_nGen_ = 38139240;
    mySample_sampleName_ = "ttjets";//"TTJets";
  }
  else if( insample==2300 ){
    mySample_xSec_ = 3*2008.4;
    mySample_nGen_ = 96658928;
    mySample_sampleName_ = "zjets";//"DYJetsToLL";
  }
  else if( insample==2310 ){
    mySample_xSec_ = 18610;
    mySample_nGen_ = 35291552;
    mySample_sampleName_ = "lowMasszjets";
  }
  else if( insample==2514 ){
    mySample_xSec_ = 35.6;  
    mySample_nGen_ = 992024; 
    mySample_sampleName_ = "singletW";
  }
  else if( insample==2515 ){
    mySample_xSec_ = 35.6;  
    mySample_nGen_ = 998276; 
    mySample_sampleName_ = "singletbarW";//"Tbar_tW_DR";
  }
  else if( insample==2600 ){
    mySample_xSec_ =  12.178; 
    mySample_nGen_ = 1999000;
    mySample_sampleName_ = "WW";
  }
  
  //*****Data*****
  else if( insample==-100 ){
    mySample_xSec_ = 1; 
    mySample_nGen_ = 1; 
    mySample_sampleName_ = "DoubleEG";
  }
  else if( insample==-200 ){
    mySample_xSec_ = 1; 
    mySample_nGen_ = 1; 
    mySample_sampleName_ = "DoubleMuon";
  }
  else if( insample==-300 ){
    mySample_xSec_ = 1; 
    mySample_nGen_ = 1; 
    mySample_sampleName_ = "MuonEG";
  }

  std::string treefilename = mySample_inputDir_ + mySample_sampleName_ + "*.root";
  std::cout << "\tInput treefilename = " << treefilename.c_str() << std::endl;

  TChain *chain = new TChain("ttHTreeMaker/worldTree");
  chain->Add(treefilename.c_str());

  EventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  //Load the PU reweighting file
  
  TFile* f_PUwgt = new TFile ("PileUPReweighting.root");
  TH1D* h_PU = (TH1D*)f_PUwgt->Get("numPVs_PUratio")->Clone();
  
 //=================================================================================================================================
 //Prepare Output .root file and Output Histograms
 //=================================================================================================================================
 
  std::cout << "===> Preparing the output .root file..." << std::endl;
  std::string s_end = "_histo_All.root"; 
  
  const char* dirprefix = "Output/";
  
  struct stat st;
  if( stat(dirprefix,&st) != 0 )  mkdir(dirprefix,0777);

  std::string histofilename = Form("%s%s_rwt_hf_%s%s", dirprefix, taggerName.c_str(), mySample_sampleName_.c_str(), s_end.c_str());
  if( !isHF ) histofilename = Form("%s%s_rwt_lf_%s%s", dirprefix, taggerName.c_str(), mySample_sampleName_.c_str(), s_end.c_str());
  
  std::cout << "\tOutput histofilename = " << histofilename.c_str() << std::endl;
 
  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();

  TH1::SetDefaultSumw2();
  double maxPt1 = 300. , maxPt2 = 200.;
  if (!isHF) {maxPt1 = 200. ; maxPt2 = 100. ;}
  TH1D* h_nJets_noSF  = new TH1D("h_nJets_noSF",";numJet", 10, 0, 10 );
  TH1D* h_nJets  = new TH1D("h_nJets",";numJet", 10, 0, 10 );

  int nBinsBTag = 102;
  double xMinBTag = -0.01;

  TString btags = taggerName;
  btags.ToUpper();

  TH1D* h_all_jet_pt = new TH1D("h_all_jet_pt","; all jet p_{T}", 100, 0, maxPt1 );
  
  TH1D* h_all_jet_csv = new TH1D("h_all_jet_csv",";all jet "+btags, nBinsBTag, xMinBTag, 1.01 );
  TH1D* h_all_jet_csv_noSF = new TH1D("h_all_jet_csv_noSF",";all jet "+btags, nBinsBTag, xMinBTag, 1.01 );
  
  TH1D* h_first_jet_csv = new TH1D("h_first_jet_csv",";first jet "+btags, nBinsBTag, xMinBTag, 1.01 );
  TH1D* h_first_jet_csv_noSF = new TH1D("h_first_jet_csv_noSF",";first jet "+btags, nBinsBTag, xMinBTag, 1.01 );

  TH1D* h_second_jet_csv = new TH1D("h_second_jet_csv",";second jet "+btags, nBinsBTag, xMinBTag, 1.01 );
  TH1D* h_second_jet_csv_noSF = new TH1D("h_second_jet_csv_noSF",";second jet "+btags, nBinsBTag, xMinBTag, 1.01 );
  
  int NumCutsHF = 10;
  TH1D* h_hf_event_selection  = new TH1D("h_hf_event_selection",";cut", NumCutsHF, 0, NumCutsHF );

  h_hf_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_hf_event_selection->GetXaxis()->SetBinLabel(2,"==2 jets");
  h_hf_event_selection->GetXaxis()->SetBinLabel(3,"Dilepton trigger");
  h_hf_event_selection->GetXaxis()->SetBinLabel(4,"==2 leptons");
  h_hf_event_selection->GetXaxis()->SetBinLabel(5,"Opposite charge");
  h_hf_event_selection->GetXaxis()->SetBinLabel(6,"#Delta R(lep,lep) > 0.2");
  h_hf_event_selection->GetXaxis()->SetBinLabel(7,"M(lep,lep) > 12");
  h_hf_event_selection->GetXaxis()->SetBinLabel(8,"Zmass window");
  h_hf_event_selection->GetXaxis()->SetBinLabel(9,"MET > 30");
  h_hf_event_selection->GetXaxis()->SetBinLabel(10,"jet passes medium b-tag");

  int NumCutsLF = 10;
  TH1D* h_lf_event_selection  = new TH1D("h_lf_event_selection",";cut", NumCutsLF, 0, NumCutsLF );

  h_lf_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_lf_event_selection->GetXaxis()->SetBinLabel(2,"==2 jets");
  h_lf_event_selection->GetXaxis()->SetBinLabel(3,"Dilepton trigger");
  h_lf_event_selection->GetXaxis()->SetBinLabel(4,"==2 leptons");
  h_lf_event_selection->GetXaxis()->SetBinLabel(5,"Opposite charge");
  h_lf_event_selection->GetXaxis()->SetBinLabel(6,"#Delta R(lep,lep) > 0.2");
  h_lf_event_selection->GetXaxis()->SetBinLabel(7,"M(lep,lep) > 12");
  h_lf_event_selection->GetXaxis()->SetBinLabel(8,"Zmass window");
  h_lf_event_selection->GetXaxis()->SetBinLabel(9,"MET < 30");
  h_lf_event_selection->GetXaxis()->SetBinLabel(10,"jet fails loose b-tag");
  
    
 //=================================================================================================================================
 //Load the .csv file containing the SFs using BTagCalibration Standalone tool === Incomplete ===
 //=================================================================================================================================
  
  
  // <Lines of code to initialize the BTagCalibration Standalone tool>
  
  
  //Working points of CSV:
  
  double Mwp = 0.8484; //Medium WP
  double Lwp = 0.5426; //Loose WP
  

 //=================================================================================================================================
 //Event Loop
 //=================================================================================================================================
  int numEvents_all=0;
  int numEvents_2jets=0;

  int numEvents_lepselection2=0;
  int numEvents_lepselection1a=0;
  int numEvents_lepselection1b=0;
  int numEvents_lepselection1c=0;

  int numEvents_exselection=0;

  int nentries = chain->GetEntries();
  std::cout << "\n\tNumber of entries = " << nentries << std::endl;
  std::cout << "\tMax number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;

  int cnt = 0;
  int nPass = 0;
  std::cout << "========  Starting Event Loop  ========" << std::endl;
  
  
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) { 	//Event loop starts
    cnt++;
    if( ievt<firstEvent ) continue;
    if( ievt==lastEvent ) break;

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%reportEvery==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					     << int(double(ievt-firstEvent)/double(NeventsPerJob)*100) << "% done" << std::endl;

    if( ievt==(maxNentries+1) && ievt!=0 ) break;

    chain->GetEntry(ievt);
    numEvents_all++;

    double Xsec = mySample_xSec_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    double lumi = ( intLumi > 0 ) ? intLumi : eve->wgt_lumi_ ;
    if (insample < 0) lumi = 1;
    double wgt_gen = ( insample==2300 || insample==2310) ? (eve->wgt_generator_/fabs(eve->wgt_generator_) ): 1;
    double wgt = wgt_gen * lumi * (Xsec/nGen);//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights


   //=================================================================================================================================
   //Trigger and PU Reweighting
   //=================================================================================================================================

    // Trigger Efficiency weights
    int TwoMuon = eve->TwoMuon_;
    int TwoElectron = eve->TwoElectron_ ;
    int MuonElectron = eve->MuonElectron_ ;

    double triggerWgt = 1;
    if(insample >= 0){
      if(TwoMuon)           triggerWgt = (isHF) ? 0.962157 : 1.08586;
      else if(TwoElectron)  triggerWgt = (isHF) ? 0.939806 : 1.05302;
      else if(MuonElectron) triggerWgt = 0.916678;
	}
    wgt *= triggerWgt;

    // PU wgt
    int numPVs = eve->numTruePV_ ; // official PU recipe, using numTruePV for MC
    double PUwgt = (insample < 0) ? 1 : h_PU->GetBinContent(h_PU->FindBin(numPVs));
    wgt *= PUwgt; 

   //=================================================================================================================================
   //CSV Reweighting
   //=================================================================================================================================

    vecTLorentzVector jet_vect_TLV_tmp = eve->jet_vect_TLV_[iSys];
    vdouble jet_CSV_tmp;
    jet_CSV_tmp = eve->jet_CSV_[iSys];    

    vdouble jet_CMVA_tmp = eve->jet_cMVA_[iSys];

    vint jet_flavour_tmp = eve->jet_flavour_[iSys];
    vint jet_partonflavour_tmp = eve->jet_partonflavour_[iSys];

    vint jet_PUID_passWPLoose_tmp = eve->jet_PUID_passWPLoose_[iSys];

    vecTLorentzVector jet_vect_TLV;
    vdouble jet_CSV;
    vdouble jet_CMVA;
    vint jet_flavour;
    vint jet_partonflavour;
    vint jet_PUID_passWPLoose;
    for( int iJet=0; iJet<int(jet_vect_TLV_tmp.size()); iJet++ ){
      TLorentzVector myJet = jet_vect_TLV_tmp[iJet];

      double myCSV = jet_CSV_tmp[iJet];
      double myCMVA = jet_CMVA_tmp[iJet];
      double myJetPt = myJet.Pt();
      int myFlavor = jet_flavour_tmp[iJet];
      int mypartonFlavor = jet_partonflavour_tmp[iJet];
      int myPUID = jet_PUID_passWPLoose_tmp[iJet];

      if(!myPUID) continue;

      jet_vect_TLV.push_back(myJet);
      jet_CSV.push_back(myCSV);
      jet_CMVA.push_back(myCMVA);
      jet_flavour.push_back(myFlavor);
      jet_partonflavour.push_back(mypartonFlavor);
      jet_PUID_passWPLoose.push_back(myPUID);

    }
	
	double newCSVwgt = 1.;
	
    // <Script to implement BTagCalibrationReader> === Incomplete ===
   
    wgt *= newCSVwgt;

    if (insample < 0) wgt = 1;
    

   //=================================================================================================================================
   //Selections
   //=================================================================================================================================
    
	bool inclusiveSelection = true;
    bool emuOnlyHF = false;   // using MuonEG events only for HF?
    bool tpj =  false;  // tag and probe
	
    h_hf_event_selection->Fill(0.5, wgt);
    h_lf_event_selection->Fill(0.5, wgt);


    //////------- exactly 2 jets -----
    int numJets = int(jet_vect_TLV.size()) ;
    if (tpj || inclusiveSelection){
      if (numJets < 2) continue;
    }
    else{
      if (numJets != 2) continue; //// loosen nJets cut
    }
    numEvents_2jets++;

    h_hf_event_selection->Fill(1.5, wgt);
    h_lf_event_selection->Fill(1.5, wgt);

    double jet1_btag = jet_CSV[0];
    double jet2_btag = jet_CSV[1];

    bool passTightBtag = false;
    bool failLooseBtag = false;
    if( jet1_btag>Mwp || jet2_btag>Mwp ) passTightBtag = true;
    if( jet1_btag<Lwp || jet2_btag<Lwp ) failLooseBtag = true;

    double MHT = eve->MHT_[iSys];
    double met_pt = eve->MET_[iSys];

    bool PassZmask = eve->PassZmask_ ;
    double mass_leplep = eve->mass_leplep_;
    double dR_leplep = eve->dR_leplep_;
    int oppositeLepCharge = eve->oppositeLepCharge_;

    vecTLorentzVector lepton_vect_TLV = eve->lepton_vect_TLV_;
    vint lepton_trkCharge = eve->lepton_trkCharge_;
    if (tpj){
      if (! (lepton_vect_TLV[0].Pt()>20 && lepton_vect_TLV[1].Pt()>20) ) continue;
      if (mass_leplep < 20 ) continue;    
    }

    // Triggers 
    bool isDoubleMuTriggerPass = eve->passDoubleMuonTriggerNew_;
    bool isDoubleElectronTriggerPass = eve->passDoubleElectronTriggerNew_;
    bool isMuEGTriggerPass = eve->passElectronMuonTriggerNew_;

    bool lepselection1a = ( TwoMuon && isDoubleMuTriggerPass && abs(mass_leplep-91)>10 && (met_pt>30) ); //Selection for TwoMuon data events
    bool lepselection1b = ( TwoElectron && isDoubleElectronTriggerPass && abs(mass_leplep-91)>10 && (met_pt>30) ); //Selection for TwoEle data events

    bool lepselection1c = ( MuonElectron && isMuEGTriggerPass ); //Selection for MuonEle data events
    if (!isHF){
      lepselection1a = ( TwoMuon && isDoubleMuTriggerPass && (PassZmask==0) && (met_pt<30) && abs(mass_leplep-91)<10 ); //Selection for TwoMuon data events
      lepselection1b = ( TwoElectron && isDoubleElectronTriggerPass && (PassZmask==0) && (met_pt<30) && abs(mass_leplep-91)<10 ); //Selection for TwoEle data events
      lepselection1c = 0; //Selection for MuonEle data events
    }
    // for MC events
    bool lepselection2 = ( lepselection1a || lepselection1b || lepselection1c ) ;
    if (tpj || emuOnlyHF ) {
      if(isHF) lepselection2 = lepselection1c; // tpj, checking different dilepton categories
      else   lepselection2 = (lepselection1a || lepselection1b);
    }
    // different lepton flavor

    if ( insample == -200 ) lepselection2 = lepselection1a;
    if ( insample == -100 ) lepselection2 = lepselection1b;
    if ( insample == -300 ) lepselection2 = lepselection1c;

    if( lepselection1a ) numEvents_lepselection1a++;
    if( lepselection1b ) numEvents_lepselection1b++;
    if( lepselection1c ) numEvents_lepselection1c++;


    bool isCleanEvent = 1;
    bool firstGoodPV = eve->GoodFirstPV_;
    bool exselection = ((firstGoodPV) && (dR_leplep > 0.2) && (mass_leplep > 12) && (isCleanEvent == 1) && (oppositeLepCharge == 1)); //General dilepton selection   


    // trigger
    if(firstGoodPV){
    if( isDoubleMuTriggerPass || isDoubleElectronTriggerPass || isMuEGTriggerPass ){
      h_hf_event_selection->Fill(2.5, wgt);
      if( TwoMuon || TwoElectron || MuonElectron ){
	h_hf_event_selection->Fill(3.5, wgt);
	if( oppositeLepCharge == 1 ){
	  h_hf_event_selection->Fill(4.5, wgt);
	  if( (dR_leplep > 0.2) ){
	    h_hf_event_selection->Fill(5.5, wgt);
	    if( mass_leplep > 12 ){
	      h_hf_event_selection->Fill(6.5, wgt);
	      if( fabs(mass_leplep-91)>10 ){
		h_hf_event_selection->Fill(7.5, wgt);
		if( met_pt>30 ){
		  h_hf_event_selection->Fill(8.5, wgt);
		  if( passTightBtag ){
		    h_hf_event_selection->Fill(9.5, wgt);
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    if( isDoubleMuTriggerPass || isDoubleElectronTriggerPass ){
      h_lf_event_selection->Fill(2.5, wgt);
      if( TwoMuon || TwoElectron ){
	h_lf_event_selection->Fill(3.5, wgt);
	if( oppositeLepCharge == 1 ){
	  h_lf_event_selection->Fill(4.5, wgt);
	  if( (dR_leplep > 0.2) ){
	    h_lf_event_selection->Fill(5.5, wgt);
	    if( mass_leplep > 12 ){
	      h_lf_event_selection->Fill(6.5, wgt);
	      if( PassZmask==0 && fabs(mass_leplep-91)<10 ){
		h_lf_event_selection->Fill(7.5, wgt);
		if( met_pt<30 ){
		  h_lf_event_selection->Fill(8.5, wgt);
		  if( failLooseBtag ){
		    h_lf_event_selection->Fill(9.5, wgt);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    }
    
    if ( !lepselection2 ) continue;
    numEvents_lepselection2++;


    //------ extra seletions -----
    if ( !exselection ) continue;
    numEvents_exselection++;

    nPass++;    
    
   //=================================================================================================================================
   //Populate Output Histograms
   //=================================================================================================================================
    
    double first_jet_pt = -9.;   
    double first_jet_eta = -9.;
    double first_jet_csv = -9.; 
    int first_jet_flavour = -9;  
    int first_jet_partonflavour = -9; 
    

    double second_jet_pt = -9.;
    double second_jet_eta = -9.;
    double second_jet_csv = -9.;
    int    second_jet_flavour = -9;
    int    second_jet_partonflavour = -9; 

    if( verbose ) std::cout << "--for event  "<< cnt << std::endl;    
    if( verbose ) std::cout << " --number of jets is "<< numJets << std::endl;  
    
    int numTags = 0;
    int numJets30 = 0;
    
    for( int iJet=0; iJet<int(jet_vect_TLV.size()); iJet++ ){
      TLorentzVector myJet = jet_vect_TLV[iJet];     

      double myCSV = jet_CSV[iJet];
      double myJetPt = myJet.Pt();
      double myJetEta = myJet.Eta();
      int myFlavor = jet_flavour[iJet];
      int mypartonFlavor = jet_partonflavour[iJet];

      //tpj
      if(myJetPt > 30) numJets30++;
      else if(tpj) continue;

      h_all_jet_pt->Fill(myJetPt, wgt);
      h_all_jet_csv->Fill(myCSV, wgt);
      h_all_jet_csv_noSF->Fill(myCSV, wgt/newCSVwgt);

      if (myCSV > Mwp) numTags++;

      if( (tpj && numJets30==1) || (!tpj && iJet==0) ){ //tpj
		first_jet_pt = myJetPt;
		first_jet_eta = myJetEta;
		first_jet_csv = myCSV;
		first_jet_flavour = myFlavor;
		first_jet_partonflavour = mypartonFlavor;

	  }
		  
 	  if( (tpj && numJets30==2) || (!tpj && iJet==1) ){ //tpj
		second_jet_pt = myJetPt;
		second_jet_eta = myJetEta;
		second_jet_csv = myCSV;
		second_jet_flavour = myFlavor;
		second_jet_partonflavour = mypartonFlavor;

      }
    } // end loop over jets
    
    if (tpj && numJets30 != 2) continue; //tpj

    if( verbose ) std::cout << "   -first jet pt, eta is " << first_jet_pt << ";  "<<first_jet_eta << std::endl;
    if( verbose ) std::cout << "   -second jet pt, eta is " << second_jet_pt << ";  "<<second_jet_eta << std::endl;

    h_nJets->Fill(numJets, wgt);
    h_nJets_noSF->Fill(numJets, wgt/newCSVwgt);

    h_first_jet_csv->Fill(first_jet_csv, wgt);
    h_first_jet_csv_noSF->Fill(first_jet_csv, wgt/newCSVwgt);

    h_second_jet_csv->Fill(second_jet_csv, wgt);
    h_second_jet_csv_noSF->Fill(second_jet_csv, wgt/newCSVwgt);

  } // end loop over events

  
  //=================================================================================================================================
  //Report
  //=================================================================================================================================
  std::cout << "\nTotal number of selected events is " << nPass << std::endl;

  std::cout << "===========================================" << std::endl;
  std::cout << "  Number of all events = " << numEvents_all << std::endl;
  std::cout << "  Number of events with == 2 jets      = " << numEvents_2jets << std::endl;
  std::cout << "  Number of events with lepselection2  = " << numEvents_lepselection2 << std::endl;
  std::cout << "  Number of events with lepselection1a = " << numEvents_lepselection1a << std::endl;
  std::cout << "  Number of events with lepselection1b = " << numEvents_lepselection1b << std::endl;
  std::cout << "  Number of events with lepselection1c = " << numEvents_lepselection1c << std::endl;
  std::cout << "  Number of events with exselection    = " << numEvents_exselection << std::endl;
  std::cout << "===========================================" << std::endl;


  std::cout << "Done!\n" << std::endl;

  histofile.Write();
  histofile.Close();

}
