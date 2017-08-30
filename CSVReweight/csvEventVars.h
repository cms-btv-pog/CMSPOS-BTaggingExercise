#ifndef csvReweightingRun2_csvTreeMaker_csvEventVars_h
#define csvReweightingRun2_csvTreeMaker_csvEventVars_h

//
// Dependencies (#includes)
//
#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

using namespace std;



typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<double> vdouble;
typedef std::vector<string> vstring;
typedef std::vector<bool> vbool;
typedef std::vector<int> vint;
typedef std::vector< TLorentzVector > vecTLorentzVector;

//
// Utility Class for Handling Event Variables
//

const Int_t rNumSys = 3;

struct EventVars{


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  explicit EventVars() { }

  Int_t  PassDIL_;

  Int_t  TwoMuon_;
  Int_t  TwoElectron_;
  Int_t  MuonElectron_;

  double mass_leplep_;
  double dR_leplep_;
  Int_t  oppositeLepCharge_;

  Int_t  PassZmask_;
  Int_t  PassZmaskMET_;

  int passDoubleElectronTrigger_;
  int passDoubleMuonTrigger_;
  int passElectronMuonTrigger_;

  int passDoubleElectronTriggerNew_;
  int passDoubleMuonTriggerNew_;
  int passElectronMuonTriggerNew_;

  //  vstring TriggerPaths_ ;
  //  vint TriggerAcceps_ ;

  int run_;
  int lumi_;
  long evt_;

  int     numTruePV_;
  int     numGenPV_;
  Int_t   numPVs_;
  
  bool    GoodFirstPV_;
  

  int numTightMuons_;
  int numTightElectrons_;
  int numLooseMuons_;
  int numLooseElectrons_;
  
  
  vecTLorentzVector lepton_vect_TLV_;
  vint lepton_trkCharge_;
  vint lepton_isMuon_;


  double  wgt_generator_;
  Float_t wgt_lumi_;
  Float_t wgt_xs_;
  Float_t wgt_nGen_;
  Float_t wgt_lepSF_;
  Float_t wgt_pu_;


  double  wgt_[rNumSys];

  Float_t numJets_float_[rNumSys];
  Float_t numTags_float_[rNumSys];

  int numJets_[rNumSys];
  int numTags_[rNumSys];


  Float_t first_jet_pt_[rNumSys];
  Float_t second_jet_pt_[rNumSys];

  Float_t MET_[rNumSys];
  Float_t MET_phi_[rNumSys];
  Float_t MHT_[rNumSys];
  /* Float_t METNoHF_[rNumSys]; */
  /* Float_t METNoHF_phi_[rNumSys]; */

  vecTLorentzVector jet_vect_TLV_[rNumSys];
  vdouble jet_CSV_[rNumSys];
  vdouble jet_cMVA_[rNumSys];
  vint jet_flavour_[rNumSys];
  vint jet_partonflavour_[rNumSys];
  vdouble jet_pt_[rNumSys];
  vdouble jet_eta_[rNumSys];

  vdouble jet_PUID_mva_[rNumSys];
  vint jet_PUID_flag_[rNumSys];
  vint jet_PUID_passWPLoose_[rNumSys];

  void initialize();

};


typedef std::vector<EventVars> vEventVars;


void EventVars::initialize(){

  PassDIL_ = 0;

  TwoMuon_ = -99;
  TwoElectron_ = -99;
  MuonElectron_ = -99;

  mass_leplep_ = -99;
  dR_leplep_ = -99;
  oppositeLepCharge_ = -99;

  PassZmask_ = -99;
  PassZmaskMET_ = -99;
  
  passDoubleElectronTrigger_ = -99;
  passDoubleMuonTrigger_     = -99;
  passElectronMuonTrigger_   = -99;

  passDoubleElectronTriggerNew_ = -99;
  passDoubleMuonTriggerNew_     = -99;
  passElectronMuonTriggerNew_   = -99;
  
  //  TriggerPaths_.clear();
  //  TriggerAcceps_.clear() ;

  run_  = -99;
  lumi_ = -99;
  evt_ = -99;

  numTruePV_ = -99;
  numGenPV_ = -99;
  numPVs_ = -99;
  
  GoodFirstPV_ = false;

  
  numTightMuons_=-99;
  numTightElectrons_=-99;
  numLooseMuons_=-99;
  numLooseElectrons_=-99;


  lepton_vect_TLV_.clear();
  lepton_trkCharge_.clear();
  lepton_isMuon_.clear();


  wgt_generator_        = -99.9;
  wgt_lumi_             = -99.9;
  wgt_xs_               = -99.9;
  wgt_nGen_             = -99.9;
  wgt_lepSF_            = -99.9;
  wgt_pu_               = -99.9;
 
  for(int iSys=0; iSys<rNumSys; iSys++){

    wgt_[iSys] = -99.9;  

    numJets_float_[iSys] = -99.9;
    numTags_float_[iSys] = -99.9;

    numJets_[iSys] = -99;
    numTags_[iSys] = -99;

    first_jet_pt_[iSys]                   = -99.9;
    second_jet_pt_[iSys]                  = -99.9;

    MET_[iSys]                            = -99.9;
    MET_phi_[iSys]                        = -99.9;
    MHT_[iSys]     = -99.9;
    /* METNoHF_[iSys]                            = -99.9; */
    /* METNoHF_phi_[iSys]                        = -99.9; */

  
    jet_vect_TLV_[iSys].clear();
    jet_CSV_[iSys].clear();
    jet_cMVA_[iSys].clear();
    jet_pt_[iSys].clear();
    jet_eta_[iSys].clear();
    jet_flavour_[iSys].clear();
    jet_partonflavour_[iSys].clear();

    jet_PUID_mva_[iSys].clear();
    jet_PUID_flag_[iSys].clear();
    jet_PUID_passWPLoose_[iSys].clear();

  }


  return;
}

  

#endif
