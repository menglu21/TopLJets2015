#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExYukawa.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/CtagScaleFactorsWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "TMath.h"

//#include <filesystem>
//namespace fs = std::filesystem;

using namespace std;
float DeltaEta(float eta1, float eta2);

//
void RunExYukawa(const TString in_fname,
                      TString outname,
                      TH1F *normH,
                      TH1F *genPU,
                      TString era,
                      Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////

  bool isLowPUrun(false);
  if(in_fname.Contains("2017H")) {
    isLowPUrun=true;
    cout << "Running with low PU run settings" << endl;
  }

  //preselection cuts to apply
  float minLeptonPt( isLowPUrun ? 20. : 20);
  size_t minJetMultiplicity(3);
  int minNum_btags = 2;

  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();

  //CORRECTIONS: LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(in_fname.Contains("Data13TeV"),era,"TIGHT");

  //CORRECTIONS: L1-prefire
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);

  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

  //CORRECTIONS: C-TAG CALIBRATION
  CtagScaleFactorsWrapper CtagSFWR(in_fname.Contains("Data13TeV"),era);


  //READ Cross section over sum(weight) for the MC sample

  TString dirin = gSystem->DirName(in_fname);
  TString baseMC = gSystem->BaseName(dirin);  cout<<"Base MC name: "<<baseMC<<endl;
/*
  TString githashfromdir = gSystem->BaseName(gSystem->DirName(dirin));
  cout<<"Git hash: "<<githashfromdir<<endl;
  TFile *norm_file=TFile::Open("$CMSSW_BASE/src/TopLJets2015/TopAnalysis/data/era2017/genweights_"+githashfromdir+".root");
  norm_file->cd();
  TKey *key=norm_file->FindKey(baseMC);
  if (key==0){
   cout<<"No such histogram"<<endl;
   throw 1;
  }
  TH1F *h_norm = (TH1F*)norm_file->Get(baseMC);
  float norm=h_norm->GetBinContent(1);
  cout<<"first bin content of the normalization histogram:  "<<norm<<endl;
*/
//  norm_file->Close();


  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname);
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
////////  TFile *f_tmva=TFile::Open("tree_tmva.root","recreate");
  fOut->cd();



  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("control_Z_mass_overall", new TH1F("control_Z_mass_overall", ";M(Z) [GeV];Events",80,50,170));
  ht.addHist("control_Z_mass",    new TH1F("control_Z_mass",    ";M(Z) [GeV];Events",80,50,170));
  ht.addHist("control_lep_pt_bc", new TH1F("control_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("control_lep_eta_bc",new TH1F("control_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));
  ht.addHist("control_m_ll_bc",		new TH1F("control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("control_pT_ll_bc",  new TH1F("control_pT_ll_bc",  ";p_{T}(l+,l-) [GeV] ; Events", 80,0,200));
  ht.addHist("control_y_ll_bc",  new TH1F("control_y_ll_bc",  ";y(l+,l-) [GeV] ; Events", 20,-5,5));
  ht.addHist("control_2_lep_pt_bc", new TH1F("control_2_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("control_2_lep_eta_bc",new TH1F("control_2_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));

  ht.addHist("isgsf_control_Z_mass_overall", new TH1F("isgsf_control_Z_mass_overall", ";M(Z) [GeV];Events",80,50,170));
  ht.addHist("isgsf_control_Z_mass",    new TH1F("isgsf_control_Z_mass",    ";M(Z) [GeV];Events",80,50,170));
  ht.addHist("isgsf_control_lep_pt_bc", new TH1F("isgsf_control_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("isgsf_control_lep_eta_bc",new TH1F("isgsf_control_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));
  ht.addHist("isgsf_control_m_ll_bc",         new TH1F("isgsf_control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("isgsf_control_pT_ll_bc",  new TH1F("isgsf_control_pT_ll_bc",  ";p_{T}(l+,l-) [GeV] ; Events", 80,0,200));
  ht.addHist("isgsf_control_y_ll_bc",  new TH1F("isgsf_control_y_ll_bc",  ";y(l+,l-) [GeV] ; Events", 20,-5,5));
  ht.addHist("isgsf_control_2_lep_pt_bc", new TH1F("isgsf_control_2_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("isgsf_control_2_lep_eta_bc",new TH1F("isgsf_control_2_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));

  ht.addHist("OS_isgsf_control_m_ll_bc",         new TH1F("OS_isgsf_control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("SS_isgsf_control_m_ll_bc",         new TH1F("SS_isgsf_control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("OS_control_m_ll_bc",         new TH1F("OS_control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("SS_control_m_ll_bc",         new TH1F("SS_control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));

  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",35,0,140));

  ht.addHist("njets_bc",        new TH1F("njets_bc",       ";Jet multiplicity;Events",12,0.5,12.5));
  ht.addHist("jet_pt_bc",       new TH1F("jet_pt_bc",      ";p_{T}(jet) [GeV];Events",15,0,300));
  ht.addHist("jet_eta_bc",       new TH1F("jet_eta_bc",      ";#eta(jet) ;Events",12,-3,3));

  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",12,0.5,12.5));
  ht.addHist("jet_pt",       new TH1F("jet_pt",      ";p_{T}(jet) [GeV];Events",15,0,300));
  ht.addHist("jet_eta",       new TH1F("jet_eta",      ";#eta(jet) ;Events",12,-3,3));
  ht.addHist("jet_pt1",       new TH1F("jet_pt1",      ";p_{T}(1st jet) [GeV];Events",15,0,300));
  ht.addHist("jet_pt2",       new TH1F("jet_pt2",      ";p_{T}(2nd jet) [GeV];Events",15,0,300));
  ht.addHist("jet_pt3",       new TH1F("jet_pt3",      ";p_{T}(3rd jet) [GeV];Events",15,0,300));

  ht.addHist("nbjets",		 new TH1F("nbjets",      ";b-jet multiplicity;Events",5,0.5,5.5));
  ht.addHist("bjet_pt1",       new TH1F("bjet_pt1",      ";p_{T}(1st jet) [GeV];Events",15,0,300));
  ht.addHist("bjet_pt2",       new TH1F("bjet_pt2",      ";p_{T}(2nd jet) [GeV];Events",15,0,300));
  ht.addHist("bjet_pt3",       new TH1F("bjet_pt3",      ";p_{T}(3rd jet) [GeV];Events",15,0,300));

  ht.addHist("HT",            new TH1F("HT",           ";H_{T} [GeV];Events",20,0,1000));
  ht.addHist("mllb",          new TH1F("mllb",         ";m(l,b) [GeV];Events",12,0,1000));
  ht.addHist("signal_mllb",          new TH1F("signal_mllb",         ";m(l,b) [GeV];Events",12,0,1000));


  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));
  ht.addHist("ratevsrun_bc",    new TH1F("ratevsrun_bc",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));

  ht.addHist("nmuons",   new TH1F("nmuons",       ";N(muons);Events",6,-0.5,5.5));
  ht.addHist("nelectrons",   new TH1F("nelectrons",       ";N(electrons);Events",6,-0.5,5.5));
  ht.addHist("nelmu",   new TH1F("nelmu",       ";N(elmu);Events",6,-0.5,5.5));
  ht.addHist("q",   new TH1F("q",       ";charge;Events",5,-2.5,2.5));

  ht.addHist("lep_pt",		 new TH1F("lep_pt",       ";p_{T}(l) [GeV]; Events", 20,0,300));
  ht.addHist("lep_pt1",            new TH1F("lep_pt1",       ";p_{T}(Leading lepton) [GeV]; Events", 20,0,300));
  ht.addHist("lep_pt2",            new TH1F("lep_pt2",       ";p_{T}(Sub-leading lepton) [GeV]; Events", 20,0,300));
  ht.addHist("lep_eta",		 new TH1F("lep_eta",      ";#eta(lepton) ; Events", 10,-2.5,2.5));

  ht.addHist("lep_2_pt",		 new TH1F("lep_2_pt",       ";p_{T}(l) [GeV]; Events", 20,15,300));
  ht.addHist("lep_2_eta",		 new TH1F("lep_2_eta",       ";#eta(l) [GeV]; Events", 10,-2.5,2.5));

  ht.addHist("lep_phi1",		 new TH1F("lep_phi1",      ";#phi(lepton1) ; Events", 10,-3.2,3.2));
  ht.addHist("lep_phi2",		 new TH1F("lep_phi2",      ";#phi(lepton2) ; Events", 10,-3.2,3.2));
  ht.addHist("dphi_dilep",        new TH1F("dphi_dilep",   ";#Delta#phi(l,l);Events", 10, -4,4));
  ht.addHist("deta_dilep",        new TH1F("deta_dilep",   ";#Delta#eta(l,l);Events", 10, -4,4));
  ht.addHist("bjet_pt",      new TH1F("bjet_pt",     ";p_{T}(b jet) [GeV]; Events", 15,0,300));
  ht.addHist("bjet_eta",	 new TH1F("bjet_eta",    ";#eta(b jet) ; Events", 12,-3,3));
  ht.addHist("met",       new TH1F("met",      ";MET [GeV]; Events", 20,0,300));
  ht.addHist("met_phi",       new TH1F("met_phi",      ";MET #phi; Events", 10, -4,4));
  ht.addHist("DR_jl",       new TH1F("DR_jl",      ";#DeltaR(j,l); Events", 12, 0.,4.));

  ht.addHist("h_dr_jq",       new TH1F("h_dr_jq",      ";#DeltaR(j,q); Events", 12, 0.,4.));


  ht.addHist("hf_csv", new TH1F("hf_csv", ";csv ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_deepcsv", new TH1F("hf_deepcsv", ";deepcsv ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probc", new TH1F("hf_probc", ";P[charm] ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probudsg", new TH1F("hf_probudsg", ";P[udsg] ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb", new TH1F("hf_probb", ";P[b] ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probbb", new TH1F("hf_probbb", ";P[bb] ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_CvsL", new TH1F("hf_CvsL", ";CvsL ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_CvsB", new TH1F("hf_CvsB", ";CvsB ; Events", 20.,-0.5,1.5));

  ht.addHist("hf_CvsL_vs_CvsB", new TH2D("hf_CvsL_vs_CvsB", ";CvsL ; CvsB", 20.,-0.5,1.5, 20.,-0.5,1.5));
  ht.addHist("hf_CvsL_vs_CvsB_gen_c", new TH2D("hf_CvsL_vs_CvsB_gen_c", ";CvsL ; CvsB", 20.,-0.5,1.5, 20.,-0.5,1.5));
  ht.addHist("hf_CvsL_vs_CvsB_gen_b", new TH2D("hf_CvsL_vs_CvsB_gen_b", ";CvsL ; CvsB", 20.,-0.5,1.5, 20.,-0.5,1.5));
  ht.addHist("hf_CvsL_vs_CvsB_gen_lightjet", new TH2D("hf_CvsL_vs_CvsB_gen_lightjet", ";CvsL ; CvsB", 20.,-0.5,1.5, 20.,-0.5,1.5));

  ht.addHist("h_CvsL1_300", new TH1F("h_CvsL1_300", ";CvsL1 300 GeV ; Events", 20.,-0.5,1.5));

  ht.addHist("h_deepjet", new TH1F("h_deepjet", ";deepjet ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_btag_loose", new TH1F("h_deepjet_btag_loose", ";deep jet btag loose ; Events", 2,-0.5,1.5));
  ht.addHist("h_deepjet_btag_medium", new TH1F("h_deepjet_btag_medium", ";deep jet btag medium ; Events", 2,-0.5,1.5));
  ht.addHist("h_deepjet_btag_tight", new TH1F("h_deepjet_btag_tight", ";deep jet btag tight ; Events", 2,-0.5,1.5));
  ht.addHist("h_deepjet_probb", new TH1F("h_deepjet_probb", ";deepjet probb ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_probbb", new TH1F("h_deepjet_probbb", ";deepjet probbb ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_problepb", new TH1F("h_deepjet_problepb", ";deepjet problepb ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_probc", new TH1F("h_deepjet_probc", ";deepjet probc ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_probuds", new TH1F("h_deepjet_probuds", ";deepjet probuds ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_probg", new TH1F("h_deepjet_probg", ";deepjet probg ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_CvsL", new TH1F("h_deepjet_CvsL", ";deepjet CvsL ; Events", 20.,-0.5,1.5));
  ht.addHist("h_deepjet_CvsB", new TH1F("h_deepjet_CvsB", ";deepjet CvsB ; Events", 20.,-0.5,1.5));



  ht.addHist("hf_probb_gen_b", new TH1F("hf_probb_b", ";P[b] w/ Matched GEN b ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_c", new TH1F("hf_probb_c", ";P[c] w/ Matched GEN c ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_lightjet", new TH1F("hf_probb_gen_lightjet", ";P[udsg] w/ Matched GEN light quarks ; Events", 20.,-0.5,1.5));

  ht.addHist("h_scan_mass", new TH1F("h_scan_mass",";particle mass (GeV) ; Events", 20., 200,1000));
  ht.addHist("h_scan_rho", new TH1F("h_scan_rho",";rho ; Events", 10., 0,1));
  ht.addHist("h_scan_coup", new TH1F("h_scan_coup",";coupling ; Events", 4, 0,4));

  ht.addHist("h_m_top_charm", new TH1F("h_m_top_charm",";m(top,charm) ; Events", 100, 0,1000));
  ht.addHist("h_m_top_charm_x", new TH1F("h_m_top_charm_x",";m(top,charm) ; Events", 100, 0,1000));

  ht.addHist("m_ll",		new TH1F("m_ll",   ";M(l+,l-) [GeV] ; Events", 50,0,600));
  ht.addHist("h_flags", new TH1F("h_flags",   ";Flag ; Events", 10,-0.5,9.5));
  ht.addHist("h_hadflav", new TH1F("h_hadflav",  ";Flavour ; Events", 6, -0.5,5.5));

//  TH1F *a_test1 = new TH1F("a_test1","a_test1",30,0,60);//for debugging

//  ht.addHist("b_matched_jet",  new TH1F("b_matched_jet", ";p_T(b matched jet) [GeV]; Events", 24,0,600));
//  ht.addHist("c_matched_jet",  new TH1F("c_matched_jet", ";p_T(c matched jet) [GeV]; Events", 24,0,600));

  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
    ht.getPlots()["ratevsrun_bc"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }


  //INPUT
  MiniEvent_t ev;

  const char* treename = "TreeInput";
  const char* tree300 = "Tree300";
  const char* tree350 = "Tree350";
  const char* tree400 = "Tree400";
  const char* tree500 = "Tree500";
  const char* tree600 = "Tree600";
  const char* tree700 = "Tree700";

    TTree t_input(treename,treename);

//  int t_event;
    float CvsL1,CvsB1;
    float CvsL2,CvsB2;
    float CvsL3,CvsB3;
    float t_m_lep_charm;
    float t_m_lep_bottom;
    float t_HT;
    float t_dphi_ll;
    float t_deepcsv;
    float t_deepjet;
    float t_pt_l1,t_pt_l2,t_eta_l1,t_eta_l2,t_phi_l1,t_phi_l2;
    float t_pt_j1,t_pt_j2,t_pt_j3,t_eta_j1,t_eta_j2,t_eta_j3;
    float t_phi_j1,t_phi_j2,t_phi_j3;
    float t_met;
    float t_weight, t_normWgt, t_norm;
    float t_scan_mass, t_scan_rho, t_scan_coup;

    int Flag_HBHENoiseFilter;
    int Flag_HBHENoiseIsoFilter;
    int Flag_EcalDeadCellTriggerPrimitiveFilter;
    int Flag_goodVertices;
    int Flag_eeBadScFilter;
    int Flag_globalSuperTightHalo2016Filter;
    //int Flag_BadChargedCandidateFilter;//Not recommended yet.
    int Flag_BadPFMuonFilter;
    int Flag_badPFMuonDzFilter;
    int Flag_ecalBadCalibFilter;
    
    
    t_input.Branch("event",&ev.event,"event/I");
    t_input.Branch("run",&ev.run,"run/i");
    t_input.Branch("lumi",&ev.lumi,"lumi/i");
    t_input.Branch("t_weight",&t_weight,"t_weight/F");
    t_input.Branch("t_normWgt",&t_normWgt,"t_normWgt/F");
    t_input.Branch("t_norm",&t_norm,"t_norm/F");

    t_input.Branch("CvsL1",&CvsL1,"CvsL1/F");
    t_input.Branch("CvsB1",&CvsB1,"CvsB1/F");
    t_input.Branch("CvsL2",&CvsL2,"CvsL2/F");
    t_input.Branch("CvsB2",&CvsB2,"CvsB2/F");
    t_input.Branch("CvsL3",&CvsL3,"CvsL3/F");
    t_input.Branch("CvsB3",&CvsB3,"CvsB3/F");
    t_input.Branch("t_m_lep_jet1",&t_m_lep_charm,"t_m_lep_jet1/F");
    t_input.Branch("t_m_lep_bottom",&t_m_lep_bottom,"t_m_lep_bottom/F");
    t_input.Branch("t_HT",&t_HT,"t_HT/F");
    t_input.Branch("t_dphi_ll",&t_dphi_ll,"t_dphi_ll/F");
    t_input.Branch("t_deepcsv",&t_deepcsv,"t_deepcsv/F");
    t_input.Branch("t_deepjet",&t_deepjet,"t_deepjet/F");

    t_input.Branch("t_pt_l1", &t_pt_l1, "t_pt_l1/F");
    t_input.Branch("t_pt_l2", &t_pt_l2, "t_pt_l2/F");
    t_input.Branch("t_eta_l1", &t_eta_l1, "t_eta_l1/F");
    t_input.Branch("t_eta_l2", &t_eta_l2, "t_eta_l2/F");
    t_input.Branch("t_phi_l1", &t_phi_l1, "t_phi_l1/F");
    t_input.Branch("t_phi_l2", &t_phi_l2, "t_phi_l2/F");
    t_input.Branch("t_pt_j1", &t_pt_j1, "t_pt_j1/F");
    t_input.Branch("t_pt_j2", &t_pt_j2, "t_pt_j2/F");
    t_input.Branch("t_pt_j3", &t_pt_j3, "t_pt_j3/F");
    t_input.Branch("t_eta_j1", &t_eta_j1, "t_eta_j1/F");
    t_input.Branch("t_eta_j2", &t_eta_j2, "t_eta_j2/F");
    t_input.Branch("t_eta_j3", &t_eta_j3, "t_eta_j3/F");
    t_input.Branch("t_phi_j1", &t_phi_j1, "t_phi_j1/F");
    t_input.Branch("t_phi_j2", &t_phi_j2, "t_phi_j2/F");
    t_input.Branch("t_phi_j3", &t_phi_j3, "t_phi_j3/F");
    t_input.Branch("t_met", &t_met, "t_met/F");

    t_input.Branch("t_scan_mass", &t_scan_mass, "t_scan_mass/F");
    t_input.Branch("t_scan_rho", &t_scan_rho, "t_scan_rho/F");
    t_input.Branch("t_scan_coup", &t_scan_coup, "t_scan_coup/F");
    
  TTree *t_300 = t_input.CloneTree();
  t_300->SetName(tree300);
  TTree *t_350 = t_input.CloneTree();
  t_350->SetName(tree350);
  TTree *t_400 = t_input.CloneTree();
  t_400->SetName(tree400);
  TTree *t_500 = t_input.CloneTree();
  t_500->SetName(tree500);
  TTree *t_600 = t_input.CloneTree();
  t_600->SetName(tree600);
  TTree *t_700 = t_input.CloneTree();
  t_700->SetName(tree700);
  //TFile *f_300 = TFile::Open(dirName+"/"+"300_"+baseName,"RECREATE");
  //  t_300->Print();

  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }

  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(10000000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;


  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);

  //EVENT LOOP
  //select 2mu+>=3 jets events triggered by a single muon trigger


  int Ntotal = 0;
  int N_after_all_selections = 0;
  int Ntotal_after_trig = 0;

  /*
  int Ntotal_lepton = 0;
  int Ntotal_after_samesign = 0;
  int Ntotal_after_jet = 0;
  int Ntotal_after_num_btags = 0;
  int Ntotal_after_met = 0;

  int Ntotal_Z = 0;
*/
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
      //trigger
      int passtrigger_ee = 0;
      int passtrigger_mm = 0;
      int passtrigger_em = 0;
      if(era.Contains("2017")) {
        if (baseMC.Contains("DoubleMuon",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits);
          int b =0;
          if(ev.run < 299330){ b = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", ev.triggerBits);}
          if(ev.run > 299330){ b = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits);}
          if (a || b) passtrigger_mm = 1;
        }
        if (baseMC.Contains("SingleMuon",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits);
          int b=0;
          if(ev.run < 299330) {b = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", ev.triggerBits);}
          if(ev.run > 299330) {b = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits);}
          int c = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int d = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int e = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int f = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int g = selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits);
          if (!(a||b) && g) passtrigger_mm = 1;
          if (!(c||d||e||f) && g) passtrigger_em = 1;
        }
        if (baseMC.Contains("DoubleEG",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int b = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ", ev.triggerBits);
          if (a||b) passtrigger_ee = 1;
        }
        if (baseMC.Contains("SingleElectron",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int b = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ", ev.triggerBits);
          int c = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int d = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int e = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int f = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int g = selector.hasTriggerBit("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v", ev.triggerBits);
          int h = selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits);
          if (!(a||b) && (g||h)) passtrigger_ee = 1;
          if (!(c||d||e||f) && (g||h)) passtrigger_em = 1;
        }
        if (baseMC.Contains("MuonEG",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int b = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int c = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int d = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          if (a || b || c || d) passtrigger_em = 1;
        }
        if (baseMC.Contains("MC",TString::kIgnoreCase)){
          int a = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits);
          int b = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits);
          int c = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", ev.triggerBits);
          int d = selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits);
          int e = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int f = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ", ev.triggerBits);
          int g = selector.hasTriggerBit("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v", ev.triggerBits);
          int h = selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits);
          int i = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int j = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          int k = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          int l = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          if (a||b||c||d) passtrigger_mm = 1;
          if (e||f||g||h) passtrigger_ee = 1;
          if (d||g||h||i||j||k||l) passtrigger_em = 1;
        }
     }

        Ntotal++;
//      cout << "Trigger = "<<hasMTrigger << endl;
//      if(ev.isData && !hasMTrigger) continue;


//      Filters
        Flag_HBHENoiseFilter=(ev.met_filterBits)&0x1;
        Flag_HBHENoiseIsoFilter=(ev.met_filterBits>>1)&0x1;
        Flag_EcalDeadCellTriggerPrimitiveFilter=(ev.met_filterBits>>2)&0x1;
        Flag_goodVertices=(ev.met_filterBits>>3)&0x1;
        Flag_eeBadScFilter=(ev.met_filterBits>>4)&0x1;//Not suggested for MC
        Flag_globalSuperTightHalo2016Filter=(ev.met_filterBits>>5)&0x1;
        Flag_BadPFMuonFilter=(ev.met_filterBits>>6)&0x1;
 //     TO BE USED IN NEXT ROUND OF NTUPLE PRODUCTION!!!
/* 
       Flag_badPFMuonDzFilter=(ev.met_filterBits>>7)&0x1;
        Flag_ecalBadCalibFilter=(ev.met_filterBits>>8)&0x1;
*/

        ht.fill("h_flags",Flag_eeBadScFilter,1);      
          
 
  /* 
      cout<<Flag_HBHENoiseFilter<<"  "<<Flag_HBHENoiseIsoFilter<<"  "
          <<Flag_EcalDeadCellTriggerPrimitiveFilter<<"  "
          <<Flag_goodVertices<<"  "<<Flag_eeBadScFilter<<"  "
          <<Flag_globalSuperTightHalo2016Filter<<"  "
          <<Flag_BadPFMuonFilter<<"  "<<Flag_badPFMuonDzFilter<<"  "
          <<Flag_ecalBadCalibFilter<<endl;
  */

        if (!(Flag_HBHENoiseFilter && 
              Flag_HBHENoiseIsoFilter && 
              Flag_EcalDeadCellTriggerPrimitiveFilter && 
              Flag_goodVertices && 
              Flag_globalSuperTightHalo2016Filter && 
              Flag_BadPFMuonFilter)) continue;
 // ----TO BE USED IN NEXT ROUND OF NTUPLE PRODUCTION!!!---             
/*
        if (!(Flag_badPFMuonDzFilter &&
              Flag_ecalBadCalibFilter)) continue;
        if (ev.isData && !Flag_eeBadScFilter) continue;
*/
//---------------------------------------------------------


      ht.fill("h_scan_mass_bc",ev.scan_mass,1);


      //if(!hasMTrigger) continue;
      //cout<<dimutrig<<dielectrig<<mueltrig<<singmutrig<<singeltrig<<endl;
     // if (!passtrigger) continue;


      //Ntotal_after_trig++;

      //select two offline muons
      std::vector<Particle> flaggedleptons = selector.flaggedLeptons(ev);
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
      //std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::MVA90,minLeptonPt,2.4);
      std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::TIGHT,minLeptonPt,2.4);



      if(leptons.size()<2) continue;
//      Ntotal_lepton++;


      sort(leptons.begin(),leptons.end(),
	            [](const Particle& a, const Particle& b){
		          return a.Pt() > b.Pt();
	           }
      );

/*
      a_test1->Fill(leptons[0].pt());
      a_test1->SetDirectory(0);//
*/

      if (leptons[0].pt() < 30.) continue;
      //if (leptons[0].pt() < 40. && leptons[0].id() == 11) continue;
      if (leptons[1].pt() < 20.) continue;
      if (leptons.size() > 2 && leptons[2].pt() > 20.) continue;

      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun_bc",runBin,lumi,"test");
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }



      int dimuon_event = 0;
      int dielectron_event = 0;
      int emu_event = 0;
      int mue_event = 0;
//      if (leptons[0].id() == 11 && leptons[1].id() == 11) dielectron_event = 1;
//      if (leptons[0].id() == 13 && leptons[1].id() == 13) dimuon_event = 1;
//      if (leptons[0].id() == 11 && leptons[1].id() == 13) emu_event = 1;
//      if (leptons[0].id() == 13 && leptons[1].id() == 11) mue_event = 1;
      if (baseMC.Contains("DoubleMuon",TString::kIgnoreCase)){
                if (!((leptons[0].id() == 13 && leptons[1].id() == 13) && passtrigger_mm)) continue;
                dimuon_event=1;
        }
        if (baseMC.Contains("DoubleEG",TString::kIgnoreCase)){
                if (!((leptons[0].id() == 11 && leptons[1].id() == 11) && passtrigger_ee)) continue;
                dielectron_event = 1;
        }
        if (baseMC.Contains("MuonEG",TString::kIgnoreCase)){
                if (!((leptons[0].id()+leptons[1].id())==24 && passtrigger_em)) continue;
                if (leptons[0].id() == 11 && leptons[1].id() == 13) emu_event = 1;
                if (leptons[0].id() == 13 && leptons[1].id() == 11) mue_event = 1;
        }
        if (baseMC.Contains("SingleElectron",TString::kIgnoreCase)){
                if (!(((leptons[0].id() == 11 && leptons[1].id() == 11) && passtrigger_ee) || ((leptons[0].id()+leptons[1].id())==24 && passtrigger_em))) continue;
                if ((leptons[0].id()+leptons[1].id())==24) emu_event = 1;
                else dielectron_event = 1;
        }
        if (baseMC.Contains("SingleMuon",TString::kIgnoreCase)){
                if (!(((leptons[0].id() == 13 && leptons[1].id() == 13) && passtrigger_mm) || ((leptons[0].id()+leptons[1].id())==24 && passtrigger_em))) continue;
                if ((leptons[0].id()+leptons[1].id())==24) emu_event = 1;
                else dimuon_event = 1;
        }
        if (baseMC.Contains("MC",TString::kIgnoreCase)){
                if ((leptons[0].id() == 11 && leptons[1].id() == 11) && passtrigger_ee) dielectron_event = 1;
                if ((leptons[0].id() == 13 && leptons[1].id() == 13) && passtrigger_mm) dimuon_event = 1;
                if (((leptons[0].id() + leptons[1].id()) == 24) && passtrigger_em) emu_event = 1;
        }

      std::vector<TString> tags2={"inc"};
      if (dimuon_event) tags2.push_back("mm");
      if (dielectron_event) tags2.push_back("ee");
      if (emu_event || mue_event) tags2.push_back("emu");

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);
      std::vector<Jet> allJets = selector.getGoodJets(ev,50.,2.4,leptons,{});
      bool passJets(allJets.size()>=minJetMultiplicity);

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);

      //event weight
      float evWgt(1.0);

      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,tags2);
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //MC specific: compute event weight
      if (!ev.isData) {
        float normWgt(normH? normH->GetBinContent(1) : 1.0);
        TString period = lumi.assignRunPeriod();
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
//        EffCorrection_t selSF(1.0,0.0);// = lepEffH.getOfflineCorrection(leptons[0], period);
        float muonSF = 1.;
        float electronSF = 1.;
        float emuSF = 1.;
        float dilepton_trig_SF = 1.;
        if (dimuon_event == 1){
          EffCorrection_t muon1SF = lepEffH.getMuonSF(leptons[0].pt(),leptons[0].eta());
          EffCorrection_t muon2SF = lepEffH.getMuonSF(leptons[1].pt(),leptons[1].eta());
          muonSF = muon1SF.first*muon2SF.first;
          EffCorrection_t dimuon_pt_trig_factor = lepEffH.getMuMuPtSF(leptons[0].pt(),leptons[1].pt());
          EffCorrection_t dimuon_eta_trig_factor = lepEffH.getMuMuEtaSF(abs(leptons[0].eta()),abs(leptons[1].eta()));
          dilepton_trig_SF = dimuon_pt_trig_factor.first*dimuon_eta_trig_factor.first;
        }
        if (dielectron_event == 1){
          EffCorrection_t electron1SF = lepEffH.getElectronSF(leptons[0].pt(),leptons[0].eta());
          EffCorrection_t electron2SF = lepEffH.getElectronSF(leptons[1].pt(),leptons[1].eta());
          electronSF = electron1SF.first*electron2SF.first;
          EffCorrection_t dielectron_pt_trig_factor = lepEffH.getEEPtSF(leptons[0].pt(),leptons[1].pt());
          EffCorrection_t dielectron_eta_trig_factor = lepEffH.getEEEtaSF(abs(leptons[0].eta()),abs(leptons[1].eta()));
          //1.25 from fit to Z->e+e- data/MC ratio
          //dilepton_trig_SF = 1.25*dielectron_pt_trig_factor.first*dielectron_eta_trig_factor.first;
          dilepton_trig_SF = dielectron_pt_trig_factor.first*dielectron_eta_trig_factor.first;
        }
        if (emu_event == 1){
          EffCorrection_t electron1SF = lepEffH.getElectronSF(leptons[0].pt(),leptons[0].eta());
          EffCorrection_t muon2SF = lepEffH.getMuonSF(leptons[1].pt(),leptons[1].eta());
          emuSF = electron1SF.first*muon2SF.first;
          EffCorrection_t emu_pt_trig_factor = lepEffH.getEMuPtSF(leptons[0].pt(),leptons[1].pt());
          EffCorrection_t emu_eta_trig_factor = lepEffH.getEMuEtaSF(abs(leptons[0].eta()),abs(leptons[1].eta()));
          dilepton_trig_SF = emu_pt_trig_factor.first*emu_eta_trig_factor.first;
        }
        if (mue_event == 1){
          EffCorrection_t muon1SF = lepEffH.getMuonSF(leptons[0].pt(),leptons[0].eta());
          EffCorrection_t electron2SF = lepEffH.getElectronSF(leptons[1].pt(),leptons[1].eta());
          emuSF = muon1SF.first*electron2SF.first;
          EffCorrection_t emu_pt_trig_factor = lepEffH.getEMuPtSF(leptons[0].pt(),leptons[1].pt());
          EffCorrection_t emu_eta_trig_factor = lepEffH.getEMuEtaSF(abs(leptons[0].eta()),abs(leptons[1].eta()));
          dilepton_trig_SF = emu_pt_trig_factor.first*emu_eta_trig_factor.first;
        }
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});
        float leptonSF = muonSF*electronSF*emuSF;

        EffCorrection_t Chargeflip_SF=lepEffH.getChargeFlipSF(in_fname,"Chi2",leptons[0].pt(),leptons[0].eta(),leptons[0].charge(),leptons[1].pt(),leptons[1].eta(),leptons[1].charge());

        evWgt  = normWgt*puWgt*l1prefireProb.first*leptonSF*dilepton_trig_SF*Chargeflip_SF.first;
        t_normWgt = normWgt;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);//generator weights
      }


      //if (evWgt < 0.) continue;

      float invariant_mass = (leptons[0]+leptons[1]).M();
      ht.fill("control_m_ll_bc",  invariant_mass,        evWgt, tags2);

      int tag_isgsf = 0;
      if (dielectron_event && (ev.l_isGsfCtfScPixChargeConsistent[1] && ev.l_isGsfCtfScPixChargeConsistent[1])) tag_isgsf = 1;
      if (emu_event && ev.l_isGsfCtfScPixChargeConsistent[0]) tag_isgsf = 1;
      if (mue_event && ev.l_isGsfCtfScPixChargeConsistent[1]) tag_isgsf = 1;

      float zmass = (leptons[0]+leptons[1]).M();
      if (leptons[0].charge()*leptons[1].charge() < 0){
        ht.fill("control_Z_mass_overall", zmass,        evWgt, tags2);
//        if (zmass > 50. && zmass < 170.){
//        if (zmass > 50. && zmass < 170. &&  leptons[0].pt()>50 && leptons[0].pt()<100 && leptons[1].pt()>50 && leptons[1].pt()<100){
//          if (allJets.size() == 1){      
          ht.fill("control_Z_mass", zmass,        evWgt, tags2);
          ht.fill("control_2_lep_pt_bc", leptons[0].pt(),evWgt, tags2);
          ht.fill("control_2_lep_pt_bc", leptons[1].pt(),evWgt, tags2);
          ht.fill("control_2_lep_eta_bc", leptons[0].eta(),evWgt, tags2);
          ht.fill("control_2_lep_eta_bc", leptons[1].eta(),evWgt, tags2);
          for (size_t in_nlep=0; in_nlep<leptons.size();in_nlep++){
            ht.fill("control_lep_pt_bc", leptons[in_nlep].pt(), evWgt, tags2);
            ht.fill("control_lep_eta_bc",leptons[in_nlep].eta(), evWgt, tags2);
          }
          ht.fill("control_pT_ll_bc",(leptons[0]+leptons[1]).Pt(),evWgt, tags2);
          ht.fill("control_Y_ll_bc",(leptons[0]+leptons[1]).Eta(),evWgt, tags2);
//        }
      }

      if (tag_isgsf == 1){ 
         ht.fill("isgsf_control_m_ll_bc",  invariant_mass,        evWgt, tags2);
         if (leptons[0].charge()*leptons[1].charge() < 0 && tag_isgsf == 1){
             ht.fill("isgsf_control_Z_mass_overall", zmass,        evWgt, tags2);
             ht.fill("isgsf_control_Z_mass", zmass,        evWgt, tags2);
             ht.fill("isgsf_control_2_lep_pt_bc", leptons[0].pt(),evWgt, tags2);
             ht.fill("isgsf_control_2_lep_pt_bc", leptons[1].pt(),evWgt, tags2);
             ht.fill("isgsf_control_2_lep_eta_bc", leptons[0].eta(),evWgt, tags2);
             ht.fill("isgsf_control_2_lep_eta_bc", leptons[1].eta(),evWgt, tags2);
             for (size_t in_nlep=0; in_nlep<leptons.size();in_nlep++){
               ht.fill("isgsf_control_lep_pt_bc", leptons[in_nlep].pt(), evWgt, tags2);
               ht.fill("isgsf_control_lep_eta_bc",leptons[in_nlep].eta(), evWgt, tags2);
             }
             ht.fill("isgsf_control_pT_ll_bc",(leptons[0]+leptons[1]).Pt(),evWgt, tags2);
             ht.fill("isgsf_control_Y_ll_bc",(leptons[0]+leptons[1]).Eta(),evWgt, tags2);           
         }
      }

     if (allJets.size() <= 3){
      if (leptons[0].charge()*leptons[1].charge() < 0){
         if (tag_isgsf == 1) ht.fill("OS_isgsf_control_m_ll_bc",  invariant_mass,        evWgt, tags2);  
         ht.fill("OS_control_m_ll_bc",  invariant_mass,        evWgt, tags2);
      }
      if (leptons[0].charge()*leptons[1].charge() > 0){
         if (tag_isgsf == 1) ht.fill("SS_isgsf_control_m_ll_bc",  invariant_mass,        evWgt, tags2);
         ht.fill("SS_control_m_ll_bc",  invariant_mass,        evWgt, tags2);
      }
     }

      int num_btags = 0;
      for(size_t ij=0; ij<allJets.size(); ij++){
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;
		      num_btags++;
        }
/*
      if(!passJets) continue;
      if(num_btags < minNum_btags) continue;
*/

//      if(ev.met_pt < 30.) continue;

      //sort jets w.r.t. pt
      /*
      sort(allJets.begin(),allJets.end(),
        [](const Jet& a, const Jet& b){
                return a.Pt() > b.Pt();
        }
      );
      */

      //select gen-level objects
/*
      std::vector<Particle> gen_leptons=selector.getGenLeptons(ev,20.,2.5);
      std::vector<Particle> gen_photons=selector.getGenPhotons(ev,50.,2.4);
      std::vector<Jet> gen_jets = selector.getGenJets(ev,25,2.5,gen_leptons,gen_photons);
      cout<<"Event: "<<iev<<endl;
      cout<<"Genjets: "<<gen_jets.size()<<endl;
      for (size_t i_gen=0;i_gen<gen_jets.size();i_gen++) cout<<i_gen<<"  "<<gen_jets[i_gen].pt()<<"  "<<gen_jets[i_gen].eta()<<"   "<<gen_jets[i_gen].Phi()<<endl;
      for (int i=0;i<ev.ng;i++){
        if (abs(ev.g_id[i])>5) continue;
        cout<<ev.g_id[i]<<"  "<<ev.g_pt[i]<<"  "<<ev.g_eta[i]<<"  "<<ev.g_phi[i]<<endl;
        TLorentzVector genjet4mom;
        genjet4mom.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
      }
*/


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    std::vector<Jet> jets;
    for(auto &j: allJets) {

      int idx=j.getJetIndex();

      //pileup jet id
      int jid=ev.j_id[idx];
      bool passLoosePu((jid>>2)&0x1);
      if(!passLoosePu) continue;

      //ECAL noise (2017)
      float emf=ev.j_emf[idx];
      if(era.Contains("2017") && fabs(j.Eta())>2.650 && fabs(j.Eta())<3.139 && emf>0.55) continue;

      //HEM 15/16 issue (2018)
      //if(is2018 && j.Eta()>-3.0 && j.Eta()<-1.3 && j.Phi()>-1.57 && j.Phi()<-0.87) continue;

      jets.push_back(j);
    }

  
        
  if (leptons[0].charge()*leptons[1].charge() < 0 && zmass > 70. && zmass < 110. && jets.size()>0) ht.fill("njets_bc",      jets.size(), evWgt, tags2);

        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (leptons[0].charge()*leptons[1].charge() < 0){
            if (zmass > 70. && zmass < 110.){
                for(size_t ij=0; ij<jets.size(); ij++){
                    ht.fill("jet_pt_bc",    jets[ij].pt(), evWgt, tags2);
                    ht.fill("jet_eta_bc",   jets[ij].eta(), evWgt, tags2);
                }
            }
        }

// note that in productions prior and including githash  29490641, the index of isGsfCtfScPixChargeConsistent is +1 w.r.t. the ones below
  if (dielectron_event && (!ev.l_isGsfCtfScPixChargeConsistent[0] || !ev.l_isGsfCtfScPixChargeConsistent[1])) continue;
  if (emu_event && !ev.l_isGsfCtfScPixChargeConsistent[0]) continue;
  if (mue_event && !ev.l_isGsfCtfScPixChargeConsistent[1]) continue;

   
        
  if (leptons[0].charge()*leptons[1].charge() < 0) continue;

  invariant_mass = (leptons[0]+leptons[1]).M();
  std::vector<TString> tags3={"inc"};
  tags3.push_back(to_string(ev.scan_mass));

  if (ev.met_pt < 30.) continue;
  if (jets.size() < 3) continue;
  if (leptons[0].DeltaR(leptons[1]) < 0.3) continue; 
  if (invariant_mass < 20.) continue;
  if (invariant_mass > 76. && invariant_mass < 106.) continue;


  float HT = 0;
  for(size_t ij=0; ij<jets.size(); ij++){
    ht.fill("jet_pt",    jets[ij].pt(), evWgt, tags2);
    ht.fill("jet_eta",   jets[ij].eta(), evWgt, tags2);
    HT += jets[ij].pt();
    for (size_t ik = 0;ik<leptons.size();ik++) ht.fill("DR_jl",sqrt(pow(leptons[ik].DeltaPhi(jets[ij]),2)+pow(DeltaEta(leptons[ik].eta(),jets[ij].eta()),2)),evWgt,tags2);
    int idx=jets[ij].getJetIndex();
    bool passBtag(ev.j_btag[idx]>0);
    if(!passBtag) continue;
    float mllb( (leptons[0]+leptons[1]+jets[ij]).M() );//M is the magnitude.
    ht.fill("mllb",mllb,evWgt,tags2);
    ht.fill("bjet_pt",jets[ij].pt(),evWgt,tags2);
    ht.fill("bjet_eta",jets[ij].eta(),evWgt,tags2);
    ht.fill("signal_mllb", mllb, evWgt, tags3);
}



float delta_phi = leptons[0].DeltaPhi(leptons[1]);

t_pt_l1=leptons[0].pt();
t_pt_l2=leptons[1].pt();
t_eta_l1=leptons[0].eta();
t_eta_l2=leptons[1].eta();
t_phi_l1=leptons[0].phi();
t_phi_l2=leptons[1].phi();
t_met=ev.met_pt;
//t_norm=norm;
t_HT = HT;
t_dphi_ll=delta_phi;
t_scan_mass=ev.scan_mass;
t_scan_rho=ev.scan_rho;
t_scan_coup=ev.scan_coup;
t_weight=evWgt;

//  bool passJets(jets.size()>=minJetMultiplicity);

 if(!passJets) continue;//not in recent previous tests. 
 if(num_btags < minNum_btags) continue;//not in recent previous tests.



 sort(jets.begin(),jets.end(),
   [](const Jet& a, const Jet& b){
           return a.Pt() > b.Pt();
   }
 );

  float Ctag_SF=1.0;
  int jet_index=0;
  for(size_t ij=0; ij<jets.size(); ij++){
    int idx=jets[ij].getJetIndex();
    bool passBtag(ev.j_btag[idx]>0);
    ht.fill("h_hadflav",abs(ev.j_hadflav[idx]),1);
    if (jet_index==0){
      CvsL1 = ev.j_CvsL[idx];
      CvsB1 = ev.j_CvsB[idx];
      float mlc(((leptons[0]+jets[ij]).M()));
      t_m_lep_charm = mlc;
      t_pt_j1=jets[ij].pt();
      t_eta_j1=jets[ij].eta();
      t_phi_j1=jets[ij].Phi();
      if (ev.scan_mass > 299. && ev.scan_mass < 301.) ht.fill("h_CvsL1_300",CvsL1,1,tags2);
      if (passBtag) t_m_lep_bottom = mlc;
      if (CvsL1 > 0.5 and CvsB1 > 0.5){
        ht.fill("h_m_top_charm",mlc,evWgt,tags3);
        ht.fill("h_m_top_charm_x",mlc,evWgt,tags2);
      }
    }
    if (jet_index==1){
      CvsL2 = ev.j_CvsL[idx];
      CvsB2 = ev.j_CvsB[idx];
      t_pt_j2=jets[ij].pt();
      t_eta_j2=jets[ij].eta();
      t_phi_j2=jets[ij].Phi();
    }
    if (jet_index==2){
      CvsL3 = ev.j_CvsL[idx];
      CvsB3 = ev.j_CvsB[idx];
      t_pt_j3=jets[ij].pt();
      t_eta_j3=jets[ij].eta();
      t_phi_j3=jets[ij].Phi();
    }
    jet_index++;
    // Just apply the C-tagging scale factor to CvsL and CvsB plots in current phase.
    Ctag_SF = 1.0;
    if(!ev.isData)   Ctag_SF=CtagSFWR.GetCtagSF("DeepCSV_alljet", "DeepCSV", jets[ij], ev, evWgt);
    ht.fill("hf_csv",ev.j_csv[idx],evWgt,tags2);
    ht.fill("hf_deepcsv",ev.j_deepcsv[idx],evWgt,tags2);
    ht.fill("hf_probc",ev.j_probc[idx],evWgt,tags2);
    ht.fill("hf_probudsg",ev.j_probudsg[idx],evWgt,tags2);
    ht.fill("hf_probb",ev.j_probb[idx],evWgt,tags2);
    ht.fill("hf_probbb",ev.j_probbb[idx],evWgt,tags2);
    ht.fill("hf_CvsL",ev.j_CvsL[idx],evWgt*Ctag_SF,tags2);
    ht.fill("hf_CvsB",ev.j_CvsB[idx],evWgt*Ctag_SF,tags2);
    ht.fill2D("hf_CvsL_vs_CvsB",ev.j_CvsL[idx],ev.j_CvsB[idx],evWgt*Ctag_SF,tags2);
    t_deepcsv=ev.j_deepcsv[idx];

    Ctag_SF = 1.0;
    if(!ev.isData)   Ctag_SF=CtagSFWR.GetCtagSF("DeepJet_alljet", "DeepJet", jets[ij], ev, evWgt);
    ht.fill("h_deepjet",ev.j_deepjet[idx],evWgt,tags2);
    ht.fill("h_deepjet_btag_loose",ev.j_deepjet_btag_loose[idx],evWgt,tags2);
    ht.fill("h_deepjet_btag_medium",ev.j_deepjet_btag_medium[idx],evWgt,tags2);
    ht.fill("h_deepjet_btag_tight",ev.j_deepjet_btag_tight[idx],evWgt,tags2);
    ht.fill("h_deepjet_probb",ev.j_deepjet_probb[idx],evWgt,tags2);
    ht.fill("h_deepjet_probbb",ev.j_deepjet_probbb[idx],evWgt,tags2);
    ht.fill("h_deepjet_problepb",ev.j_deepjet_problepb[idx],evWgt,tags2);
    ht.fill("h_deepjet_probc",ev.j_deepjet_probc[idx],evWgt,tags2);
    ht.fill("h_deepjet_probuds",ev.j_deepjet_probuds[idx],evWgt,tags2);
    ht.fill("h_deepjet_probg",ev.j_deepjet_probg[idx],evWgt,tags2);
    ht.fill("h_deepjet_CvsL",ev.j_deepjet_CvsL[idx],evWgt*Ctag_SF,tags2);
    ht.fill("h_deepjet_CvsB",ev.j_deepjet_CvsB[idx],evWgt*Ctag_SF,tags2);
    t_deepjet=ev.j_deepjet[idx];






 //if(!passJets) continue;

// if(num_btags < minNum_btags) continue;


     for (int i=0;i<ev.ng;i++){
       if (abs(ev.g_id[i]) < 6 || abs(ev.g_id[i]) == 21){
         TLorentzVector genjet4mom;
         genjet4mom.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
         ht.fill("h_dr_jq",jets[ij].DeltaR(genjet4mom),evWgt,tags2);
         if (jets[ij].DeltaR(genjet4mom) < 0.4){
           if (abs(ev.g_id[i]) == 4){
             ht.fill("hf_probc_gen_c",ev.j_probc[idx],evWgt,tags2);
             ht.fill("hf_probb_gen_c",ev.j_probc[idx],evWgt,tags2);
             ht.fill2D("hf_CvsL_vs_CvsB_gen_c",ev.j_CvsL[idx],ev.j_CvsB[idx],evWgt,tags2);
           }
           if (abs(ev.g_id[i]) == 5){
             ht.fill("hf_probb_gen_b",ev.j_probb[idx],evWgt,tags2);
             ht.fill2D("hf_CvsL_vs_CvsB_gen_b",ev.j_CvsL[idx],ev.j_CvsB[idx],evWgt,tags2);
           }
           if (abs(ev.g_id[i]) < 4 || abs(ev.g_id[i]) == 21){
             ht.fill("hf_probb_gen_lightjet",ev.j_probudsg[idx],evWgt,tags2);
             ht.fill2D("hf_CvsL_vs_CvsB_gen_lightjet",ev.j_CvsL[idx],ev.j_CvsB[idx],evWgt,tags2);
           }
         }
       }
     }
   }
   N_after_all_selections++;



//   t_event = iev;
   t_input.Fill();
   lepEffH.recordToChargeflipNormRatio(evWgt); //This can actually apply to any step that you want to normalize. But be careful that only apply once in the analysis.
   if (baseMC.Contains("scan",TString::kIgnoreCase)){
      if (ev.scan_mass > 299. && ev.scan_mass < 301.) t_300->Fill();
      if (ev.scan_mass > 349. && ev.scan_mass < 351.) t_350->Fill();
      if (ev.scan_mass > 399. && ev.scan_mass < 401.) t_400->Fill();
      if (ev.scan_mass > 499. && ev.scan_mass < 501.) t_500->Fill();
      if (ev.scan_mass > 599. && ev.scan_mass < 601.) t_600->Fill();
      if (ev.scan_mass > 699. && ev.scan_mass < 701.) t_700->Fill();
   }

/*
     sort(jets.begin(),jets.end(),
       [](const Jet& a, const Jet& b){
               return a.Pt() > b.Pt();
       }
     );
*/
//     float delta_phi = leptons[0].DeltaPhi(leptons[1]);

      ht.fill("nvtx",       ev.nvtx,        evWgt, tags2);
      ht.fill("nvtx_uw",    ev.nvtx,        1., tags2);
      ht.fill("njets",      jets.size(), evWgt, tags2);
      ht.fill("jet_pt1",   jets[0].pt(), evWgt, tags2);
      ht.fill("jet_pt2",   jets[1].pt(), evWgt, tags2);
      ht.fill("jet_pt3",   jets[2].pt(), evWgt, tags2);
      ht.fill("HT", HT, evWgt, tags2);

      for (size_t in_nlep=0; in_nlep<leptons.size();in_nlep++){
      	ht.fill("lep_pt",          leptons[in_nlep].pt(), evWgt, tags2);
      	ht.fill("lep_eta",     leptons[in_nlep].eta(), evWgt, tags2);
      }

      ht.fill("lep_2_pt", leptons[0].pt(), evWgt, tags2);
      ht.fill("lep_2_pt", leptons[1].pt(), evWgt, tags2);
      ht.fill("lep_2_eta", leptons[0].eta(), evWgt, tags2);
      ht.fill("lep_2_eta", leptons[1].eta(), evWgt, tags2);

      ht.fill("lep_pt1",          leptons[0].pt(), evWgt, tags2);
      ht.fill("lep_pt2",          leptons[1].pt(), evWgt, tags2);
      ht.fill("dphi_dilep", delta_phi,evWgt,tags2);
      ht.fill("lep_phi1",leptons[0].phi(),evWgt,tags2);
      ht.fill("lep_phi2",leptons[1].phi(),evWgt,tags2);
      ht.fill("deta_dilep", DeltaEta(leptons[0].eta(),leptons[1].eta()),evWgt,tags2);
      ht.fill("nbjets",num_btags,evWgt,tags2);
      ht.fill("met",ev.met_pt,evWgt,tags2);
      ht.fill("met_phi",ev.met_phi,evWgt,tags2);
      ht.fill("h_scan_mass",ev.scan_mass,1);
      ht.fill("h_scan_rho",ev.scan_rho,evWgt);
      ht.fill("h_scan_coup",ev.scan_coup,evWgt);
      ht.fill("m_ll",(leptons[0]+leptons[1]).M(),evWgt,tags2);
/*
      t_HT = HT;
      t_dphi_ll=delta_phi;
      t_scan_mass=ev.scan_mass;
      t_scan_rho=ev.scan_rho;
      t_scan_coup=ev.scan_coup;
      t_weight=evWgt;
*/
    }
    //close input file

    cout<<endl;
     cout<<"Notal: "<<Ntotal<<endl;
     cout<<"Ntotal_after_trigger: "<<Ntotal_after_trig<<"   "<<100.*Ntotal_after_trig/Ntotal<<" %"<<endl;
     cout<<"Ntotal_after_all_selections: "<<N_after_all_selections<<"   "<<100.*N_after_all_selections/Ntotal<<" %"<<endl;

     float r_chargeflip = lepEffH.getChargeflipNormRatio(in_fname.Contains("Data13TeV"));
     float r_DeepCSV = CtagSFWR.GetEvwgtRatio(in_fname.Contains("Data13TeV"),"DeepCSV_alljet");
     float r_DeepJet = CtagSFWR.GetEvwgtRatio(in_fname.Contains("Data13TeV"),"DeepJet_alljet");
     cout<<"Chargeflip ratio r: "<<r_chargeflip<<endl;
     cout<<"DeepCSV ratio r: "<<r_DeepCSV<<endl; // Can be used to scale the plots to corresponding ctagging plot.
     cout<<"DeepJet ratio r: "<<r_DeepJet<<endl;


    f->Close();

  //save histos to file
  fOut->cd();
  for (auto& it : ht.getPlots())  {
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write();
  }
  for (auto& it : ht.get2dPlots())  {
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write();
  }

  //save chargeflip normalized infomation
    TBranch* b_chargeflip_evWgt = t_input.Branch("norm_evWgt_chargeflip", &r_chargeflip, "norm_evWgt_chargeflip/F");
    b_chargeflip_evWgt->Fill();
  

  t_input.Write();
  if (baseMC.Contains("scan",TString::kIgnoreCase)){
     t_300->Write();
     t_350->Write();
     t_400->Write();
     t_500->Write();
     t_600->Write();
     t_700->Write();
  }



  fOut->Close();

}

float DeltaEta(float eta1, float eta2)
{
  float deta = eta2 - eta1;
  return deta;
}
