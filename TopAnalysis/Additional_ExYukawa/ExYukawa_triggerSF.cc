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

  //READ Cross section over sum(weight) for the MC sample
  TString dirin = gSystem->DirName(in_fname);
  TString baseMC = gSystem->BaseName(dirin);  cout<<"Base MC name: "<<baseMC<<endl;
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
//  norm_file->Close();


  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname);
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
////////  TFile *f_tmva=TFile::Open("tree_tmva.root","recreate");
  fOut->cd();



  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("control_Z_mass",    new TH1F("control_Z_mass",    ";M(Z) [GeV];Events",40,70,110));
  ht.addHist("control_lep_pt_bc", new TH1F("control_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("control_lep_eta_bc",new TH1F("control_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));
  ht.addHist("control_m_ll_bc",		new TH1F("control_m_ll_bc",   ";M(l+,l-) [GeV] ; Events", 20,12,600));
  ht.addHist("control_pT_ll_bc",  new TH1F("control_pT_ll_bc",  ";p_{T}(l+,l-) [GeV] ; Events", 80,0,200));
  ht.addHist("control_y_ll_bc",  new TH1F("control_y_ll_bc",  ";y(l+,l-) [GeV] ; Events", 20,-5,5));

  ht.addHist("control_2_lep_pt_bc", new TH1F("control_2_lep_pt_bc", ";p_{T}(l) [GeV]; Events", 30,0,600));
  ht.addHist("control_2_lep_eta_bc",new TH1F("control_2_lep_eta_bc",";#eta(lepton) ; Events", 10,-2.5,2.5));

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

  ht.addHist("hf_probb_gen_b", new TH1F("hf_probb_b", ";P[b] w/ Matched GEN b ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_c", new TH1F("hf_probb_c", ";P[c] w/ Matched GEN c ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_lightjet", new TH1F("hf_probb_gen_lightjet", ";P[udsg] w/ Matched GEN light quarks ; Events", 20.,-0.5,1.5));

  ht.addHist("h_scan_mass", new TH1F("h_scan_mass",";particle mass (GeV) ; Events", 20., 200,1000));
  ht.addHist("h_scan_rho", new TH1F("h_scan_rho",";rho ; Events", 10., 0,1));
  ht.addHist("h_scan_coup", new TH1F("h_scan_coup",";coupling ; Events", 4, 0,4));

  ht.addHist("h_m_top_charm", new TH1F("h_m_top_charm",";m(top,charm) ; Events", 100, 0,1000));
  ht.addHist("h_m_top_charm_x", new TH1F("h_m_top_charm_x",";m(top,charm) ; Events", 100, 0,1000));


//  TH1F *a_test1 = new TH1F("a_test1","a_test1",30,0,60);//for debugging

//  ht.addHist("b_matched_jet",  new TH1F("b_matched_jet", ";p_T(b matched jet) [GeV]; Events", 24,0,600));
//  ht.addHist("c_matched_jet",  new TH1F("c_matched_jet", ";p_T(c matched jet) [GeV]; Events", 24,0,600));

  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
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
  float t_pt_l1,t_pt_l2,t_eta_l1,t_eta_l2,t_phi_l1,t_phi_l2;
  int t_id_l1, t_id_l2;
  float t_pt_j1,t_pt_j2,t_pt_j3,t_eta_j1,t_eta_j2,t_eta_j3;
  float t_phi_j1,t_phi_j2,t_phi_j3;
  float t_met;
  float t_weight, t_normWgt, t_norm;
  float t_scan_mass, t_scan_rho, t_scan_coup;
  int passtrigger_met1;
  int passtrigger_met2;
  int passtrigger_met3;
  int passtrigger_met4;
  int passtrigger_met5;
  int passtrigger_dimu1;
  int passtrigger_dimu2;
  int passtrigger_sinmu1;
  int passtrigger_diele1;
  int passtrigger_diele2;
  int passtrigger_sinele1;
  int passtrigger_sinele2;
  int passtrigger_emu1;
  int passtrigger_emu2;
  int passtrigger_emu3;
  int passtrigger_emu4;
  int Njets;
  int Flag_HBHENoiseFilter;
  int Flag_HBHENoiseIsoFilter;
  int Flag_EcalDeadCellTriggerPrimitiveFilter;
  int Flag_goodVertices;
  int Flag_eeBadScFilter;
  int Flag_globalTightHalo2016Filter;

  t_input.Branch("event",&ev.event,"event/I");
  t_input.Branch("run",&ev.run,"run/i");
  t_input.Branch("lumi",&ev.lumi,"lumi/i");
  t_input.Branch("t_weight",&t_weight,"t_weight/F");
  t_input.Branch("t_normWgt",&t_normWgt,"t_normWgt/F");
  t_input.Branch("t_norm",&t_norm,"t_norm/F");

  t_input.Branch("passtrigger_met1",&passtrigger_met1,"passtrigger_met1/I");
  t_input.Branch("passtrigger_met2",&passtrigger_met2,"passtrigger_met2/I");
  t_input.Branch("passtrigger_met3",&passtrigger_met3,"passtrigger_met3/I");
  t_input.Branch("passtrigger_met4",&passtrigger_met4,"passtrigger_met4/I");
  t_input.Branch("passtrigger_met5",&passtrigger_met5,"passtrigger_met5/I");
  t_input.Branch("passtrigger_dimu1",&passtrigger_dimu1,"passtrigger_dimu1/I");
  t_input.Branch("passtrigger_dimu2",&passtrigger_dimu2,"passtrigger_dimu2/I");
  t_input.Branch("passtrigger_sinmu1",&passtrigger_sinmu1,"passtrigger_sinmu1/I");
  t_input.Branch("passtrigger_diele1",&passtrigger_diele1,"passtrigger_diele1/I");
  t_input.Branch("passtrigger_diele2",&passtrigger_diele2,"passtrigger_diele2/I");
  t_input.Branch("passtrigger_sinele1",&passtrigger_sinele1,"passtrigger_sinele1/I");
  t_input.Branch("passtrigger_sinele2",&passtrigger_sinele2,"passtrigger_sinele2/I");
  t_input.Branch("passtrigger_emu1",&passtrigger_emu1,"passtrigger_emu1/I");
  t_input.Branch("passtrigger_emu2",&passtrigger_emu2,"passtrigger_emu2/I");
  t_input.Branch("passtrigger_emu3",&passtrigger_emu3,"passtrigger_emu3/I");
  t_input.Branch("passtrigger_emu4",&passtrigger_emu4,"passtrigger_emu4/I");
  t_input.Branch("t_id_l1",&t_id_l1,"t_id_l1/I");
  t_input.Branch("t_id_l2",&t_id_l2,"t_id_l2/I");
  t_input.Branch("Njets",&Njets,"Njets/I");
  t_input.Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/I");
  t_input.Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/I");
  t_input.Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/I");
  t_input.Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/I");
  t_input.Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter,"Flag_eeBadScFilter/I");
  t_input.Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/I");

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

  /*
  MET trigger list:
met1: HLT_PFMET120_PFMHT120_IDTight_v*
met2: HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*
met3: HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*
met4: HLT_PFMET250_HBHECleaned_v*
met5: HLT_PFHT500_PFMET100_PFMHT100_IDTight_v* 

  di-Mu trigger list:
dimu1: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v
dimu2: (run<299330)HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v, HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v

  single-Mu trigger list:
sinmu1: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v

  di-Ele trigger list:
diele1: HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
diele2: HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v

  single-ele trigger list:
sinele1: HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
sinele2: HLT_Ele35_WPTight_Gsf_v

  emu trigger list:
emu1: HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
emu2: HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
emu3: HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
emu4: HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
  */
  int Ntotal = 0;
  int Ntotal1 = 0;
  int Ntotal2 = 0;
  int Ntotal3 = 0;
  int Ntotal4 = 0;
  int Ntotal5 = 0;
  int N_after_all_selections = 0;
  int Ntotal_after_trig = 0;

  for (Int_t iev=0;iev<nentries;iev++)
//  for (Int_t iev=0;iev<10;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
      //trigger
      passtrigger_met1 = 0;
      passtrigger_met2 = 0;
      passtrigger_met3 = 0;
      passtrigger_met4 = 0;
      passtrigger_met5 = 0;
      passtrigger_dimu1 = 0;
      passtrigger_dimu2 = 0;
      int passtrigger_dimu2_temp1 = 0;
      int passtrigger_dimu2_temp2 = 0;
      passtrigger_sinmu1 = 0;
      passtrigger_diele1 = 0;
      passtrigger_diele2 = 0;
      passtrigger_sinele1 = 0;
      passtrigger_sinele2 = 0;
      passtrigger_emu1 = 0;
      passtrigger_emu2 = 0;
      passtrigger_emu3 = 0;
      passtrigger_emu4 = 0;

      if(era.Contains("2017")) {
        if (baseMC.Contains("MET",TString::kIgnoreCase)){
	  passtrigger_met1 = selector.hasTriggerBit("HLT_PFMET120_PFMHT120_IDTight_v", ev.addTriggerBits);
	  passtrigger_met2 = selector.hasTriggerBit("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v", ev.addTriggerBits);
	  passtrigger_met3 = selector.hasTriggerBit("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v", ev.addTriggerBits);
	  passtrigger_met4 = selector.hasTriggerBit("HLT_PFMET250_HBHECleaned_v", ev.addTriggerBits);
//	  passtrigger_met5 = selector.hasTriggerBit("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v", ev.addTriggerBits);
	  if (ev.addTriggerBits<0) passtrigger_met5=1;
	  passtrigger_dimu1 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits);
	  if(ev.run < 299330){ passtrigger_dimu2 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", ev.triggerBits);}
          if(ev.run > 299330){ passtrigger_dimu2 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits);}  
	  passtrigger_sinmu1 = selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits);
	  passtrigger_diele1 = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
	  passtrigger_diele2 = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
	  passtrigger_sinele1 = selector.hasTriggerBit("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v", ev.triggerBits);
	  passtrigger_sinele2 = selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits);
          passtrigger_emu1 = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          passtrigger_emu2 = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          passtrigger_emu3 = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          passtrigger_emu4 = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
        }
        if (baseMC.Contains("MC",TString::kIgnoreCase)){
	  passtrigger_met1 = selector.hasTriggerBit("HLT_PFMET120_PFMHT120_IDTight_v", ev.addTriggerBits);
	  passtrigger_met2 = selector.hasTriggerBit("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v", ev.addTriggerBits);
	  passtrigger_met3 = selector.hasTriggerBit("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v", ev.addTriggerBits);
	  passtrigger_met4 = selector.hasTriggerBit("HLT_PFMET250_HBHECleaned_v", ev.addTriggerBits);
//	  passtrigger_met5 = selector.hasTriggerBit("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v", ev.addTriggerBits);
	  if (ev.addTriggerBits<0) passtrigger_met5=1;
	  passtrigger_dimu1 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits);
	  passtrigger_dimu2_temp1 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", ev.triggerBits);
          passtrigger_dimu2_temp2 = selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits); 
 	  if (passtrigger_dimu2_temp1 || passtrigger_dimu2_temp2) passtrigger_dimu2=1;
	  passtrigger_sinmu1 = selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits);
	  passtrigger_diele1 = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
	  passtrigger_diele2 = selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
	  passtrigger_sinele1 = selector.hasTriggerBit("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v", ev.triggerBits);
	  passtrigger_sinele2 = selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits);
	  
          passtrigger_emu1 = selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          passtrigger_emu2 = selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
          passtrigger_emu3 = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
          passtrigger_emu4 = selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
        }
     }

//	if (!(passtrigger_met1 || passtrigger_met2 || passtrigger_met3 || passtrigger_met4 || passtrigger_met5)) continue;

     Ntotal++;

      ht.fill("h_scan_mass_bc",ev.scan_mass,1);

      //select two offline muons
      std::vector<Particle> flaggedleptons = selector.flaggedLeptons(ev);
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
      //std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::MVA90,minLeptonPt,2.4);
      std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::TIGHT,minLeptonPt,2.4);


      if(leptons.size()<2) continue;
     Ntotal1++;

      sort(leptons.begin(),leptons.end(),
	            [](const Particle& a, const Particle& b){
		          return a.Pt() > b.Pt();
	           }
      );

      if (leptons[0].pt() < 30.) continue;
     Ntotal2++;
      //if (leptons[0].pt() < 40. && leptons[0].id() == 11) continue;
      if (leptons[1].pt() < 20.) continue;
     Ntotal3++;
      if (leptons.size() > 2 && leptons[2].pt() > 20.) continue;
     Ntotal4++;

      t_id_l1=leptons[0].id();
      t_id_l2=leptons[1].id();


      int dimuon_event = 0;
      int dielectron_event = 0;
      int emu_event = 0;
      int mue_event = 0;
      if (leptons[0].id() == 11 && leptons[1].id() == 11) dielectron_event = 1;
      if (leptons[0].id() == 13 && leptons[1].id() == 13) dimuon_event = 1;
      if (leptons[0].id() == 11 && leptons[1].id() == 13) emu_event = 1;
      if (leptons[0].id() == 13 && leptons[1].id() == 11) mue_event = 1;

      std::vector<TString> tags2={"inc"};
      if (dimuon_event) tags2.push_back("mm");
      if (dielectron_event) tags2.push_back("ee");
      if (emu_event || mue_event) tags2.push_back("emu");

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);
      std::vector<Jet> allJets = selector.getGoodJets(ev,20.,2.4,leptons,{});
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

        evWgt  = normWgt*puWgt*l1prefireProb.first*leptonSF*dilepton_trig_SF;
        t_normWgt = normWgt;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);//generator weights
      }


      ht.fill("met",ev.met_pt,evWgt,"inc");
      //if (evWgt < 0.) continue;
	

      float invariant_mass = (leptons[0]+leptons[1]).M();
      ht.fill("control_m_ll_bc",  invariant_mass,        evWgt, tags2);

      float zmass = (leptons[0]+leptons[1]).M();
      if (leptons[0].charge()*leptons[1].charge() < 0){
        ht.fill("control_Z_mass", zmass,        evWgt, tags2);
        if (zmass > 70. && zmass < 110.){
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
        }
      }

      int num_btags = 0;
      for(size_t ij=0; ij<allJets.size(); ij++){
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;
		      num_btags++;
        }

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
      Njets=jets.size();
      Flag_HBHENoiseFilter=(ev.met_filterBits)&0x1;
      Flag_HBHENoiseIsoFilter=(ev.met_filterBits>>1)&0x1;
      Flag_EcalDeadCellTriggerPrimitiveFilter=(ev.met_filterBits>>2)&0x1;
      Flag_goodVertices=(ev.met_filterBits>>3)&0x1;
      Flag_eeBadScFilter=(ev.met_filterBits>>4)&0x1;
      Flag_globalTightHalo2016Filter=(ev.met_filterBits>>5)&0x1;

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

  //if (leptons[0].charge()*leptons[1].charge() < 0) continue;
  if (leptons[0].charge()*leptons[1].charge() < 0) Ntotal5++;

  invariant_mass = (leptons[0]+leptons[1]).M();
  std::vector<TString> tags3={"inc"};
  tags3.push_back(to_string(ev.scan_mass));



 // if(ev.met_pt < 30.) continue;

//  if (jets.size() < 3) continue;

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
t_norm=norm;
t_HT = HT;
t_dphi_ll=delta_phi;
t_scan_mass=ev.scan_mass;
t_scan_rho=ev.scan_rho;
t_scan_coup=ev.scan_coup;
t_weight=evWgt;

//  bool passJets(jets.size()>=minJetMultiplicity);
 // if(!passJets) continue;

 // if(num_btags < minNum_btags) continue;

 sort(jets.begin(),jets.end(),
   [](const Jet& a, const Jet& b){
           return a.Pt() > b.Pt();
   }
 );


  int jet_index=0;
  for(size_t ij=0; ij<jets.size(); ij++){
    int idx=jets[ij].getJetIndex();
    bool passBtag(ev.j_btag[idx]>0);
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
    ht.fill("hf_csv",ev.j_csv[idx],evWgt,tags2);
    ht.fill("hf_deepcsv",ev.j_deepcsv[idx],evWgt,tags2);
    ht.fill("hf_probc",ev.j_probc[idx],evWgt,tags2);
    ht.fill("hf_probudsg",ev.j_probudsg[idx],evWgt,tags2);
    ht.fill("hf_probb",ev.j_probb[idx],evWgt,tags2);
    ht.fill("hf_probbb",ev.j_probbb[idx],evWgt,tags2);
    ht.fill("hf_CvsL",ev.j_CvsL[idx],evWgt,tags2);
    ht.fill("hf_CvsB",ev.j_CvsB[idx],evWgt,tags2);
    ht.fill2D("hf_CvsL_vs_CvsB",ev.j_CvsL[idx],ev.j_CvsB[idx],evWgt,tags2);
    t_deepcsv=ev.j_deepcsv[idx];





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
//      ht.fill("jet_pt1",   jets[0].pt(), evWgt, tags2);
//      ht.fill("jet_pt2",   jets[1].pt(), evWgt, tags2);
//      ht.fill("jet_pt3",   jets[2].pt(), evWgt, tags2);
//      ht.fill("HT", HT, evWgt, tags2);

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
//      ht.fill("nbjets",num_btags,evWgt,tags2);
//      ht.fill("met",ev.met_pt,evWgt,tags2);
      ht.fill("met_phi",ev.met_phi,evWgt,tags2);
      ht.fill("h_scan_mass",ev.scan_mass,1);
      ht.fill("h_scan_rho",ev.scan_rho,evWgt);
      ht.fill("h_scan_coup",ev.scan_coup,evWgt);

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
//    cout<<Ntotal<<" "<<Ntotal1<<" "<<Ntotal2<<" "<<Ntotal3<<" "<<Ntotal4<<" "<<Ntotal5<<endl;

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
