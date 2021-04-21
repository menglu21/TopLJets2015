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
#include "DataFormats/MuonReco/interface/Muon.h"

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
enum QualityFlags {VETO, LOOSE, MEDIUM, TIGHT, CONTROL, QCDTEMP, RELAXEDTIGHT, MVA80, MVA90,LOOSEIDONLY,MEDIUMIDONLY,TIGHTIDONLY,MVANONISOWPLOOSE,HIGHPT,HIGHPTIDONLY,TIGHTIDNOSIHIH};

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
  float t_reliso_l1, mt;
  int t_id_l1, t_id_l2;
  float t_pt_j1,t_pt_j2,t_pt_j3,t_eta_j1,t_eta_j2,t_eta_j3;
  float t_phi_j1,t_phi_j2,t_phi_j3;
  float t_met;
  float t_weight, t_normWgt, t_norm;
  float t_scan_mass, t_scan_rho, t_scan_coup;
  int Njets;

  t_input.Branch("event",&ev.event,"event/I");
  t_input.Branch("run",&ev.run,"run/i");
  t_input.Branch("lumi",&ev.lumi,"lumi/i");
  t_input.Branch("t_weight",&t_weight,"t_weight/F");
  t_input.Branch("t_normWgt",&t_normWgt,"t_normWgt/F");
  t_input.Branch("t_norm",&t_norm,"t_norm/F");

  t_input.Branch("t_id_l1",&t_id_l1,"t_id_l1/I");
  t_input.Branch("t_id_l2",&t_id_l2,"t_id_l2/I");
  t_input.Branch("Njets",&Njets,"Njets/I");

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


  t_input.Branch("t_reliso_l1", &t_reliso_l1, "t_reliso_l1/F");
  t_input.Branch("mt", &mt, "mt/F");
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
  int passtrigger_sinmu1 =0;
  for (Int_t iev=0;iev<nentries;iev++)
//  for (Int_t iev=0;iev<10;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
      //trigger

      passtrigger_sinmu1=0;
      if(era.Contains("2017")) {
        if (baseMC.Contains("SingleMuon",TString::kIgnoreCase)){
	  passtrigger_sinmu1 = selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits);
        }
     }

     if (!passtrigger_sinmu1) continue;

      ht.fill("h_scan_mass_bc",ev.scan_mass,1);

//      //select two offline muons
//      std::vector<Particle> flaggedleptons = selector.flaggedLeptons(ev);
//      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
//      //std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::MVA90,minLeptonPt,2.4);
//      std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::TIGHT,minLeptonPt,2.4);

	// select tight ID muon but with iso<0.15 or 0.2<iso<0.4
	std::vector<Particle> leptons;
	std::vector<Particle> veto_muons;
	std::vector<Particle> veto_electrons;
	std::vector<float> rel_isos;
	for (int il=0; il<ev.nl; il++) {
		TLorentzVector lp4_temp;
		float pt_temp(ev.l_pt[il]);
		float eta_temp(fabs(ev.l_eta[il]));
		unsigned int pid_temp(ev.l_pid[il]);
		float relIso_temp(ev.l_relIso[il]);
		int qualityFlagsWord_temp(0);
		Float_t unc_temp(0.);
		if (abs(ev.l_id[il])==13 && pt_temp>20 && eta_temp<2.5)
		{
		  if( (pid_temp&reco::Muon::Selector::CutBasedIdTight)==reco::Muon::Selector::CutBasedIdTight ){
			if ((relIso_temp>0.2 && relIso_temp<0.4) ||relIso_temp<0.15)
			{qualityFlagsWord_temp |= (0x1 << TIGHT);
			 rel_isos.push_back(relIso_temp);
			 lp4_temp.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
			 leptons.push_back(Particle(lp4_temp, ev.l_charge[il], ev.l_id[il], qualityFlagsWord_temp, il, 1.0, unc_temp));}
			if (relIso_temp>0.15 &&relIso_temp<0.2)
			{
			 qualityFlagsWord_temp |= (0x1 << TIGHT);
			 lp4_temp.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
                         veto_muons.push_back(Particle(lp4_temp, ev.l_charge[il], ev.l_id[il], qualityFlagsWord_temp, il, 1.0, unc_temp));
			}
		}}
		if (abs(ev.l_id[il])==11 && pt_temp>20 && eta_temp<2.4)
		{
		  if((pid_temp>>1)&0x1)
			{
			  qualityFlagsWord_temp |= (0x1 << VETO);
			  lp4_temp.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
			  veto_electrons.push_back(Particle(lp4_temp, ev.l_charge[il], ev.l_id[il], qualityFlagsWord_temp, il, 1.0, unc_temp));
			}
		}
	}


      if(leptons.size()!=1) continue;
      if(veto_electrons.size()>0) continue;
      if(veto_muons.size()>0) continue;

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,4.7,leptons,{});
      if (allJets.size()!=1)continue;

      if (ev.met_pt>30) continue;
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);

      //lepton transverse
      TLorentzVector lep_T(0,0,0,0);
      lep_T.SetPxPyPzE(leptons[0].Px(),leptons[0].Py(),0,leptons[0].E()*leptons[0].Pt()/leptons[0].P());

      mt=(met+lep_T).M();
      if(mt>20) continue;

      //event weight
      float evWgt(1.0);

      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }


      ht.fill("met",ev.met_pt,evWgt,"inc");
      //if (evWgt < 0.) continue;
	

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    std::vector<Jet> jets;
    for(auto &j: allJets) {

      int idx=j.getJetIndex();

      jets.push_back(j);
    }
      Njets=jets.size();

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  float HT = 0;

t_pt_l1=leptons[0].pt();
t_eta_l1=leptons[0].eta();
t_phi_l1=leptons[0].phi();
t_reliso_l1=rel_isos[0];
t_met=ev.met_pt;
t_norm=norm;
t_HT = HT;
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
