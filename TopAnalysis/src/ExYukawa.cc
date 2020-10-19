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
  EfficiencyScaleFactorsWrapper lepEffH(in_fname.Contains("Data13TeV"),era);

  //CORRECTIONS: L1-prefire
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);

  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

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

  ht.addHist("hf_probb_gen_b", new TH1F("hf_probb_b", ";P[b] w/ Matched GEN b ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_c", new TH1F("hf_probb_c", ";P[c] w/ Matched GEN c ; Events", 20.,-0.5,1.5));
  ht.addHist("hf_probb_gen_lightjet", new TH1F("hf_probb_gen_lightjet", ";P[udsg] w/ Matched GEN light quarks ; Events", 20.,-0.5,1.5));

  ht.addHist("h_scan_mass", new TH1F("h_scan_mass",";particle mass (GeV) ; Events", 20., 200,1000));
  ht.addHist("h_scan_rho", new TH1F("h_scan_rho",";rho ; Events", 10., 0,1));
  ht.addHist("h_scan_coup", new TH1F("h_scan_coup",";coupling ; Events", 4, 0,4));

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

  TTree t_input("TreeInput","TreeInput");
//  int t_event;
  float CvsL1,CvsB1;
  float CvsL2,CvsB2;
  float CvsL3,CvsB3;
  bool btag1,btag2,btag3;


  t_input.Branch("event",&ev.event,"event/I");
  t_input.Branch("run",&ev.run,"run/i");
  t_input.Branch("lumi",&ev.lumi,"lumi/i");
  t_input.Branch("CvsL1",&CvsL1,"CvsL1/F");
  t_input.Branch("CvsB1",&CvsB1,"CvsB1/F");
  t_input.Branch("CvsL2",&CvsL2,"CvsL2/F");
  t_input.Branch("CvsB2",&CvsB2,"CvsB2/F");
  t_input.Branch("CvsL3",&CvsL3,"CvsL3/F");
  t_input.Branch("CvsB3",&CvsB3,"CvsB3/F");
  t_input.Branch("btag1",&btag1,"btag1/O");
  t_input.Branch("btag2",&btag2,"btag2/O");
  t_input.Branch("btag3",&btag3,"btag3/O");


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
      bool hasMTrigger(false);

      if(era.Contains("2016")) hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) );
      if(era.Contains("2017")) {
	hasMTrigger=(
    // previous random selection
    //    selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) ||
    //			  selector.hasTriggerBit("HLT_Mu50_v", ev.triggerBits) ||
    //			  selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits) ||
    //			  selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits) ||
    //    selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits) ||
    //    selector.hasTriggerBit("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v", ev.triggerBits) ||
    //    selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits) || // OK
    //    selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits)

    //From AN2019_140_v3
        // emu triggers
        //selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) || //??
        //selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits) ||
        //selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf", ev.triggerBits) ||
        //selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) ||
        // ee triggers
        //selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits) ||
        selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits) ||
        // mumu triggers
        //selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits) ||
        selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits)
      );
      }

     Ntotal++;
//      cout << "Trigger = "<<hasMTrigger << endl;
//      if(ev.isData && !hasMTrigger) continue;

      ht.fill("h_scan_mass_bc",ev.scan_mass,1);


      if(!hasMTrigger) continue;


      Ntotal_after_trig++;

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
      if (leptons[1].pt() < 20.) continue;
      if (leptons.size() > 2 && leptons[2].pt() > 20.) continue;


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
//      bool passJets(allJets.size()>=minJetMultiplicity);

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

        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);//generator weights
      }




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
      for(size_t ij=0; ij<allJets.size(); ij++)
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;
		      num_btags++;
        }

//      if(num_btags < minNum_btags) continue;

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

  if (leptons[0].charge()*leptons[1].charge() < 0 && zmass > 70. && zmass < 110.) ht.fill("njets_bc",      jets.size(), evWgt, tags2);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (leptons[0].charge()*leptons[1].charge() < 0){
    if (zmass > 70. && zmass < 110.){
      for(size_t ij=0; ij<jets.size(); ij++){
          ht.fill("jet_pt_bc",    jets[ij].pt(), evWgt, tags2);
          ht.fill("jet_eta_bc",   jets[ij].eta(), evWgt, tags2);
      }
    }
  }

  if (leptons[0].charge()*leptons[1].charge() < 0) continue;


  invariant_mass = (leptons[0]+leptons[1]).M();
  std::vector<TString> tags3={"inc"};
  tags3.push_back(to_string(ev.scan_mass));

  if(ev.met_pt < 30.) continue;

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

  bool passJets(jets.size()>=minJetMultiplicity);
  if(!passJets) continue;

  if(num_btags < minNum_btags) continue;

  int jet_index=0;
  btag1=false;
  btag2=false;
  btag3=false;
  for(size_t ij=0; ij<jets.size(); ij++){
    int idx=jets[ij].getJetIndex();
    bool passBtag(ev.j_btag[idx]>0);
    if (jet_index==0){
      CvsL1 = ev.j_CvsL[idx];
      CvsB1 = ev.j_CvsB[idx];
      btag1=passBtag;
    }
    if (jet_index==1){
      CvsL2 = ev.j_CvsL[idx];
      CvsB2 = ev.j_CvsB[idx];
      btag2=passBtag;
    }
    if (jet_index==2){
      CvsL3 = ev.j_CvsL[idx];
      CvsB3 = ev.j_CvsB[idx];
      btag3=passBtag;
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

     sort(jets.begin(),jets.end(),
       [](const Jet& a, const Jet& b){
               return a.Pt() > b.Pt();
       }
     );

     float delta_phi = leptons[0].DeltaPhi(leptons[1]);

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


    }
    //close input file

    cout<<endl;
     cout<<"Notal: "<<Ntotal<<endl;
     cout<<"Ntotal_after_trigger: "<<Ntotal_after_trig<<"   "<<100.*Ntotal_after_trig/Ntotal<<" %"<<endl;
     cout<<"Ntotal_after_all_selections: "<<N_after_all_selections<<"   "<<100.*N_after_all_selections/Ntotal<<" %"<<endl;

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

  fOut->Close();
}

float DeltaEta(float eta1, float eta2)
{
  float deta = eta2 - eta1;
  return deta;
}
