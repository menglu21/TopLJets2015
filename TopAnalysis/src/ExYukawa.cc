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
  float minLeptonPt( isLowPUrun ? 20. : 27);
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
  fOut->cd();
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("h_Z_mass", new TH1F("h_Z_mass",    ";M(Z) [GeV];Events",62,75,106));
  ht.addHist("nvtx_before_cuts_w", new TH1F("nvtx_before_cuts_w",        ";Vertex multiplicity;Events",35,0,140));
  ht.addHist("nvtx_before_cuts_uw", new TH1F("nvtx_before_cuts_uw",        ";Vertex multiplicity;Events",35,0,140));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",35,0,140));
  ht.addHist("nvtx_w",      new TH1F("nvtx_w",     ";Vertex multiplicity w pileupWeight",35,0,140));
  ht.addHist("nvtx_uw",      new TH1F("nvtx_uw",     ";Vertex multiplicity w/o pileupWeight",35,0,140));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",12,0.5,12.5));
  ht.addHist("njets_bc",        new TH1F("njets_bc",       ";Jet multiplicity;Events",12,0.5,12.5));

  ht.addHist("jet_pt",       new TH1F("jet_pt",      ";p{_T}(jet) [GeV];Events",24,0,600));
  ht.addHist("jet_pt_bc",       new TH1F("jet_pt_bc",      ";p{_T}(jet) [GeV];Events",24,0,600));

  ht.addHist("jet_pt1",       new TH1F("jet_pt1",      ";p{_T}(1st jet) [GeV];Events",24,0,600));
  ht.addHist("jet_pt2",       new TH1F("jet_pt2",      ";p{_T}(2nd jet) [GeV];Events",24,0,600));
  ht.addHist("jet_pt3",       new TH1F("jet_pt3",      ";p{_T}(3rd jet) [GeV];Events",24,0,600));
  ht.addHist("bjet_pt1",       new TH1F("bjet_pt1",      ";p{_T}(1st jet) [GeV];Events",24,0,600));
  ht.addHist("bjet_pt2",       new TH1F("bjet_pt2",      ";p{_T}(2nd jet) [GeV];Events",24,0,600));
  ht.addHist("bjet_pt3",       new TH1F("bjet_pt3",      ";p{_T}(3rd jet) [GeV];Events",24,0,600));

  ht.addHist("HT",            new TH1F("HT",           ";H{_T} [GeV];Events",50,0,1000));
  ht.addHist("HT_bc",            new TH1F("HT_bc",           ";H{_T} [GeV];Events",50,0,1000));

  ht.addHist("nbjets",		 new TH1F("nbjets",      ";b-jet multiplicity;Events",5,0.5,5.5));

  ht.addHist("mllb",          new TH1F("mllb",         ";m(l,b) [GeV];Events",50,0,1000));
  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));
  ht.addHist("nmuons",   new TH1F("nmuons",       ";N(muons);Events",6,-0.5,5.5));
  ht.addHist("mu_q",   new TH1F("mu_q",       ";Muon charge;Events",5,-2.5,2.5));

  ht.addHist("mu_pt",		 new TH1F("mu_pt",       ";p_T(#mu) [GeV]; Events", 24,0,600));
  ht.addHist("mu_pt_bc",		 new TH1F("mu_pt_bc",       ";p_T(#mu) [GeV]; Events", 24,0,600));

  ht.addHist("mu_pt1",            new TH1F("mu_pt1",       ";p_T(Leading #mu) [GeV]; Events", 24,0,600));
  ht.addHist("mu_pt2",            new TH1F("mu_pt2",       ";p_T(Sub-leadeing #mu) [GeV]; Events", 24,0,600));

  ht.addHist("mu_eta",		 new TH1F("mu_eta",      ";#eta(#mu) ; Events", 25,-2.5,2.5));
  ht.addHist("mu_eta_bc",		 new TH1F("mu_eta_bc",      ";#eta(#mu) ; Events", 25,-2.5,2.5));

  ht.addHist("mu_phi1",		 new TH1F("mu_phi1",      ";#phi(#mu1) ; Events", 25,-3.2,3.2));
  ht.addHist("mu_phi2",		 new TH1F("mu_phi2",      ";#phi(#mu2) ; Events", 25,-3.2,3.2));
  ht.addHist("dphi_mumu",        new TH1F("dphi_mumu",   ";#Delta#phi(#mu,#mu);Events", 25, -6,6));
  ht.addHist("deta_mumu",        new TH1F("deta_mumu",   ";#Delta#eta(#mu,#mu);Events", 25, -6,6));
  ht.addHist("bjet_pt",      new TH1F("bjet_pt",     ";p_T(b jet) [GeV]; Events", 24,0,600));
  ht.addHist("bjet_eta",	 new TH1F("bjet_eta",    ";#eta(b jet) ; Events", 25,-2.5,2.5));
  ht.addHist("met",       new TH1F("met",      ";MET [GeV]; Events", 24,0,600));
  ht.addHist("met_phi",       new TH1F("met_phi",      ";MET #phi; Events", 25, -4,4));
  ht.addHist("b_matched_jet",  new TH1F("b_matched_jet", ";p_T(b matched jet) [GeV]; Events", 24,0,600));
  ht.addHist("c_matched_jet",  new TH1F("c_matched_jet", ";p_T(c matched jet) [GeV]; Events", 24,0,600));
  ht.addHist("DR_reco_gen_jet", new TH1F("DR_reco_gen_jet", ";#Delta R; Events", 30,0,1.5));

  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }


  //INPUT
  MiniEvent_t ev;
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }

  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(300000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);

  //EVENT LOOP
  //select 2mu+>=3 jets events triggered by a single muon trigger


  int Ntotal = 0;
  int Ntotal_lepton = 0;
  int Ntotal_after_trig = 0;
  int Ntotal_after_samesign = 0;
  int Ntotal_after_jet = 0;
  int Ntotal_after_num_btags = 0;
  int Ntotal_after_met = 0;
  int N_after_all_selections = 0;

  int Ntotal_Z = 0;
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      //trigger
      bool hasMTrigger(false);
      if(era.Contains("2016")) hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) );
      if(era.Contains("2017")) {
	hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) ||
			selector.hasTriggerBit("HLT_Mu50_v", ev.triggerBits) ||
			selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", ev.triggerBits) ||
			selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits));
      }
      Ntotal++;




//      cout << "Trigger = "<<hasMTrigger << endl;
      if(ev.isData && !hasMTrigger) continue;
      Ntotal_after_trig++;

      //select two offline muons
      std::vector<Particle> flaggedleptons = selector.flaggedLeptons(ev);
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
      std::vector<Particle> leptons = selector.selLeptons(flaggedleptons,muId,SelectionTool::MVA90,minLeptonPt,2.4);
      if(leptons.size()<2) continue;
      Ntotal_lepton++;

      sort(leptons.begin(),leptons.end(),
	[](const Particle& a, const Particle& b){
		return a.Pt() > b.Pt();
	}
      );

//      if (!ev.isData && (leptons[0].id() != 13 || leptons[1].id() != 13)) continue;
      if (leptons[0].id() != 13 || leptons[1].id() != 13) continue;

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,2.4,leptons,{});
      bool passJets(allJets.size()>=minJetMultiplicity);

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);

      //event weight
      float evWgt(1.0);
      float evWgt_test(1.0);

      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,"inc");
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //MC specific: compute event weight

      if (!ev.isData) {

        float normWgt(normH? normH->GetBinContent(1) : 1.0);
        TString period = lumi.assignRunPeriod();
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        EffCorrection_t selSF(1.0,0.0);// = lepEffH.getOfflineCorrection(leptons[0], period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        evWgt  = normWgt*puWgt*selSF.first*l1prefireProb.first;
        evWgt_test  = normWgt*puWgt;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);//generator weights
        evWgt_test *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
      }

      ht.fill("nvtx_before_cuts_w",       ev.nvtx,        evWgt_test, "inc");
      ht.fill("nvtx_before_cuts_uw",       ev.nvtx,        1., "inc");

      for (size_t in_nmu=0; in_nmu<leptons.size();in_nmu++){
        ht.fill("mu_pt_bc",          leptons[in_nmu].pt(), evWgt, "inc");
        ht.fill("mu_eta_bc",     leptons[in_nmu].eta(), evWgt, "inc");
      }

      if (leptons[0].charge()*leptons[1].charge() < 0){
        float zmass = (leptons[0]+leptons[1]).M();
        ht.fill("h_Z_mass", zmass,        evWgt, "inc");
        Ntotal_Z++;
      }

      if (leptons[0].charge()*leptons[1].charge() < 0) continue;
      Ntotal_after_samesign++;

      ht.fill("njets_bc",      allJets.size(), evWgt, "inc");
      float HT_bc = 0;
      for(size_t ij=0; ij<allJets.size(); ij++){
          ht.fill("jet_pt_bc",    allJets[ij].pt(), evWgt, "inc");
          //cout<<"Jet flavor"<<"  "<<allJets[ij].flavor()<<"  "<<allJets[ij].pt()<<"  "<<allJets[ij].eta()<<"  "<<allJets[ij].Phi()<<endl;
          HT_bc += allJets[ij].pt();
      }
      ht.fill("HT_bc", HT_bc, evWgt, "inc");
      if(!passJets) continue;
      Ntotal_after_jet++;

      //lepton-b systems
      int num_btags = 0;
      for(size_t ij=0; ij<allJets.size(); ij++)
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;

		      num_btags++;
        }

       if(num_btags < minNum_btags) continue;
       Ntotal_after_num_btags++;

      if(ev.met_pt < 30.) continue;
      Ntotal_after_met++;
      N_after_all_selections++;

      //sort jets w.r.t. pt
      sort(allJets.begin(),allJets.end(),
        [](const Jet& a, const Jet& b){
                return a.Pt() > b.Pt();
        }
      );


      //select gen-level objects
      std::vector<Particle> gen_leptons=selector.getGenLeptons(ev,20.,2.5);
      std::vector<Particle> gen_photons=selector.getGenPhotons(ev,50.,2.4);
      std::vector<Jet> gen_jets = selector.getGenJets(ev,20,2.5,gen_leptons,gen_photons);

/*
            cout<<"Event: "<<iev<<endl;
            cout<<"Genjets: "<<gen_jets.size()<<endl;

            for (size_t i_gen=0;i_gen<gen_jets.size();i_gen++){
              cout<<i_gen<<"  "<<gen_jets[i_gen].pt()<<"  "<<gen_jets[i_gen].eta()<<"   "<<gen_jets[i_gen].Phi()<<endl;
            }
            */

/*
            for (int i=0;i<ev.ng;i++){
              if (abs(ev.g_id[i])>5) continue;
              cout<<ev.g_id[i]<<"  "<<ev.g_pt[i]<<"  "<<ev.g_eta[i]<<"  "<<ev.g_phi[i]<<endl;
              TLorentzVector genjet4mom;
              genjet4mom.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
            }
*/


     float delta_phi = leptons[0].DeltaPhi(leptons[1]);

     float HT = 0;
     for(size_t ij=0; ij<allJets.size(); ij++){
         ht.fill("jet_pt",    allJets[ij].pt(), evWgt, "inc");
         //cout<<"Jet flavor"<<"  "<<allJets[ij].flavor()<<"  "<<allJets[ij].pt()<<"  "<<allJets[ij].eta()<<"  "<<allJets[ij].Phi()<<endl;
	       HT += allJets[ij].pt();
         int idx=allJets[ij].getJetIndex();
         bool passBtag(ev.j_btag[idx]>0);
         if(!passBtag) continue;
         float mllb( (leptons[0]+leptons[1]+allJets[ij]).M() );//M is the magnitude.
         std::vector<TString> tags={"inc",leptons[0].charge()>0 ? "plus" : "minus"};
         ht.fill("mllb",mllb,evWgt,tags);
         ht.fill("bjet_pt",allJets[ij].pt(),evWgt,"inc");
         ht.fill("bjet_eta",allJets[ij].eta(),evWgt,"inc");
     }

     for(size_t ij=0; ij<allJets.size(); ij++){
       for (int i=0;i<ev.ng;i++){
         if (abs(ev.g_id[i]>5)) continue;
         TLorentzVector genjet4mom;
         genjet4mom.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
//         if (DeltaR(allJets[ij].eta(),ev.g_eta[i],allJets[ij].Phi(),ev.g_phi[i])<0.4 && (allJets[ij].pt()/ev.g_pt[i]) > 0.6 ){
         ht.fill("DR_reco_gen_jet",allJets[ij].DeltaR(genjet4mom),evWgt,"inc");
         if (allJets[ij].DeltaR(genjet4mom) < 0.4 && (allJets[ij].pt()/ev.g_pt[i]) > 0.6 ){
           if (ev.g_id[i] == 4) ht.fill("c_matched_jet",allJets[ij].pt(),evWgt,"inc");
           if (ev.g_id[i] == 5) ht.fill("b_matched_jet",allJets[ij].pt(),evWgt,"inc");
         }
       }
     }

     int num_mu = 0;
//     if (!ev.isData){
        for (size_t ind_mu = 0;ind_mu<leptons.size();ind_mu++){
          if(abs(leptons[ind_mu].id())==13){
            num_mu++;
            ht.fill("mu_q",     leptons[ind_mu].charge(), evWgt, "inc");
          }
        }
//      }



      ht.fill("nvtx",       ev.nvtx,        evWgt, "inc");
      ht.fill("nvtx_w",     ev.nvtx,        evWgt_test, "inc");
      ht.fill("nvtx_uw",    ev.nvtx,        1., "inc");
      ht.fill("njets",      allJets.size(), evWgt, "inc");
      ht.fill("jet_pt1",   allJets[0].pt(), evWgt, "inc");
      ht.fill("jet_pt2",   allJets[1].pt(), evWgt, "inc");
      ht.fill("jet_pt3",   allJets[2].pt(), evWgt, "inc");
      ht.fill("HT", HT, evWgt, "inc");
      for (size_t in_nmu=0; in_nmu<leptons.size();in_nmu++){
      	ht.fill("mu_pt",          leptons[in_nmu].pt(), evWgt, "inc");
      	ht.fill("mu_eta",     leptons[in_nmu].eta(), evWgt, "inc");
      }
      ht.fill("mu_pt1",          leptons[0].pt(), evWgt, "inc");
      ht.fill("mu_pt2",          leptons[1].pt(), evWgt, "inc");
      ht.fill("nmuons",   num_mu, evWgt, "inc");
      ht.fill("dphi_mumu", delta_phi,evWgt,"inc");
      ht.fill("mu_phi1",leptons[0].phi(),evWgt,"inc");
      ht.fill("mu_phi2",leptons[1].phi(),evWgt,"inc");
      ht.fill("deta_mumu", DeltaEta(leptons[0].eta(),leptons[1].eta()),evWgt,"inc");
      ht.fill("nbjets",num_btags,evWgt,"inc");
      ht.fill("met",ev.met_pt,evWgt,"inc");
      ht.fill("met_phi",ev.met_phi,evWgt,"inc");

    }


  int nbins_label = 8;
  TH1I *h_yields = new TH1I("h_yields","h_yields",nbins_label,-0.5,nbins_label-0.5);
  TH1I *h_yields_Z = new TH1I("h_yields_Z","h_yields_Z",nbins_label,-0.5,nbins_label-0.5);


  cout<<endl;
  cout<<"Notal: "<<Ntotal<<endl;
  cout<<"Ntotal_after_trigger: "<<Ntotal_after_trig<<"   "<<100.*Ntotal_after_trig/Ntotal<<" %"<<endl;
  cout<<"Ntotal_after_all_selections: "<<N_after_all_selections<<"   "<<100.*N_after_all_selections/Ntotal<<" %"<<endl;

  const char *labels[nbins_label]  = {"Total","Trig","Lepton","Same Sign","jet","b-tag","MET","All"};
  const char *labels_Z[nbins_label]  = {"Total","Opp. Sign Lep.","Null","Null","Null","Null","Null","Null"};

  h_yields->SetBinContent(1,Ntotal);
  h_yields->SetBinContent(2,Ntotal_after_trig);
  h_yields->SetBinContent(3,Ntotal_lepton);
  h_yields->SetBinContent(4,Ntotal_after_samesign);
  h_yields->SetBinContent(5,Ntotal_after_jet);
  h_yields->SetBinContent(6,Ntotal_after_num_btags);
  h_yields->SetBinContent(7,Ntotal_after_met);
  h_yields->SetBinContent(8,N_after_all_selections);
  for (int i=1;i<=nbins_label;i++) h_yields->GetXaxis()->SetBinLabel(i,labels[i-1]);
  h_yields->SetDirectory(0);

  h_yields_Z->SetBinContent(1,Ntotal);
  h_yields_Z->SetBinContent(2,Ntotal_Z);
  for (int i=3;i<=nbins_label;i++)  h_yields_Z->SetBinContent(i,0.);
  for (int i=1;i<=nbins_label;i++) h_yields_Z->GetXaxis()->SetBinLabel(i,labels_Z[i-1]);
  h_yields_Z->SetDirectory(0);

  //close input file
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

  h_yields->Write();
  h_yields_Z->Write();
  fOut->Close();
}

float DeltaEta(float eta1, float eta2)
{
  float deta = eta2 - eta1;
  return deta;
}