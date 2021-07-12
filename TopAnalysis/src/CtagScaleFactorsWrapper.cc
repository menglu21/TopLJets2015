#include "TopLJets2015/TopAnalysis/interface/CtagScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"


#include "TFile.h"
#include "TKey.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
CtagScaleFactorsWrapper::CtagScaleFactorsWrapper(bool isData,TString era)
{
  if(isData) return;
  init(era);
}

void CtagScaleFactorsWrapper::init(TString era)
{

  cout << "[CtagScaleFactorsWrapper] with efficiencies for " << era << endl;

  TString url(era+"/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2.root");
  gSystem->ExpandPathName(url);
  TFile *fIn=TFile::Open(url);
  TString key_c("SFc_hist");
  TString key_b("SFb_hist");
  TString key_l("SFl_hist");
  ScaleFactorH_["DeepCSV_c"]=(TH2F *) fIn->Get(key_c);
  ScaleFactorH_["DeepCSV_c"]->SetDirectory(0);
  ScaleFactorH_["DeepCSV_b"]=(TH2F *) fIn->Get(key_b);
  ScaleFactorH_["DeepCSV_b"]->SetDirectory(0);
  ScaleFactorH_["DeepCSV_l"]=(TH2F *) fIn->Get(key_l);
  ScaleFactorH_["DeepCSV_l"]->SetDirectory(0);
  fIn->Close();

  TString url_deepjet(era+"/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2.root");
  gSystem->ExpandPathName(url_deepjet);
  TFile *fIn_deepjet=TFile::Open(url_deepjet);
  ScaleFactorH_["DeepJet_c"]=(TH2F *) fIn_deepjet->Get(key_c);
  ScaleFactorH_["DeepJet_c"]->SetDirectory(0);
  ScaleFactorH_["DeepJet_b"]=(TH2F *) fIn_deepjet->Get(key_b);
  ScaleFactorH_["DeepJet_b"]->SetDirectory(0);
  ScaleFactorH_["DeepJet_l"]=(TH2F *) fIn_deepjet->Get(key_l);
  ScaleFactorH_["DeepJet_l"]->SetDirectory(0);
  fIn_deepjet->Close();

}

float CtagScaleFactorsWrapper::GetCtagSF(TString algorithm_type, std::vector<Jet> &jets, MiniEvent_t ev)
{
  if (ev.isData) return 1.0;
  float ctagWt = 1.;
  for (size_t ij=0; ij<jets.size();ij++) {
    float SF = GetCtagSF(algorithm_type, jets[ij],ev);
    ctagWt*=SF;
  }
  return ctagWt;

}

float CtagScaleFactorsWrapper::GetCtagSF(TString algorithm_type, Jet jet, MiniEvent_t ev)
{
  if(ev.isData) return 1.0;
  float ctagWt = 1.;
  int flav, xbin, ybin;
  float CvsLval,CvsBval;
  TH2F *wtHist;

  if(algorithm_type.Contains("DeepCSV")){
    int idx = jet.getJetIndex();
    flav = int(ev.j_hadflav[idx]);
    CvsLval = ev.j_CvsL[idx];
    CvsBval = ev.j_CvsB[idx];
    if (flav == 4) {
            wtHist = ScaleFactorH_["DeepCSV_c"];
    }
    else if (flav == 5) {
            wtHist = ScaleFactorH_["DeepCSV_b"];
    }
    else {
            wtHist = ScaleFactorH_["DeepCSV_l"];
    }
    xbin = wtHist->GetXaxis()->FindBin(CvsLval);
    ybin = wtHist->GetYaxis()->FindBin(CvsBval);
    ctagWt *= wtHist->GetBinContent(xbin,ybin);
  }
  if(algorithm_type.Contains("DeepJet")){
    int idx = jet.getJetIndex();
    flav = int(ev.j_hadflav[idx]);
    CvsLval = ev.j_deepjet_CvsL[idx];
    CvsBval = ev.j_deepjet_CvsB[idx];
    if (flav == 4) {
            wtHist = ScaleFactorH_["DeepJet_c"];
    }
    else if (flav == 5) {
            wtHist = ScaleFactorH_["DeepJet_b"];
    }
    else {
            wtHist = ScaleFactorH_["DeepJet_l"];
    }
    xbin = wtHist->GetXaxis()->FindBin(CvsLval);
    ybin = wtHist->GetYaxis()->FindBin(CvsBval);
    ctagWt *= wtHist->GetBinContent(xbin,ybin);
  }
  return ctagWt;
}


CtagScaleFactorsWrapper::~CtagScaleFactorsWrapper()
{
}
