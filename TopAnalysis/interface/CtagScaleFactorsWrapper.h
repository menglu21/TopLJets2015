#ifndef _CtagScaleFactorsWrapper_h_
#define _CtagScaleFactorsWrapper_h_

#include "TH2.h"
#include "TString.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

class CtagScaleFactorsWrapper 
{
 public:
  CtagScaleFactorsWrapper(bool isData,TString era);
  float GetCtagSF(TString name, TString algorithm_type, std::vector<Jet> &jets, MiniEvent_t ev, float evWgt);
  float GetCtagSF(TString name, TString algorithm_type, Jet jet, MiniEvent_t ev, float evWgt);
  float GetEvwgtRatio(bool isData, TString name);
  ~CtagScaleFactorsWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2F *> ScaleFactorH_;
  std::map<TString,float > evWgt_woSF_H_;
  std::map<TString,float > evWgt_SF_H_;
};

#endif
