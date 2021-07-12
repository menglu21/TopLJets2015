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
  float GetCtagSF(TString algorithm_type, std::vector<Jet> &jets, MiniEvent_t ev);
  float GetCtagSF(TString algorithm_type, Jet jet, MiniEvent_t ev);
  ~CtagScaleFactorsWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2F *> ScaleFactorH_;
};

#endif
