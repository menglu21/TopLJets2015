#ifndef _EfficiencyScaleFactorsWrapper_h_
#define _EfficiencyScaleFactorsWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

typedef std::pair<float,float> EffCorrection_t;

class EfficiencyScaleFactorsWrapper
{
 public:
  EfficiencyScaleFactorsWrapper(bool isData,TString era);
  EffCorrection_t getMuonSF(float pt,float eta);
  EffCorrection_t getElectronSF(float pt,float eta);
  EffCorrection_t getPhotonSF(float pt,float eta);
  EffCorrection_t getMuMuPtSF(float pt1, float pt2);
  EffCorrection_t getMuMuEtaSF(float eta1, float eta2);
  EffCorrection_t getEMuPtSF(float pt1, float pt2);
  EffCorrection_t getEMuEtaSF(float eta1, float eta2);
  EffCorrection_t getEEPtSF(float pt1, float pt2);
  EffCorrection_t getEEEtaSF(float eta1, float eta2);
  ~EfficiencyScaleFactorsWrapper();

 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2 *> scaleFactorsH_;
  std::map<TString,TH1 *> scaleFactors1DH_;
  std::map<TString,TGraphAsymmErrors *> scaleFactorsGr_;
};

#endif
