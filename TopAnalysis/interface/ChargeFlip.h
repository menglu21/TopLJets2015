#ifndef _ChargeFlip_h_
#define _ChargeFlip_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunChargeFlip(const TString filename,
                      TString outname,
                      TH1F *normH,
                      TH1F *genPU,
                      TString era,
                      Bool_t debug=false);

#endif
