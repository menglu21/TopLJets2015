# Charge Flip Rate Study

## Produce Analysis Files

To produce the root files used to analysis the charge flip rate, please run the following code.

```
cd $CMSSW_BASE/src/TopLJets2015/TopAnalysis
scram b -j 8
sh test/analysis/ExYukawa/ChargeFlipStudy/steerAnalysis_ChargeFlip.sh -o SEL

# You can check the main production code in $CMSSW_BASE/src/TopLJets2015/TopAnalysis/src/ChargeFlip.cc

cd test/analysis/ExYukawa/ChargeFlipStudy

```

## Analysis charge flip rate

There are several method to analysis. Please edit the plot directory in the run\_all\_CF.sh. And also edit the plotdir in following command.

```
# Analysis Z mass distribution
python Zmassdistribution.py [plotdir]
# Analysis CF rate with detailed kinematic region
python SB_Detail_Region.py [plotdir] [1 or 0] # Consider Covariance : 1. Otherwise : 0.
# Analysis CF rate with combine kinematic region
python SB_Combine_Region.py [plotdir] [1 or 0]
# Analysis CF rate with independent model (P = f(Pt)g(eta))
python SB_seperateModel.py [plotdir] [1 or 0]
# If you just want to run all of them
sh run_all_CF.sh
```

