# TopLJets2015

## Installation instructions

These installation instructions correspond to the 2017 data/MC production.
To install execute the following in your work area.
Notice: if you are not creating the ntuples, you can skip the part of the instructions
marked with the `##OPTIONAL/##END OPTIONAL` markers.
If compilation fails for some reason repeat the scram b...

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_16
cd CMSSW_10_6_16/src
cmsenv
git cms-init

#EGM id
git cms-merge-topic jainshilpi:ULV1_backport10616_forUsers
git clone https://github.com/jainshilpi/EgammaPostRecoTools.git -b ULV0
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b UL2017SSV2 EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
scram b -j 8

#B-fragmentation analyzer
mkdir TopQuarkAnalysis
cd TopQuarkAnalysis
git clone -b 94x https://gitlab.cern.ch/psilva/BFragmentationAnalyzer.git
scram b -j 8
cd -

#This package
cd $CMSSW_BASE/src
git clone https://github.com/efeyazgan/TopLJets2015.git -b 106_protonreco
cd TopLJets2015

scram b -j 8
```

## Running ntuple creation and checking the selection

The ntuplizer is steered with test/runMiniAnalyzer_cfg.py.
It takes several options from command line (see cfg for details).
To run locally the ntuplizer, for testing purposes do something like:

```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False era=era2017 outFilename=MC13TeV_TTJets.root noParticleLevel=True
cmsRun test/runMiniAnalyzer_cfg.py runOnData=True  era=era2017 outFilename=Data13TeV_DoubleEG_UL.root
cmsRun test/runL1PrefireAna_cfg.py runOnData=True  era=era2017 outFilename=Data13TeV_SinglePhoton_l1prefire.root
```

To submit the ntuplizer to the grid start by setting the environment for crab3.
More details can be found in [CRAB3CheatSheet](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
The following script helps submitting a list of files described in a json file.
Partial submission can be made adding "-o csv_list" as an option.
Adding "-s" will trigger the submission to the grid (otherwise the script only writes down the crab cfg files)

If you receive something like the message below:
```
ERROR: CMSSW_10_6_3 on slc7_amd64_gcc820 is not among supported releases; Use config.JobType.allowUndistributedCMSSW = True if you are sure of what you are doing
```
Do
```
 export SCRAM_ARCH= slc7_amd64_gcc700
```
before cmsenv and compile before submission.


```
python scripts/submitToGrid.py -j test/analysis/ExYukawa/samples_UL2017.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --era era2017 -s
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py -w grid_2016 --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt --era era2016
python scripts/submitToGrid.py -j test/analysis/ExYukawa/samples_2017_ttxy.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only MC -s
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only Data -s
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only 2017H --addParents --rawParents --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt  -s
python scripts/submitToGrid.py -j data/era2018/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only Data -w grid_2018 --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --era era2018
python scripts/submitToGrid.py -j data/era2018/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only MC   -w grid_2018 --era era2018 -s
```

You can get the hyperlink to monitoring by typing
```
crab status
```

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script

```
python scripts/mergeGridOutputs.py -i /store/cmst3/group/top/psilva/ab05162/ -o /store/cmst3/group/top/RunIIReReco/ab05162/
python scripts/mergeGridOutputs.py -i /store/cmst3/group/top/grid_2016/113427a -o /store/cmst3/group/top/RunIIReReco/113427a_2016
```
If condor has a problem, first, do:
```
export _CONDOR_SCHEDD_HOST=bigbird15.cern.ch
export _CONDOR_CREDD_HOST=bigbird15.cern.ch
```

## Luminosity

After ntuples are processed, you can create the list of runs/lumi sections processed using crab as:
```
a=(`find grid/ -maxdepth 1 | grep crab_Data `)
for i in ${a[@]}; do
    crab kill ${i};
    crab status ${i};
    crab report ${i};
done
```
In case of failed jobs the missing lumis can be processed with the following script to wrap the tedious process of
updating the cfg with a finer grain luminosity per job and the missing lumis json
```
for i in ${a[@]}; do
    python scripts/runMissingLumiSecs.py ${i}
done
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run
(see http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html for more details).

```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
brilcalc lumi -b "STABLE BEAMS" -u /pb -i processedLumis.json --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json
```

The following script runs brilcalc inclusively and per trigger path, and stores the results in a ROOT file with the total integrated lumi per run.
It takes a bit to run, depending on the number of triggers configured to use in the analysis

```
python scripts/convertLumiTable.py --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt -y 2018 -o data/era2018
python scripts/convertLumiTable.py --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json
python scripts/convertLumiTable.py -o data/era2016/ -y 2016 --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
```

## Running ntuplizer on condor

You can use this option also, although there won't be a lumi report in the end to run brilcalc on.
Check the options with
```
python scripts/submitLocalNtuplizer.py -h

```

After jobs run, check the integrity by doing
```
python scripts/checkLocalNtuplizerInteg.py condor_file
```


## Preparing the analysis

Correction and uncertainty files are stored under data by era directories (e.g data/era2017) in order no to mix different periods.

* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --out data/era2017/pileupWgts.root
python scripts/runPileupEstimation.py --out data/era2016/pileupWgts.root \
       --json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt \
       --puJson /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /eos/cms/store/group/phys_top/efe/ntuples_a02ce4df_all/MC13TeV_2017_TTTo2L2Nu   -o data/era2017/expectedBtagEff.root;
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/RunIIReReco/ab05162      -o data/era2017/genweights_ab05162.root
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/RunIIReReco/2016/0c522df -o data/era2016/genweights_0c522df.root
```
The lepton/photon trigger/id/iso efficiencies should also be placed under data/era2017.
The src/EfficiencyScaleFactorsWrapper.cc  should then be updated to handle the reading of the ROOT files and the application of the scale factors
event by event.

## Updating the code

Commit your changes regularly with
```
git commit -a -m 'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository click on the green button "Compare,review and create a pull request" to create the PR to merge with your colleagues.

# Local analyses

The ROOT trees created by MiniAnalyzer.cc can be analyzed with a simple executable.
See some examples under `src`.
The new executable should be included in `bin/analysisWrapper.cc` so that it can be used with the runLocalAnalysis.py script
which allows to run over single files or full directories.
See examples under test/ in ```steer*Analysis.sh```.
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis/   -j data/era2017/samples.json  -l 12870
```

main analysis code: ExYukawa.cc (ExYukawa.h) included in ``bin/analysisWrapper.cc`` so that it can be used with ``scripts/runLocalAnalysis.py`` that allows to run over single files or full directories.

Test code:
```
python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py -i /store/cmst3/group/top/RunIIReReco/ab05162/MC13TeV_2017_TTTo2L2Nu_psweights/Chunk_71_ext0.root --tag MC13TeV_2017_dilepton -o ttbar_dilepton_71.root --njobs 1 -q local --debug --era era2017 -m RunExYukawa
```

Submitting the analysis to condor:
```
source test/analysis/ExYukawa/steerAnalysis.sh -o SEL
```
(note the two groups of ntuples that needs to be commented/uncommented exclusively:
1) https://github.com/efeyazgan/TopLJets2015/blob/e745fad072c0ecd1c66b343a87691daaf4dedfe4/TopAnalysis/test/analysis/ExYukawa/steerAnalysis.sh#L39-L40
2) https://github.com/efeyazgan/TopLJets2015/blob/e745fad072c0ecd1c66b343a87691daaf4dedfe4/TopAnalysis/test/analysis/ExYukawa/steerAnalysis.sh#L41-L42
  )

Making plots:
```
python scripts/plotter.py -i /eos/user/e/efe/DataAnalysis/ntuples/ -l 41500 -j test/analysis/ExYukawa/samples_2017.json  -o whatever.root
```
(note that the code uses the cross sections entered in the json files)

To separately the signal from the stack histogram:
```
python scripts/plotter.py -i /eos/user/e/efe/DataAnalysis/ntuples/ -l 41500    -j test/analysis/ExYukawa/samples_2017.json  -o final_plotter.root --signalJson test/analysis/ExYukawa/samples_2017_signal.json
```

Once you make plots you can copy it to a web accessible area and there you should also copy https://github.com/efeyazgan/TopLJets2015/blob/106_protonreco/TopAnalysis/test/index.php

Legend sizes are controlled in:
https://github.com/efeyazgan/TopLJets2015/blob/e745fad072c0ecd1c66b343a87691daaf4dedfe4/TopAnalysis/python/Plot.py#L259
https://github.com/efeyazgan/TopLJets2015/blob/e745fad072c0ecd1c66b343a87691daaf4dedfe4/TopAnalysis/python/Plot.py#L290-L291
