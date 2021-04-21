#!/bin/bash

#parse input arguments
ERA=2017
while getopts "o:y:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done


#check an operation has been given
if [ -z "$WHAT" ]; then
    echo "steerAnalysis_ChargeFlip.sh -o <TEST/SEL/MERGE/PLOT> [ -y 2016/2017/2017lowpu/2018 ] ";
    echo "   TEST          - test locally the code on a single file";
    echo "   SEL           - launches selection jobs to the batch, output will contain summary trees and control plots";
    echo "   MERGE         - merge output"
    echo "   PLOT          - make plots"
    exit 1;
fi

#configuration parameters
queue=tomorrow
samples=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/ExYukawa/ChargeFlipStudy/samples_${ERA}_chargeflip.json
samples_signal=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/ExYukawa/samples_${ERA}_signal.json
outdir=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/ExYukawa/ChargeFlipStudy/analysis_${ERA}/
#outdir=/eos/cms/store/group/phys_top/efe/ExYukawa/analysis_${ERA}
if [[ ${ERA} = "2016" ]]; then
    githash=0c522df
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
    dataeosdir=${eosdir}
    lumi=35882
    lumiUnc=0.025
    testtag=Data13TeV_2016B_SingleMuon
    testfile=${eosdir}/${testtag}/Chunk_0_ext0.root
elif [[ ${ERA} = "2017" ]]; then
###    githash=848840ab
###    githash=ab05162
#    githash=6bfa3f2e
#    eosdir=/store/cmst3/group/top/RunIIUL/2017/${githash}
#    githash=a02ce4df
#    eosdir=/store/group/phys_top/efe/ntuples_${githash}_all
    githash=39fa2880
    eosdir=/store/group/phys_top/efe/ntuples_${githash}
#    githash=ae6e08e
#    eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
#    dataeosdir=${eosdir}
    lumi=41367
    lumiUnc=0.025
    testtag=MC13TeV_2017_TAToTTQ_MA200_G2HDM_rtc04
    #testtag=MC13TeV_2017_TTWH
    #testfile=/store/cmst3/group/top/RunIIReReco/2017/ae6e08e/MC13TeV_2017_TTWH/Chunk_0_ext0.root
    #testfile=/store/cmst3/group/top/RunIIUL/2017/848840ab/MC13TeV_2017_TAToTTQ_MA200_G2HDM_rtc04/Chunk_0_ext0.root
    #testfile=/store/cmst3/group/hintt/psilva/6bfa3f2e/TAToTTQ_MA-200_TuneCP5_13TeV_G2HDM-rtc04-madgraphMLM-pythia8/crab_MC13TeV_2017_TAToTTQ_MA200_G2HDM_rtc04/200915_073154/0000/MiniEvents_6.root
    testfile=${eosdir}/MC13TeV_2017_DY50toInf_fxfx/Chunk_0_ext0.root
#    testfile=${eosdir}/Data13TeV_2017B_DoubleEG/Chunk_0_ext0.root
    #testfile=/store/cmst3/group/top/RunIIReReco/2017/ae6e08e/MC13TeV_2017_TAToTTQ_MA200_G2HDM_rtc04/Chunk_0_ext0.root
    #testfile=/store/cmst3/group/top/RunIIReReco/2017/ae6e08e/MC13TeV_2017_TTZZ/Chunk_0_ext0.root
#  elif [[ ${ERA} = "2017_add_bkgs" ]]; then
#      githash=ae6e08e
#      eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
#      dataeosdir=${eosdir}
#      lumi=41367
#      lumiUnc=0.025
#      testtag=MC13TeV_2017_TTHH
#      testfile=/eos/cms/store/cmst3/group/top/RunIIReReco/2017/ae6e08e/MC13TeV_2017_TTHH/Chunk_0_ext0.root
elif [[ ${ERA} = "2017lowpu" ]]; then
    ERA=2017
    githash=ab05162
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    dataeosdir=/store/cmst3/group/top/RunIIReReco/2017/newproton_calib/
    lumi=220
    lumiUnc=0.025
    testtag=Data13TeV_2017H_SingleMuon_v2
    testfile=${dataeosdir}/${testtag}/Chunk_0_ext0.root
fi

#run the operation required
RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TEST )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${testfile} --tag ${testtag} \
            -o testsel_${ERA}.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era${ERA} -m RunChargeFlip;
        ;;

    SEL )
        baseOpt="--genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era${ERA} -m RunChargeFlip"
        baseOpt="${baseOpt} --only ${samples}";
  #echo ${samples};
  #return;
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} -i ${eosdir};
#        if [ "${eosdir}" != "{dataeosdir}" ]; then
#	    python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} -i ${dataeosdir} --farmappendix data;
#
#        fi
	;;

    CHECKSELINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARM${githash}/
        ;;

    MERGE )
	mergeOutputs.py ${outdir} True;
	;;

    PLOT )
  plotinputdir="/eos/user/e/efe/DataAnalysis/ntuples_and_plots/1Mars/"
  plotoutdir="/eos/user/e/efe/www/exyukawa/1mars"
  mkdir -p ${plotoutdir}
  wget https://raw.githubusercontent.com/efeyazgan/TopLJets2015/106_protonreco/TopAnalysis/test/index.php -P ${plotoutdir}
  commonOpts="-i ${outdir} -l ${lumi} --mcUnc ${lumiUnc}"
  #plotOpts="-i ${plotinputdir} -l ${lumi} -j ${samples}  --signalJson ${samples_signal} -O ${plotoutdir}"
  plotOpts="-i ${plotinputdir} -l ${lumi} -j ${samples} -O ${plotoutdir}"
  echo ${plotOpts}
  python scripts/plotter.py ${plotOpts}
	;;

    PLOTANAPERERA )
  #plotinputdir="/eos/user/e/efe/DataAnalysis/ntuples_and_plots/2oct/"
  plotinputdir="/eos/user/e/efe/DataAnalysis/ntuples_and_plots/Production_5Oct/"

  plots=""
  for evcat in ee mm emu inc; do
      plist=(control_Z_mass)
      for p in ${plist[@]}; do
          plots="${evcat}_${p},${plots}"
      done
  done
  echo ${plots}

  for era in B C D E; do
      eralumi=${lumi}
      if [ "${era}" = "B" ]; then
          eralumi=`echo ${eralumi}*0.115 | bc`
      elif [ "${era}" = "C" ]; then
          eralumi=`echo ${eralumi}*0.233 | bc`
      elif [ "${era}" = "D" ]; then
          eralumi=`echo ${eralumi}*0.103 | bc`
      elif [ "${era}" = "E" ]; then
          eralumi=`echo ${eralumi}*0.22 | bc`
#      elif [ "${era}" = "F" ]; then
#          eralumi=`echo ${eralumi}*0.329 | bc`
      fi
      plotoutdir="/eos/user/e/efe/www/exyukawa/6oct_test2"
      mkdir -p ${plotoutdir}
      wget https://raw.githubusercontent.com/efeyazgan/TopLJets2015/106_protonreco/TopAnalysis/test/index.php -P ${plotoutdir}
      plotoutdir=${plotoutdir}/${era}
      mkdir -p ${plotoutdir}
      echo ${plotoutdir}
      wget https://raw.githubusercontent.com/efeyazgan/TopLJets2015/106_protonreco/TopAnalysis/test/index.php -P ${plotoutdir}
      era_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/ExYukawa/samples_2017_${era}.json;
      echo ${era_json}
      plotOpts="-i ${plotinputdir} -l ${eralumi} -j ${era_json}  --signalJson ${samples_signal} -O ${plotoutdir}"
      python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${plotOpts} --only ${plots} --strictOnly;
  done
  ;;

esac
