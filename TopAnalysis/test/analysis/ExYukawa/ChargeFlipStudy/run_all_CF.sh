plotdir=/eos/user/t/tihsu/DataAnalysis/ChargeFlipStudy/
python Zmassdistribution.py ${plotdir}
python SB_seperateModel.py ${plotdir} 1
python SB_seperateModel.py ${plotdir} 0
python SB_Detail_Region.py ${plotdir} 1
python SB_Detail_Region.py ${plotdir} 0
python SB_Combine_Region.py ${plotdir} 1
python SB_Combine_Region.py ${plotdir} 0
