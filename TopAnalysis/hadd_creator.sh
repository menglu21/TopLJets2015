#!/bin/bash
if [ -z "$1" ]
  then
    echo "No argument given. Exiting..."
    return
fi
rm -f hadd.sh
rm -f hadd1.sh
rm -f hadd2.sh
foldername="/eos/user/e/efe/DataAnalysis/ntuples_and_plots/"$1
foldername2="test/analysis/ExYukawa/analysis_2017/Chunks/"
echo $foldername2
if [ ! -d $foldername ]; then
  mkdir $foldername
else
  echo 'Folder exists. Do you want to erase it and a create a new one? (y/n)'
  read -r
  response=$REPLY
fi
if [ $response = "y" ]; then
  rm -rf $foldername
  mkdir $foldername
else
  return
fi
find /eos/cms/store/cmst3/group/top/RunIIUL/2017/6bfa3f2e/. -type d -exec basename {} \; > hadd.sh
find /eos/cms/store/cmst3/group/top/RunIIUL/2017/6bfa3f2e/. -type d -exec basename {} \; > hadd1.sh
tail -n +2 hadd.sh > hadd_tmp.sh
mv hadd_tmp.sh hadd.sh
tail -n +2 hadd1.sh > hadd_tmp.sh
mv hadd_tmp.sh hadd1.sh
sed -i "s|^|hadd $foldername/|g" hadd.sh
sed -i "s|$|.root $foldername2/|" hadd.sh
paste -d'\0' hadd.sh hadd1.sh > hadd2.sh
sed -i "s|$|_*.root|" hadd2.sh
mv hadd2.sh hadd.sh
rm -f hadd1.sh
