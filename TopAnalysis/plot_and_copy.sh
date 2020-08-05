mkdir -p $2
python scripts/plotter.py -i $1 -l 41500    -j test/analysis/ExYukawa/samples_2017.json  -o final_plotter.root --signalJson test/analysis/ExYukawa/samples_2017_signal.json -O $2
cp test/index.php $2
