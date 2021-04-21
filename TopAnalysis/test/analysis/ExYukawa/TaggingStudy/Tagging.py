import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np

#Basic Setting
Range = 10

#Read files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)
f = open('tagging.txt','w')

#Read json file
samples = dict()
jsonfile = open(os.path.join(workdir+"/../samples_2017.json"))
samplesList = json.load(jsonfile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
for s,desc in samplesList:
  samples[s] = [desc[0]]
jsonfile.close()

#Setting reading target
channel = ["ee","emu","mm","inc"]
name = ["_compare_btag", "_compare_ctag", "_compare_2b1ctag"]
histo = []
for i in range(4):
  for j in range(3):
    histo.append(channel[i]+name[j])

#Setting Counting list 
result = dict()
for t in histo:
  result[t] = [0. for i in range(Range)]
Title = ["1","2","3.MC","4.Tagging","5.Both","6","7","8","9","10"]

#Reading root file
for s,desc in samplesList:
  weight = desc[0]*41.5*1000
  subresult = dict()
  for t in histo:
    subresult[t]=[0. for i in range(Range)]
  f.write("------------"+str(s)+"--------------\n")

  for ff in Allfile:
    if ((s in ff) and ("MC" in ff)) :
      inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
      for t in histo:
        h = inf.Get(t)
        for i in range(Range):
          try:
            subresult[t][i] += h.GetBinContent(i)*weight
            result[t][i] += h.GetBinContent(i)*weight
          except:
            subresult[t][i] += 0
      inf.Close()
  for t in histo:
    f.write("* "+t)
    for i in range(Range):
      f.write(Title[i]+": "+ str(subresult[t][i])+" ")
    try:
      f.write("Efficiency: "+str(round(subresult[t][4]/subresult[t][2],3)*100.)+"%\n") 
    except:
      f.write("Efficiency: N/A \n")
f.write("------------Total------------\n")
for t in histo:
  f.write("* "+t)
  for i in range(Range):
    f.write(Title[i]+": "+ str(result[t][i])+" ")
  try:
    f.write("Efficiency: "+str(round(result[t][4]/result[t][2],3)*100.)+"%\n")   
  except:
    f.write("Efficiency: N\A \n")

  


