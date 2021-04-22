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
f = open('list.txt','w')

#Read json file
samples = dict()
jsonfile = open(os.path.join(workdir+"/../samples_2017.json"))
samplesList = json.load(jsonfile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
for s,desc in samplesList:
  samples[s] = [desc[0]]
jsonfile.close()

#Setting reading target
channel = ["ee","emu","mm","inc"]
name = "_n_after_selection"
histo  = []
for i in range(4):
  histo.append(channel[i]+name)

#Setting Counting list
result = dict()
for t in histo:
  result["MC_"+t] = [0. for i in range(Range)]
  result["Data_"+t]= [0. for i in range(Range)]

Selection = ["Total","Trigger","$N(\ell) \geq 2$","$p_T(\ell1) > 30$ GeV", "$p_T(\ell2) > 20$ GeV", " $p_T(\ell3) < 20$ GeV ", "Opposite sign leptons", "N(jets)$ \geq 3$","N(btags)$ \geq 2$","pass 2b1c tag",""]
Data = [0. for i in range(Range)]
MC = [0. for i in range(Range)]

#Reading root file
for s,desc in samplesList:
  weight = desc[0]*41.5*1000  
  subresult = dict()
  for t in histo:
    subresult[t]=[0. for i in range(Range)]
  
  f.write("---------------"+str(s)+"----------------\n")  

  for ff in Allfile:
    if s in ff:
      th = dict()
      for t in histo:
        th[t] = [0. for i in range(Range)]
      inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
      for t in histo:
        h = inf.Get(t)
        for i in range(Range):
          try:
            subresult[t][i] += h.GetBinContent(i)*weight
            th[t][i] = h.GetBinContent(i)
          except:
            subresult[t][i] += 0.
      if "MC" in ff:
        for t in histo:
          for i in range(Range):
            result["MC_"+t][i] +=th[t][i]*weight
            MC[i]+=th[t][i]*weight
      if "Data" in ff:
        for t in histo:
          for i in range(Range):
            result["Data_"+t][i] +=th[t][i]
            Data[i]+=th[t][i]
      inf.Close()
#  for t in histo:
 #   f.write("-----"+t+"-----\n")  
  #  for i in range(9):
   #   f.write(str(i)+" : "+str(subresult[t][i])+"    \n")
f.write("----------MC/Data-------------\n")
for t in histo:
  f.write("-----"+t+"------\n")
  for i in range(Range):
    try:
      f.write(str(i)+".  Data: "+str(result["Data_"+t][i])+" MC: "+str(result["MC_"+t][i])+"  MC/Data: "+str(round(result["MC_"+t][i]/result["Data_"+t][i],3)*100)+"%\n")
    except:
      f.write(str(i)+".  Data: "+str(result["Data_"+t][i])+" MC: "+str(result["MC_"+t][i])+"  MC/Data: 0.0 %\n")
f.write("----------Overleaf----------\n")
for i in range(Range):
  f.write("\hline\n")
  f.write(Selection[i])
  for t in histo:
    try:
      f.write(" & "+str(round(result["Data_"+t][i]/result["MC_"+t][i],3)*100)+"\%")
    except:
      f.write(" & 0.0 \%")
  f.write("\n")
 # try: 
   # f.write(" & "+str(round(MC[i]/Data[i],3)*100)+"\%\n")
 # except:
  #  f.write(" & 0.0\%\n")
f.write("---------Data Number------------\n")
for i in range(Range):
  f.write("\hline\n")
  f.write(Selection[i])
  for t in histo:
    #print(t)
    f.write(" & "+str(result["Data_"+t][i]))
  #f.write(" & "+ str(Data[i]) +"\n")
  f.write("\n")
f.write("--------MC Number----------------\n")
for i in range(Range):
  f.write("\hline\n")
  f.write(Selection[i])
  for t in histo:
    #print(t)
    f.write(" & "+str(round(result["MC_"+t][i],1)))
  f.write("\n")
  #f.write(" & "+ str(MC[i]) +"\n")
f.write("-------Efficiency---------------\n")
for i in range(Range):
  f.write("\hline\n")
  f.write(Selection[i])
  for t in histo:
    #print(t)
    f.write(" & "+str(round(result["Data_"+t][i]/result["Data_"+t][2],6)*100)+"\%")
  #f.write(" & "+ str(Data[i]) +"\n")
  f.write("\n")


f.close()
