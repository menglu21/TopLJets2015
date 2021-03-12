import argparse
import ROOT
import json

#Setting Parse
parser = argparse.ArgumentParser(description="Process the root file")
parser.add_argument('fin', type=str)
args = parser.parse_args()

#Read File
inf = ROOT.TFile.Open(args.fin)
t = inf.Get("analysis/data")
print(t)

#Main code
while(1):
  entry = input('Entry :')
  t.GetEntry(int(entry))
  while(1):
    try:
      n = input('Object : ')
      try:
        for i in range (len(n)):
          print(str(i)+" : " + str(n[i]))
      except:
        print(n)
    except:
      break
        
