import os
import ROOT
import json
import sys
from collections import OrderedDict
# Read all files
workdir =  os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)
f = open('out.txt','w')
#Read json file
samples = dict()
jsonfile = open(os.path.join(workdir+"/../samples_2017.json"),'r')
samplesList = json.load(jsonfile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
for s,desc in samplesList:
  samples[s] = [desc[0]]
jsonfile.close()

Data = [0,0,0,0,0,0,0,0,0]
MC = [0,0,0,0,0,0,0,0,0]
Control_data = [0,0,0,0]
Control_MC = [0,0,0,0]
Selection = ["","","Total","Trigger","$N(\ell) \geq 2$","$p_T(\ell1) > 30$ GeV", "$p_T(\ell2) > 20$ GeV", " $p_T(\ell3) < 20$ GeV ", "Opposite sign leptons", "N(jets)$ \geq 3$"]
for s, desc in samplesList:  

  Data_Ntotal = 0
  Data_N_after_trig = 0
  Data_N_after_lep = 0
  Data_N_after_lep1 = 0
  Data_N_after_lep2 = 0
  Data_N_after_lep3 = 0
  Data_N_after_ss = 0
  Data_N_after_jet = 0
  Data_N_after_all = 0

  control_n_lep = 0
  control_n_lep1 = 0
  control_n_lep2 = 0
  control_n_lep3 = 0

  for tt in Allfile :
    if s in tt :
      inf = ROOT.TFile.Open(os.path.join(targetdir+tt))
      t = inf.Get("TreeInput")
      t.GetEntry(0)
      Data_Ntotal = Data_Ntotal + t.Ntotal
      Data_N_after_trig = Data_N_after_trig + t.N_after_trig
      Data_N_after_lep = Data_N_after_lep + t.N_after_lep_selec
      Data_N_after_lep1 += t.N_after_lep1_selec
      Data_N_after_lep2 += t.N_after_lep2_selec
      Data_N_after_lep3 += t.N_after_lep3_selec
      Data_N_after_ss += t.N_after_ss
      Data_N_after_jet += t.N_after_jet_selec
      Data_N_after_all += t.N_after_all
      control_n_lep += t.control_n_lep
      control_n_lep1 += t.control_n_lep1
      control_n_lep2 += t.control_n_lep2
      control_n_lep3 += t.control_n_lep3
      inf.Close()
  f.write("---------------------------"+str(s)+"---------------------------\n")
  if "MC" in s:  
    f.write("Cross section = "+ str(desc[0])+"\n")
    f.write("Expected number = " + str(desc[0]*1000*41.5)+"\n")
    weight = desc[0]*1000*41.5
    f.write("weight = "+str(weight)+"\n")
  #f.write("0: Ntotal = " + str(Data_Ntotal)+"\n")
  if Data_Ntotal>0:
   # f.write("1: N_after_trig = " + str(Data_N_after_trig) +"  "+ str(round(float(Data_N_after_trig)/float(Data_Ntotal)*100.0,3))+"%"+"\n")
    f.write("2: N_after_lep = " + str(Data_N_after_lep) +"  "+ str(round(float(Data_N_after_lep)/float(Data_Ntotal)*100.0,3))+"%"+"\n")
    f.write("3: N_after_lep1 = " + str(Data_N_after_lep1) +"  "+ str(round(float(Data_N_after_lep1)/float(Data_Ntotal)*100.0,3))+"%\n")
    f.write("4: N_after_lep2 = " + str(Data_N_after_lep2) +"  "+ str(round(float(Data_N_after_lep2)/float(Data_Ntotal)*100.0,3))+"%\n")
    f.write("5: N_after_lep3 = " + str(Data_N_after_lep3) +"  "+ str(round(float(Data_N_after_lep3)/float(Data_Ntotal)*100.0,3))+"%\n")
    f.write("6: N_after_samesign = " + str(Data_N_after_ss) +"  "+ str(round(float(Data_N_after_ss)/float(Data_Ntotal)*100.0,3))+"%\n")
    f.write("7: N_after_jet = " + str(Data_N_after_jet) +"  "+ str(round(float(Data_N_after_jet)/float(Data_Ntotal)*100.0,3))+"%\n")
    f.write("8: N_after_all = " + str(Data_N_after_all) +"  "+ str(round(float(Data_N_after_all)/float(Data_Ntotal)*100.0,3))+"%\n")
  f.write("------------Control----------\n")
  if(control_n_lep>0):
    f.write("2: control_n_after_lep = " + str(control_n_lep) + " "+str(round(control_n_lep/control_n_lep*100.0,3))+"%\n")
    f.write("3: control_n_after_lep1 = " + str(control_n_lep1) + " "+str(round(control_n_lep1/control_n_lep*100.0,3))+"%\n")
    f.write("4: control_n_after_lep2 = " + str(control_n_lep2) + " "+str(round(control_n_lep2/control_n_lep*100.0,3))+"%\n")
    f.write("5: control_n_after_lep3 = " + str(control_n_lep3) + " "+str(round(control_n_lep3/control_n_lep*100.0,3))+"%\n")

  if "Data" in s:
    Data[0]+=Data_Ntotal
    Data[1]+=Data_N_after_trig
    Data[2]+=Data_N_after_lep
    Data[3]+=Data_N_after_lep1
    Data[4]+=Data_N_after_lep2
    Data[5]+=Data_N_after_lep3
    Data[6]+=Data_N_after_ss
    Data[7]+=Data_N_after_jet
    Data[8]+=Data_N_after_all
    Control_data[0]+= control_n_lep
    Control_data[1]+= control_n_lep1
    Control_data[2]+= control_n_lep2
    Control_data[3]+= control_n_lep3
  if"MC" in s:
    MC[0]+=Data_Ntotal*weight
    MC[1]+=Data_N_after_trig*weight
    MC[2]+=Data_N_after_lep*weight
    MC[3]+=Data_N_after_lep1*weight
    MC[4]+=Data_N_after_lep2*weight
    MC[5]+=Data_N_after_lep3*weight
    MC[6]+=Data_N_after_ss*weight
    MC[7]+=Data_N_after_jet*weight
    MC[8]+=Data_N_after_all*weight
    Control_MC[0]+= control_n_lep*weight
    Control_MC[1]+= control_n_lep1*weight
    Control_MC[2]+= control_n_lep2*weight
    Control_MC[3]+= control_n_lep3*weight


f.write("--------Data/MC-------\n")
for i in range (2,9):
  f.write(str(i)+":  Data: "+str(Data[i])+"    MC: "+str(MC[i])+"  MC/Data =  "+str(float(MC[i])/float(Data[i])*100.)+"%\n")  
f.write("---------Data/MC(Ratio compared to Number after lep selec)------\n")
for i in range (2,9):
  f.write(str(i)+":  Data: "+str(float(Data[i])/float(Data[2])*100.)+"%     MC: "+str(float(MC[i])/float(MC[2])*100.)+"%  MC/Data =  "+str((float(MC[i])/float(MC[2]))/(float(Data[i])/float(Data[2]))*100.)+"%\n")

f.write("---------Control--------\n")
for i in range (4):
  f.write(str(i)+":  Data: "+str(Control_data[i])+"    MC: "+str(Control_MC[i])+"  MC/Data =  "+str(float(Control_MC[i])/float(Control_data[i])*100.)+"%\n")

f.close()


