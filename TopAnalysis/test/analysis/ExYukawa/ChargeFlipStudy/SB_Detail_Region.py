import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import math
import argparse


##########
# Parser #
##########

parser = argparse.ArgumentParser(description = "Charge Flip Rate")
parser.add_argument('plotdir',type=str)
parser.add_argument('Covariance',type=int)
args = parser.parse_args()

#################
# Basic Setting #
#################

plotdir = args.plotdir # Change to your own directory
Covariance =  (args.Covariance==1)   # Consider the covariance of the uncertainty of Pij or not.

# Root file directory
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/analysis_2017/Chunks/")
Allfile = os.listdir(targetdir)

# Json Files
jsonfile = open(os.path.join(workdir+"/samples_2017_chargeflip.json"))
samplesList = json.load(jsonfile,encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Basic Set Up
pt_region = np.array(( 20.0, 50.0, 100.0, 200., 300.)) # It will actually include overflow events ( Already in the root file )
eta_region = np.array((0.0, 0.8, 1.479, 2.4))
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
channel = ['MC_oc','MC_ss','data_oc','data_ss']
Number = dict()
Pij = dict()
Vij = dict()
for ch in channel:
  Number[ch] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]

#############
# Read Data #
#############

for s ,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in Allfile:
    if s in ff:
      print(ff)
      try:
        inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
        t = inf.Get("TreeInput")
        n_entry = t.GetEntries()
        t_Nss = t.t_Nss
        t_Noc = t.t_Noc
        if "MC" in s:
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  Number['MC_ss'][i][j][ii][jj] += t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
                  Number['MC_oc'][i][j][ii][jj] += t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]*weight
        elif "Data" in s:
          t.GetEntry(0)
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  if not (type(t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])==type(22.2)): # Use to Debug
                    print(ff)
                    print(type(t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]))
                    print(t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                    print(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)
                  if not (math.isnan(t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])): # I have got nan in data root files. But it disappeared later. Still don't know the reason.
                    Number['data_ss'][i][j][ii][jj] += float(t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  else:
                    print(ff+"ss nan happen")
                  if not (math.isnan(t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])):
                    Number['data_oc'][i][j][ii][jj] += float(t_Noc[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i])
                  else:
                    print(ff+"oc nan happen")
#                    print("******* "+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+str(t.t_Nss[jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i]))
        inf.Close()
      except:
        print(ff+"-->Exception")

#####################
# Calculate CF Rate #
#####################

# Record Pij and Vij ( The Uncertainty of Pij )
DvMC = ['data','MC']
for h in DvMC:
  Pij[h] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
  Vij[h] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          print("l1_pt = " + str(pt_region[i]) + " l1_eta = " + str(eta_region[j]) + " l2_pt = " + str(pt_region[ii]) + " l2_eta = " + str(eta_region[jj]))
          print('ss : '+str(Number[h+'_ss'][i][j][ii][jj])+' oc : '+str(Number[h+'_oc'][i][j][ii][jj]))
          try:
            N_ss = Number[h+'_ss'][i][j][ii][jj]
            N_T = N_ss + Number[h+'_oc'][i][j][ii][jj]
            pij = N_ss/N_T
            if Covariance: 
              Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3)-2*(N_ss**1.5)/(N_T**2.5)
            else:  
              Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3) # without consider covariance 
            if(pij<0.0): 
              Pij[h][i][j][ii][jj] = 0.0 # Avoid the background > signal region.
              Vij[h][i][j][ii][jj] = 0.0
            else:
              Pij[h][i][j][ii][jj] = pij
          except:
            pass


# Use chi2 to fit ( The code is dirty. Need to be more elegant. )
for h in DvMC:
  w = ROOT.RooWorkspace(h+'w')
  chi2_string = ''
  para_string = ''
  sum_string = ''
  for i in range(pt_bins):
    for j in range(eta_bins):
      if Number[h+'_oc'][i][j][i][j]>10000: # High statistics region is used to give the initial value of Pij in fitting.
        para_string += ',P'+str(i)+str(j)+'['+str(1.-(1.-Pij[h][i][j][i][j])**0.5)+',0,0.1]' 
      else:
        para_string += ',P'+str(i)+str(j)+'[0,0.1]'
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          if not Vij[h][i][j][ii][jj]==0 and Number[h+'_oc'][i][j][ii][jj]>100 and Number[h+'_ss'][i][j][ii][jj]>1: # Avoid low statistics region.
            P1 = 'P'+str(i)+str(j)
            P2 = 'P'+str(ii)+str(jj)
            chi2_string = '('+str(Pij[h][i][j][ii][jj])+'-('+P1+'+'+P2+'-'+P1+'*'+P2+'))**2'+'/('+str(Vij[h][i][j][ii][jj])+')'
            w.factory('expr::chi'+str(i)+str(j)+str(ii)+str(jj)+'(\''+chi2_string+'\''+para_string+')')
          else:
            w.factory('expr::chi'+str(i)+str(j)+str(ii)+str(jj)+'(\'0.0'+'\''+para_string+')') # Low statistics region doesn't join the fitting.
          if i==0 and j ==0 and ii==0 and jj==0:
            pass
          elif i==0 and j==0 and ii==0 and jj==1:
            w.factory('sum::chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+'(chi0000,chi0001)')
          else:
            w.factory('sum:chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+'(chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i-1)+',chi'+str(i)+str(j)+str(ii)+str(jj)+')')
  w.Print()
# Minimize the Chi2 Model
  mini = ROOT.RooMinimizer(w.function('chi2_'+str(eta_bins*pt_bins*eta_bins*pt_bins-1)))
  mini.minimize('Minuit2')
  r = mini.save()

###############
# Result Plot #
###############

  c = ROOT.TCanvas(h+'c','c',600,600)
  ROOT.gStyle.SetOptStat("kFALSE")
  ROOT.gStyle.SetPaintTextFormat(".2e");
  ROOT.gStyle.SetPalette(69);
  r.correlationHist(h).Draw('colz')
  c.Update()
  if Covariance:
    c.SaveAs(plotdir+h+'_CorrelationHist_SB_detail_Cov.png') # Can change the file name.
  else:
    c.SaveAs(plotdir+h+'_CorrelationHist_SB_detail_woCov.png')
  
  h_chargeflip = ROOT.TH2D(h+"CFRate",";P_{T}[GeV] ; \eta;"+h+"_ChargeFlip rate",pt_bins,pt_region,eta_bins,eta_region)
  for i in range(pt_bins):
    for j in range(eta_bins):
      h_chargeflip.SetBinContent(i+1,j+1,w.var('P'+str(i)+str(j)).getValV())
      h_chargeflip.SetBinError(i+1,j+1,w.var('P'+str(i)+str(j)).getError())
  if h=='data':
    h_data = h_chargeflip.Clone()
  else:
    h_MC = h_chargeflip.Clone()
  c.SetLogx()
  h_chargeflip.Draw('COLZTEXT e')
  c.Update()
  if Covariance:
    c.SaveAs(plotdir+h+'_CFRate_SB_detail_Cov.png') # Can change the file name
  else:
    c.SaveAs(plotdir+h+'_CFRate_SB_detail_woCov.png')

################
# Scale Factor #
################

h_data.Divide(h_data,h_MC)
h_data.Draw('COLZTEXT e')
c.Update()
if Covariance:
  c.SaveAs(plotdir+'CFRate_SF_SB_detail_Cov.png') # Can chage the file name
else:
  c.SaveAs(plotdir+'CFRate_SF_SB_detail_woCov.png')
