import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import math
# Read Files
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/Chunks/")
Allfile = os.listdir(targetdir)

# Read json Files
jsonfile = open(os.path.join(workdir+"/../samples_2017_chargeflip.json"))
samplesList = json.load(jsonfile,encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Basic Setting
pt_region = np.array(( 20.0, 50.0, 100.0, 200., 300.))
eta_region = np.array((0.0, 0.8, 1.5, 2.5))
pt_bins = len(pt_region)-1
eta_bins = len(eta_region)-1
channel = ['MC_oc','MC_ss','data_oc','data_ss']
Number = dict()
Pij = dict()
Vij = dict()
for ch in channel:
  Number[ch] = [[[[0. for i in range(eta_bins)] for j in range(pt_bins)] for ii in range(eta_bins)] for jj in range(pt_bins)]
Mean = [90.620,89.991,90.218,89.113]
Gamma = [4.267,5.230,4.691,6.350]
Zmass = dict()
Zwidth = dict()
for i in range(len(channel)):
  Zmass[channel[i]] = Mean[i]
  Zwidth[channel[i]] = Gamma[i]
# S+B model
w = ROOT.RooWorkspace()
w.factory('BreitWigner::BW(x[50,130],mean[70,110],gamma[0,40])')
w.factory('Exponential::e(x,alpha[-10,0])')
w.factory('CBShape::CB(x,m[0,150],s[0,150],a[0,150],n[-5,5])')
w.factory('FCONV::signal(x,BW,CB)')
w.factory('SUM::pdf(sig[0,100000]*signal,bg[100,0,1000]*e)')


pdf = w.pdf('pdf')
signal_pdf = w.pdf('BW')
x = w.var('x')
gamma = w.var('gamma')
mean = w.var('mean')
sig = w.var('sig')
bg = w.var('bg')
xframe = x.frame()

# Plot Setting
plotdir = "/eos/user/t/tihsu/DataAnalysis/Charge_Flip/"
fit_plotdir = plotdir+"Fit_Result/"

Result_hist = dict()
NTotal = dict()
for i in range(pt_bins):
  for j in range(eta_bins):
    for ii in range(pt_bins):
      for jj in range(eta_bins):
        MC_hss_name = "Fit_MC_ss_"+str(i)+str(j)+str(ii)+str(jj)
        MC_hoc_name = "Fit_MC_oc_"+str(i)+str(j)+str(ii)+str(jj)
        data_hss_name = "Fit_data_ss_"+str(i)+str(j)+str(ii)+str(jj)
        data_hoc_name = "Fit_data_oc_"+str(i)+str(j)+str(ii)+str(jj)
        Result_hist[MC_hss_name] = ROOT.TH1F(MC_hss_name,"; M(Z) [GeV] ; Events",80,50,130)
        Result_hist[MC_hoc_name] = ROOT.TH1F(MC_hoc_name,"; M(Z) [GeV] ; Events",80,50,130)
        Result_hist[data_hss_name] = ROOT.TH1F(data_hss_name,"; M(Z) [GeV] ; Events",80,50,130)
        Result_hist[data_hoc_name] = ROOT.TH1F(data_hoc_name,"; M(Z) [GeV] ; Events",80,50,130)

# Run Through All files
for s ,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in Allfile:
    if s in ff:
      print(ff)
      try:
        inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
        if "MC" in s:
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  hss_name = "Fit_MC_ss_"+str(i)+str(j)+str(ii)+str(jj)
                  hoc_name = "Fit_MC_oc_"+str(i)+str(j)+str(ii)+str(jj)
                  hh_ss = inf.Get("ee_CF"+str(i)+str(j)+str(ii)+str(jj)+"_ss_Zmass")
                  hh_oc = inf.Get("ee_CF"+str(i)+str(j)+str(ii)+str(jj)+"_oc_Zmass")
                  if(type(hh_ss)==type(Result_hist[hss_name])):
                    Result_hist[hss_name].Add(Result_hist[hss_name],hh_ss,1.0,weight)
                  if(type(hh_oc)==type(Result_hist[hss_name])):
                    Result_hist[hoc_name].Add(Result_hist[hoc_name],hh_oc,1.0,weight)
        if "Data" in s and "Muon" not in s:
          for i in range(pt_bins):
            for j in range(eta_bins):
              for ii in range(pt_bins):
                for jj in range(eta_bins):
                  hss_name = "Fit_data_ss_"+str(i)+str(j)+str(ii)+str(jj)
                  hoc_name = "Fit_data_oc_"+str(i)+str(j)+str(ii)+str(jj)
                  hh_ss = inf.Get("ee_CF"+str(i)+str(j)+str(ii)+str(jj)+"_ss_Zmass")
                  hh_oc = inf.Get("ee_CF"+str(i)+str(j)+str(ii)+str(jj)+"_oc_Zmass")
                  if(type(hh_ss)==type(Result_hist[hss_name])):
                    Result_hist[hss_name].Add(Result_hist[hss_name],hh_ss,1.0,1.0)
                  if(type(hh_oc)==type(Result_hist[hss_name])):
                    Result_hist[hoc_name].Add(Result_hist[hoc_name],hh_oc,1.0,1.0)
        inf.Close()
      except:
        print(ff+"-->Trigger exception")
c = ROOT.TCanvas('c1','c1',600,600)
ROOT.gStyle.SetOptStat("kFALSE")
ROOT.gStyle.SetPaintTextFormat(".2e");
ROOT.gStyle.SetPalette(57);

for i in range(pt_bins):
  for j in range(eta_bins):
    for ii in range(pt_bins):
      for jj in range(eta_bins):
        MC_hss_name = "Fit_MC_ss_"+str(i)+str(j)+str(ii)+str(jj)
        MC_hoc_name = "Fit_MC_oc_"+str(i)+str(j)+str(ii)+str(jj)
        data_hss_name = "Fit_data_ss_"+str(i)+str(j)+str(ii)+str(jj)
        data_hoc_name = "Fit_data_oc_"+str(i)+str(j)+str(ii)+str(jj)
        dh_ss = ROOT.RooDataHist(MC_hss_name,MC_hss_name,ROOT.RooArgList(x),Result_hist[MC_hss_name])
        dh_oc = ROOT.RooDataHist(MC_hoc_name,MC_hoc_name,ROOT.RooArgList(x),Result_hist[MC_hoc_name])

        xframe = x.frame()
        pdf.fitTo(dh_ss)
        dh_ss.plotOn(xframe)
        pdf.plotOn(xframe)
        pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
        xframe.Draw()
        c.SaveAs(fit_plotdir+MC_hss_name+'.png')
        Result_hist[MC_hss_name].Draw()
        c.Update()
        c.SaveAs(fit_plotdir+MC_hss_name+'_original.png')
        Number['MC_ss'][i][j][ii][jj] += sig.getValV()
        print("sig"+str(sig.getValV()))

        xframe = x.frame()
        pdf.fitTo(dh_oc)
        dh_oc.plotOn(xframe)
        pdf.plotOn(xframe)
        pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
        xframe.Draw()
        c.SaveAs(fit_plotdir+MC_hoc_name+'.png')
        Number['MC_oc'][i][j][ii][jj] += sig.getValV()
     
        dh_ss = ROOT.RooDataHist(data_hss_name,data_hss_name,ROOT.RooArgList(x),Result_hist[data_hss_name])
        dh_oc = ROOT.RooDataHist(data_hoc_name,data_hoc_name,ROOT.RooArgList(x),Result_hist[data_hoc_name])
        
        xframe = x.frame()
        pdf.fitTo(dh_ss)
        dh_ss.plotOn(xframe)
        pdf.plotOn(xframe)
        pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
        xframe.Draw()
        c.SaveAs(fit_plotdir+data_hss_name+'.png')
        Number['data_ss'][i][j][ii][jj] += sig.getValV()
        
        xframe = x.frame()
        pdf.fitTo(dh_oc)
        dh_oc.plotOn(xframe)
        pdf.plotOn(xframe)
        pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
        xframe.Draw() 
        c.SaveAs(fit_plotdir+data_hoc_name+'.png')
        Number['data_oc'][i][j][ii][jj] += sig.getValV()


       
# Calculate Probability
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
            Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3)-2*(N_ss**1.5)/(N_T**2.5)
            #Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3) # without consider covariance 
            if(pij<0.0):
              Pij[h][i][j][ii][jj] = 0.0
            else:
              Pij[h][i][j][ii][jj] = pij
          except:
            pass
print(Pij['data'])
print(Pij['MC'])

# Use chi2 to fit
for h in DvMC:
  w = ROOT.RooWorkspace(h+'w')
  chi2_string = ''
  para_string = ''
  sum_string = ''
  for i in range(pt_bins):
    for j in range(eta_bins):
      if Number[h+'_oc'][i][j][i][j]>10000:
        para_string += ',P'+str(i)+str(j)+'['+str(1.-(1.-Pij[h][i][j][i][j])**0.5)+',0,0.1]'
      else:
        para_string += ',P'+str(i)+str(j)+'[0,0.1]'
  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          if not Vij[h][i][j][ii][jj]==0 and Number[h+'_oc'][i][j][ii][jj]>100 and Number[h+'_ss'][i][j][ii][jj]>1:
            P1 = 'P'+str(i)+str(j)
            P2 = 'P'+str(ii)+str(jj)
            chi2_string = '('+str(Pij[h][i][j][ii][jj])+'-('+P1+'+'+P2+'-'+P1+'*'+P2+'))**2'+'/('+str(Vij[h][i][j][ii][jj])+')'
            w.factory('expr::chi'+str(i)+str(j)+str(ii)+str(jj)+'(\''+chi2_string+'\''+para_string+')')
          else:
            w.factory('expr::chi'+str(i)+str(j)+str(ii)+str(jj)+'(\'0.0'+'\''+para_string+')')
          if i==0 and j ==0 and ii==0 and jj==0:
            pass
          elif i==0 and j==0 and ii==0 and jj==1:
            w.factory('sum::chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+'(chi0000,chi0001)')
          else:
            w.factory('sum:chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i)+'(chi2_'+str(jj+eta_bins*ii+eta_bins*pt_bins*j+eta_bins*pt_bins*eta_bins*i-1)+',chi'+str(i)+str(j)+str(ii)+str(jj)+')')
  w.Print()
  mini = ROOT.RooMinimizer(w.function('chi2_'+str(eta_bins*pt_bins*eta_bins*pt_bins-1)))
  mini.minimize('Minuit2')
  r = mini.save()
# Result Plot
  c = ROOT.TCanvas(h+'c','c',600,600)
  ROOT.gStyle.SetOptStat("kFALSE")
  ROOT.gStyle.SetPaintTextFormat(".2e");
  ROOT.gStyle.SetPalette(57);
  r.correlationHist(h).Draw('colz')
  c.Update()
  c.SaveAs(plotdir+h+'_CorrelationHist_SB_detail_Cov_FitMethod.png')
  
  h_chargeflip = ROOT.TH2D(h+"CFRate",";P_{T}[GeV] ; \eta;"+h+"_ChargeFlip rate",pt_bins,pt_region,eta_bins,eta_region)
  for i in range(pt_bins):
    for j in range(eta_bins):
      h_chargeflip.SetBinContent(i+1,j+1,w.var('P'+str(i)+str(j)).getValV())
      h_chargeflip.SetBinError(i+1,j+1,w.var('P'+str(i)+str(j)).getError())
  c.SetLogx()
  h_chargeflip.Draw('COLZTEXT e')
  c.Update()
  c.SaveAs(plotdir+h+'_CFRate_SB_detail_Cov_FitMethod.png')
