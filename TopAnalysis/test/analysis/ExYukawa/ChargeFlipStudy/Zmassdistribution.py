import os
import ROOT
import json
import sys
from collections import OrderedDict
import numpy as np
import argparse

##########
# Parser #
##########

parser = argparse.ArgumentParser(description = "Charge Flip Rate")
parser.add_argument('plotdir',type=str)
args = parser.parse_args()

#################
# Basic Setting #
#################

plotdir = args.plotdir # Change to your own directory

# Root Files Directory
workdir = os.getcwd()
targetdir = os.path.join(workdir+"/analysis_2017/Chunks/")
Allfile = os.listdir(targetdir)

# Json File
jsonfile = open(os.path.join(workdir+"/samples_2017_chargeflip.json"))
samplesList = json.load(jsonfile, encoding = 'utf-8', object_pairs_hook=OrderedDict).items()
jsonfile.close()

# Set up histograms
h_oc_Zmass_MC = ROOT.TH1F('oc_Zmass_MC',"; M(l,l) [GeV] ; Events",150,0,150)
h_oc_Zmass_data = ROOT.TH1F('oc_Zmass_data',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_MC = ROOT.TH1F('ss_Zmass_MC',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_data = ROOT.TH1F('ss_Zmass_data',"; M(l,l) [GeV] ; Events",150,0,150)
h_ss_Zmass_DMinC = ROOT.TH1F('ss_Zmass_DMinC',"; M(l,l) [GeV] ; Events",150,0,150)

diagrams = ['ee_ss_Zmass', 'ee_oc_Zmass']
Zmass_result = dict()

#############
# Read Data #
#############

for s,desc in samplesList:
  weight = desc[0]*41.5*1000
  for ff in Allfile:
    if s in ff:
      inf = ROOT.TFile.Open(os.path.join(targetdir+ff))
      try:
        if "MC" in ff:
          h_oc = inf.Get('ee_oc_Zmass')
          h_ss = inf.Get('ee_ss_Zmass')
          h_oc_Zmass_MC.Add(h_oc_Zmass_MC,h_oc,1.0,weight)
          h_ss_Zmass_MC.Add(h_ss_Zmass_MC,h_ss,1.0,weight)
        else:
          h_oc = inf.Get('ee_oc_Zmass')
          h_ss = inf.Get('ee_ss_Zmass')
          h_oc_Zmass_data.Add(h_oc_Zmass_data,h_oc,1.0,1.0)
          h_ss_Zmass_data.Add(h_ss_Zmass_data,h_ss,1.0,1.0)
#        print(s)
      except:
        pass
      inf.Close()

###############
# Z Mass Plot #
###############

c = ROOT.TCanvas('c','c',600,400)
#ROOT.gStyle.SetOptStat("kFALSE")

h_ss_Zmass_DMinC.Add(h_ss_Zmass_data,h_ss_Zmass_MC,1.0,-1.0)
h_ss_Zmass_DMinC.Draw()
h_ss_Zmass_DMinC.SetMarkerStyle(20)
c.Update()
c.SaveAs(plotdir+'SS_DminC.png')

h_oc_Zmass_MC.GetXaxis().SetRange(70,110)
h_oc_Zmass_MC.Draw("HIST")
c.Update()

h_oc_Zmass_data.GetXaxis().SetRange(70,110)
h_oc_Zmass_data.SetMarkerStyle(20)
h_oc_Zmass_data.Draw("E2 SAME")
c.Update()

#print(h_ss_Zmass_MC.GetMaximum())

rightmax = 1.1 * float(h_ss_Zmass_MC.GetMaximum())
scale = ROOT.gPad.GetUymax()/rightmax

h_ss_Zmass_MC.GetXaxis().SetRange(70,110)
h_ss_Zmass_MC.SetLineColor(2)
h_ss_Zmass_MC.Scale(scale)
h_ss_Zmass_MC.Draw("HIST SAME")
c.Update()

h_ss_Zmass_data.GetXaxis().SetRange(70,110)
h_ss_Zmass_data.SetLineColor(2)
h_ss_Zmass_data.Scale(scale)
h_ss_Zmass_data.SetMarkerStyle(20)
h_ss_Zmass_data.SetMarkerColor(2)
h_ss_Zmass_data.Draw("E2 SAME")
c.Update()
legend = ROOT.TLegend(0.1,0.7,0.44448,0.9)
legend.AddEntry(h_ss_Zmass_data,"Data(SS)","p")
legend.AddEntry(h_ss_Zmass_MC,"MC(SS)","l")
legend.AddEntry(h_oc_Zmass_data,"Data(OC)","p")
legend.AddEntry(h_oc_Zmass_MC,"MC(OC)","l")
legend.Draw("SAME")

axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(),ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymax(),0,rightmax,20,"+L")
axis.SetLineColor(2)
axis.SetLabelColor(2)
axis.Draw()
c.Update()

c.SaveAs(plotdir+'Zmass.png')

###########
# Fitting #
###########

w = ROOT.RooWorkspace()
#w.factory('Gaussian::pdf(x[0,150],mean[90,50,130],gamma[0,20])')
w.factory('BreitWigner::BW(x[0,150],mean[90,70,110],gamma[0,20])')
w.factory('Exponential::e(x,alpha[-5,-10,0.1])')
w.factory('SUM::pdf(a[9000,0,10000000]*BW,b[10,0,10000]*e)')
w.factory('Gaussian::g(x,M[90,70,110],sigma[0,20])')
#w.factory('CBShape::CB(x,m[0,150],s[0,150],a[0,150],n[-5,5])')
#w.factory('FCONV::pdf(x,BW,CB)')

pdf = w.pdf('pdf')
x = w.var('x')
gamma = w.var('gamma')
mean = w.var('mean')
BW = w.pdf('BW')
g = w.pdf('g')

# Data minus MC Plot in Same Sign channel ( Fit with Gaussain )
dh_ss_DMinC = ROOT.RooDataHist('dh_ss_DmC','dh_ss_DmC',ROOT.RooArgList(x),h_ss_Zmass_DMinC)
print(type(dh_ss_DMinC))
g.fitTo(dh_ss_DMinC)
xframe = x.frame()
dh_ss_DMinC.plotOn(xframe)
g.plotOn(xframe)
xframe.Draw()
c.SaveAs(plotdir+'fit_ss_DminC.png')


meanlist = [0,0,0,0]
gammalist = [0,0,0,0]
widthlist = [0,0,0,0]
wglist = [0,0,0,0]
goallist = [0,0,0,0]

# Same Sign Data
dh_ss_data = ROOT.RooDataHist('dh_ss_data','dh_ss_data',ROOT.RooArgList(x),h_ss_Zmass_data)
pdf.fitTo(dh_ss_data)
xframe = x.frame()
dh_ss_data.plotOn(xframe)
pdf.plotOn(xframe)
pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
xframe.Draw()
mean.Print()
gamma.Print()
c.SaveAs(plotdir+'fit_ss_data.png')

meanlist[0] = mean.getValV()
gammalist[0] = gamma.getValV()
ig = BW.createIntegral(ROOT.RooArgSet(x))
goal = 0.
width = 0.
while(goal<0.95):
  width+=0.1
  x.setRange(mean.getValV()-width,mean.getValV()+width)
  goal = ig.getVal()
  x.setRange(0,150)
  goal = goal/ig.getVal()
  print(str(width)+" : "+str(goal)+" "+ str(width/gamma.getValV()))
widthlist[0] = width
goallist[0] = goal
wglist[0] = width/gamma.getValV()

# Same Sign MC
dh_ss_MC = ROOT.RooDataHist('dh_ss_MC','dh_ss_MC',ROOT.RooArgList(x),h_ss_Zmass_MC)
pdf.fitTo(dh_ss_MC)
xframe = x.frame()
dh_ss_MC.plotOn(xframe)
pdf.plotOn(xframe)
pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
xframe.Draw()
mean.Print()
gamma.Print()
c.SaveAs(plotdir+'fit_ss_MC.png')

meanlist[1] = mean.getValV()
gammalist[1] = gamma.getValV()
ig = BW.createIntegral(ROOT.RooArgSet(x))
goal = 0.
width = 0.
while(goal<0.95):
  width+=0.1
  x.setRange(mean.getValV()-width,mean.getValV()+width)
  goal = ig.getVal()
  x.setRange(0,150)
  goal = goal/ig.getVal()
  print(str(width)+" : "+str(goal)+" "+ str(width/gamma.getValV()))
widthlist[1] = width
goallist[1] = goal
wglist[1] = width/gamma.getValV()

# Opposite Sign Charge Data
dh_oc_data = ROOT.RooDataHist('dh_oc_data','dh_oc_data',ROOT.RooArgList(x),h_oc_Zmass_data)
pdf.fitTo(dh_oc_data)
xframe = x.frame()
dh_oc_data.plotOn(xframe)
pdf.plotOn(xframe)
pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineStyle(2))
xframe.Draw()
mean.Print()
gamma.Print()
c.SaveAs(plotdir+'fit_oc_data.png')

meanlist[2] = mean.getValV()
gammalist[2] = gamma.getValV()
ig = BW.createIntegral(ROOT.RooArgSet(x))
goal = 0.
width = 0.
while(goal<0.95):
  width+=0.1
  x.setRange(mean.getValV()-width,mean.getValV()+width)
  goal = ig.getVal()
  x.setRange(0,150)
  goal = goal/ig.getVal()
  print(str(width)+" : "+str(goal)+" "+ str(width/gamma.getValV()))
widthlist[2] = width
goallist[2] = goal
wglist[2] = width/gamma.getValV()

# Opposite Sign Charge MC
dh_oc_MC = ROOT.RooDataHist('dh_oc_MC','dh_oc_MC',ROOT.RooArgList(x),h_oc_Zmass_MC)
pdf.fitTo(dh_oc_MC)
xframe = x.frame()
dh_oc_MC.plotOn(xframe)
pdf.plotOn(xframe)
pdf.plotOn(xframe,ROOT.RooFit.Components("e"),ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
xframe.Draw()
mean.Print()
gamma.Print()
c.SaveAs(plotdir+'fit_oc_MC.png')
xframe = x.frame()


meanlist[3] = mean.getValV()
gammalist[3] = gamma.getValV()

###################
# Determine Width #
###################

# Choose the width that contain 95% of the sign event.
ig = BW.createIntegral(ROOT.RooArgSet(x))
goal = 0.
width = 0.
while(goal<0.95):
  width+=0.1
  x.setRange(mean.getValV()-width,mean.getValV()+width)
  goal = ig.getVal()
  x.setRange(0,150)
  goal = goal/ig.getVal()
  print(str(width)+" : "+str(goal)+" "+ str(width/gamma.getValV()))
widthlist[3] = width
goallist[3] = goal
wglist[3] = width/gamma.getValV()
#ig.PlotOn(xframe)
#c.Update()
#c.SaveAs(plotdir+"CDF.png")

print(meanlist)
print(gammalist)
print(widthlist)
print(goallist)
print(wglist)
