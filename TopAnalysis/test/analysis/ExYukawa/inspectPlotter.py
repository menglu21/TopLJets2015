import ROOT
import sys

#open root file
fIn=ROOT.TFile.Open(sys.argv[1])

#draw these plots
pname='csi'
pList=[('inc_csi','inclusive',1),
       ('rp123_csi','RP=123',ROOT.kOrange+1),
       ('rp23_csi','RP=23',ROOT.kMagenta+2)]

#pname='mlb'
#pList=[('inc_mlb','inclusive',1),
#       ('plus_mlb','#mu^{+}',ROOT.kOrange+1),
#       ('minus_mlb','#mu^{-}',ROOT.kMagenta+2)]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
histos=[]
for p,title,color in pList:
    histos.append( fIn.Get(p) )
    histos[-1].SetLineWidth(2)
    histos[-1].SetLineColor(color)
    histos[-1].SetTitle(title)
    histos[-1].Draw('same' if len(histos)>1 else 'hist')
c.BuildLegend(0.65,0.95,0.95,0.8)
c.Modified()
c.Update()
c.SaveAs(pname+'.png')

#all done
fIn.Close()
