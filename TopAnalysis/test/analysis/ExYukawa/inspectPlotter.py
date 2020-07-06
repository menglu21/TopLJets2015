import ROOT
import sys
import argparse


parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
		Enter the MC and data files'''))
args = parser.parse_args()

#open root file
f1=ROOT.TFile.Open(sys.argv[1])
f2=ROOT.TFile.Open(sys.argv[2])

#draw these plots
#pname='csi'
#pList=[('inc_csi','inclusive',1),
#       ('rp123_csi','RP=123',ROOT.kOrange+1),
#       ('rp23_csi','RP=23',ROOT.kMagenta+2)]

pname='mlb'
#pList=[('inc_mlb','inclusive',1),
#       ('plus_mlb','#mu^{+}',ROOT.kOrange+1),
#       ('minus_mlb','#mu^{-}',ROOT.kMagenta+2)]
pList=[('inc_mlb','MC',1)]
pList2=[('inc_mlb','Data',2)]


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
histos=[]
for p,title,color in pList:
    histos.append( f1.Get(p) )
    histos[-1].SetLineWidth(2)
    histos[-1].SetLineColor(color)
    histos[-1].SetTitle(title)
    histos[-1].DrawNormalized('same' if len(histos)>1 else 'hist')



for p,title,color in pList2:
    histos.append( f2.Get(p) )
    histos[-1].SetLineWidth(2)
    histos[-1].SetLineColor(color)
    histos[-1].SetTitle(title)
    histos[-1].DrawNormalized('Esame' if len(histos)>1 else 'hist')


c.BuildLegend(0.65,0.95,0.95,0.8)
c.Modified()
c.Update()
c.SaveAs(pname+'.png')

#all done
f1.Close()
f2.Close()
