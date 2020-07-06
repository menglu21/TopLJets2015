import ROOT 
import sys
#import argparse
#import textwrap

#parser = argparse.ArgumentParser(
#    formatter_class=argparse.RawDescriptionHelpFormatter,
#    description=textwrap.dedent('''\
#		Enter the MC and data files'''))
#args = parser.parse_args()

def plot_mc_data(outname,list1,list2,list3):
	pname=outname
	pList=list1
	pList2=list2
	pList3=list3

	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)
	c=ROOT.TCanvas('c','c',500,500)

	histos=[]
	for p,title,color in pList:
		histos.append( f1.Get(p) )
		histos[-1].SetLineWidth(2)
		histos[-1].SetLineColor(color)
		histos[-1].SetTitle(title)
		histos[-1].DrawNormalized('Esame' if len(histos)>1 else 'hist')
	for p,title,color in pList2:
		histos.append( f2.Get(p) )
		histos[-1].SetLineWidth(2)
		histos[-1].SetLineColor(color)
		histos[-1].SetTitle(title)
		histos[-1].DrawNormalized('same' if len(histos)>1 else 'hist')
	for p,title,color in pList3:
		histos.append( f3.Get(p) )
		histos[-1].SetLineWidth(2)
		histos[-1].SetLineColor(color+1)
		histos[-1].SetTitle(title)
		histos[-1].DrawNormalized('same' if len(histos)>1 else 'hist')	

	c.BuildLegend(0.65,0.95,0.95,0.8)
	c.Modified()
	c.Update()
	c.SaveAs(pname+'.png')


#open root file
f1=ROOT.TFile.Open(sys.argv[1])
f2=ROOT.TFile.Open(sys.argv[2])
f3=ROOT.TFile.Open(sys.argv[3])

pname='mlb'
pList=[('inc_mlb','MC',1)]
pList2=[('inc_mlb','Data',2)]
pList3=[('inc_mlb','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='mu_pt'
pList=[('inc_mu_pt','MC',1)]
pList2=[('inc_mu_pt','Data',2)]
pList3=[('inc_mu_pt','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='mu_eta'
pList=[('inc_mu_eta','MC',1)]
pList2=[('inc_mu_eta','Data',2)]
pList3=[('inc_mu_eta','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='njets'
pList=[('inc_njets','MC',1)]
pList2=[('inc_njets','Data',2)]
pList3=[('inc_njets','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='nbjets'
pList=[('inc_nbjets','MC',1)]
pList2=[('inc_nbjets','Data',2)]
pList3=[('inc_nbjets','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='nvtx'
pList2=[('inc_nvtx','MC',1)]
pList=[('inc_nvtx','Data',2)]
pList3=[('inc_nvtx','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='bjet_eta'
pList=[('inc_bjet_eta','MC',1)]
pList2=[('inc_bjet_eta','Data',2)]
pList3=[('inc_bjet_eta','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='bjet_pt'
pList=[('inc_bjet_pt','MC',1)]
pList2=[('inc_bjet_pt','Data',2)]
pList3=[('inc_bjet_pt','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

pname='met'
pList=[('inc_met_pt','MC',1)]
pList2=[('inc_met_pt','Data',2)]
pList3=[('inc_met_pt','CGToBHpm_MH-300-rtc04-rtt06',3)]
plot_mc_data(pname,pList,pList2,pList3)

#all done
f1.Close()
f2.Close()
f3.Close()
