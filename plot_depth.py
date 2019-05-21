import os,sys,socket,argparse,time
import re
import ROOT
import math
import numpy as np
from array import array
from tdrStyle import *
setTDRStyle()  

f_in = ROOT.TFile("jetmetNtuples.root","READ")
t = f_in.Get("events")
_pt = np.array([25,50,80,200])
_eta = np.array([1.3,2.1,2.5,3.0,5.0])
_rho = np.array([0,30,40,70])
_depth = np.array(range(8))
_ratio = np.arange(1,5,0.4) #try to plot the distribution of layer2/layer1

ROOT.gStyle.SetOptStat(0)


#init depth distribution histograms with different Pt
h_pt = []
r_pt = []
for ipt,pt in enumerate(_pt[:-1]):
	rangename="pt_"+str(pt)+"_"+str(_pt[ipt+1])
	h_pt.append([ROOT.TH1D("LV_"+rangename,"LV_"+rangename,len(_depth)-1,array('d',_depth)),
		ROOT.TH1D("PU_"+rangename,"PU_"+rangename,len(_depth)-1,array('d',_depth))])
	r_pt.append([ROOT.TH1D("r_LV_"+rangename,"LV_"+rangename,len(_ratio)-1,array('d',_ratio)),
		ROOT.TH1D("r_PU_"+rangename,"PU_"+rangename,len(_ratio)-1,array('d',_ratio))])
#init depth distribution histograms with different Pt
h_eta = []
r_eta = []
for ieta,eta in enumerate(_eta[:-1]):
	rangename="eta_"+str(eta)+"_"+str(_eta[ieta+1])
	h_eta.append([ROOT.TH1D("LV_"+rangename,"LV_"+rangename,len(_depth)-1,array('d',_depth)),
		ROOT.TH1D("PU_"+rangename,"PU_"+rangename,len(_depth)-1,array('d',_depth))])
	r_eta.append([ROOT.TH1D("r_LV_"+rangename,"LV_"+rangename,len(_ratio)-1,array('d',_ratio)),
		ROOT.TH1D("r_PU_"+rangename,"PU_"+rangename,len(_ratio)-1,array('d',_ratio))])
#init depth distribution histograms with different Pt
h_rho = []
r_rho = []
for irho,rho in enumerate(_rho[:-1]):
	rangename="rho_"+str(rho)+"_"+str(_rho[irho+1])
	h_rho.append([ROOT.TH1D("LV_"+rangename,"LV_"+rangename,len(_depth)-1,array('d',_depth)),
		ROOT.TH1D("PU_"+rangename,"PU_"+rangename,len(_depth)-1,array('d',_depth))])
	r_rho.append([ROOT.TH1D("r_LV_"+rangename,"LV_"+rangename,len(_ratio)-1,array('d',_ratio)),
		ROOT.TH1D("r_PU_"+rangename,"PU_"+rangename,len(_ratio)-1,array('d',_ratio))])

#init depth distribution histograms with different PU and Pt
h_pt_rho = []
r_pt_rho = []
for ipt, pt in enumerate(_pt[:-1]):
	h_rho_tmp = []
	r_rho_tmp = []
	for irho,rho in enumerate(_rho[:-1]):
		rangename="rho_"+str(rho)+"_"+str(_rho[irho+1])+"_pt_"+str(pt)+"_"+str(_pt[ipt+1])
		h_rho_tmp.append([ROOT.TH1D("LV_"+rangename,"LV_"+rangename,len(_depth)-1,array('d',_depth)),
			ROOT.TH1D("PU_"+rangename,"PU_"+rangename,len(_depth)-1,array('d',_depth))])
		r_rho_tmp.append([ROOT.TH1D("r_LV_"+rangename,"LV_"+rangename,len(_ratio)-1,array('d',_ratio)),
			ROOT.TH1D("r_PU_"+rangename,"PU_"+rangename,len(_ratio)-1,array('d',_ratio))])
	h_pt_rho.append(h_rho_tmp)
	r_pt_rho.append(r_rho_tmp)


#Load branches
from load_jetmet_tree import *
declare_branches(t)

#in order to print out the progress
def print_same_line(s):
	sys.stdout.write(s)                     # just print
	sys.stdout.flush()                      # needed for flush when using \x08
	sys.stdout.write((b'\x08' * len(s)).decode())# back n chars    	
	#time.sleep(0.2)

ROOT.TH1.SetDefaultSumw2()
#Loop through tree
nentries = t.GetEntriesFast()
nbytes,nb=0,0
for jentry in range(nentries):
	print_same_line('Now loading ... '+str(round(100.*jentry/nentries,2))+'%')
	ientry = t.LoadTree(jentry)
	if (ientry < 0): break
	nb = t.GetEntry(jentry)
	nbytes+=nb
	#print njet,npjet
	for ipjet in range(npjet[0]):
		ipt=_pt.searchsorted(pjet_pt[ipjet])-1
		ieta=_eta.searchsorted(pjet_eta[ipjet])-1
		irho=_rho.searchsorted(rhoall[0])-1
		isPU = [0,1][genpjet_pt[ipjet]<0] #if the puppi jet doesn't have a genjet, it will be -999
		if pjet_eta[ipjet]<1.3: continue #ignore barrel region for now
		if pjet_depth[ipjet][0]==0: continue
		if ipt>=0 and ipt<len(h_pt): 
			for idepth in range(len(_depth)-1):
				h_pt[ipt][isPU].Fill(idepth+0.5, pjet_depth[ipjet][idepth])
			r_pt[ipt][isPU].Fill(pjet_depth[ipjet][1]/pjet_depth[ipjet][0],1)
			if irho>=0 and irho<len(h_rho): 
				for idepth in range(len(_depth)-1):
					h_pt_rho[ipt][irho][isPU].Fill(idepth+0.5, pjet_depth[ipjet][idepth])
				r_pt_rho[ipt][irho][isPU].Fill(pjet_depth[ipjet][1]/pjet_depth[ipjet][0],1)
		if ieta>=0 and ieta<len(h_eta): 
			for idepth in range(len(_depth)-1):
				h_eta[ieta][isPU].Fill(idepth+0.5, pjet_depth[ipjet][idepth])
			r_eta[ieta][isPU].Fill(pjet_depth[ipjet][1]/pjet_depth[ipjet][0],1)
		if irho>=0 and irho<len(h_rho): 
			for idepth in range(len(_depth)-1):
				h_rho[irho][isPU].Fill(idepth+0.5, pjet_depth[ipjet][idepth])
			r_rho[irho][isPU].Fill(pjet_depth[ipjet][1]/pjet_depth[ipjet][0],1)

f_out = ROOT.TFile("hists.root","recreate")

#Plot the histograms
def plot_hists(hists,filename):
	c = ROOT.TCanvas("c","c",800,800) #who cares about your name?
	colors=[ROOT.kCyan+1,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kRed+1,ROOT.kOrange,ROOT.kYellow+1,ROOT.kGreen+1,ROOT.kGray] 
	firstDraw=True
	legend=ROOT.TLegend(0.60,0.65,0.95,1)
	for ihist,hist in enumerate(hists):
		hist[0].SetLineColor(colors[ihist])
		hist[1].SetLineColor(colors[ihist])
		hist[1].SetLineStyle(8) #dash line
		if firstDraw:
			firstDraw=False
			hist[0].GetXaxis().SetTitle("layer")
			hist[0].GetXaxis().SetTitleOffset(1.2)
			hist[0].GetYaxis().SetTitle("Percentage of distributed energy")
			hist[0].GetYaxis().SetTitleOffset(1.3)
			hist[0].SetTitle("")
			if hist[0].GetName()[0]=='r':
				hist[0].GetXaxis().SetTitle("layer2/layer1")
				hist[0].GetYaxis().SetTitle("Distribution")
			htmp=hist[0].DrawNormalized()
			htmp.GetYaxis().SetRangeUser(0,1.05)
			htmp.Draw()
		else:
			hist[0].DrawNormalized("same")
		hist[1].DrawNormalized("same")
		legend.AddEntry(hist[0], hist[0].GetName(),"lp")
		legend.AddEntry(hist[1], hist[1].GetName(),"lp")
		hist[0].Write()
		hist[1].Write()
	legend.Draw("same")
	c.Print("result/"+filename)

plot_hists(h_pt, "depth_pt_distribution.png")
plot_hists(h_eta, "depth_eta_distribution.png")
plot_hists(h_rho, "depth_rho_distribution.png")
plot_hists(r_pt, "r21_pt_distribution.png")
plot_hists(r_eta, "r21_eta_distribution.png")
plot_hists(r_rho, "r21_rho_distribution.png")
for ihists,hists in enumerate(h_pt_rho):
	plot_hists(hists, "depth_pt_"+str(_pt[ihists])+"_rho_distribution.png")
for ihists,hists in enumerate(r_pt_rho):
	plot_hists(hists, "r21_pt_"+str(_pt[ihists])+"_rho_distribution.png")
