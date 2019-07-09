
import os,sys,socket,argparse
import re
import ROOT
import math
import numpy as np
from math import sqrt
from array import array
from tdrStyle import *
setTDRStyle()  
ROOT.gROOT.SetBatch(True)

# RooFit
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
ROOT.gROOT.SetStyle("Plain")
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )

from helper_functions import *

f_in = ROOT.TFile("jetmetNtuples.root","READ")
t = f_in.Get("events")
_pt = [0,10,20,30,40,50,60,80,90,100,200]
_eta = np.arange(-5,5,0.4)
_npv = np.arange(20,64,4)
folder = "result"
os.system("mkdir -p "+folder)

#ROOT.gStyle.SetPalette(ROOT.kOcean)
h_jet_tot=ROOT.TH1F("h_jet_tot","h_jet_tot",len(_pt)-1,array('d',_pt))
h_jet_charged=ROOT.TH1F("h_jet_charged","charged",len(_pt)-1,array('d',_pt))
h_jet_neutral=ROOT.TH1F("h_jet_neutral","neutral",len(_pt)-1,array('d',_pt))
h_jet_photon=ROOT.TH1F("h_jet_photon","photon",len(_pt)-1,array('d',_pt))
h_jet_muon=ROOT.TH1F("h_jet_muon","muon",len(_pt)-1,array('d',_pt))
h_jet_electron=ROOT.TH1F("h_jet_electron","electron",len(_pt)-1,array('d',_pt))
h_jet_hhf=ROOT.TH1F("h_jet_hhf","hhf",len(_pt)-1,array('d',_pt))
h_jet_ehf=ROOT.TH1F("h_jet_ehf","ehf",len(_pt)-1,array('d',_pt))
jet_components=[h_jet_charged,h_jet_neutral,h_jet_photon,h_jet_muon,h_jet_electron,h_jet_hhf,h_jet_ehf]

#Load branches
from load_jetmet_tree import *
declare_branches(t)

#in order to print out the progress
def print_same_line(s):
	sys.stdout.write(s)                     # just print
	sys.stdout.flush()                      # needed for flush when using \x08
	sys.stdout.write((b'\x08' * len(s)).decode())# back n chars    	
	#time.sleep(0.2)

nentries = t.GetEntriesFast()
nbytes,nb=0,0
for jentry in range(nentries):
	print_same_line('Now loading ... '+str(round(100.*jentry/nentries,2))+'%')
	ientry = t.LoadTree(jentry)
	if (ientry < 0): break
	nb = t.GetEntry(jentry)
	nbytes+=nb
	#print njet,npjet
	for ijet in range(njet[0]):
		if abs(rawjet_eta[ijet])<2.5 or rawjet_eta[ijet]>3.0: continue
		h_jet_tot.Fill(rawjet_pt[ijet],rawjet_energy[ijet])
		h_jet_charged.Fill(rawjet_pt[ijet],charged_e[ijet])
		h_jet_neutral.Fill(rawjet_pt[ijet],neutral_e[ijet])
		h_jet_photon.Fill(rawjet_pt[ijet],photon_e[ijet])
		h_jet_muon.Fill(rawjet_pt[ijet],muon_e[ijet])
		h_jet_electron.Fill(rawjet_pt[ijet],electron_e[ijet])
		h_jet_hhf.Fill(rawjet_pt[ijet],hhf_e[ijet])
		h_jet_ehf.Fill(rawjet_pt[ijet],ehf_e[ijet])

hstack=ROOT.THStack("hstack","jet energy fraction")
legend=ROOT.TLegend(0.8,0.6,1,0.95)
for ihist,hist in enumerate(jet_components):
	hist.Divide(h_jet_tot)
	hist.SetLineWidth(0)
	hstack.Add(hist)
	legend.AddEntry(hist,hist.GetTitle())

canv=ROOT.TCanvas("canv","canv")
hstack.Draw("hist PFC")
legend.Draw("same")
canv.SetLogx()
canv.Print(folder+"/jet_composition_2p5_3p0.png")
