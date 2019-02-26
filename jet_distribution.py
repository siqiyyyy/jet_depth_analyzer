import os,sys,socket,argparse
import re
import ROOT
import math
import numpy as np
from array import array
from tdrStyle import *
setTDRStyle()  

f_in = ROOT.TFile("jetmetNtuples.root","READ")
t_in = f_in.Get("events")
_pt = np.array([20,30,40,50,100,200,500,2000])
_eta = np.arange(0,5,0.2)
_rho = np.arange(0,100,10)

c1 = ROOT.TCanvas("c1","c1")
ROOT.gStyle.SetOptStat(0)


#Eta distribution
h_eta = [ROOT.TH1F("h_eta_pf","h_eta_pf",len(_eta)-1,array('d',_eta)),
	ROOT.TH1F("h_eta_puppi","h_eta_puppi",len(_eta)-1,array('d',_eta))]
t_in.Draw("jet_eta>>h_eta_pf","","goff")
t_in.Draw("pjet_eta>>h_eta_puppi","","goff")
h_eta[0].SetLineColor(ROOT.kLake)
h_eta[1].SetLineColor(ROOT.kRose)
h_eta[0].GetXaxis().SetTitle("jet #eta [GeV]")
h_eta[0].GetYaxis().SetTitle("Events")
h_eta[0].SetTitle("")
h_eta[0].GetXaxis().SetTitleOffset(1.2)
h_eta[0].GetYaxis().SetTitleOffset(1.3)
h_eta[0].SetMinimum(0)
h_eta[0].Draw()
h_eta[1].Draw("same")
legend=ROOT.TLegend(0.60,0.65,0.95,.85)
legend.AddEntry(h_eta[0],"pf+CHS","lp")
legend.AddEntry(h_eta[1],"pf+PUPPI","lp")
legend.Draw("same")
c1.Print("result/sample_eta_distribution.png")

#rho distribution
h_rho = ROOT.TH1F("h_rho","h_rho",len(_rho)-1,array('d',_rho))
t_in.Draw("rhoall>>h_rho","","goff")
h_rho.GetXaxis().SetTitle("event #rho [GeV]")
h_rho.GetYaxis().SetTitle("Events")
h_rho.SetTitle("")
h_rho.GetXaxis().SetTitleOffset(1.2)
h_rho.GetYaxis().SetTitleOffset(1.3)
h_rho.Draw()
c1.Print("result/sample_rho_distribution.png")

#PT distribution
c1.SetLogx()
h_pt = [ROOT.TH1F("h_pt_pf","h_pt_pf",len(_pt)-1,array('d',_pt)),
	ROOT.TH1F("h_pt_puppi","h_pt_puppi",len(_pt)-1,array('d',_pt))]
t_in.Draw("jet_pt>>h_pt_pf","","goff")
t_in.Draw("pjet_pt>>h_pt_puppi","","goff")
h_pt[0].SetLineColor(ROOT.kLake)
h_pt[1].SetLineColor(ROOT.kRose)
h_pt[0].GetXaxis().SetTitle("jet p_{T} [GeV]")
h_pt[0].GetYaxis().SetTitle("Events")
h_pt[0].SetTitle("")
h_pt[0].GetXaxis().SetTitleOffset(1.2)
h_pt[0].GetYaxis().SetTitleOffset(1.3)
h_pt[0].Draw()
h_pt[1].Draw("same")
legend=ROOT.TLegend(0.60,0.65,0.95,.85)
legend.AddEntry(h_pt[0],"pf+CHS","lp")
legend.AddEntry(h_pt[1],"pf+PUPPI","lp")
legend.Draw("same")
c1.Print("result/sample_pt_distribution.png")
