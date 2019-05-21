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
t_in = f_in.Get("events")
_pt = [0,10,20,25,30,40,50,60,80,100]
_eta = np.arange(-5,5,0.4)
_npv = np.arange(20,64,4)
folder = "result"
os.system("mkdir -p "+folder)
os.system("mkdir -p "+folder+"/fit")

def get_jet_params(rg, fit=None):
	params={}
	N_rg=t_in.GetEntries(rg)
	print "getting params in range: " + rg
	print "number of events in range: " + str(N_rg)
	if N_rg<3: return None
	_shape=np.arange(0.1, 5., 0.1)
	shape_pf = ROOT.TH1F("shape_pf","shape_pf",len(_shape)-1,array('d',_shape))
	shape_puppi = ROOT.TH1F("shape_puppi","shape_puppi",len(_shape)-1,array('d',_shape))
	t_in.Draw("rawjet_pt/genjet_pt>>shape_pf",rg,"goff")
	t_in.Draw("rawpjet_pt/genpjet_pt>>shape_puppi",rg,"goff")
	params["response_pf"]=abs(shape_pf.GetMean())
	params["response_pf_error"]=abs(shape_pf.GetMeanError())
	params["response_puppi"]=abs(shape_puppi.GetMean())
	params["response_puppi_error"]=abs(shape_puppi.GetMeanError())
	params["sigma_pf"]=abs(shape_pf.GetStdDev())
	params["sigma_pf_error"]=abs(shape_pf.GetStdDevError())
	params["sigma_puppi"]=abs(shape_puppi.GetStdDev())
	params["sigma_puppi_error"]=abs(shape_puppi.GetStdDevError())
	if params["response_pf"]>0: 
		params["resolution_pf"]=params["sigma_pf"]/params["response_pf"]
		params["resolution_pf_error"]=sqrt((params["sigma_pf_error"]/params["sigma_pf"])**2+(params["response_pf_error"]/params["response_pf"])**2)*params["resolution_pf"]
	else: 
		params["resolution_pf"]=0
		params["resolution_pf_error"]=0
	if params["response_puppi"]>0: 
		params["resolution_puppi"]=params["sigma_puppi"]/params["response_puppi"]
		params["resolution_puppi_error"]=sqrt((params["sigma_puppi_error"]/params["sigma_puppi"])**2+(params["response_puppi_error"]/params["response_puppi"])**2)*params["resolution_puppi"]
	else: 
		params["resolution_puppi"]=0
		params["resolution_puppi_error"]=0
	plot_hists([shape_pf,shape_puppi], legend_title_list=["PF jet","PUPPI jet"], x_title="Jet p_{T} response", y_title="Events", plot_name="fit/jet_response_"+rg,text_description=rg)
	return params


#create histograms to store response, sigma, resolution, in different phase space
h_response_pt       =  [ ROOT.TH1F("h_response_pt_pf"          , "h_response_pt_pf"          , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_response_pt_puppi"          , "h_response_pt_puppi"          , len(_pt)-1 , array('d' , _pt))]
h_sigma_pt      =  [ ROOT.TH1F("h_sigma_pt_pf"      , "h_sigma_pt_pf"      , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_sigma_pt_puppi"      , "h_sigma_pt_puppi"      , len(_pt)-1 , array('d' , _pt))]
h_resolution_pt =  [ ROOT.TH1F("h_resolution_pt_pf" , "h_resolution_pt_pf" , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_resolution_pt_puppi" , "h_resolution_pt_puppi" , len(_pt)-1 , array('d' , _pt))]
h_response_eta       =  [ ROOT.TH1F("h_response_eta_pf"          , "h_response_eta_pf"          , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_response_eta_puppi"          , "h_response_eta_puppi"          , len(_eta)-1 , array('d' , _eta))]
h_sigma_eta      =  [ ROOT.TH1F("h_sigma_eta_pf"      , "h_sigma_eta_pf"      , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_sigma_eta_puppi"      , "h_sigma_eta_puppi"      , len(_eta)-1 , array('d' , _eta))]
h_resolution_eta =  [ ROOT.TH1F("h_resolution_eta_pf" , "h_resolution_eta_pf" , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_resolution_eta_puppi" , "h_resolution_eta_puppi" , len(_eta)-1 , array('d' , _eta))]
h_response_npv       =  [ ROOT.TH1F("h_response_npv_pf"          , "h_response_npv_pf"          , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_response_npv_puppi"          , "h_response_npv_puppi"          , len(_npv)-1 , array('d' , _npv))]
h_sigma_npv      =  [ ROOT.TH1F("h_sigma_npv_pf"      , "h_sigma_npv_pf"      , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_sigma_npv_puppi"      , "h_sigma_npv_puppi"      , len(_npv)-1 , array('d' , _npv))]
h_resolution_npv =  [ ROOT.TH1F("h_resolution_npv_pf" , "h_resolution_npv_pf" , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_resolution_npv_puppi" , "h_resolution_npv_puppi" , len(_npv)-1 , array('d' , _npv))]


#Make pf vs. puppi plots
ROOT.gStyle.SetOptStat(0) #tdrStyle not in effect?
def plot_hists(hist_list, title="", legend_title_list=None, x_title="", y_title="", text_description=None, plot_name=None, limits=None):
	colors=[ROOT.kCyan+1,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kRed+1,ROOT.kOrange,ROOT.kYellow+1,ROOT.kGreen+1,ROOT.kGray]
	canv = ROOT.TCanvas("canv","canv")
	mg=ROOT.TMultiGraph() #Use a multiGraph to auto adjust the frame size:
	mg.SetTitle(title+"; "+x_title+"; "+y_title+";") 
	if len(hist_list)>0:
		for ihist,hist in enumerate(hist_list):
			hist.SetLineWidth(2)
			hist.SetMarkerStyle(1)
			hist.SetMarkerSize(0.8)
			hist.SetLineColor(colors[ihist])
			mg.Add(ROOT.TGraphErrors(hist))
	if limits != None:
		mg.SetMinimum(limits[0])
		mg.SetMaximum(limits[1])
	mg.Draw("AP")
	if(text_description):
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextSize(0.4*canv.GetTopMargin())
		latex.SetTextFont(42)
		latex.SetTextAlign(31) # align right
		latex.DrawLatex(0.90, 0.93,text_description)
		latex.Draw("same")
	if legend_title_list==None:
		legend_title_list = [i.GetTitle() for i in hist_list]
	legend=ROOT.TLegend(0.60,0.80,0.95,.95)
	for ihist,hist in enumerate(hist_list):
		legend.AddEntry(hist,legend_title_list[ihist],"lp")
	legend.Draw("same")
	if plot_name==None:
		plot_name=base_hist.GetTitle().replace("_pf","").replace("_puppi","")
	canv.Print(folder+"/"+plot_name+".png")


#call function get_jet_params and set the histograms' bin content
for ipt,pt in enumerate(_pt[:-1]):
	params=get_jet_params("qt>"+str(_pt[ipt])+" && qt<"+str(_pt[ipt+1]))
	if params == None: continue
	h_response_pt[0].SetBinContent(ipt+1,params["response_pf"]);
	h_response_pt[0].SetBinError(ipt+1,params["response_pf_error"]);
	h_sigma_pt[0].SetBinContent(ipt+1,params["sigma_pf"]);
	h_sigma_pt[0].SetBinError(ipt+1,params["sigma_pf_error"]);
	h_resolution_pt[0].SetBinContent(ipt+1,params["resolution_pf"]);
	h_resolution_pt[0].SetBinError(ipt+1,params["resolution_pf_error"]);
	h_response_pt[1].SetBinContent(ipt+1,params["response_puppi"]);
	h_response_pt[1].SetBinError(ipt+1,params["response_puppi_error"]);
	h_sigma_pt[1].SetBinContent(ipt+1,params["sigma_puppi"]);
	h_sigma_pt[1].SetBinError(ipt+1,params["sigma_puppi_error"]);
	h_resolution_pt[1].SetBinContent(ipt+1,params["resolution_puppi"]);
	h_resolution_pt[1].SetBinError(ipt+1,params["resolution_puppi_error"]);
for ieta,eta in enumerate(_eta[:-1]):
	params=get_jet_params("q_eta>"+ "%.1f"%_eta[ieta] +" && q_eta<"+"%.1f"%_eta[ieta+1])
	if params == None: continue
	h_response_eta[0].SetBinContent(ieta+1,params["response_pf"]);
	h_response_eta[0].SetBinError(ieta+1,params["response_pf_error"]);
	h_sigma_eta[0].SetBinContent(ieta+1,params["sigma_pf"]);
	h_sigma_eta[0].SetBinError(ieta+1,params["sigma_pf_error"]);
	h_resolution_eta[0].SetBinContent(ieta+1,params["resolution_pf"]);
	h_resolution_eta[0].SetBinError(ieta+1,params["resolution_pf_error"]);
	h_response_eta[1].SetBinContent(ieta+1,params["response_puppi"]);
	h_response_eta[1].SetBinError(ieta+1,params["response_puppi_error"]);
	h_sigma_eta[1].SetBinContent(ieta+1,params["sigma_puppi"]);
	h_sigma_eta[1].SetBinError(ieta+1,params["sigma_puppi_error"]);
	h_resolution_eta[1].SetBinContent(ieta+1,params["resolution_puppi"]);
	h_resolution_eta[1].SetBinError(ieta+1,params["resolution_puppi_error"]);
for inpv,npv in enumerate(_npv[:-1]):
	params=get_jet_params("npv>"+str(_npv[inpv])+" && npv<"+str(_npv[inpv+1]))
	if params == None: continue
	h_response_npv[0].SetBinContent(inpv+1,params["response_pf"]);
	h_response_npv[0].SetBinError(inpv+1,params["response_pf_error"]);
	h_sigma_npv[0].SetBinContent(inpv+1,params["sigma_pf"]);
	h_sigma_npv[0].SetBinError(inpv+1,params["sigma_pf_error"]);
	h_resolution_npv[0].SetBinContent(inpv+1,params["resolution_pf"]);
	h_resolution_npv[0].SetBinError(inpv+1,params["resolution_pf_error"]);
	h_response_npv[1].SetBinContent(inpv+1,params["response_puppi"]);
	h_response_npv[1].SetBinError(inpv+1,params["response_puppi_error"]);
	h_sigma_npv[1].SetBinContent(inpv+1,params["sigma_puppi"]);
	h_sigma_npv[1].SetBinError(inpv+1,params["sigma_puppi_error"]);
	h_resolution_npv[1].SetBinContent(inpv+1,params["resolution_puppi"]);
	h_resolution_npv[1].SetBinError(inpv+1,params["resolution_puppi_error"]);



plot_hists(h_response_pt, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="jet p_{T} [GeV]", y_title="#frac{jet p_{T}}{genjet p_{T}}", plot_name="jet_response_pt", limits=[0,1.5])
plot_hists(h_sigma_pt, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="jet p_{T} [GeV]", y_title="#sigma", plot_name="jet_sigma_pt")
plot_hists(h_resolution_pt, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="jet p_{T} [GeV]", y_title="#corrected resolution", plot_name="jet_resolution_pt", limits=[0,1])
plot_hists(h_response_eta, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="#eta", y_title="#frac{jet p_{T}}{genjet p_{T}}", plot_name="jet_response_eta", limits=[0,1.5])
plot_hists(h_sigma_eta, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="#eta", y_title="#sigma", plot_name="jet_sigma_eta")
plot_hists(h_resolution_eta, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="#eta", y_title="#corrected resolution", plot_name="jet_resolution_eta", limits=[0,1])
plot_hists(h_response_npv, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="Number of vertices", y_title="#frac{jet p_{T}}{genjet p_{T}}", plot_name="jet_response_npv", limits=[0,1.5])
plot_hists(h_sigma_npv, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="Number of vertices", y_title="#sigma", plot_name="jet_sigma_npv")
plot_hists(h_resolution_npv, legend_title_list=["PF+CHS jet", "PUPPI jet"], x_title="Number of vertices", y_title="corrected resolution", plot_name="jet_resolution_npv", limits=[0,1])
