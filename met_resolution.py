import os,sys,socket,argparse
import re
import ROOT
import math
from math import sqrt
from array import array
from tdrStyle import *
import numpy as np
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
_pt = [0,5,10,15,20,25,30,40,50,60,80,100]
_eta = np.arange(-5,5,0.2)
_npv = np.arange(0,70,2)
folder = "result"
os.system("mkdir -p "+folder)
os.system("mkdir -p "+folder+"/fit")



def get_met_params(rg): #given a phase space range, return response, sigma u_pll, sigma uprp, and their errors
	params={} #this dict is to be returned 
	#get the responses
	#mean of upll devided by the mean of qt
	#do the same thing for both pf and puppi
	print "getting params in range: " + rg
	print "number of events in range: " + str(t_in.GetEntries(rg))
	_shape=np.arange(-200,200,5)
	shape_qt = ROOT.TH1F("shape_qt","shape_qt",len(_shape)-1,array('d',_shape))
	shape_pf = ROOT.TH1F("shape_pf","shape_pf",len(_shape)-1,array('d',_shape)) #this can be used to fill several different variables
	shape_puppi = ROOT.TH1F("shape_puppi","shape_puppi",len(_shape)-1,array('d',_shape)) #this can be used to fill several different variables
	t_in.Draw("qt>>shape_qt", rg)
	t_in.Draw("u_pll_pf>>shape_pf", rg)
	if shape_qt.GetMean()==0 or shape_pf.GetMean()==0: return None
	params["response_pf"]=abs(shape_pf.GetMean()/shape_qt.GetMean())
	params["response_pf_error"]=sqrt((shape_pf.GetMeanError()/shape_pf.GetMean())**2+(shape_qt.GetMeanError()/shape_qt.GetMean())**2) * params["response_pf"]
	t_in.Draw("u_pll_puppi>>shape_puppi",rg)
	if shape_puppi.GetMean()==0: return None
	params["response_puppi"]=abs(shape_puppi.GetMean()/shape_qt.GetMean())
	params["response_puppi_error"]=sqrt((shape_puppi.GetMeanError()/shape_puppi.GetMean())**2+(shape_qt.GetMeanError()/shape_qt.GetMean())**2) * params["response_puppi"]
	plot_hists([shape_qt,shape_pf,shape_puppi], legend_title_list=["q_{T}", "u_{||}(PF)","u_{||}(PUPPI)"], x_title="p_{T}(GeV)", y_title="Events/5GeV", plot_name="fit/qt_upll"+rg,text_description=rg)
	t_in.Draw("u_pll_pf>>shape_pf",rg)
	if shape_pf.GetStdDev()<=0: return None
	params["sigma_pll_pf"]=shape_pf.GetStdDev()
	params["sigma_pll_pf_error"]=shape_pf.GetStdDevError()
	params["resolution_pll_pf"]=params["sigma_pll_pf"]/params["response_pf"]
	params["resolution_pll_pf_error"]=sqrt((params["sigma_pll_pf_error"]/params["sigma_pll_pf"])**2+(params["response_pf_error"]/params["response_pf"])**2)*params["resolution_pll_pf"]
	t_in.Draw("u_pll_puppi>>shape_puppi",rg)
	if shape_puppi.GetStdDev()<=0: return None
	params["sigma_pll_puppi"]=shape_puppi.GetStdDev()
	params["sigma_pll_puppi_error"]=shape_puppi.GetStdDevError()
	params["resolution_pll_puppi"]=params["sigma_pll_puppi"]/params["response_puppi"]
	params["resolution_pll_puppi_error"]=sqrt((params["sigma_pll_puppi_error"]/params["sigma_pll_puppi"])**2+(params["response_puppi_error"]/params["response_puppi"])**2)*params["resolution_pll_puppi"]
	t_in.Draw("u_prp_pf>>shape_pf",rg)
	if shape_pf.GetStdDev()<=0: return None
	params["sigma_prp_pf"]=shape_pf.GetStdDev()
	params["sigma_prp_pf_error"]=shape_pf.GetStdDevError()
	params["resolution_prp_pf"]=params["sigma_prp_pf"]/params["response_pf"]
	params["resolution_prp_pf_error"]=sqrt((params["sigma_prp_pf_error"]/params["sigma_prp_pf"])**2+(params["response_pf_error"]/params["response_pf"])**2)*params["resolution_prp_pf"]
	t_in.Draw("u_prp_puppi>>shape_puppi",rg)
	if shape_puppi.GetStdDev()<=0: return None
	params["sigma_prp_puppi"]=shape_puppi.GetStdDev()
	params["sigma_prp_puppi_error"]=shape_puppi.GetStdDevError()
	params["resolution_prp_puppi"]=params["sigma_prp_puppi"]/params["response_puppi"]
	params["resolution_prp_puppi_error"]=sqrt((params["sigma_prp_puppi_error"]/params["sigma_prp_puppi"])**2+(params["response_puppi_error"]/params["response_puppi"])**2)*params["resolution_prp_puppi"]
	plot_hists([shape_pf,shape_puppi], legend_title_list=["PF met","PUPPI MET"], x_title="u_{#perp}(GeV)", y_title="Events/5GeV", plot_name="fit/uprp"+rg,text_description=rg)
	return params
	

#create histograms to store response, sigma, resolution, in different phase space
h_response_pt       =  [ ROOT.TH1F("h_response_pt_pf"          , "h_response_pt_pf"          , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_response_pt_puppi"          , "h_response_pt_puppi"          , len(_pt)-1 , array('d' , _pt))]
h_sigma_pll_pt      =  [ ROOT.TH1F("h_sigma_pll_pt_pf"      , "h_sigma_pll_pt_pf"      , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_sigma_pll_pt_puppi"      , "h_sigma_pll_pt_puppi"      , len(_pt)-1 , array('d' , _pt))]
h_sigma_prp_pt      =  [ ROOT.TH1F("h_sigma_prp_pt_pf"      , "h_sigma_prp_pt_pf"      , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_sigma_prp_pt_puppi"      , "h_sigma_prp_pt_puppi"      , len(_pt)-1 , array('d' , _pt))]
h_resolution_pll_pt =  [ ROOT.TH1F("h_resolution_pll_pt_pf" , "h_resolution_pll_pt_pf" , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_resolution_pll_pt_puppi" , "h_resolution_pll_pt_puppi" , len(_pt)-1 , array('d' , _pt))]
h_resolution_prp_pt =  [ ROOT.TH1F("h_resolution_prp_pt_pf" , "h_resolution_prp_pt_pf" , len(_pt)-1 , array('d' , _pt)) , ROOT.TH1F("h_resolution_prp_pt_puppi" , "h_resolution_prp_pt_puppi" , len(_pt)-1 , array('d' , _pt))]
h_response_eta       =  [ ROOT.TH1F("h_response_eta_pf"          , "h_response_eta_pf"          , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_response_eta_puppi"          , "h_response_eta_puppi"          , len(_eta)-1 , array('d' , _eta))]
h_sigma_pll_eta      =  [ ROOT.TH1F("h_sigma_pll_eta_pf"      , "h_sigma_pll_eta_pf"      , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_sigma_pll_eta_puppi"      , "h_sigma_pll_eta_puppi"      , len(_eta)-1 , array('d' , _eta))]
h_sigma_prp_eta      =  [ ROOT.TH1F("h_sigma_prp_eta_pf"      , "h_sigma_prp_eta_pf"      , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_sigma_prp_eta_puppi"      , "h_sigma_prp_eta_puppi"      , len(_eta)-1 , array('d' , _eta))]
h_resolution_pll_eta =  [ ROOT.TH1F("h_resolution_pll_eta_pf" , "h_resolution_pll_eta_pf" , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_resolution_pll_eta_puppi" , "h_resolution_pll_eta_puppi" , len(_eta)-1 , array('d' , _eta))]
h_resolution_prp_eta =  [ ROOT.TH1F("h_resolution_prp_eta_pf" , "h_resolution_prp_eta_pf" , len(_eta)-1 , array('d' , _eta)) , ROOT.TH1F("h_resolution_prp_eta_puppi" , "h_resolution_prp_eta_puppi" , len(_eta)-1 , array('d' , _eta))]
h_response_npv       =  [ ROOT.TH1F("h_response_npv_pf"          , "h_response_npv_pf"          , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_response_npv_puppi"          , "h_response_npv_puppi"          , len(_npv)-1 , array('d' , _npv))]
h_sigma_pll_npv      =  [ ROOT.TH1F("h_sigma_pll_npv_pf"      , "h_sigma_pll_npv_pf"      , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_sigma_pll_npv_puppi"      , "h_sigma_pll_npv_puppi"      , len(_npv)-1 , array('d' , _npv))]
h_sigma_prp_npv      =  [ ROOT.TH1F("h_sigma_prp_npv_pf"      , "h_sigma_prp_npv_pf"      , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_sigma_prp_npv_puppi"      , "h_sigma_prp_npv_puppi"      , len(_npv)-1 , array('d' , _npv))]
h_resolution_pll_npv =  [ ROOT.TH1F("h_resolution_pll_npv_pf" , "h_resolution_pll_npv_pf" , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_resolution_pll_npv_puppi" , "h_resolution_pll_npv_puppi" , len(_npv)-1 , array('d' , _npv))]
h_resolution_prp_npv =  [ ROOT.TH1F("h_resolution_prp_npv_pf" , "h_resolution_prp_npv_pf" , len(_npv)-1 , array('d' , _npv)) , ROOT.TH1F("h_resolution_prp_npv_puppi" , "h_resolution_prp_npv_puppi" , len(_npv)-1 , array('d' , _npv))]




#Make pf vs. puppi plots
ROOT.gStyle.SetOptStat(0) #tdrStyle not in effect?
def plot_hists(hist_list, title="", legend_title_list=None, x_title="", y_title="", text_description=None, plot_name=None, limits=None):
	colors=[ROOT.kCyan+1,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kRed+1,ROOT.kOrange,ROOT.kYellow+1,ROOT.kGreen+1,ROOT.kGray]
	canv = ROOT.TCanvas("canv","canv")
	mg=ROOT.TMultiGraph() #Use a multiGraph to auto adjust the frame size:
	base_hist = hist_list[0]
	base_hist.SetTitle(title)
	base_hist.GetXaxis().SetTitle(x_title)
	base_hist.GetYaxis().SetTitle(y_title)
	base_hist.SetLineWidth(2)
	base_hist.SetLineColor(colors[0])
	base_hist.SetMarkerStyle(1)
	base_hist.SetMarkerSize(0.8)
	mg.Add(ROOT.TGraphErrors(base_hist))
	if len(hist_list)>1:
		for ihist,hist in enumerate(hist_list[1:]):
			hist.SetMarkerStyle(1)
			hist.SetMarkerSize(0.8)
			hist.SetLineColor(colors[ihist+1])
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
		



#call function get_met_params and set the histograms' bin content
for ipt,pt in enumerate(_pt[:-1]):
	params=get_met_params("qt>"+str(_pt[ipt])+" && qt<"+str(_pt[ipt+1]))
	if params == None: continue
	h_response_pt[0].SetBinContent(ipt+1,params["response_pf"]);
	h_response_pt[0].SetBinError(ipt+1,params["response_pf_error"]);
	h_sigma_pll_pt[0].SetBinContent(ipt+1,params["sigma_pll_pf"]);
	h_sigma_pll_pt[0].SetBinError(ipt+1,params["sigma_pll_pf_error"]);
	h_resolution_pll_pt[0].SetBinContent(ipt+1,params["resolution_pll_pf"]);
	h_resolution_pll_pt[0].SetBinError(ipt+1,params["resolution_pll_pf_error"]);
	h_sigma_prp_pt[0].SetBinContent(ipt+1,params["sigma_prp_pf"]);
	h_sigma_prp_pt[0].SetBinError(ipt+1,params["sigma_prp_pf_error"]);
	h_resolution_prp_pt[0].SetBinContent(ipt+1,params["resolution_prp_pf"]);
	h_resolution_prp_pt[0].SetBinError(ipt+1,params["resolution_prp_pf_error"]);
	h_response_pt[1].SetBinContent(ipt+1,params["response_puppi"]);
	h_response_pt[1].SetBinError(ipt+1,params["response_puppi_error"]);
	h_sigma_pll_pt[1].SetBinContent(ipt+1,params["sigma_pll_puppi"]);
	h_sigma_pll_pt[1].SetBinError(ipt+1,params["sigma_pll_puppi_error"]);
	h_resolution_pll_pt[1].SetBinContent(ipt+1,params["resolution_pll_puppi"]);
	h_resolution_pll_pt[1].SetBinError(ipt+1,params["resolution_pll_puppi_error"]);
	h_sigma_prp_pt[1].SetBinContent(ipt+1,params["sigma_prp_puppi"]);
	h_sigma_prp_pt[1].SetBinError(ipt+1,params["sigma_prp_puppi_error"]);
	h_resolution_prp_pt[1].SetBinContent(ipt+1,params["resolution_prp_puppi"]);
	h_resolution_prp_pt[1].SetBinError(ipt+1,params["resolution_prp_puppi_error"]);
for ieta,eta in enumerate(_eta[:-1]):
	params=get_met_params("q_eta>"+ "%.1f"%_eta[ieta] +" && q_eta<"+"%.1f"%_eta[ieta+1])
	if params == None: continue
	h_response_eta[0].SetBinContent(ieta+1,params["response_pf"]);
	h_response_eta[0].SetBinError(ieta+1,params["response_pf_error"]);
	h_sigma_pll_eta[0].SetBinContent(ieta+1,params["sigma_pll_pf"]);
	h_sigma_pll_eta[0].SetBinError(ieta+1,params["sigma_pll_pf_error"]);
	h_resolution_pll_eta[0].SetBinContent(ieta+1,params["resolution_pll_pf"]);
	h_resolution_pll_eta[0].SetBinError(ieta+1,params["resolution_pll_pf_error"]);
	h_sigma_prp_eta[0].SetBinContent(ieta+1,params["sigma_prp_pf"]);
	h_sigma_prp_eta[0].SetBinError(ieta+1,params["sigma_prp_pf_error"]);
	h_resolution_prp_eta[0].SetBinContent(ieta+1,params["resolution_prp_pf"]);
	h_resolution_prp_eta[0].SetBinError(ieta+1,params["resolution_prp_pf_error"]);
	h_response_eta[1].SetBinContent(ieta+1,params["response_puppi"]);
	h_response_eta[1].SetBinError(ieta+1,params["response_puppi_error"]);
	h_sigma_pll_eta[1].SetBinContent(ieta+1,params["sigma_pll_puppi"]);
	h_sigma_pll_eta[1].SetBinError(ieta+1,params["sigma_pll_puppi_error"]);
	h_resolution_pll_eta[1].SetBinContent(ieta+1,params["resolution_pll_puppi"]);
	h_resolution_pll_eta[1].SetBinError(ieta+1,params["resolution_pll_puppi_error"]);
	h_sigma_prp_eta[1].SetBinContent(ieta+1,params["sigma_prp_puppi"]);
	h_sigma_prp_eta[1].SetBinError(ieta+1,params["sigma_prp_puppi_error"]);
	h_resolution_prp_eta[1].SetBinContent(ieta+1,params["resolution_prp_puppi"]);
	h_resolution_prp_eta[1].SetBinError(ieta+1,params["resolution_prp_puppi_error"]);
for inpv,npv in enumerate(_npv[:-1]):
	params=get_met_params("npv>"+str(_npv[inpv])+" && npv<"+str(_npv[inpv+1]))
	if params == None: continue
	h_response_npv[0].SetBinContent(inpv+1,params["response_pf"]);
	h_response_npv[0].SetBinError(inpv+1,params["response_pf_error"]);
	h_sigma_pll_npv[0].SetBinContent(inpv+1,params["sigma_pll_pf"]);
	h_sigma_pll_npv[0].SetBinError(inpv+1,params["sigma_pll_pf_error"]);
	h_resolution_pll_npv[0].SetBinContent(inpv+1,params["resolution_pll_pf"]);
	h_resolution_pll_npv[0].SetBinError(inpv+1,params["resolution_pll_pf_error"]);
	h_sigma_prp_npv[0].SetBinContent(inpv+1,params["sigma_prp_pf"]);
	h_sigma_prp_npv[0].SetBinError(inpv+1,params["sigma_prp_pf_error"]);
	h_resolution_prp_npv[0].SetBinContent(inpv+1,params["resolution_prp_pf"]);
	h_resolution_prp_npv[0].SetBinError(inpv+1,params["resolution_prp_pf_error"]);
	h_response_npv[1].SetBinContent(inpv+1,params["response_puppi"]);
	h_response_npv[1].SetBinError(inpv+1,params["response_puppi_error"]);
	h_sigma_pll_npv[1].SetBinContent(inpv+1,params["sigma_pll_puppi"]);
	h_sigma_pll_npv[1].SetBinError(inpv+1,params["sigma_pll_puppi_error"]);
	h_resolution_pll_npv[1].SetBinContent(inpv+1,params["resolution_pll_puppi"]);
	h_resolution_pll_npv[1].SetBinError(inpv+1,params["resolution_pll_puppi_error"]);
	h_sigma_prp_npv[1].SetBinContent(inpv+1,params["sigma_prp_puppi"]);
	h_sigma_prp_npv[1].SetBinError(inpv+1,params["sigma_prp_puppi_error"]);
	h_resolution_prp_npv[1].SetBinContent(inpv+1,params["resolution_prp_puppi"]);
	h_resolution_prp_npv[1].SetBinError(inpv+1,params["resolution_prp_puppi_error"]);




plot_hists(h_response_pt, legend_title_list=["PF MET", "PUPPI MET"], x_title="q_{T} [GeV]", y_title="-<u_{||}>/<q_{T}>", plot_name="met_response_pt", limits=[0,2])
plot_hists(h_sigma_pll_pt, legend_title_list=["PF MET", "PUPPI MET"], x_title="q_{T} [GeV]", y_title="#sigma (u_{||}) [GeV]", plot_name="met_sigma_pll_pt")
plot_hists(h_sigma_prp_pt, legend_title_list=["PF MET", "PUPPI MET"], x_title="q_{T} [GeV]", y_title="#sigma (u_{#perp }) [GeV]", plot_name="met_sigma_prp_pt")
plot_hists(h_resolution_pll_pt, legend_title_list=["PF MET", "PUPPI MET"], x_title="q_{T} [GeV]", y_title="#sigma (u_{||})/response [GeV]", plot_name="met_resolution_pll_pt", limits=[0,70])
plot_hists(h_resolution_prp_pt, legend_title_list=["PF MET", "PUPPI MET"], x_title="q_{T} [GeV]", y_title="#sigma (u_{#perp })/response [GeV]", plot_name="met_resolution_prp_pt", limits=[0,70])
plot_hists(h_response_eta, legend_title_list=["PF MET", "PUPPI MET"], x_title="#eta", y_title="-<u_{||}>/<q_{T}>", plot_name="met_response_eta", limits=[0,2])
plot_hists(h_sigma_pll_eta, legend_title_list=["PF MET", "PUPPI MET"], x_title="#eta", y_title="#sigma (u_{||}) [GeV]", plot_name="met_sigma_pll_eta")
plot_hists(h_sigma_prp_eta, legend_title_list=["PF MET", "PUPPI MET"], x_title="#eta", y_title="#sigma (u_{#perp }) [GeV]", plot_name="met_sigma_prp_eta")
plot_hists(h_resolution_pll_eta, legend_title_list=["PF MET", "PUPPI MET"], x_title="#eta", y_title="#sigma (u_{||})/response [GeV]", plot_name="met_resolution_pll_eta", limits=[0,500])
plot_hists(h_resolution_prp_eta, legend_title_list=["PF MET", "PUPPI MET"], x_title="#eta", y_title="#sigma (u_{#perp })/response [GeV]", plot_name="met_resolution_prp_eta", limits=[0,600])
plot_hists(h_response_npv, legend_title_list=["PF MET", "PUPPI MET"], x_title="Number of vertices", y_title="-<u_{||}>/<q_{T}>", plot_name="met_response_npv", limits=[0,2])
plot_hists(h_sigma_pll_npv, legend_title_list=["PF MET", "PUPPI MET"], x_title="Number of vertices", y_title="#sigma (u_{||}) [GeV]", plot_name="met_sigma_pll_npv")
plot_hists(h_sigma_prp_npv, legend_title_list=["PF MET", "PUPPI MET"], x_title="Number of vertices", y_title="#sigma (u_{#perp }) [GeV]", plot_name="met_sigma_prp_npv")
plot_hists(h_resolution_pll_npv, legend_title_list=["PF MET", "PUPPI MET"], x_title="Number of vertices", y_title="#sigma (u_{||})/response [GeV]", plot_name="met_resolution_pll_npv", limits=[0,130])
plot_hists(h_resolution_prp_npv, legend_title_list=["PF MET", "PUPPI MET"], x_title="Number of vertices", y_title="#sigma (u_{#perp })/response [GeV]", plot_name="met_resolution_prp_npv", limits=[0,130])

#Make some plots directly from TTree:
c1=ROOT.TCanvas("c1","c1")
#Make a pt plot of two muons
t_in.Draw("mu1_pt")
c1.GetPrimitive("htemp").SetLineColor(ROOT.kRed)
t_in.Draw("mu2_pt","","same")
c1.Print("result/two_muon_pt.png")
#Make a eta plot of two muons
t_in.Draw("mu1_eta")
c1.GetPrimitive("htemp").SetLineColor(ROOT.kRed)
t_in.Draw("mu2_eta","","same")
c1.Print("result/two_muon_eta.png")
#Make a phi plot of two muons
t_in.Draw("mu1_phi")
c1.GetPrimitive("htemp").SetLineColor(ROOT.kRed)
t_in.Draw("mu2_phi","","same")
c1.Print("result/two_muon_phi.png")





