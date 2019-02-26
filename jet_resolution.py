import os,sys,socket,argparse
import re
import ROOT
import math
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


parser = argparse.ArgumentParser(description='Plot resolution and response')
parser.add_argument("-i", "--input" ,dest="input", help="input file name", type=str)
parser.add_argument("-o", "--output" ,dest="output" , help="output folder name", type=str)
parser.add_argument("-p", "--putype" ,dest="pu" ,help="define the sample, options are 'flatpu','nopu'", type=str)
parser.add_argument("-e", "--eta" ,dest="eta" ,help="define the eta range, options are '0_1p3','1p3_2p1','2p1_2p5','2p5_3p0','3p0_5p0'",type=str)
parser.add_argument("-r", "--rho" ,dest="rho", help="define the rho range", type=str)
args = parser.parse_args()

print args.pu, args.eta
f_in = ROOT.TFile(str(args.input),"READ")
t_in = f_in.Get("events")
_pt = [20,50,100,200,500,2000]
_eta = [float(args.eta.split("_")[0].replace("p",".")), float(args.eta.split("_")[1].replace("p","."))]
_rho = [float(args.rho.split("_")[0]), float(args.rho.split("_")[1])]
folder = args.output+"/"+args.pu+"/"
os.system("mkdir -p "+folder)
os.system("mkdir -p "+folder+"/fit")

h_mean = [ROOT.TH1F("h_mean_pf","h_mean_pf",len(_pt)-1,array('d',_pt)),
	ROOT.TH1F("h_mean_puppi","h_mean_puppi",len(_pt)-1,array('d',_pt))]
h_sigma = [ROOT.TH1F("h_sigma_pf","h_sigma_pf",len(_pt)-1,array('d',_pt)),
	ROOT.TH1F("h_sigma_puppi","h_sigma_puppi",len(_pt)-1,array('d',_pt))] #Absolute sigma
h_sigmal1 = [ROOT.TH1F("h_sigmal1_pf","h_sigmal1_pf",len(_pt)-1,array('d',_pt)),
	ROOT.TH1F("h_sigmal1_puppi","h_sigmal1_puppi",len(_pt)-1,array('d',_pt))] #sigma/response

def get_params(numerator,denominator,ipt,index=0):
	assert index in range(len(h_mean))
	shape = ROOT.TH1F("shape","shape",50,0.25,2.5)
	cut_arg=denominator+">"+str(_pt[ipt])+" && "+denominator+"<"+str(_pt[ipt+1])+" && abs(genjet_eta)>="+str(_eta[0])+" && abs(genjet_eta)<"+str(_eta[1])+rho_cut
	print("cut="+cut_arg)
	t_in.Draw(numerator+"/"+denominator+">>shape",cut_arg,"goff")
	mean,mean_error,sigma,sigma_error = ConvFit(shape ,False,"ratio",numerator,folder+"/fit","FIT_"+numerator+"_"+str(_pt[ipt])+"_"+str(_pt[ipt+1])+"_eta_"+str(_eta[0])+"_"+str(_eta[1])+"_rho_"+str(_rho[0])+"_"+str(_rho[1]))
	h_mean[index].SetBinContent(ipt+1,mean)
	h_mean[index].SetBinError(ipt+1,mean_error)
	#scale correct the resolution
	if mean>0:
		h_sigma[index].SetBinContent(ipt+1,max(0,sigma/mean))
	else:
		h_sigma[index].SetBinContent(ipt+1,0)
	h_sigma[index].SetBinError(ipt+1,sigma_error)
	h_sigmal1[index].SetBinContent(ipt+1,sigma_error)
	h_sigmal1[index].SetBinError(ipt+1,sigma_error)


rho_cut = " && rhoall>="+str(_rho[0])+" && rhoall<"+str(_rho[1])
for i in range(len(_pt)-1):
	get_params("rawjet_pt","genjet_pt",i,0) #pf jets
	get_params("rawpjet_pt","genpjet_pt",i,1) #puppi jets


#Make plot:
c = ROOT.TCanvas("sigma","sigma")
c.cd()
c.SetLogx()
h_sigma[0].GetXaxis().SetTitle("Gen jet p_{T} [GeV]")
h_sigma[0].GetYaxis().SetTitle("Resolution / Response")
h_sigma[0].SetTitle("")
h_sigma[0].GetXaxis().SetTitleOffset(1.2)
h_sigma[0].GetYaxis().SetTitleOffset(1.3)
h_sigma[0].SetLineWidth(2)
h_sigma[0].SetLineColor(ROOT.kLake)
h_sigma[0].SetMaximum(0.5)
h_sigma[0].SetMinimum(0)
ROOT.gStyle.SetOptStat(0) #tdrStyle not in effect?
h_sigma[0].Draw()
h_sigma[0].SetMarkerStyle(1)
h_sigma[0].SetMarkerSize(0.8)
h_sigma[1].Draw("same")
h_sigma[1].SetMarkerStyle(1)
h_sigma[1].SetMarkerSize(0.8)
h_sigma[1].SetLineColor(ROOT.kRose)

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.4*c.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right
latex2.DrawLatex(0.90, 0.93,str(_eta[0])+" < #eta <"+ str(_eta[1]) + ", "+ str(_rho[0])+" < #rho < " + str(_rho[1]))
latex2.Draw("same")

legend=ROOT.TLegend(0.60,0.65,0.95,.85)
legend.AddEntry(h_sigma[0],"pf+CHS","lp")
legend.AddEntry(h_sigma[1],"puppi","lp")
legend.Draw("same")
c.SaveAs(folder+"FIT_sigma_eta_"+args.eta+"_rho_"+args.rho+".png")
c.SaveAs(folder+"FIT_sigma_eta_"+args.eta+"_rho_"+args.rho+".pdf")
