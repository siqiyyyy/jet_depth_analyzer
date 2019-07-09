# import ROOT in batch mode
import sys
import time
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv
from array import array

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

#default options
options.inputFiles="../step3_inMINIAODSIM.root"
options.outputFile="jetmetNtuples.root"
options.maxEvents=-1

#overwrite if given any command line arguments
options.parseArguments()
#in case of txt input file, read the information from txt
li_f=[]
for iF,F in enumerate(options.inputFiles):
	print F
	if F.split('.')[-1] == "txt":
		options.inputFiles.pop(iF)
		with open(F) as f:
			li_f+=f.read().splitlines()
options.inputFiles=li_f
print "analyzing files:"
for F in options.inputFiles: print F

# define deltaR
from math import hypot, pi, sqrt, fabs
import numpy as n

#from functions import *

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

pfs, pfLabel		= Handle("std::vector<pat::PackedCandidate>"), "packedPFCandidates"
jets, jetLabel	  = Handle("std::vector<pat::Jet>"), "slimmedJets"
pjets, pjetLabel	  = Handle("std::vector<pat::Jet>"), "slimmedJetsPuppi"
genjets, genjetLabel	  = Handle("std::vector<reco::GenJet>"), "slimmedGenJets"
mets, metLabel	  = Handle("std::vector<pat::MET>"), "slimmedMETs"
pmets, pmetLabel	  = Handle("std::vector<pat::MET>"), "slimmedMETsPuppi"
vertex, vertexLabel = Handle("std::vector<reco::Vertex>"),"offlineSlimmedPrimaryVertices"

rhoall_, rhoallLabel		 = Handle("double"), "fixedGridRhoFastjetAll"
rhocentral_, rhocentralLabel = Handle("double"), "fixedGridRhoFastjetCentral"
rhoneutral_, rhoneutralLabel = Handle("double"), "fixedGridRhoFastjetCentralNeutral"
rhochargedpileup_, rhochargedpileupLabel = Handle("double"), "fixedGridRhoFastjetCentralChargedPileUp"

#in order to print out the progress
def print_same_line(s):
	sys.stdout.write(s)					 # just print
	sys.stdout.flush()					  # needed for flush when using \x08
	sys.stdout.write((b'\x08' * len(s)).decode())# back n chars		
	#time.sleep(0.2)

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events("file:/eos/cms/store/relval/CMSSW_9_4_0_pre3/RelValTTbar_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_PixFailScenario_Run305081_FIXED_HS_AVE50-v1/10000/02B605A1-86C2-E711-A445-4C79BA28012B.root")
events = Events(options)
nevents = int(events.size())
print "total events: ", events.size()

outfile=ROOT.TFile("hists.root","recreate")

ROOT.TH1.SetDefaultSumw2()
#Init histograms in npv,eta,pt
_eta = n.arange(-5,5,0.2)
eff_eta = [] #efficiency
eff_eta.append(ROOT.TH1F("eff_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta)))
eff_eta.append(ROOT.TH1F("eff_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta)))
prt_eta = [] #Purity
prt_eta.append(ROOT.TH1F("prt_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta)))
prt_eta.append(ROOT.TH1F("prt_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta)))
_npv = n.arange(0,70,2)
eff_npv = [] #efficiency
eff_npv.append(ROOT.TH1F("eff_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv)))
eff_npv.append(ROOT.TH1F("eff_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv)))
prt_npv = [] #Purity
prt_npv.append(ROOT.TH1F("prt_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv)))
prt_npv.append(ROOT.TH1F("prt_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv)))
_pt = n.arange(30,250,5)
eff_pt = [] #efficiency
eff_pt.append(ROOT.TH1F("eff_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt)))
eff_pt.append(ROOT.TH1F("eff_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt)))
prt_pt = [] #Purity
prt_pt.append(ROOT.TH1F("prt_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt)))
prt_pt.append(ROOT.TH1F("prt_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt)))

h_chsjet_pt=ROOT.TH1F("recojets_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt))
h_chsjet_matched_pt=ROOT.TH1F("matchedrecojets_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt))
h_chsjet_eta=ROOT.TH1F("recojets_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta))
h_chsjet_matched_eta=ROOT.TH1F("matchedrecojets_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta))
h_chsjet_npv=ROOT.TH1F("recojets_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv))
h_chsjet_matched_npv=ROOT.TH1F("matchedrecojets_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv))
h_puppijet_pt=ROOT.TH1F("recojets_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt))
h_puppijet_matched_pt=ROOT.TH1F("matchedrecojets_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt))
h_puppijet_eta=ROOT.TH1F("recojets_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta))
h_puppijet_matched_eta=ROOT.TH1F("matchedrecojets_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta))
h_puppijet_npv=ROOT.TH1F("recojets_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv))
h_puppijet_matched_npv=ROOT.TH1F("matchedrecojets_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv))
h_gen_chsjet_pt=ROOT.TH1F("genjets_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt))
h_gen_chsjet_matched_pt=ROOT.TH1F("matchedgenjets_pt_CHS","pf+CHS",len(_pt)-1,array('d',_pt))
h_gen_chsjet_eta=ROOT.TH1F("genjets_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta))
h_gen_chsjet_matched_eta=ROOT.TH1F("matchedgenjets_eta_CHS","pf+CHS",len(_eta)-1,array('d',_eta))
h_gen_chsjet_npv=ROOT.TH1F("genjets_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv))
h_gen_chsjet_matched_npv=ROOT.TH1F("matchedgenjets_npv_CHS","pf+CHS",len(_npv)-1,array('d',_npv))
h_gen_puppijet_pt=ROOT.TH1F("genjets_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt))
h_gen_puppijet_matched_pt=ROOT.TH1F("matchedgenjets_pt_PUPPI","PUPPI",len(_pt)-1,array('d',_pt))
h_gen_puppijet_eta=ROOT.TH1F("genjets_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta))
h_gen_puppijet_matched_eta=ROOT.TH1F("matchedgenjets_eta_PUPPI","PUPPI",len(_eta)-1,array('d',_eta))
h_gen_puppijet_npv=ROOT.TH1F("genjets_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv))
h_gen_puppijet_matched_npv=ROOT.TH1F("matchedgenjets_npv_PUPPI","PUPPI",len(_npv)-1,array('d',_npv))

maxjet=1000
for ievent,event in enumerate(events):

	if options.maxEvents is not -1:
		if ievent > options.maxEvents: continue
	
	event.getByLabel(pfLabel, pfs)
	event.getByLabel(jetLabel, jets)
	event.getByLabel(pjetLabel, pjets)
	event.getByLabel(genjetLabel, genjets)
	event.getByLabel(metLabel, mets)
	event.getByLabel(pmetLabel, pmets)
	event.getByLabel(vertexLabel, vertex)

	event.getByLabel(rhoallLabel,rhoall_)
	event.getByLabel(rhocentralLabel,rhocentral_)
	event.getByLabel(rhoneutralLabel,rhoneutral_)
	event.getByLabel(rhochargedpileupLabel,rhochargedpileup_)

	#print "\nEvent %d: run %6d, lumi %4d, event %12d" % (ievent,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
	print_same_line(str(round(100.*ievent/nevents,2))+'%')

	###CHS JETS
	matched_gen_jets=[]
	matched_rec_jets=[]
	for i,j in enumerate(jets.product()):
		if i>=maxjet: break
		if j.pt()<30: continue
		h_chsjet_pt.Fill(j.pt())
		h_chsjet_eta.Fill(j.eta())
		h_chsjet_npv.Fill(vertex.product().size())
		if not (j.genJet() == None):
			matched_gen_jets.append(j.genJet())
			matched_rec_jets.append(j)
			if j.genJet().pt()<15: continue
			h_chsjet_matched_pt.Fill(j.pt())
			h_chsjet_matched_eta.Fill(j.eta())
			h_chsjet_matched_npv.Fill(vertex.product().size())
	for i,j in enumerate(genjets.product()):
		if i>=maxjet: break
		if j.pt()<30: continue
		h_gen_chsjet_pt.Fill(j.pt())
		h_gen_chsjet_eta.Fill(j.eta())
		h_gen_chsjet_npv.Fill(vertex.product().size())
		if j in matched_gen_jets:
			if matched_rec_jets[matched_gen_jets.index(j)].pt()<20: continue
			h_gen_chsjet_matched_pt.Fill(j.pt())
			h_gen_chsjet_matched_eta.Fill(j.eta())
			h_gen_chsjet_matched_npv.Fill(vertex.product().size())
		else:
			continue


	###Puppi JETS
	matched_gen_jets=[]
	matched_rec_jets=[]
	for i,j in enumerate(pjets.product()):
		if i>=maxjet: break
		if j.pt()<30: continue
		h_puppijet_pt.Fill(j.pt())
		h_puppijet_eta.Fill(j.eta())
		h_puppijet_npv.Fill(vertex.product().size())
		if not (j.genJet() == None):
			matched_gen_jets.append(j.genJet())
			matched_rec_jets.append(j)
			if j.genJet().pt()<15: continue
			h_puppijet_matched_pt.Fill(j.pt())
			h_puppijet_matched_eta.Fill(j.eta())
			h_puppijet_matched_npv.Fill(vertex.product().size())
	for i,j in enumerate(genjets.product()):
		if i>=maxjet: break
		if j.pt()<30: continue
		h_gen_puppijet_pt.Fill(j.pt())
		h_gen_puppijet_eta.Fill(j.eta())
		h_gen_puppijet_npv.Fill(vertex.product().size())
		if j in matched_gen_jets:
			if matched_rec_jets[matched_gen_jets.index(j)].pt()<20: continue
			h_gen_puppijet_matched_pt.Fill(j.pt())
			h_gen_puppijet_matched_eta.Fill(j.eta())
			h_gen_puppijet_matched_npv.Fill(vertex.product().size())
		else:
			continue

prt_pt[0].Divide(h_chsjet_matched_pt, h_chsjet_pt)
prt_eta[0].Divide(h_chsjet_matched_eta, h_chsjet_eta)
prt_npv[0].Divide(h_chsjet_matched_npv, h_chsjet_npv)
prt_pt[1].Divide(h_puppijet_matched_pt, h_puppijet_pt)
prt_eta[1].Divide(h_puppijet_matched_eta, h_puppijet_eta)
prt_npv[1].Divide(h_puppijet_matched_npv, h_puppijet_npv)
eff_pt[0].Divide(h_gen_chsjet_matched_pt, h_gen_chsjet_pt)
eff_eta[0].Divide(h_gen_chsjet_matched_eta, h_gen_chsjet_eta)
eff_npv[0].Divide(h_gen_chsjet_matched_npv, h_gen_chsjet_npv)
eff_pt[1].Divide(h_gen_puppijet_matched_pt, h_gen_puppijet_pt)
eff_eta[1].Divide(h_gen_puppijet_matched_eta, h_gen_puppijet_eta)
eff_npv[1].Divide(h_gen_puppijet_matched_npv, h_gen_puppijet_npv)


ROOT.gStyle.SetOptStat(0)
#Plot the histograms
def plot(hists,filename,xvarname,yvarname):
	c = ROOT.TCanvas("c","c",800,800) #who cares about your name?
	colors=[ROOT.kCyan+1,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kRed+1,ROOT.kOrange,ROOT.kYellow+1,ROOT.kGreen+1,ROOT.kGray] 
	legend=ROOT.TLegend(0.70,0.85,1,0.93)
	binedge=array('d',n.zeros(hists[0].GetNbinsX()))
	hists[0].GetXaxis().GetLowEdge(binedge)
	hframe=c.DrawFrame(binedge[0]-hists[0].GetBinWidth(0),0,binedge[-1]+2*hists[0].GetBinWidth(hists[0].GetNbinsX()),1.1,"")
	hframe.GetXaxis().SetTitle(xvarname)
	hframe.GetXaxis().SetTitleOffset(1.2)
	hframe.GetYaxis().SetTitle("Jet "+yvarname)
	hframe.GetYaxis().SetTitleOffset(1.3)
	for ihist,hist in enumerate(hists):
		hist.SetLineColor(colors[ihist])
		hist.Draw("same")
		legend.AddEntry(hist, hist.GetTitle(),"lp")
	legend.Draw("same")
	c.Print("result/"+filename)

outfile.cd()
prt_pt[0].Write()
prt_pt[1].Write()
prt_eta[0].Write()
prt_eta[1].Write()
prt_npv[0].Write()
prt_npv[1].Write()
eff_pt[0].Write()
eff_pt[1].Write()
eff_eta[0].Write()
eff_eta[1].Write()
eff_npv[0].Write()
eff_npv[1].Write()
plot(prt_pt, "purity_pt.png", "p_T/GeV","purity")
plot(prt_eta, "purity_eta.png", "#eta","purity")
plot(prt_npv, "purity_npv.png", "npv","purity")
plot(eff_pt  , "efficiency_pt.png"  , "p_T/GeV" , "efficiency")
plot(eff_eta , "efficiency_eta.png" , "#eta"    , "efficiency")
plot(eff_npv , "efficiency_npv.png" , "npv"     , "efficiency")
