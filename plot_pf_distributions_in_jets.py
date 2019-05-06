# import ROOT in batch mode
import sys
import time
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

# RooFit
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
ROOT.gROOT.SetStyle("Plain")
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )

#default options
#options.inputFiles="../step3_inMINIAODSIM.root"
options.inputFiles="files.txt"
options.outputFile="result/pf_eta_in_jets.root"
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
options.inputFiles+=li_f
print "analyzing files:"
for F in options.inputFiles: print F

# define deltaR
from math import hypot, pi, sqrt, fabs
import numpy as n

from jetmet_tree import *
from functions import *

# create an oput tree.

fout = ROOT.TFile (options.outputFile,"recreate")
#t    = ROOT.TTree ("events","events")

#declare_branches(t)

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

pfs, pfLabel        = Handle("std::vector<pat::PackedCandidate>"), "packedPFCandidates"
jets, jetLabel      = Handle("std::vector<pat::Jet>"), "slimmedJets"
pjets, pjetLabel      = Handle("std::vector<pat::Jet>"), "slimmedJetsPuppi"
mets, metLabel      = Handle("std::vector<pat::MET>"), "slimmedMETs"
pmets, pmetLabel      = Handle("std::vector<pat::MET>"), "slimmedMETsPuppi"
vertex, vertexLabel = Handle("std::vector<reco::Vertex>"),"offlineSlimmedPrimaryVertices"

rhoall_, rhoallLabel         = Handle("double"), "fixedGridRhoFastjetAll"
rhocentral_, rhocentralLabel = Handle("double"), "fixedGridRhoFastjetCentral"
rhoneutral_, rhoneutralLabel = Handle("double"), "fixedGridRhoFastjetCentralNeutral"
rhochargedpileup_, rhochargedpileupLabel = Handle("double"), "fixedGridRhoFastjetCentralChargedPileUp"

#in order to print out the progress
def print_same_line(s):
	sys.stdout.write(s)                     # just print
	sys.stdout.flush()                      # needed for flush when using \x08
	sys.stdout.write((b'\x08' * len(s)).decode())# back n chars    	
	#time.sleep(0.2)

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events("file:/eos/cms/store/relval/CMSSW_9_4_0_pre3/RelValTTbar_13/MINIAODSIM/PU25ns_94X_mc2017_realistic_PixFailScenario_Run305081_FIXED_HS_AVE50-v1/10000/02B605A1-86C2-E711-A445-4C79BA28012B.root")
events = Events(options)
nevents = int(events.size())
print "total events: ", events.size()
num_jet=0
max_num_jet=5
v_eta = ROOT.RooRealVar("eta","eta",0,5.0)

fout.cd()
for ievent,event in enumerate(events):

    if options.maxEvents is not -1:
        if ievent > options.maxEvents: continue
    
    event.getByLabel(pfLabel, pfs)
    event.getByLabel(jetLabel, jets)
    event.getByLabel(pjetLabel, pjets)
    event.getByLabel(metLabel, mets)
    event.getByLabel(pmetLabel, pmets)
    event.getByLabel(vertexLabel, vertex)

    event.getByLabel(rhoallLabel,rhoall_)
    event.getByLabel(rhocentralLabel,rhocentral_)
    event.getByLabel(rhoneutralLabel,rhoneutral_)
    event.getByLabel(rhochargedpileupLabel,rhochargedpileup_)

    rhoall[0]     = rhoall_.product()[0]
    rhocentral[0] = rhocentral_.product()[0]
    rhoneutral[0] = rhoneutral_.product()[0]
    rhochargedpileup[0] = rhochargedpileup_.product()[0]

    #print "\nEvent %d: run %6d, lumi %4d, event %12d" % (ievent,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    print_same_line(str(round(100.*num_jet/max_num_jet,2))+'%')
    num_jet+=1
    if num_jet>max_num_jet: break

    nrun[0]   = event.eventAuxiliary().run()
    nlumi[0]  = event.eventAuxiliary().luminosityBlock()
    nevent[0] = event.eventAuxiliary().event()
    npv[0]    = vertex.product().size()

    njet[0]   = jets.product().size()
    npjet[0]   = pjets.product().size()

    ###CHS JETS
    for i,j in enumerate(jets.product()):

        if i>=maxjet: break
        if num_jet>max_num_jet: break

        sourceCandidate = set()

        # now get a list of the PF candidates used to build this jet
        for isource in xrange(j.numberOfSourceCandidatePtrs()):
            sourceCandidate.add(j.sourceCandidatePtr(isource).key()) # the key is the index in the pf collection

        tmp_eta_set=ROOT.RooDataSet("jet_eta_"+str(j.eta()),"jet_eta_"+str(j.eta()),ROOT.RooArgSet(v_eta))
        for ipf,pf in enumerate(pfs.product()):            
            if ipf in sourceCandidate:
                v_eta.setVal(pf.eta())
                tmp_eta_set.add(ROOT.RooArgSet(v_eta))
        tmp_eta_set.Write()


    ###Puppi JETS
    num_jet=0
    for i,j in enumerate(pjets.product()):
	

        if i>=maxjet: break
        if num_jet>max_num_jet: break

        sourceCandidate = set()

        # now get a list of the PF candidates used to build this jet
        for isource in xrange(j.numberOfSourceCandidatePtrs()):
            sourceCandidate.add(j.sourceCandidatePtr(isource).key()) # the key is the index in the pf collection

        tmp_eta_set=ROOT.RooDataSet("pjet_eta_"+str(j.eta()),"pjet_eta_"+str(j.eta()),ROOT.RooArgSet(v_eta))
        for ipf,pf in enumerate(pfs.product()):            
            if ipf in sourceCandidate:
                v_eta.setVal(pf.eta())
                tmp_eta_set.add(ROOT.RooArgSet(v_eta),pf.puppiWeight())
        tmp_eta_set.Write()




fout.Close()
