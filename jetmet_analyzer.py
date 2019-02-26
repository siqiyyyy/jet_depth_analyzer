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
options.inputFiles+=li_f
print "analyzing files:"
for F in options.inputFiles: print F

# define deltaR
from math import hypot, pi, sqrt, fabs
import numpy as n

from jetmet_tree import *
from functions import *

# create an oput tree.

#fout = ROOT.TFile ("jetmet.root","recreate")
fout = ROOT.TFile (options.outputFile,"recreate")
t    = ROOT.TTree ("events","events")

declare_branches(t)

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
    print_same_line(str(round(100.*ievent/nevents,2))+'%')

    nrun[0]   = event.eventAuxiliary().run()
    nlumi[0]  = event.eventAuxiliary().luminosityBlock()
    nevent[0] = event.eventAuxiliary().event()
    npv[0]    = vertex.product().size()

    njet[0]   = jets.product().size()
    npjet[0]   = pjets.product().size()

    ###CHS JETS
    for i,j in enumerate(jets.product()):

        if i>=maxjet: break

        jet_pt[i]  = j.pt()
        jet_eta[i] = j.eta()
        jet_phi[i] = j.phi()
        jet_energy[i]  = j.energy()

        rawjet_pt[i]  = j.pt()*j.jecFactor("Uncorrected")
        rawjet_eta[i] = j.eta()
        rawjet_phi[i] = j.phi()
        rawjet_energy[i] = j.energy()*j.jecFactor("Uncorrected")

        NHF[i]                 = j.neutralHadronEnergyFraction()
        NEMF[i]                = j.neutralEmEnergyFraction();
        CHF[i]                 = j.chargedHadronEnergyFraction();
        MUF[i]                 = j.muonEnergyFraction();
        CEMF[i]                = j.chargedEmEnergyFraction();
        NumConst[i]            = j.chargedMultiplicity()+j.neutralMultiplicity();
        NumNeutralParticle[i]  = j.neutralMultiplicity();
        CHM[i]                 = j.chargedMultiplicity();
        eta                    = j.eta();
        
        jet_loose[i] = (NHF[i]<0.99 and NEMF[i]<0.99 and NumConst[i]>1) and ((abs(eta)<=2.4 and CHF[i]>0 and CHM[i]>0 and CEMF[i]<0.99) or abs(eta)>2.4) and abs(eta)<=2.7
        jet_loose[i] = jet_loose[i] or (NHF[i]<0.98 and NEMF[i]>0.01 and NumNeutralParticle[i]>2 and abs(eta)>2.7 and abs(eta)<=3.0 )
        jet_loose[i] = jet_loose[i] or (NEMF[i]<0.90 and NumNeutralParticle[i]>10 and abs(eta)>3.0)

        if not (j.genJet() == None):
            genjet_pt[i]  = j.genJet().pt()
            genjet_eta[i] = j.genJet().eta()
            genjet_phi[i] = j.genJet().phi()
            genjet_energy[i] = j.genJet().energy()

        else:
            genjet_pt[i] = -999.
            genjet_eta[i] = -999.
            genjet_phi[i] = -999.
            genjet_energy[i] = -999.


        sourceCandidate = set()

        # now get a list of the PF candidates used to build this jet
        for isource in xrange(j.numberOfSourceCandidatePtrs()):
            sourceCandidate.add(j.sourceCandidatePtr(isource).key()) # the key is the index in the pf collection

        neutral[i], charged[i], photon[i], muon[i], electron[i], hhf[i], ehf[i], other[i] = 0, 0, 0, 0, 0, 0, 0, 0
        neutral_e[i], charged_e[i], photon_e[i],muon_e[i],electron_e[i],hhf_e[i], ehf_e[i], other_e[i] = 0, 0, 0, 0, 0, 0, 0, 0
        neutral_n[i], charged_n[i], photon_n[i],muon_n[i],electron_n[i],hhf_n[i], ehf_n[i], other_n[i] = 0, 0, 0, 0, 0, 0, 0, 0

        for ipf,pf in enumerate(pfs.product()):            

            # if pf candidate is part of the jet, lets categorize
            if ipf in sourceCandidate:

                if ( abs(pf.pdgId()) == 211 ):
                    charged[i]   += pf.pt();
                    charged_e[i] += pf.energy();
                    charged_n[i] += 1;
                    
                elif abs(pf.pdgId()) == 130:
                    neutral[i] += pf.pt();
                    neutral_e[i] += pf.energy();
                    neutral_n[i] += 1;

                elif abs(pf.pdgId()) == 22:
                    photon[i] += pf.pt();
                    photon_e[i] += pf.energy();
                    photon_n[i] += 1;

                elif abs(pf.pdgId()) == 13:
                    muon[i] += pf.pt();
                    muon_e[i] += pf.energy();
                    muon_n[i] += 1;

                elif abs(pf.pdgId()) == 11:
                    electron[i] += pf.pt();
                    electron_e[i] += pf.energy();
                    electron_n[i] += 1;

                elif abs(pf.pdgId()) == 1:
                    hhf[i] += pf.pt();
                    hhf_e[i] += pf.energy();
                    hhf_n[i] += 1;

                elif abs(pf.pdgId()) == 2:
                    ehf[i] += pf.pt();
                    ehf_e[i] += pf.energy();
                    ehf_n[i] += 1;
                
                else:
                    other[i] += pf.pt();
                    other_e[i] += pf.energy();
                    other_n[i] += 1;

            	# store the hcal depth information
                for depth in range(1,8):
					if abs(pf.pdgId())==130 and abs(pf.eta())>1.3 and abs(pf.eta())<5.0:
						jet_depth[i][depth-1] += pf.hcalDepthEnergyFraction(depth)*pf.energy()*pf.hcalFraction()
						jet_depth_uncorrected[i][depth-1] += pf.hcalDepthEnergyFraction(depth)*pf.energy()

    ###Puppi JETS
    for i,j in enumerate(pjets.product()):
	

        if i>=maxjet: break

        pjet_pt[i]  = j.pt()
        pjet_eta[i] = j.eta()
        pjet_phi[i] = j.phi()
        pjet_energy[i]  = j.energy()

        rawpjet_pt[i]  = j.pt()*j.jecFactor("Uncorrected")
        rawpjet_eta[i] = j.eta()
        rawpjet_phi[i] = j.phi()
        rawpjet_energy[i] = j.energy()*j.jecFactor("Uncorrected")

        NHF[i]                 = j.neutralHadronEnergyFraction()
        NEMF[i]                = j.neutralEmEnergyFraction();
        CHF[i]                 = j.chargedHadronEnergyFraction();
        MUF[i]                 = j.muonEnergyFraction();
        CEMF[i]                = j.chargedEmEnergyFraction();
        NumConst[i]            = j.chargedMultiplicity()+j.neutralMultiplicity();
        NumNeutralParticle[i]  = j.neutralMultiplicity();
        CHM[i]                 = j.chargedMultiplicity();
        eta                    = j.eta();
        
        pjet_loose[i] = (NHF[i]<0.99 and NEMF[i]<0.99 and NumConst[i]>1) and ((abs(eta)<=2.4 and CHF[i]>0 and CHM[i]>0 and CEMF[i]<0.99) or abs(eta)>2.4) and abs(eta)<=2.7
        pjet_loose[i] = jet_loose[i] or (NHF[i]<0.98 and NEMF[i]>0.01 and NumNeutralParticle[i]>2 and abs(eta)>2.7 and abs(eta)<=3.0 )
        pjet_loose[i] = jet_loose[i] or (NEMF[i]<0.90 and NumNeutralParticle[i]>10 and abs(eta)>3.0)

        if not (j.genJet() == None):
            genpjet_pt[i]  = j.genJet().pt()
            genpjet_eta[i] = j.genJet().eta()
            genpjet_phi[i] = j.genJet().phi()
            genpjet_energy[i] = j.genJet().energy()

        else:
            genpjet_pt[i] = -999.
            genpjet_eta[i] = -999.
            genpjet_phi[i] = -999.
            genpjet_energy[i] = -999.

        sourceCandidate = set()

        # now get a list of the PF candidates used to build this jet
        for isource in xrange(j.numberOfSourceCandidatePtrs()):
            sourceCandidate.add(j.sourceCandidatePtr(isource).key()) # the key is the index in the pf collection

        for ipf,pf in enumerate(pfs.product()):            
            if ipf in sourceCandidate:

                if ( abs(pf.pdgId()) == 211 ):
                    charged_pjet[i]   += pf.pt()*pf.puppiWeight();
                    charged_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    charged_n_pjet[i] += 1*pf.puppiWeight();
                    
                elif abs(pf.pdgId()) == 130:
                    neutral_pjet[i] += pf.pt()*pf.puppiWeight();
                    neutral_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    neutral_n_pjet[i] += 1*pf.puppiWeight();

                elif abs(pf.pdgId()) == 22:
                    photon_pjet[i] += pf.pt()*pf.puppiWeight();
                    photon_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    photon_n_pjet[i] += 1*pf.puppiWeight();

                elif abs(pf.pdgId()) == 13:
                    muon_pjet[i] += pf.pt()*pf.puppiWeight();
                    muon_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    muon_n_pjet[i] += 1*pf.puppiWeight();

                elif abs(pf.pdgId()) == 11:
                    electron_pjet[i] += pf.pt()*pf.puppiWeight();
                    electron_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    electron_n_pjet[i] += 1*pf.puppiWeight();

                elif abs(pf.pdgId()) == 1:
                    hhf_pjet[i] += pf.pt()*pf.puppiWeight();
                    hhf_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    hhf_n_pjet[i] += 1*pf.puppiWeight();

                elif abs(pf.pdgId()) == 2:
                    ehf_pjet[i] += pf.pt()*pf.puppiWeight();
                    ehf_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    ehf_n_pjet[i] += 1*pf.puppiWeight();
                
                else:
                    other_pjet[i] += pf.pt()*pf.puppiWeight();
                    other_e_pjet[i] += pf.energy()*pf.puppiWeight();
                    other_n_pjet[i] += 1*pf.puppiWeight();

            	# store the hcal depth information
                for depth in range(1,8):
					if abs(pf.pdgId())==130 and abs(pf.eta())>1.3 and abs(pf.eta())<5.0:
						pjet_depth[i][depth-1] += pf.hcalDepthEnergyFraction(depth)*pf.energy()*pf.hcalFraction()*pf.puppiWeight()
						pjet_depth_uncorrected[i][depth-1] += pf.hcalDepthEnergyFraction(depth)*pf.energy()*pf.puppiWeight()


    #met information
    metprod   = mets.product().front()
    met[0]    = metprod.pt()                
    genmet[0] = metprod.genMET().pt()                
    rawmet[0] = metprod.uncorPt()
    
    for ipf,pf in enumerate(pfs.product()):            
        if ( abs(pf.pdgId()) == 211 ):
            charged_met[0]   += pf.pt();
            #charged_n_met[0] += 1;
                    
        elif abs(pf.pdgId()) == 130:
            neutral_met[0] += pf.pt();
            #neutral_n_met[0] += 1;

        elif abs(pf.pdgId()) == 22:
            photon_met[0] += pf.pt();
            #photon_n_met[0] += 1;

        elif abs(pf.pdgId()) == 13:
            muon_met[0] += pf.pt();
            #muon_n_met[0] += 1;
            
        elif abs(pf.pdgId()) == 11:
            electron_met[0] += pf.pt();
            #electron_n_met[0] += 1;

        elif abs(pf.pdgId()) == 1:
            hhf_met[0] += pf.pt();            
            #hhf_n_met[0] += 1;
            
        elif abs(pf.pdgId()) == 2:
            ehf_met[0] += pf.pt();
            #ehf_n_met[0] += 1;
            
        else:
            other_met[0] += pf.pt();
            #other_n_met[0] += 1;
            
    #puppi met information
    pmetprod   = pmets.product().front()
    pmet[0]    = pmetprod.pt()                
    genpmet[0] = pmetprod.genMET().pt()                
    rawpmet[0] = pmetprod.uncorPt()
    
    t.Fill()


fout.Write()
fout.Close()
