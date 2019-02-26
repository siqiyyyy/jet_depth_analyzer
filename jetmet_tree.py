import numpy as n

rhoall     = n.zeros(1,dtype=float)
maxjet = 10000
    
    #event information
nrun       = n.zeros(1,dtype=int)
nlumi      = n.zeros(1,dtype=int)
nevent     = n.zeros(1,dtype=float)
npv        = n.zeros(1,dtype=int)
rhoall     = n.zeros(1,dtype=float)
rhocentral = n.zeros(1,dtype=float)
rhoneutral = n.zeros(1,dtype=float)
rhochargedpileup = n.zeros(1,dtype=float)
dphipfmet = n.zeros(1,dtype=float)

    #jet information
njet       = n.zeros(1,dtype=int)
jet_pt     = n.zeros(maxjet,dtype=float)
jet_energy = n.zeros(maxjet,dtype=float)
jet_eta    = n.zeros(maxjet,dtype=float)
jet_phi    = n.zeros(maxjet,dtype=float)
genjet_pt  = n.zeros(maxjet,dtype=float)
genjet_energy = n.zeros(maxjet,dtype=float)
genjet_eta = n.zeros(maxjet,dtype=float)
genjet_phi = n.zeros(maxjet,dtype=float)
rawjet_pt  = n.zeros(maxjet,dtype=float)
rawjet_energy = n.zeros(maxjet,dtype=float)
rawjet_eta = n.zeros(maxjet,dtype=float)
rawjet_phi = n.zeros(maxjet,dtype=float)
jet_loose  = n.zeros(maxjet,dtype=int)
jet_depth = n.zeros([maxjet,7],dtype=float) #hcal energy distributed in each layer
jet_depth_uncorrected = n.zeros([maxjet,7],dtype=float) #similar as above but used total energy rather than hcal energy
    #puppi jet information
npjet       = n.zeros(1,dtype=int)
pjet_pt     = n.zeros(maxjet,dtype=float)
pjet_energy = n.zeros(maxjet,dtype=float)
pjet_eta    = n.zeros(maxjet,dtype=float)
pjet_phi    = n.zeros(maxjet,dtype=float)
genpjet_pt  = n.zeros(maxjet,dtype=float)
genpjet_energy = n.zeros(maxjet,dtype=float)
genpjet_eta = n.zeros(maxjet,dtype=float)
genpjet_phi = n.zeros(maxjet,dtype=float)
rawpjet_pt  = n.zeros(maxjet,dtype=float)
rawpjet_energy = n.zeros(maxjet,dtype=float)
rawpjet_eta = n.zeros(maxjet,dtype=float)
rawpjet_phi = n.zeros(maxjet,dtype=float)
pjet_loose  = n.zeros(maxjet,dtype=int)
pjet_depth = n.zeros([maxjet,7],dtype=float) #hcal energy distributed in each layer
pjet_depth_uncorrected = n.zeros([maxjet,7],dtype=float) #similar as above but used total energy rather than hcal energy


#jet energy fraction
NHF = n.zeros(maxjet,dtype=float)
NEMF = n.zeros(maxjet,dtype=float)
CHF = n.zeros(maxjet,dtype=float)
MUF = n.zeros(maxjet,dtype=float)
CEMF = n.zeros(maxjet,dtype=float)
NumConst = n.zeros(maxjet,dtype=int)
NumNeutralParticle = n.zeros(maxjet,dtype=int)
CHM = n.zeros(maxjet,dtype=int)


    #jet fraction information for each jet
charged    = n.zeros(maxjet,dtype=float)
neutral    = n.zeros(maxjet,dtype=float)
photon     = n.zeros(maxjet,dtype=float)
muon       = n.zeros(maxjet,dtype=float)
electron   = n.zeros(maxjet,dtype=float)
hhf        = n.zeros(maxjet,dtype=float)
ehf        = n.zeros(maxjet,dtype=float)
other      = n.zeros(maxjet,dtype=float)

charged_e  = n.zeros(maxjet,dtype=float)
neutral_e  = n.zeros(maxjet,dtype=float)
photon_e   = n.zeros(maxjet,dtype=float)
muon_e     = n.zeros(maxjet,dtype=float)
electron_e = n.zeros(maxjet,dtype=float)
hhf_e      = n.zeros(maxjet,dtype=float)
ehf_e      = n.zeros(maxjet,dtype=float)
other_e    = n.zeros(maxjet,dtype=float)

charged_n  = n.zeros(maxjet,dtype=int)
neutral_n  = n.zeros(maxjet,dtype=int)
photon_n   = n.zeros(maxjet,dtype=int)
muon_n     = n.zeros(maxjet,dtype=int)
electron_n = n.zeros(maxjet,dtype=int)
hhf_n      = n.zeros(maxjet,dtype=int)
ehf_n      = n.zeros(maxjet,dtype=int)
other_n    = n.zeros(maxjet,dtype=int)

    #puppi jet fraction information for each jet
charged_pjet    = n.zeros(maxjet,dtype=float)
neutral_pjet    = n.zeros(maxjet,dtype=float)
photon_pjet     = n.zeros(maxjet,dtype=float)
muon_pjet       = n.zeros(maxjet,dtype=float)
electron_pjet   = n.zeros(maxjet,dtype=float)
hhf_pjet        = n.zeros(maxjet,dtype=float)
ehf_pjet        = n.zeros(maxjet,dtype=float)
other_pjet      = n.zeros(maxjet,dtype=float)

charged_e_pjet  = n.zeros(maxjet,dtype=float)
neutral_e_pjet  = n.zeros(maxjet,dtype=float)
photon_e_pjet   = n.zeros(maxjet,dtype=float)
muon_e_pjet     = n.zeros(maxjet,dtype=float)
electron_e_pjet = n.zeros(maxjet,dtype=float)
hhf_e_pjet      = n.zeros(maxjet,dtype=float)
ehf_e_pjet      = n.zeros(maxjet,dtype=float)
other_e_pjet    = n.zeros(maxjet,dtype=float)

charged_n_pjet  = n.zeros(maxjet,dtype=int)
neutral_n_pjet  = n.zeros(maxjet,dtype=int)
photon_n_pjet   = n.zeros(maxjet,dtype=int)
muon_n_pjet     = n.zeros(maxjet,dtype=int)
electron_n_pjet = n.zeros(maxjet,dtype=int)
hhf_n_pjet      = n.zeros(maxjet,dtype=int)
ehf_n_pjet      = n.zeros(maxjet,dtype=int)
other_n_pjet    = n.zeros(maxjet,dtype=int)

    #met information
met    = n.zeros(1,dtype=float)
mex    = n.zeros(1,dtype=float)
mey    = n.zeros(1,dtype=float)
met_phi= n.zeros(1,dtype=float)
genmet = n.zeros(1,dtype=float)
rawmet = n.zeros(1,dtype=float)
charged_met = n.zeros(1,dtype=float)
neutral_met = n.zeros(1,dtype=float)
photon_met = n.zeros(1,dtype=float)
muon_met = n.zeros(1,dtype=float)
electron_met = n.zeros(1,dtype=float)
hhf_met = n.zeros(1,dtype=float)
ehf_met = n.zeros(1,dtype=float)
other_met = n.zeros(1,dtype=float)
chsmet = n.zeros(1,dtype=float)
trkmet = n.zeros(1,dtype=float)
phomet = n.zeros(1,dtype=float)
neumet = n.zeros(1,dtype=float)
    #pmet information
pmet    = n.zeros(1,dtype=float)
pmex    = n.zeros(1,dtype=float)
pmey    = n.zeros(1,dtype=float)
pmet_phi= n.zeros(1,dtype=float)
genpmet = n.zeros(1,dtype=float)
rawpmet = n.zeros(1,dtype=float)
    
def declare_branches(t):

    t.Branch("run", nrun, 'run/I')
    t.Branch("lumi", nlumi, 'lumi/I')
    t.Branch("event", nevent, 'event/D')
    t.Branch("npv", npv, 'npv/I')

    t.Branch("dphipfmet", dphipfmet, 'dphipfmet/D')

    t.Branch("rhoall", rhoall, 'rhoall/D')
    t.Branch("rhocentral", rhocentral, 'rhocentral/D')
    t.Branch("rhoneutral", rhoneutral, 'rhoneutral/D')
    t.Branch("rhochargedpileup", rhochargedpileup, 'rhochargedpileup/D')

    t.Branch("njet", njet, 'njet/I')
    t.Branch("npjet", npjet, 'npjet/I')

    t.Branch("jet_pt",jet_pt,'jet_pt[njet]/D')
    t.Branch("jet_energy",jet_energy,'jet_energy[njet]/D')
    t.Branch("jet_eta",jet_eta,'jet_eta[njet]/D')
    t.Branch("jet_phi",jet_phi,'jet_phi[njet]/D')

    t.Branch("genjet_pt",genjet_pt,'genjet_pt[njet]/D')
    t.Branch("genjet_energy",genjet_energy,'genjet_energy[njet]/D')
    t.Branch("genjet_eta",genjet_eta,'genjet_eta[njet]/D')
    t.Branch("genjet_phi",genjet_phi,'genjet_phi[njet]/D')

    t.Branch("rawjet_pt",rawjet_pt,'rawjet_pt[njet]/D')
    t.Branch("rawjet_energy",rawjet_energy,'rawjet_energy[njet]/D')
    t.Branch("rawjet_eta",rawjet_eta,'rawjet_eta[njet]/D')
    t.Branch("rawjet_phi",rawjet_phi,'rawjet_phi[njet]/D')

    t.Branch("pjet_pt",pjet_pt,'pjet_pt[npjet]/D')
    t.Branch("pjet_energy",pjet_energy,'pjet_energy[npjet]/D')
    t.Branch("pjet_eta",pjet_eta,'pjet_eta[npjet]/D')
    t.Branch("pjet_phi",pjet_phi,'pjet_phi[npjet]/D')

    t.Branch("genpjet_pt",genpjet_pt,'genpjet_pt[npjet]/D')
    t.Branch("genpjet_energy",genpjet_energy,'genpjet_energy[npjet]/D')
    t.Branch("genpjet_eta",genpjet_eta,'genpjet_eta[npjet]/D')
    t.Branch("genpjet_phi",genpjet_phi,'genpjet_phi[npjet]/D')

    t.Branch("rawpjet_pt",rawpjet_pt,'rawpjet_pt[npjet]/D')
    t.Branch("rawpjet_energy",rawpjet_energy,'rawpjet_energy[npjet]/D')
    t.Branch("rawpjet_eta",rawpjet_eta,'rawpjet_eta[npjet]/D')
    t.Branch("rawpjet_phi",rawpjet_phi,'rawpjet_phi[npjet]/D')

    t.Branch("NHF",NHF,'NHF[njet]/D')
    t.Branch("NEMF",NEMF,'NEMF[njet]/D')
    t.Branch("CHF",CHF,'CHF[njet]/D')
    t.Branch("MUF",MUF,'MUF[njet]/D')
    t.Branch("CEMF",CEMF,'CEMF[njet]/D')
    t.Branch("NumConst",NumConst,'NumConst[njet]/I')
    t.Branch("NumNeutralParticle",NumNeutralParticle,'NumNeutralParticle[njet]/I')
    t.Branch("CHM",CHM,'CHM[njet]/I')

    t.Branch("jet_loose",jet_loose,'jet_loose[njet]/I')
    t.Branch("pjet_loose",pjet_loose,'pjet_loose[npjet]/I')
    t.Branch("jet_depth",jet_depth,'jet_depth[njet][7]/D')
    t.Branch("jet_depth_uncorrected",jet_depth_uncorrected,'jet_depth_uncorrected[njet][7]/D')
    t.Branch("pjet_depth",pjet_depth,'pjet_depth[njet][7]/D')
    t.Branch("pjet_depth_uncorrected",pjet_depth_uncorrected,'pjet_depth_uncorrected[njet][7]/D')

    t.Branch("charged", charged, 'charged[njet]/D')
    t.Branch("neutral", neutral, 'neutral[njet]/D')
    t.Branch("photon", photon, 'photon[njet]/D')
    t.Branch("muon", muon, 'muon[njet]/D')
    t.Branch("electron", electron, 'electron[njet]/D')
    t.Branch("hhf", hhf, 'hhf[njet]/D')
    t.Branch("ehf", ehf, 'ehf[njet]/D')
    t.Branch("other", other, 'other[njet]/D')

    t.Branch("charged_e", charged_e, 'charged_e[njet]/D')
    t.Branch("neutral_e", neutral_e, 'neutral_e[njet]/D')
    t.Branch("photon_e", photon_e, 'photon_e[njet]/D')
    t.Branch("muon_e", muon_e, 'muon_e[njet]/D')
    t.Branch("electron_e", electron_e, 'electron_e[njet]/D')
    t.Branch("hhf_e", hhf_e, 'hhf_e[njet]/D')
    t.Branch("ehf_e", ehf_e, 'ehf_e[njet]/D')
    t.Branch("other_e", other_e, 'other_e[njet]/D')

    t.Branch("charged_n", charged_n, 'charged_n[njet]/I')
    t.Branch("neutral_n", neutral_n, 'neutral_n[njet]/I')
    t.Branch("photon_n", photon_n, 'photon_n[njet]/I')
    t.Branch("muon_n", muon_n, 'muon_n[njet]/I')
    t.Branch("electron_n", electron_n, 'electron_n[njet]/I')
    t.Branch("hhf_n", hhf_n, 'hhf_n[njet]/I')
    t.Branch("ehf_n", ehf_n, 'ehf_n[njet]/I')
    t.Branch("other_n", other_n, 'other_n[njet]/I')

    t.Branch("charged_pjet", charged_pjet, 'charged[njet]/D')
    t.Branch("neutral_pjet", neutral_pjet, 'neutral[njet]/D')
    t.Branch("photon_pjet", photon_pjet, 'photon[njet]/D')
    t.Branch("muon_pjet", muon_pjet, 'muon[njet]/D')
    t.Branch("electron_pjet", electron_pjet, 'electron[njet]/D')
    t.Branch("hhf_pjet", hhf_pjet, 'hhf[njet]/D')
    t.Branch("ehf_pjet", ehf_pjet, 'ehf[njet]/D')
    t.Branch("other_pjet", other_pjet, 'other[njet]/D')

    t.Branch("charged_e_pjet", charged_e_pjet, 'charged_e_pjet[njet]/D')
    t.Branch("neutral_e_pjet", neutral_e_pjet, 'neutral_e_pjet[njet]/D')
    t.Branch("photon_e_pjet", photon_e_pjet, 'photon_e_pjet[njet]/D')
    t.Branch("muon_e_pjet", muon_e_pjet, 'muon_e_pjet[njet]/D')
    t.Branch("electron_e_pjet", electron_e_pjet, 'electron_e_pjet[njet]/D')
    t.Branch("hhf_e_pjet", hhf_e_pjet, 'hhf_e_pjet[njet]/D')
    t.Branch("ehf_e_pjet", ehf_e_pjet, 'ehf_e_pjet[njet]/D')
    t.Branch("other_e_pjet", other_e_pjet, 'other_e_pjet[njet]/D')

    t.Branch("charged_n_pjet", charged_n_pjet, 'charged_n_pjet[njet]/I')
    t.Branch("neutral_n_pjet", neutral_n_pjet, 'neutral_n_pjet[njet]/I')
    t.Branch("photon_n_pjet", photon_n_pjet, 'photon_n_pjet[njet]/I')
    t.Branch("muon_n_pjet", muon_n_pjet, 'muon_n_pjet[njet]/I')
    t.Branch("electron_n_pjet", electron_n_pjet, 'electron_n_pjet[njet]/I')
    t.Branch("hhf_n_pjet", hhf_n_pjet, 'hhf_n_pjet[njet]/I')
    t.Branch("ehf_n_pjet", ehf_n_pjet, 'ehf_n_pjet[njet]/I')
    t.Branch("other_n_pjet", other_n_pjet, 'other_n_pjet[njet]/I')

    t.Branch("met", met, 'met/D')
    t.Branch("mex", mex, 'mex/D')
    t.Branch("mey", mey, 'mey/D')
    t.Branch("met_phi", met_phi, 'met_phi/D')
    t.Branch("genmet", genmet, 'genmet/D')
    t.Branch("rawmet", rawmet, 'rawmet/D')
    t.Branch("charged_met", charged_met, 'charged_met/D')
    t.Branch("neutral_met", neutral_met, 'neutral_met/D')
    t.Branch("photon_met", photon_met, 'photon_met/D')
    t.Branch("muon_met", muon_met, 'muon_met/D')
    t.Branch("electron_met", electron_met, 'electron_met/D')
    t.Branch("hhf_met", hhf_met, 'hhf_met/D')
    t.Branch("ehf_met", ehf_met, 'ehf_met/D')
    t.Branch("other_met", other_met, 'other_met/D')

    t.Branch("pmet", pmet, 'pmet/D')
    t.Branch("pmex", pmex, 'pmex/D')
    t.Branch("pmey", pmey, 'pmey/D')
    t.Branch("pmet_phi", pmet_phi, 'pmet_phi/D')
    t.Branch("genpmet", genpmet, 'genpmet/D')
    t.Branch("rawpmet", rawpmet, 'rawpmet/D')


    print "All branches configured"