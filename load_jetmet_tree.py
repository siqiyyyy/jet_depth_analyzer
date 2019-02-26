import numpy as n
import ROOT

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

#Init branches
b_rhoall                 = ROOT.TBranch()
b_maxjet                 = ROOT.TBranch()
b_nrun                   = ROOT.TBranch()
b_nlumi                  = ROOT.TBranch()
b_nevent                 = ROOT.TBranch()
b_npv                    = ROOT.TBranch()
b_rhoall                 = ROOT.TBranch()
b_rhocentral             = ROOT.TBranch()
b_rhoneutral             = ROOT.TBranch()
b_rhochargedpileup       = ROOT.TBranch()
b_dphipfmet              = ROOT.TBranch()
b_njet                   = ROOT.TBranch()
b_jet_pt                 = ROOT.TBranch()
b_jet_energy             = ROOT.TBranch()
b_jet_eta                = ROOT.TBranch()
b_jet_phi                = ROOT.TBranch()
b_genjet_pt              = ROOT.TBranch()
b_genjet_energy          = ROOT.TBranch()
b_genjet_eta             = ROOT.TBranch()
b_genjet_phi             = ROOT.TBranch()
b_rawjet_pt              = ROOT.TBranch()
b_rawjet_energy          = ROOT.TBranch()
b_rawjet_eta             = ROOT.TBranch()
b_rawjet_phi             = ROOT.TBranch()
b_jet_loose              = ROOT.TBranch()
b_jet_depth              = ROOT.TBranch()
b_jet_depth_uncorrected  = ROOT.TBranch()
b_npjet                  = ROOT.TBranch()
b_pjet_pt                = ROOT.TBranch()
b_pjet_energy            = ROOT.TBranch()
b_pjet_eta               = ROOT.TBranch()
b_pjet_phi               = ROOT.TBranch()
b_genpjet_pt             = ROOT.TBranch()
b_genpjet_energy         = ROOT.TBranch()
b_genpjet_eta            = ROOT.TBranch()
b_genpjet_phi            = ROOT.TBranch()
b_rawpjet_pt             = ROOT.TBranch()
b_rawpjet_energy         = ROOT.TBranch()
b_rawpjet_eta            = ROOT.TBranch()
b_rawpjet_phi            = ROOT.TBranch()
b_pjet_loose             = ROOT.TBranch()
b_pjet_depth             = ROOT.TBranch()
b_pjet_depth_uncorrected = ROOT.TBranch()
b_NHF                    = ROOT.TBranch()
b_NEMF                   = ROOT.TBranch()
b_CHF                    = ROOT.TBranch()
b_MUF                    = ROOT.TBranch()
b_CEMF                   = ROOT.TBranch()
b_NumConst               = ROOT.TBranch()
b_NumNeutralParticle     = ROOT.TBranch()
b_CHM                    = ROOT.TBranch()
b_charged                = ROOT.TBranch()
b_neutral                = ROOT.TBranch()
b_photon                 = ROOT.TBranch()
b_muon                   = ROOT.TBranch()
b_electron               = ROOT.TBranch()
b_hhf                    = ROOT.TBranch()
b_ehf                    = ROOT.TBranch()
b_other                  = ROOT.TBranch()
b_charged_e              = ROOT.TBranch()
b_neutral_e              = ROOT.TBranch()
b_photon_e               = ROOT.TBranch()
b_muon_e                 = ROOT.TBranch()
b_electron_e             = ROOT.TBranch()
b_hhf_e                  = ROOT.TBranch()
b_ehf_e                  = ROOT.TBranch()
b_other_e                = ROOT.TBranch()
b_charged_n              = ROOT.TBranch()
b_neutral_n              = ROOT.TBranch()
b_photon_n               = ROOT.TBranch()
b_muon_n                 = ROOT.TBranch()
b_electron_n             = ROOT.TBranch()
b_hhf_n                  = ROOT.TBranch()
b_ehf_n                  = ROOT.TBranch()
b_other_n                = ROOT.TBranch()
b_charged_pjet           = ROOT.TBranch()
b_neutral_pjet           = ROOT.TBranch()
b_photon_pjet            = ROOT.TBranch()
b_muon_pjet              = ROOT.TBranch()
b_electron_pjet          = ROOT.TBranch()
b_hhf_pjet               = ROOT.TBranch()
b_ehf_pjet               = ROOT.TBranch()
b_other_pjet             = ROOT.TBranch()
b_charged_e_pjet         = ROOT.TBranch()
b_neutral_e_pjet         = ROOT.TBranch()
b_photon_e_pjet          = ROOT.TBranch()
b_muon_e_pjet            = ROOT.TBranch()
b_electron_e_pjet        = ROOT.TBranch()
b_hhf_e_pjet             = ROOT.TBranch()
b_ehf_e_pjet             = ROOT.TBranch()
b_other_e_pjet           = ROOT.TBranch()
b_charged_n_pjet         = ROOT.TBranch()
b_neutral_n_pjet         = ROOT.TBranch()
b_photon_n_pjet          = ROOT.TBranch()
b_muon_n_pjet            = ROOT.TBranch()
b_electron_n_pjet        = ROOT.TBranch()
b_hhf_n_pjet             = ROOT.TBranch()
b_ehf_n_pjet             = ROOT.TBranch()
b_other_n_pjet           = ROOT.TBranch()
b_met                    = ROOT.TBranch()
b_mex                    = ROOT.TBranch()
b_mey                    = ROOT.TBranch()
b_met_phi                = ROOT.TBranch()
b_genmet                 = ROOT.TBranch()
b_rawmet                 = ROOT.TBranch()
b_charged_met            = ROOT.TBranch()
b_neutral_met            = ROOT.TBranch()
b_photon_met             = ROOT.TBranch()
b_muon_met               = ROOT.TBranch()
b_electron_met           = ROOT.TBranch()
b_hhf_met                = ROOT.TBranch()
b_ehf_met                = ROOT.TBranch()
b_other_met              = ROOT.TBranch()
b_chsmet                 = ROOT.TBranch()
b_trkmet                 = ROOT.TBranch()
b_phomet                 = ROOT.TBranch()
b_neumet                 = ROOT.TBranch()
b_pmet                   = ROOT.TBranch()
b_pmex                   = ROOT.TBranch()
b_pmey                   = ROOT.TBranch()
b_pmet_phi               = ROOT.TBranch()
b_genpmet                = ROOT.TBranch()
b_rawpmet                = ROOT.TBranch()
    
def declare_branches(t):

    t.SetBranchAddress("run"                    , nrun                   , b_nrun                   ) 
    t.SetBranchAddress("lumi"                   , nlumi                  , b_nlumi                  ) 
    t.SetBranchAddress("event"                  , nevent                 , b_nevent                 ) 
    t.SetBranchAddress("npv"                    , npv                    , b_npv                    ) 
    t.SetBranchAddress("dphipfmet"              , dphipfmet              , b_dphipfmet              ) 
    t.SetBranchAddress("rhoall"                 , rhoall                 , b_rhoall                 ) 
    t.SetBranchAddress("rhocentral"             , rhocentral             , b_rhocentral             ) 
    t.SetBranchAddress("rhoneutral"             , rhoneutral             , b_rhoneutral             ) 
    t.SetBranchAddress("rhochargedpileup"       , rhochargedpileup       , b_rhochargedpileup       ) 
    t.SetBranchAddress("njet"                   , njet                   , b_njet                   ) 
    t.SetBranchAddress("npjet"                  , npjet                  , b_npjet                  ) 
    t.SetBranchAddress("jet_pt"                 , jet_pt                 , b_jet_pt                 ) 
    t.SetBranchAddress("jet_energy"             , jet_energy             , b_jet_energy             ) 
    t.SetBranchAddress("jet_eta"                , jet_eta                , b_jet_eta                ) 
    t.SetBranchAddress("jet_phi"                , jet_phi                , b_jet_phi                ) 
    t.SetBranchAddress("genjet_pt"              , genjet_pt              , b_genjet_pt              ) 
    t.SetBranchAddress("genjet_energy"          , genjet_energy          , b_genjet_energy          ) 
    t.SetBranchAddress("genjet_eta"             , genjet_eta             , b_genjet_eta             ) 
    t.SetBranchAddress("genjet_phi"             , genjet_phi             , b_genjet_phi             ) 
    t.SetBranchAddress("rawjet_pt"              , rawjet_pt              , b_rawjet_pt              ) 
    t.SetBranchAddress("rawjet_energy"          , rawjet_energy          , b_rawjet_energy          ) 
    t.SetBranchAddress("rawjet_eta"             , rawjet_eta             , b_rawjet_eta             ) 
    t.SetBranchAddress("rawjet_phi"             , rawjet_phi             , b_rawjet_phi             ) 
    t.SetBranchAddress("pjet_pt"                , pjet_pt                , b_pjet_pt                ) 
    t.SetBranchAddress("pjet_energy"            , pjet_energy            , b_pjet_energy            ) 
    t.SetBranchAddress("pjet_eta"               , pjet_eta               , b_pjet_eta               ) 
    t.SetBranchAddress("pjet_phi"               , pjet_phi               , b_pjet_phi               ) 
    t.SetBranchAddress("genpjet_pt"             , genpjet_pt             , b_genpjet_pt             ) 
    t.SetBranchAddress("genpjet_energy"         , genpjet_energy         , b_genpjet_energy         ) 
    t.SetBranchAddress("genpjet_eta"            , genpjet_eta            , b_genpjet_eta            ) 
    t.SetBranchAddress("genpjet_phi"            , genpjet_phi            , b_genpjet_phi            ) 
    t.SetBranchAddress("rawpjet_pt"             , rawpjet_pt             , b_rawpjet_pt             ) 
    t.SetBranchAddress("rawpjet_energy"         , rawpjet_energy         , b_rawpjet_energy         ) 
    t.SetBranchAddress("rawpjet_eta"            , rawpjet_eta            , b_rawpjet_eta            ) 
    t.SetBranchAddress("rawpjet_phi"            , rawpjet_phi            , b_rawpjet_phi            ) 
    t.SetBranchAddress("NHF"                    , NHF                    , b_NHF                    ) 
    t.SetBranchAddress("NEMF"                   , NEMF                   , b_NEMF                   ) 
    t.SetBranchAddress("CHF"                    , CHF                    , b_CHF                    ) 
    t.SetBranchAddress("MUF"                    , MUF                    , b_MUF                    ) 
    t.SetBranchAddress("CEMF"                   , CEMF                   , b_CEMF                   ) 
    t.SetBranchAddress("NumConst"               , NumConst               , b_NumConst               ) 
    t.SetBranchAddress("NumNeutralParticle"     , NumNeutralParticle     , b_NumNeutralParticle     ) 
    t.SetBranchAddress("CHM"                    , CHM                    , b_CHM                    ) 
    t.SetBranchAddress("jet_loose"              , jet_loose              , b_jet_loose              ) 
    t.SetBranchAddress("pjet_loose"             , pjet_loose             , b_pjet_loose             ) 
    t.SetBranchAddress("jet_depth"              , jet_depth              , b_jet_depth              ) 
    t.SetBranchAddress("jet_depth_uncorrected"  , jet_depth_uncorrected  , b_jet_depth_uncorrected  ) 
    t.SetBranchAddress("pjet_depth"             , pjet_depth             , b_pjet_depth             ) 
    t.SetBranchAddress("pjet_depth_uncorrected" , pjet_depth_uncorrected , b_pjet_depth_uncorrected ) 
    t.SetBranchAddress("charged"                , charged                , b_charged                ) 
    t.SetBranchAddress("neutral"                , neutral                , b_neutral                ) 
    t.SetBranchAddress("photon"                 , photon                 , b_photon                 ) 
    t.SetBranchAddress("muon"                   , muon                   , b_muon                   ) 
    t.SetBranchAddress("electron"               , electron               , b_electron               ) 
    t.SetBranchAddress("hhf"                    , hhf                    , b_hhf                    ) 
    t.SetBranchAddress("ehf"                    , ehf                    , b_ehf                    ) 
    t.SetBranchAddress("other"                  , other                  , b_other                  ) 
    t.SetBranchAddress("charged_e"              , charged_e              , b_charged_e              ) 
    t.SetBranchAddress("neutral_e"              , neutral_e              , b_neutral_e              ) 
    t.SetBranchAddress("photon_e"               , photon_e               , b_photon_e               ) 
    t.SetBranchAddress("muon_e"                 , muon_e                 , b_muon_e                 ) 
    t.SetBranchAddress("electron_e"             , electron_e             , b_electron_e             ) 
    t.SetBranchAddress("hhf_e"                  , hhf_e                  , b_hhf_e                  ) 
    t.SetBranchAddress("ehf_e"                  , ehf_e                  , b_ehf_e                  ) 
    t.SetBranchAddress("other_e"                , other_e                , b_other_e                ) 
    t.SetBranchAddress("charged_n"              , charged_n              , b_charged_n              ) 
    t.SetBranchAddress("neutral_n"              , neutral_n              , b_neutral_n              ) 
    t.SetBranchAddress("photon_n"               , photon_n               , b_photon_n               ) 
    t.SetBranchAddress("muon_n"                 , muon_n                 , b_muon_n                 ) 
    t.SetBranchAddress("electron_n"             , electron_n             , b_electron_n             ) 
    t.SetBranchAddress("hhf_n"                  , hhf_n                  , b_hhf_n                  ) 
    t.SetBranchAddress("ehf_n"                  , ehf_n                  , b_ehf_n                  ) 
    t.SetBranchAddress("other_n"                , other_n                , b_other_n                ) 
    t.SetBranchAddress("charged_pjet"           , charged_pjet           , b_charged_pjet           ) 
    t.SetBranchAddress("neutral_pjet"           , neutral_pjet           , b_neutral_pjet           ) 
    t.SetBranchAddress("photon_pjet"            , photon_pjet            , b_photon_pjet            ) 
    t.SetBranchAddress("muon_pjet"              , muon_pjet              , b_muon_pjet              ) 
    t.SetBranchAddress("electron_pjet"          , electron_pjet          , b_electron_pjet          ) 
    t.SetBranchAddress("hhf_pjet"               , hhf_pjet               , b_hhf_pjet               ) 
    t.SetBranchAddress("ehf_pjet"               , ehf_pjet               , b_ehf_pjet               ) 
    t.SetBranchAddress("other_pjet"             , other_pjet             , b_other_pjet             ) 
    t.SetBranchAddress("charged_e_pjet"         , charged_e_pjet         , b_charged_e_pjet         ) 
    t.SetBranchAddress("neutral_e_pjet"         , neutral_e_pjet         , b_neutral_e_pjet         ) 
    t.SetBranchAddress("photon_e_pjet"          , photon_e_pjet          , b_photon_e_pjet          ) 
    t.SetBranchAddress("muon_e_pjet"            , muon_e_pjet            , b_muon_e_pjet            ) 
    t.SetBranchAddress("electron_e_pjet"        , electron_e_pjet        , b_electron_e_pjet        ) 
    t.SetBranchAddress("hhf_e_pjet"             , hhf_e_pjet             , b_hhf_e_pjet             ) 
    t.SetBranchAddress("ehf_e_pjet"             , ehf_e_pjet             , b_ehf_e_pjet             ) 
    t.SetBranchAddress("other_e_pjet"           , other_e_pjet           , b_other_e_pjet           ) 
    t.SetBranchAddress("charged_n_pjet"         , charged_n_pjet         , b_charged_n_pjet         ) 
    t.SetBranchAddress("neutral_n_pjet"         , neutral_n_pjet         , b_neutral_n_pjet         ) 
    t.SetBranchAddress("photon_n_pjet"          , photon_n_pjet          , b_photon_n_pjet          ) 
    t.SetBranchAddress("muon_n_pjet"            , muon_n_pjet            , b_muon_n_pjet            ) 
    t.SetBranchAddress("electron_n_pjet"        , electron_n_pjet        , b_electron_n_pjet        ) 
    t.SetBranchAddress("hhf_n_pjet"             , hhf_n_pjet             , b_hhf_n_pjet             ) 
    t.SetBranchAddress("ehf_n_pjet"             , ehf_n_pjet             , b_ehf_n_pjet             ) 
    t.SetBranchAddress("other_n_pjet"           , other_n_pjet           , b_other_n_pjet           ) 
    t.SetBranchAddress("met"                    , met                    , b_met                    ) 
    t.SetBranchAddress("mex"                    , mex                    , b_mex                    ) 
    t.SetBranchAddress("mey"                    , mey                    , b_mey                    ) 
    t.SetBranchAddress("met_phi"                , met_phi                , b_met_phi                ) 
    t.SetBranchAddress("genmet"                 , genmet                 , b_genmet                 ) 
    t.SetBranchAddress("rawmet"                 , rawmet                 , b_rawmet                 ) 
    t.SetBranchAddress("charged_met"            , charged_met            , b_charged_met            ) 
    t.SetBranchAddress("neutral_met"            , neutral_met            , b_neutral_met            ) 
    t.SetBranchAddress("photon_met"             , photon_met             , b_photon_met             ) 
    t.SetBranchAddress("muon_met"               , muon_met               , b_muon_met               ) 
    t.SetBranchAddress("electron_met"           , electron_met           , b_electron_met           ) 
    t.SetBranchAddress("hhf_met"                , hhf_met                , b_hhf_met                ) 
    t.SetBranchAddress("ehf_met"                , ehf_met                , b_ehf_met                ) 
    t.SetBranchAddress("other_met"              , other_met              , b_other_met              ) 
    t.SetBranchAddress("pmet"                   , pmet                   , b_pmet                   ) 
    t.SetBranchAddress("pmex"                   , pmex                   , b_pmex                   ) 
    t.SetBranchAddress("pmey"                   , pmey                   , b_pmey                   ) 
    t.SetBranchAddress("pmet_phi"               , pmet_phi               , b_pmet_phi               ) 
    t.SetBranchAddress("genpmet"                , genpmet                , b_genpmet                ) 
    t.SetBranchAddress("rawpmet"                , rawpmet                , b_rawpmet                ) 


    print "All branches configured"