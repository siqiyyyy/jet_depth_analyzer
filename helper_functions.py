# Standard importts
import os,sys,socket,argparse
import os
import ROOT
import math
from array import array
import numpy as n

ROOT.gROOT.SetBatch(True)

# RooFit
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
ROOT.gROOT.SetStyle("Plain") # Not sure this is needed
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )

def GaussFit( shape, isData, var_name, label, fit_plot_directory, fit_filename = None):

    print "Performing a fit using gaus x landau to get the mean and the width" 
    # declare the observable mean, and import the histogram to a RooDataHist

    tmp_mean = shape.GetMean()
    tmp_sigma = shape.GetRMS()

    tmp_mean_error = shape.GetMeanError()
    tmp_sigma_error = shape.GetRMSError()

    #asymmetry   = ROOT.RooRealVar(var_name,label,tmp_mean-5*tmp_sigma,tmp_mean+6*tmp_sigma) ;
    asymmetry   = ROOT.RooRealVar(var_name,label,tmp_mean-4*tmp_sigma,tmp_mean+4*tmp_sigma) ;
    dh          = ROOT.RooDataHist("datahistshape","datahistshape",ROOT.RooArgList(asymmetry),ROOT.RooFit.Import(shape)) ;
    
    # plot the data hist with error from sum of weighted events
    frame       = asymmetry.frame(ROOT.RooFit.Title(var_name))
    if isData:
        dh.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    else:
        dh.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)) ;


    # create a simple gaussian pdf
    gauss_mean  = ROOT.RooRealVar("mean","mean", tmp_mean, tmp_mean-0.5*tmp_sigma, tmp_mean+0.5*tmp_sigma)
    gauss_sigma = ROOT.RooRealVar("sigma","sigma gauss",tmp_sigma, 0, 2.0)
    gauss       = ROOT.RooGaussian("gauss","gauss",asymmetry, gauss_mean, gauss_sigma) 

    # now do the fit and extract the parameters with the correct error
    gauss.fitTo(dh,ROOT.RooFit.Save(),ROOT.RooFit.Range(dh.mean(asymmetry)-dh.sigma(asymmetry),dh.mean(asymmetry)+dh.sigma(asymmetry)))
    gauss.plotOn(frame)
    argset_fit = ROOT.RooArgSet(gauss_mean,gauss_sigma)
    gauss.paramOn(frame,ROOT.RooFit.Format("NELU",ROOT.RooFit.AutoPrecision(1)),ROOT.RooFit.Layout(0.55)) 
    frame.SetMaximum(frame.GetMaximum()*1.2)

    # add chi2 info
    chi2_text = ROOT.TPaveText(0.3,0.8,0.4,0.9,"BRNDC")
    chi2_text.AddText("#chi^{2} fit = %s" %round(frame.chiSquare(6),2))
    chi2_text.SetTextSize(0.04)
    chi2_text.SetTextColor(2)
    chi2_text.SetShadowColor(0)
    chi2_text.SetFillColor(0)
    chi2_text.SetLineColor(0)
    frame.addObject(chi2_text)

    if fit_filename is not None:
        c = ROOT.TCanvas("cfit","cfit",600,700)
        frame.Draw()
        if not os.path.exists(fit_plot_directory): os.makedirs(fit_plot_directory)
        c.SaveAs(os.path.join( fit_plot_directory, fit_filename+".pdf"))
        c.SaveAs(os.path.join( fit_plot_directory, fit_filename+".png"))
        del c

    mean_asymmetry        = gauss_mean.getVal()
    mean_asymmetry_error  = gauss_mean.getError()

    rms_asymmetry        = gauss_sigma.getVal()
    rms_asymmetry_error  = gauss_sigma.getError()

    #if the chi2 is going crazy default to the mean and rms of the histogram
    if frame.chiSquare(6)  > -1:
        print "Warning: Crazy chi2, setting mean and rms to that of the histogram"
    mean_asymmetry = tmp_mean
    mean_asymmetry_error = tmp_mean_error
    rms_asymmetry  = tmp_sigma
    rms_asymmetry_error  = tmp_sigma_error

    return mean_asymmetry, mean_asymmetry_error, rms_asymmetry, rms_asymmetry_error

def ConvFit( shape, isData, var_name, label, fit_plot_directory, fit_filename = None):

    print "Performing a fit using gaus x landau to get the mean and the width" 
    # declare the observable mean, and import the histogram to a RooDataHist

    tmp_mean = shape.GetMean()
    tmp_sigma = shape.GetRMS()

    tmp_mean_error = shape.GetMeanError()
    tmp_sigma_error = shape.GetRMSError()

    #asymmetry   = ROOT.RooRealVar(var_name,label,tmp_mean-5*tmp_sigma,tmp_mean+6*tmp_sigma) ;
    asymmetry   = ROOT.RooRealVar(var_name,label,tmp_mean-4*tmp_sigma,tmp_mean+4*tmp_sigma) ;
    dh          = ROOT.RooDataHist("datahistshape","datahistshape",ROOT.RooArgList(asymmetry),ROOT.RooFit.Import(shape)) ;
    
    # plot the data hist with error from sum of weighted events
    frame       = asymmetry.frame(ROOT.RooFit.Title(var_name))
    if isData:
        dh.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    else:
        dh.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)) ;


    # create a simple gaussian pdf
    gauss_mean  = ROOT.RooRealVar("mean","mean", 0 )
    gauss_sigma = ROOT.RooRealVar("sigma","sigma gauss",tmp_sigma,0,2.0)
    gauss       = ROOT.RooGaussian("gauss","gauss",asymmetry,gauss_mean,gauss_sigma) 

    landau_mean  = ROOT.RooRealVar("meanl","mean landau",1,0.70,1.5)
    landau_sigma = ROOT.RooRealVar("sigmal","sigma landau",tmp_sigma,0,2.0)
    landau       = ROOT.RooLandau("landau","landau",asymmetry,landau_mean,landau_sigma)

    lxg = ROOT.RooFFTConvPdf("lxg","landau x gauss",asymmetry,landau,gauss)

    # now do the fit and extract the parameters with the correct error
    if True:#isData: 
        gauss.fitTo(dh,ROOT.RooFit.Save(),ROOT.RooFit.Range(dh.mean(asymmetry)-2*dh.sigma(asymmetry),dh.mean(asymmetry)+2*dh.sigma(asymmetry)))
        gauss.plotOn(frame)
    else:
        #lxg.fitTo(dh,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Range(dh.mean(asymmetry)-3*dh.sigma(asymmetry),dh.mean(asymmetry)+4*dh.sigma(asymmetry)))
        #lxg.fitTo(dh,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Range(dh.mean(asymmetry)-2*dh.sigma(asymmetry),dh.mean(asymmetry)+2*dh.sigma(asymmetry)))
        lxg.fitTo(dh,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(True))
        lxg.plotOn(frame)

    argset_fit = ROOT.RooArgSet(gauss_mean,gauss_sigma)
    lxg.paramOn(frame,ROOT.RooFit.Format("NELU",ROOT.RooFit.AutoPrecision(1)),ROOT.RooFit.Layout(0.55)) 
    frame.SetMaximum(frame.GetMaximum()*1.2)

    # add chi2 info
    chi2_text = ROOT.TPaveText(0.3,0.8,0.4,0.9,"BRNDC")
    chi2_text.AddText("#chi^{2} fit = %s" %round(frame.chiSquare(6),2))
    chi2_text.SetTextSize(0.04)
    chi2_text.SetTextColor(2)
    chi2_text.SetShadowColor(0)
    chi2_text.SetFillColor(0)
    chi2_text.SetLineColor(0)
    frame.addObject(chi2_text)

    if fit_filename is not None:
        c = ROOT.TCanvas("cfit","cfit",600,700)
        frame.Draw()
        if not os.path.exists(fit_plot_directory): os.makedirs(fit_plot_directory)
        c.SaveAs(os.path.join( fit_plot_directory, fit_filename+".pdf"))
        c.SaveAs(os.path.join( fit_plot_directory, fit_filename+".png"))
        del c

    mean_asymmetry        = landau_mean.getVal()
    mean_asymmetry_error  = landau_mean.getError()

    rms_asymmetry        = gauss_sigma.getVal()
    rms_asymmetry_error  = gauss_sigma.getError()

    #if the chi2 is going crazy default to the mean and rms of the histogram
    if frame.chiSquare(6)  > -1:
        mean_asymmetry = tmp_mean
        mean_asymmetry_error = tmp_mean_error
        rms_asymmetry  = tmp_sigma
        rms_asymmetry_error  = tmp_sigma_error

    return mean_asymmetry, mean_asymmetry_error, rms_asymmetry, rms_asymmetry_error

def fraction2(t_in,charged,neutral,photon,electron,muon,hhf,ehf,pt1,pt2,eta1,eta2,phi1, phi2, jet):

    cut_gen = "("+jet+"genjet_pt>0)"
    cut = "genjet_pt>"+str(pt1)+" && ""genjet_pt<"+str(pt2)+" && "+cut_gen+" && abs(genjet_eta)>="+str(eta1)+" && abs(genjet_eta)<"+str(eta2) + "&& genjet_phi >="+str(phi1) + " && genjet_phi<"+str(phi2)

    #phi_range = [float(phi1), float(phi2)]

    #h_sum      = ROOT.TProfile("hsum","hsum",len(phi_range),phi_range)
    #h_charged  = ROOT.TProfile("hch","hch",len(phi_range),phi_range)
    
    h_sum      = ROOT.TH1F("hsum","hsum",1,0,1)
    h_charged  = ROOT.TH1F("hch" ,"hch",1,0,1)

    h_neutral  = ROOT.TH1F("hnh" ,"hnh",1,0,1)
    h_photon   = ROOT.TH1F("hph" ,"hph",1,0,1)
    h_muon     = ROOT.TH1F("hmu" ,"hmu",1,0,1)
    h_electron = ROOT.TH1F("hel" ,"hel",1,0,1)
    h_hhf      = ROOT.TH1F("hhhf","hhhf",1,0,1)
    h_ehf      = ROOT.TH1F("hehf","hehf",1,0,1)


    sum_jet = "("+charged +"+"+ neutral +"+"+ photon +"+"+ electron +"+"+ muon +"+"+ hhf +"+"+ ehf +")"

    t_in.Draw("0.5>>hsum",
              "("+sum_jet+")*("+cut+")","goff")    
    t_in.Draw("0.5>>hch",
              "(charged_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hnh",
              "(neutral_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hph",
              "(photon_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hmu",
              "(muon_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hel",
              "(electron_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hhhf",
              "(hhf_e)*("+cut+")","goff")
    t_in.Draw("0.5>>hehf",
              "(ehf_e)*("+cut+")","goff")

    print h_charged.Integral(), h_charged.Integral()/h_sum.Integral(), h_sum.Integral()
    return h_charged.GetBinContent(1)/h_sum.GetBinContent(1), h_neutral.GetBinContent(1)/h_sum.GetBinContent(1), h_photon.GetBinContent(1)/h_sum.GetBinContent(1), h_muon.GetBinContent(1)/h_sum.GetBinContent(1), h_electron.GetBinContent(1)/h_sum.GetBinContent(1), h_hhf.GetBinContent(1)/h_sum.GetBinContent(1), h_ehf.GetBinContent(1)/h_sum.GetBinContent(1)
    



def fraction(t_in,charged,neutral,photon,electron,muon,hhf,ehf,pt1,pt2,eta1,eta2,phi1, phi2, jet):
    
    h_charged  = ROOT.TH1F("hch" ,"hch",100,0,1)
    h_neutral  = ROOT.TH1F("hnh" ,"hnh",100,0,1)
    h_photon   = ROOT.TH1F("hph" ,"hph",100,0,1)
    h_muon     = ROOT.TH1F("hmu" ,"hmu",100,0,1)
    h_electron = ROOT.TH1F("hel" ,"hel",100,0,1)
    h_hhf      = ROOT.TH1F("hhhf","hhhf",100,0,1)
    h_ehf      = ROOT.TH1F("hehf","hehf",100,0,1)

    sum_jet = "("+charged +"+"+ neutral +"+"+ photon +"+"+ electron +"+"+ muon +"+"+ hhf +"+"+ ehf +")"

    cut_gen = "("+jet+"genjet_pt>0)"
    cut = "genjet_pt>"+str(pt1)+" && ""genjet_pt<"+str(pt2)+" && "+cut_gen+" && abs(genjet_eta)>="+str(eta1)+" && abs(genjet_eta)<"+str(eta2) + "&& genjet_phi >="+str(phi1) + " && genjet_phi<"+str(phi2)

    t_in.Draw(charged+"/"+sum_jet+">>hch",
              cut,"goff")

    t_in.Draw(neutral+"/"+sum_jet+">>hnh",
              "genjet_pt>"+str(pt1)+" && ""genjet_pt<"+str(pt2)+" && "+cut_gen+ "&& abs(genjet_eta)>="+str(eta1)+" && abs(genjet_eta)<"+str(eta2),"goff")

    t_in.Draw(photon+"/"+sum_jet+">>hph",
              cut,"goff")

    t_in.Draw(muon+"/"+sum_jet+">>hmu",
              cut,"goff")

    t_in.Draw(electron+"/"+sum_jet+">>hel",
              cut,"goff")

    t_in.Draw(hhf+"/"+sum_jet+">>hhhf",
              cut,"goff")

    t_in.Draw(ehf+"/"+sum_jet+">>hehf",
              cut,"goff")

    fout = ROOT.TFile("fractions.root","recreate")
    fout.cd()
    h_charged.Write()
    h_neutral.Write()
    h_photon.Write()
    fout.Close()

    median   = n.zeros(1,dtype=float)
    quantile = n.zeros(1,dtype=float)        
    quantile[0] = 0.5
    
    h_electron.GetQuantiles(1, median, quantile)    
    print median[0]

    print h_charged.GetMean(), h_neutral.GetMean(), h_photon.GetMean(), h_muon.GetMean(), h_electron.GetMean(), h_hhf.GetMean(), h_ehf.GetMean(), h_charged.GetMean()+ h_neutral.GetMean()+ h_photon.GetMean()+ h_muon.GetMean()+ h_electron.GetMean()+ h_hhf.GetMean()+ h_ehf.GetMean()

    return h_charged.GetMean(), h_neutral.GetMean(), h_photon.GetMean(), h_muon.GetMean(), h_electron.GetMean(), h_hhf.GetMean(), h_ehf.GetMean()
    
