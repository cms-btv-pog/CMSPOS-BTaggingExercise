import ROOT

f = ROOT.TFile.Open("qcd.root","read")
hData = f.Get("h_pt100t140_pretag__DATA")
hB = f.Get("h_pt100t140_pretag__b")
hC = f.Get("h_pt100t140_pretag__c")
hL = f.Get("h_pt100t140_pretag__udsg")

hB.SetMarkerColor(ROOT.kRed)
hB.SetLineColor(ROOT.kRed)
hB.SetFillColor(ROOT.kRed)
hB.SetMarkerSize(0)

hC.SetMarkerColor(ROOT.kGreen)
hC.SetLineColor(ROOT.kGreen)
hC.SetFillColor(ROOT.kGreen)
hC.SetMarkerSize(0)

hL.SetMarkerColor(ROOT.kBlue)
hL.SetLineColor(ROOT.kBlue)
hL.SetFillColor(ROOT.kBlue)
hL.SetMarkerSize(0)

mc = ROOT.TObjArray(3)
mc.Add(hL)
mc.Add(hC)
mc.Add(hB)
fit = ROOT.TFractionFitter(hData,mc)

fit.Fit()

mcL_frac = ROOT.Double()
mcL_frac_err = ROOT.Double()
mcC_frac = ROOT.Double()
mcC_frac_err = ROOT.Double()
mcB_frac = ROOT.Double()
mcB_frac_err = ROOT.Double()

fit.GetResult(0,mcL_frac,mcL_frac_err)
fit.GetResult(1,mcC_frac,mcC_frac_err)
fit.GetResult(2,mcB_frac,mcB_frac_err)

hL.Scale(mcL_frac*hData.Integral()/hL.Integral())
hC.Scale(mcC_frac*hData.Integral()/hC.Integral())
hB.Scale(mcB_frac*hData.Integral()/hB.Integral())

c1 = ROOT.TCanvas("c1","c1")
hMC = ROOT.THStack("MC","pretag")
hMC.Add(hL)
hMC.Add(hC)
hMC.Add(hB)
hMC.Draw("hist e1")
hData.Draw("e1 same")
hMC.SetMaximum(hData.GetMaximum()*1.2);
hMC.GetYaxis().SetTitle("Events")
hMC.GetXaxis().SetTitle("Jet Probability discriminator")
leg = ROOT.TLegend(0.6,0.6,0.8,0.8)
leg.AddEntry(hData,"Data","lp")
leg.AddEntry(hB,"b","f")
leg.AddEntry(hC,"c","f")
leg.AddEntry(hL,"udsg","f")
leg.Draw()

c1.Print("TFractionFitter.eps","")
c1.Print("TFractionFitter.png","")
