import ROOT

p_b_pretag = 1 # put here the pretag fit parameter result for b jets from CFIT
p_b_tag = 1 # put here the tag fit parameter result for b jets from CFIT

f = ROOT.TFile.Open("qcd.root","read")
hB_pretag = f.Get("h_pt100t140_pretag__b")
hB_tag = f.Get("h_pt100t140_tag__b")

N_b_pretag = hB_pretag.Integral()
N_b_tag = hB_tag.Integral()

Effb_MC = N_b_tag/N_b_pretag
Effb_Data = N_b_tag*p_b_tag/(N_b_pretag*p_b_tag)

SFb = Effb_Data/Effb_MC

print "Effb_MC = ", Effb_MC
print "Effb_Data = ", Effb_Data
print "Sfb = ", SFb

