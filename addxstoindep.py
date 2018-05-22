from ROOT import *
from array import *
import sys
import math
import json

def getXsec(mass,process):
  
  #kp = float(model[1]+"."+model[2])*float(model[1]+"."+model[2])
  #brn = float((model[6]+"."+model[7]))

  #sf = kp*(1-brn)
#  print "ewk singlet SF = ",sf

  my_list = ['200','210','230' ,'250', '270']

  if process=="ggH":


   if mass  in my_list:
    rootFile = TFile("/eos/cms/store/group/phys_higgs/cmshww/amassiro/Full2016_Apr17/Apr2017_summer16/lepSel__MCWeights__bSFLpTEffMulti__cleanTauMC__l2loose__hadd__l2tightOR__LepTrgFix__dorochester__formulasMC__wwSel/latino_GluGluHToWWTo2L2Nu_M"+mass+".root")

   else:

    rootFile = TFile("/eos/cms/store/group/phys_higgs/cmshww/amassiro/Full2016_Apr17/Apr2017_summer16/lepSel__MCWeights__bSFLpTEffMulti__cleanTauMC__l2loose__hadd__l2tightOR__LepTrgFix__dorochester__formulasMC__wwSel/latino_GluGluHToWWTo2L2Nu_JHUGen698_M"+mass+".root")

   latino = rootFile.Get("latino")
   latino.GetEntry(0)
   xsec = latino.Xsec

   print "M = ", mass, " ggH xsec = ", xsec

#    print "xsec*sf = ", xsec*sf
#    return xsec*sf ############################ FIXME

   return xsec

  elif process=="VBF":




   if mass  in my_list:
     rootFile = TFile("/eos/cms/store/group/phys_higgs/cmshww/amassiro/Full2016_Apr17/Apr2017_summer16/lepSel__MCWeights__bSFLpTEffMulti__cleanTauMC__l2loose__hadd__l2tightOR__LepTrgFix__dorochester__formulasMC__wwSel/latino_VBFHToWWTo2L2Nu_M"+mass+".root")

   else:
    rootFile = TFile("/eos/cms/store/group/phys_higgs/cmshww/amassiro/Full2016_Apr17/Apr2017_summer16/lepSel__MCWeights__bSFLpTEffMulti__cleanTauMC__l2loose__hadd__l2tightOR__LepTrgFix__dorochester__formulasMC__wwSel/latino_VBFHToWWTo2L2Nu_JHUGen698_M"+mass+".root")

   latino = rootFile.Get("latino")
   latino.GetEntry(0)
   xsec = latino.Xsec

   print "M = ", mass, " VBF xsec = ", xsec
#    print "xsec*sf = ", xsec*sf
   
#    return xsec*sf
   return xsec



file = open('/afs/cern.ch/work/d/dmroy/CMSSW_7_4_7/src/CombineHarvester/HIG16006/output/indep.json')
json_data = json.load(file)
file.close()

thedict = {'3000':7.7999997302e-05+0.000920000020415, '2500':0.000159999995958+0.00150000001304, '650':0.0800471380353+0.0160982366651, '800':0.0278153866529+0.0100459931418, '1000':0.00794733408839+0.0056282337755, '230':1.05834913254+0.139826774597, '1500':0.000959999975748+0.00289999996312, '200':1.40569996834+0.177031323314, '270':0.813925147057+0.108498997986, '450':0.388092696667+0.0346079431474, '2000':0.000410000007832+0.00230000005104, '350':0.725226700306+0.0683906450868, '210':1.26020014286+0.161820292473, '300':0.712544858456+0.0911082476377, '550':0.17240922153+0.0227834656835, '250':0.918380439281+0.122818663716, '900':0.0145896067843+0.00748184928671, '700':0.0556466244161+0.013708290644, '400':0.575396120548+0.0458333604038, '500':0.257722169161+0.0276690702885, '750':0.0390892401338+0.0116998897865}

for mass in json_data:
  m = mass[:-2]
  #print mass, ":", m
  #if m in ['600','850','950']: continue
  if m not in thedict:
    xsec=getXsec(m,"ggH")+getXsec(m,"VBF")
  else:
    xsec=thedict[m]
  for key in json_data[mass].keys():
    #print key
    #print json_data[mass][key]
    json_data[mass][key] = json_data[mass][key] * xsec
#with ('/afs/cern.ch/work/d/dmroy/CMSSW_7_4_7/src/CombineHarvester/HIG16006/output/indep.json', 'w') as file:
file = open('/afs/cern.ch/work/d/dmroy/CMSSW_7_4_7/src/CombineHarvester/HIG16006/output/indep.json','w')
json.dump(json_data, file)

print "Done!"
exit(0)
