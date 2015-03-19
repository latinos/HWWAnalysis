import ROOT
from array import array
import os, re
from math import *

def file2mapold(file):
    map = {}
    print file
    for line in open(file,"r"):
        fields = [x.strip() for x in line.split()]
        if len(fields) > 1 and re.match("\d+", fields[0]):
            map[float(fields[0])] = [float(y) for y in fields[1:] if re.match(r"-?\d+.*",y) ]
    return map

def file2map(x):
        ret = {}; headers = []
        for x in open(x,"r"):
            cols = x.split()
            if len(cols) < 2: continue
            if "mH" in x:
                headers = [i.strip() for i in cols[1:]]
            else:
                fields = [ float(i) for i in cols ]
                #ret[fields[0]] = dict(zip(headers,fields[1:]))
                ret[fields[0]] = [float(y) for y in fields[1:]]
        return ret

def juerg2map(file):
    map = {}
    for line in open(file,"r"):
        if line[0] == '#': continue
        fields = [x.strip() for x in line.split()]
        map[fields[0]] = [float(y) for y in fields[1:]]
    return map

YR_ggH         = {}
YR_vbfH        = {}
YR_wzttH       = {}
YR_wH          = {}
YR_zH          = {}
YR_ttH         = {}

ggH_pdfErrYR   = {}
ggH_scaErrYR   = {}
vbfH_pdfErrYR  = {}
vbfH_scaErrYR  = {}
wzttH_pdfErrYR = {}
wzttH_scaErrYR = {}

wH_pdfErrYR    = {}
wH_scaErrYR    = {}

zH_pdfErrYR    = {}
zH_scaErrYR    = {}

ttH_pdfErrYR   = {}
ttH_scaErrYR   = {}

HWW_BR         = {}
HWW_BR_vals    = {}

SYST_PATH      = os.getenv('CMSSW_BASE')+'/src/HWWAnalysis/ShapeAnalysis/data/nuisances/'
ggH_jets = dict([(m, dict(zip(['f0','f1','f2','k1','k2'], vals))) for m,vals in file2map(SYST_PATH+"ggH_jetBins.txt").items()]) 
ggH_jets2 = dict([(m, dict(zip(['0','1in0','1in1','2in1','2in2'], vals))) for m,vals in file2map(SYST_PATH+"ggH_jetBinNuisances.txt").items()]) 
ggH_UEPS = dict([(m, dict(zip(['u0','u1','u2'], vals))) for m,vals in file2map(SYST_PATH+"ggH_UEPS.txt").items()])
ggH_intf = dict([(m, dict(zip(['intf'], vals))) for m,vals in file2map(SYST_PATH+"ggH_interference.txt").items()])

#                             BR uncertainty
#HWW_BR = dict([(m, dict(zip(['brunc'], vals))) for m,vals in file2map(SYST_PATH+"HWW_BR.txt").items()])




def loadYRSyst(YRVersion=3,Energy='8TeV') :

    global YR_ggH   
    global YR_vbfH  
    global YR_wzttH  
    global YR_wH      
    global YR_zH       
    global YR_ttH       

    global HWW_BR
    global HWW_BR_vals

    global ggH_pdfErrYR  
    global ggH_scaErrYR   
    global vbfH_pdfErrYR  
    global vbfH_scaErrYR  
    global wzttH_pdfErrYR 
    global wzttH_scaErrYR 

    global wH_pdfErrYR    
    global wH_scaErrYR

    global zH_pdfErrYR 
    global zH_scaErrYR  

    global ttH_pdfErrYR  
    global ttH_scaErrYR   

    if YRVersion == 1 :
      print '\033[0;31m FIXME:  Old Yellow Report in use, 7TeV only!!!!!  \033[m'
      YR_PATH      = os.getenv('CMSSW_BASE')+'/src/HWWAnalysis/ShapeAnalysis/data/nuisances/'
      YR_ggH         = file2map(YR_PATH+'YR-ggH.txt')
      YR_vbfH        = file2map(YR_PATH+'YR-vbfH.txt')
      YR_wzttH       = file2map(YR_PATH+'YR-wzttH.txt')
      YR_wH          = file2map(YR_PATH+'YR-wH.txt')
      YR_zH          = file2map(YR_PATH+'YR-zH.txt')
      YR_ttH         = file2map(YR_PATH+'YR-ttH.txt')
    else:
      YR_PATH      = os.getenv('CMSSW_BASE')+'/src/HWWAnalysis/ShapeAnalysis/data/lhc-hxswg-YR'+str(YRVersion)+'/sm/xs/'+Energy+'/'
      YR_ggH         = file2map(YR_PATH+Energy+'-ggH.txt')
      YR_vbfH        = file2map(YR_PATH+Energy+'-vbfH.txt')
#     YR_wzttH       = file2map(YR_PATH+Energy+'-wzttH.txt')
      YR_wH          = file2map(YR_PATH+Energy+'-WH.txt')
      YR_zH          = file2map(YR_PATH+Energy+'-ZH.txt')
      YR_ttH         = file2map(YR_PATH+Energy+'-ttH.txt')

    # Only available in YR3
    HWW_BR         = file2map(os.getenv('CMSSW_BASE')+'/src/HWWAnalysis/ShapeAnalysis/data/lhc-hxswg-YR3/sm/br/BR2bosons.txt')
    HWW_BR_vals    = dict([(m, (abs(Err_Hi_WW)+abs(Err_Lo_WW))/2.) for m,(H_gg,Err_Hi_gg,Err_Lo_gg,H_gamgam,Err_Hi_gamgam,Err_Lo_gamgam,H_Zgam,Err_Hi_Zgam,Err_Lo_Zgam,H_WW,Err_Hi_WW,Err_Lo_WW,H_ZZ,Err_Hi_ZZ,Err_Lo_ZZ,Total_Width_GeV,Err_Hi,Err_Lo) in HWW_BR.items()] )

    ggH_pdfErrYR   = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )
    ggH_scaErrYR   = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )
    vbfH_pdfErrYR  = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_vbfH.items()] )
    vbfH_scaErrYR  = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_vbfH.items()] )
    wzttH_pdfErrYR = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_wzttH.items()] )
    wzttH_scaErrYR = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_wzttH.items()] )
    
    wH_pdfErrYR = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_wH.items()] )
    wH_scaErrYR = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_wH.items()] )
    
    zH_pdfErrYR = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_zH.items()] )
    zH_scaErrYR = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_zH.items()] )
    
    ttH_pdfErrYR = dict([(m, sqrt((1+0.01*pdf_hi)/(1+0.01*pdf_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ttH.items()] )
    ttH_scaErrYR = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ttH.items()] )
    
    if YRVersion == 1 :
      for X in 450, 550:
        ggH_pdfErrYR[X]  = 0.5*(ggH_pdfErrYR[X-10] +ggH_pdfErrYR[X+10])
        ggH_scaErrYR[X]  = 0.5*(ggH_scaErrYR[X-10] +ggH_scaErrYR[X+10])
        vbfH_pdfErrYR[X] = 0.5*(vbfH_pdfErrYR[X-10]+vbfH_pdfErrYR[X+10])
        vbfH_scaErrYR[X] = 0.5*(vbfH_scaErrYR[X-10]+vbfH_scaErrYR[X+10])

# ----------------------------------------------------- Add point if needed to YR ------------------------

def GetYRVal(YRDic,iMass):

    if iMass in YRDic :
       #print iMass,YRDic[iMass][Key]
       return YRDic[iMass]
    else:
       n=len(YRDic.keys())
       x=[]
       y=[]
       for jMass in sorted(YRDic.keys()):
         x.append(jMass)
         y.append(YRDic[jMass])
       gr = ROOT.TGraph(n,array('f',x),array('f',y));
       sp = ROOT.TSpline3("YR",gr);
       #print iMass,sp.Eval(iMass)
       return sp.Eval(iMass)


def getCommonSysts(mass,channel,jets,qqWWfromData,shape,options,suffix,isssactive,Energy,newInterf=False,YRVersion=3,mh_SM=125.,mh_SM2=125.,ewksinglet=False ):

    loadYRSyst(YRVersion,Energy)

    nuisances = {} 
    #MCPROC = ['ggH', 'vbfH', 'DTT', 'ggWW', 'VV', 'Vg' ]; 
    MCPROC = ['ggH', 'qqH', 'wzttH', 'WH', 'ZH', 'ttH', 'DYTT', 'VV', 'VgS', 'Vg', 'Other', 'VVV', 'WWewk', 'ggH_SM', 'qqH_SM', 'wzttH_SM' , 'WH_SM', 'ZH_SM', 'ttH_SM' ];
    MCPROC+=['Top']
    MCPROC+=['ggH_sbi', 'ggH_b', 'ggH_s'] # for Higgs width
    MCPROC+=['qqH_sbi', 'qqH_b', 'qqH_s'] # for Higgs width
    if channel == 'elmu' or channel == 'muel': MCPROC+=['DYMM','DYEE']
    if channel == 'of': MCPROC += ['DYLL']
    if not qqWWfromData: MCPROC+=['WW','ggWW']
    # -- Luminosity ---------------------
    lumiunc = 1.026
    if '7TeV' in suffix: lumiunc = 1.022
    nuisances['lumi'+suffix] = [ ['lnN'], dict([(p,lumiunc) for p in MCPROC if p!='DYTT' and p!='Top'])]
    # -- PDF ---------------------
    #nuisances['pdf_gg']    = [ ['lnN'], { 'ggH':ggH_pdfErrYR[mass], 'ggWW':(1.00 if qqWWfromData else 1.04) }]
    nuisances['pdf_gg']    = [ ['lnN'], { 'ggH'    : GetYRVal(ggH_pdfErrYR,mass), 
                                          'ggH_SM' : GetYRVal(ggH_pdfErrYR,mh_SM),
                                          'ggWW'   : 1.04 ,  # 4% uncertainty by Xavier and Guillelmo studies
                                          # for Higgsw width
                                          'ggH_sbi'   : 1.04 ,
                                          'ggH_b'     : 1.04 ,
                                          'ggH_s'     : 1.04 ,
                                        }]

    nuisances['pdf_qqbar'] = [ ['lnN'], { #'wzttH' :(1.0 if mass>300 else wzttH_pdfErrYR[mass]),  
                                          'WH'    :(1.0 if mass>300 else GetYRVal(wH_pdfErrYR,mass)),  
                                          'ZH'    :(1.0 if mass>300 else GetYRVal(zH_pdfErrYR,mass)), 
                                          'ttH'   :(1.0 if mass>300 else GetYRVal(ttH_pdfErrYR,mass)), 
                                          'qqH'   :GetYRVal(vbfH_pdfErrYR,mass), 
                                          #'wzttH_SM' :(1.0 if mh_SM>300 else wzttH_pdfErrYR[mh_SM]),  
                                          'WH_SM'    :(1.0 if mh_SM>300 else GetYRVal(wH_pdfErrYR,mh_SM)),  
                                          'ZH_SM'    :(1.0 if mh_SM>300 else GetYRVal(zH_pdfErrYR,mh_SM)), 
                                          'ttH_SM'   :(1.0 if mh_SM>300 else GetYRVal(ttH_pdfErrYR,mh_SM)), 
                                          'qqH_SM'   :GetYRVal(vbfH_pdfErrYR,mh_SM), 
                                          'VV':1.04, 
                                          'WW':(1.0 if qqWWfromData else 1.04),  # 4% uncertainty by Xavier and Guillelmo studies
                                          # for Higgsw width
                                          'qqH_sbi'   : 1.04 ,
                                          'qqH_b'     : 1.04 ,
                                          'qqH_s'     : 1.04 ,

                                        }]

    if options.WWxsec :
      nuisances['pdf_qqbar'][1]['WW'] = 1.013
      nuisances['pdf_gg'][1]['ggWW'] = 1.007

    #nuisances['pdf_qqbar_ACCEPT'] = [ ['lnN'], {
                                          #'WW':    1.013,
    #}

    #nuisances['pdf_gg_ACCEPT'] = [ ['lnN'], {
                                          #'WW':    1.008
    #}


    # -- Theory ---------------------
    #nuisances['QCDscale_ggH_offshell']    = [  ['lnN'], { 'ggH_sbi':1.15,  'ggH_b':1.15,  'ggH_s':1.15 }]
    if options.HWidth :
      if jets == 0:
        #nuisances['QCDscale_ggWW_0jet'] = [ ['lnN'], {'ggH_sbi': 1.11, 'ggH_b': 1.11, 'ggH_s': 1.11}]
        #nuisances['QCDscale_ggWW_offshell']    = [ ['lnN'], {'ggH_sbi': 1.15, 'ggH_b': 1.15, 'ggH_s': 1.15}]
        nuisances['QCDscale_ggWW_offshell']    = [ ['lnN'], {'ggH_sbi': 1.75158,  'ggH_b': 1.75158,  'ggH_s': 1.75158}]  # using NLO uncertainty for ST recipe
        nuisances['QCDscale_ggWW_offshell1in'] = [ ['lnN'], {'ggH_sbi': 0.796502, 'ggH_b': 0.796502, 'ggH_s': 0.796502}]
        #nuisances['QCDscale_qqWW_offshell0j']  = [ ['lnN'], {'qqH_sbi': 1.412499,  'qqH_b': 1.412499,  'qqH_s': 1.412499}]  # comparison phantom vs powheg
        nuisances['QCDscale_qqWW_offshell']  = [ ['lnN'], {'qqH_sbi': 1.412499,  'qqH_b': 1.412499,  'qqH_s': 1.412499}]  # comparison phantom vs powheg
      if jets == 1:
        #nuisances['QCDscale_ggWW_1jet'] = [ ['lnN'], {'ggH_sbi': 1.11, 'ggH_b': 1.11, 'ggH_s': 1.11}]
        nuisances['QCDscale_ggWW_offshell1in'] = [ ['lnN'], {'ggH_sbi': 0.676083, 'ggH_b': 0.676083, 'ggH_s': 0.676083}]
        nuisances['QCDscale_ggWW_offshell2in'] = [ ['lnN'], {'ggH_sbi': 0.996433, 'ggH_b': 0.996433, 'ggH_s': 0.996433}]
        #nuisances['QCDscale_qqWW_offshell1j']  = [ ['lnN'], {'qqH_sbi': 1.035659,  'qqH_b': 1.035659,  'qqH_s': 1.035659}]  # comparison phantom vs powheg
        nuisances['QCDscale_qqWW_offshell']  = [ ['lnN'], {'qqH_sbi': 1.035659,  'qqH_b': 1.035659,  'qqH_s': 1.035659}]  # comparison phantom vs powheg
      if jets == 2:
        #nuisances['QCDscale_ggWW_2jet'] = [ ['lnN'], {'ggH_sbi': 1.21, 'ggH_b': 1.21, 'ggH_s': 1.21}]
        nuisances['QCDscale_ggWW_offshell2in'] = [ ['lnN'], {'ggH_sbi': 1.51302, 'ggH_b': 1.51302, 'ggH_s': 1.51302}]
        nuisances['QCDscale_qqWW_offshell']  = [ ['lnN'], {'qqH_sbi': 0.9020033,  'qqH_b': 0.9020033,  'qqH_s': 0.9020033}]  # comparison phantom vs powheg (1-0.9020033 + 1)
        #nuisances['QCDscale_qqWW_offshell2j']  = [ ['lnN'], {'qqH_sbi': 1.097997,  'qqH_b': 1.097997,  'qqH_s': 1.097997}]  # comparison phantom vs powheg (1-0.9020033 + 1)

      #nuisances['QCDscale_qqWW_offshell']    = [ ['lnN'], {'qqH_sbi': 1.12,  'qqH_b': 1.12,  'qqH_s': 1.12}]  # using comparison phantom <-> powheg for central jet veto

      #QCDscale_ggH     1.75158    -          -      << with "YR" uncertainty --> NNLO ~ 15%
      #QCDscale_ggH1in  0.796502   0.676083   -
      #QCDscale_ggH2in  -          0.996433   1.51302


    if jets == 0:
        # appendix D of https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=0&confId=135333
        k0 = pow(GetYRVal(ggH_scaErrYR,mass),     1/ggH_jets[mass]['f0'])
        k1 = pow(ggH_jets[mass]['k1'], 1-1/ggH_jets[mass]['f0']) # -f1-f2=f0-1
        # ... and _SM:
        k0_SM = pow(GetYRVal(ggH_scaErrYR,mh_SM),    1/ggH_jets[mh_SM2]['f0'])
        k1_SM = pow(ggH_jets[mh_SM2]['k1'], 1-1/ggH_jets[mh_SM2]['f0']) # -f1-f2=f0-1
        #nuisances['QCDscale_ggH']    = [  ['lnN'], { 'ggH':k0, 'ggH_SM':k0_SM }]
        #nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':k1, 'ggH_SM':k1_SM }]
        nuisances['QCDscale_ggH']    = [  ['lnN'], { 'ggH':ggH_jets2[mass]['0'],    'ggH_SM':ggH_jets2[mh_SM2]['0'] }]
        nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':ggH_jets2[mass]['1in0'], 'ggH_SM':ggH_jets2[mh_SM2]['1in0'] }]
        if not qqWWfromData:
            nuisances['QCDscale_WW']    = [ ['lnN'], {'WW': 1.035, 'ggWW': 1.035 }]
            nuisances['QCDscale_WW1in'] = [ ['lnN'], {'WW': 0.987, 'ggWW': 0.987 }]
            #nuisances['QCDscale_WW']    = [ ['lnN'], {'WW': 1.042, 'ggWW': 1.042 }]  --> new values above from WW paper
            #nuisances['QCDscale_WW1in'] = [ ['lnN'], {'WW': 0.974, 'ggWW': 0.974 }]
            #nuisances['QCDscale_WW1in'] = [ ['lnN'], {'WW': 0.978, 'ggWW': 0.978 }] --> new value!
    elif jets == 1:
        k1 = pow(ggH_jets[mass]['k1'], 1+ggH_jets[mass]['f2']/ggH_jets[mass]['f1'])
        k2 = pow(ggH_jets[mass]['k2'],  -ggH_jets[mass]['f2']/ggH_jets[mass]['f1'])
        # ... and _SM:
        k1_SM = pow(ggH_jets[mh_SM2]['k1'], 1+ggH_jets[mh_SM2]['f2']/ggH_jets[mh_SM2]['f1'])
        k2_SM = pow(ggH_jets[mh_SM2]['k2'],  -ggH_jets[mh_SM2]['f2']/ggH_jets[mh_SM2]['f1'])
        #nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':k1, 'ggH_SM':k1_SM }]
        #nuisances['QCDscale_ggH2in'] = [  ['lnN'], { 'ggH':k2, 'ggH_SM':k2_SM }]
        nuisances['QCDscale_ggH1in'] = [  ['lnN'], { 'ggH':ggH_jets2[mass]['1in1'], 'ggH_SM':ggH_jets2[mh_SM2]['1in1'] }]
        nuisances['QCDscale_ggH2in'] = [  ['lnN'], { 'ggH':ggH_jets2[mass]['2in1'], 'ggH_SM':ggH_jets2[mh_SM2]['2in1'] }]
        if not qqWWfromData:
            nuisances['QCDscale_WW1in'] = [ ['lnN'], {'WW': 1.076, 'ggWW': 1.076 }]
            nuisances['QCDscale_WW2in'] = [ ['lnN'], {'WW': 0.914, 'ggWW': 0.914 }]
    elif jets == 2:
        if options.VH:
            nuisances['QCDscale_ggH2in_vh'] = [  ['lnN'], { 'ggH':1.30, 'ggH_SM':1.30 }]
            #nuisances['QCDscale_ggH2in_vh'] = [  ['lnN'], { 'ggH':1.70 }]
        else :
            nuisances['QCDscale_ggH2in_vbf'] = [  ['lnN'], { 'ggH':1.35, 'ggH125':1.35, 'ggH_SM':1.35}]
            #nuisances['QCDscale_ggH2in'] = [  ['lnN'], { 'ggH':ggH_jets[mass]['k2'], 'ggH_SM':ggH_jets[mh_SM]['k2'] }]
        if not qqWWfromData:
            nuisances['QCDscale_WW2in'] = [ ['lnN'], {'WW': 1.210, 'ggWW': 1.210 }] # reduce by 1/2 because not applicable to vbf
            if not options.VH: 
               #--> now we have Phantom WW+2jets ewk sample, no need this nuisance
               #nuisances['QCDscale_WWvbf'] = [ ['lnN'], {'WW': 1.500 }]
               nuisances['QCDscale_WWewk']     = [ ['lnN'], {'WWewk':1.20 }]

    nuisances['QCDscale_ggWW'] = [ ['lnN'], {'ggWW': 1.30}]
    nuisances['QCDscale_qqH']  = [ ['lnN'], { 'qqH':GetYRVal(vbfH_scaErrYR,mass) , 'qqH_SM':GetYRVal(vbfH_scaErrYR,mh_SM) }]

    nuisances['QCDscale_VH']  = [ ['lnN'] , {} ]
    nuisances['QCDscale_wH']  = [ ['lnN'] , {} ]
    nuisances['QCDscale_zH']  = [ ['lnN'] , {} ]
    nuisances['QCDscale_ttH'] = [ ['lnN'] , {} ]
    if mass in wzttH_scaErrYR: nuisances['QCDscale_VH']  = [ ['lnN'], { 'wzttH':wzttH_scaErrYR[mass]}]
    if mass in wH_scaErrYR :   nuisances['QCDscale_wH']  = [ ['lnN'], { 'WH'   :wH_scaErrYR[mass]}]
    if mass in zH_scaErrYR :   nuisances['QCDscale_zH']  = [ ['lnN'], { 'ZH'   :wH_scaErrYR[mass]}]
    if mass in ttH_scaErrYR:   nuisances['QCDscale_ttH'] = [ ['lnN'], {'ttH'   :wH_scaErrYR[mass]}]
    # ... and _SM:
    if mh_SM in wzttH_scaErrYR: nuisances['QCDscale_VH'][1] ['wzttH_SM']= wzttH_scaErrYR[mh_SM]
    if mh_SM in wH_scaErrYR :   nuisances['QCDscale_wH'][1] ['WH_SM']   = wH_scaErrYR[mh_SM]
    if mh_SM in zH_scaErrYR :   nuisances['QCDscale_zH'][1] ['ZH_SM']   = wH_scaErrYR[mh_SM]
    if mh_SM in ttH_scaErrYR:   nuisances['QCDscale_ttH'][1]['ttH_SM']  = wH_scaErrYR[mh_SM]
    if nuisances['QCDscale_VH'][1]  == {} : nuisances.pop('QCDscale_VH')
    if nuisances['QCDscale_wH'][1]  == {} : nuisances.pop('QCDscale_wH')
    if nuisances['QCDscale_zH'][1]  == {} : nuisances.pop('QCDscale_zH')
    if nuisances['QCDscale_ttH'][1] == {} : nuisances.pop('QCDscale_ttH')

    nuisances['QCDscale_VV']     = [ ['lnN'], { 'VV':1.03 }]
    nuisances['QCDscale_VgS']    = [ ['lnN'], {'VgS':1.30 }]
    nuisances['QCDscale_VVV']    = [ ['lnN'], {'VVV':1.50 }]

    if isssactive == True :
        # -- extrapolation from same sign (ss) to opposite sign (os) region ---- 10% ? -> to be checked on data: how much is the charge misidentification difference data/MC by CMS?
        # that it may be bin dependent, ... right? One nuisance for each jet bin
        if jets == 0 :
            nuisances['ssos_extrap_VgS_0jet']    = [ ['lnN'], {'VgS':1.10 }]
            nuisances['ssos_extrap_Vg_0jet']     = [ ['lnN'], {'Vg':1.10 }]

        if jets == 1 :
            nuisances['ssos_extrap_VgS_1jet']    = [ ['lnN'], {'VgS':1.10 }]
            nuisances['ssos_extrap_Vg_1jet']     = [ ['lnN'], {'Vg':1.10 }]

    # uncertainty on H125
    #nuisances['CMS_hww_SMH125']  = [ ['lnN'], {'ggH125':1.30, 'vbfH125':1.30, 'wzttH125':1.30} ]


    # -- Theoretical ---------------------
    nuisances['QCDscale_ggH_ACCEPT'] = [ ['lnN'], {'ggH':1.02,  'ggH_SM':1.02}]
    nuisances['QCDscale_qqH_ACCEPT'] = [ ['lnN'], {'qqH':1.02,  'qqH_SM':1.02}]
    nuisances['QCDscale_wH_ACCEPT']  = [ ['lnN'], {'WH':1.02 ,  'WH_SM':1.02}]
    nuisances['QCDscale_zH_ACCEPT']  = [ ['lnN'], {'ZH':1.02 ,  'ZH_SM':1.02}]
    nuisances['QCDscale_ttH_ACCEPT'] = [ ['lnN'], {'ttH':1.02,  'ttH_SM':1.02}]
    if   jets == 0: nuisances['UEPS'] = [ ['lnN'], {'ggH':ggH_UEPS[mass]['u0'], 'ggH_SM':ggH_UEPS[mh_SM2]['u0']}]
    elif jets == 1: nuisances['UEPS'] = [ ['lnN'], {'ggH':ggH_UEPS[mass]['u1'], 'ggH_SM':ggH_UEPS[mh_SM2]['u1']}]
    #elif jets == 2: nuisances['UEPS'] = [ ['lnN'], {'ggH':ggH_UEPS[mass]['u2'], 'ggH125':ggH_UEPS[mass]['u2']}]
    elif jets == 2:
      if options.VH:
          nuisances['UEPS_hww_vh']  = [ ['lnN'], {'ggH'   :1.05, 'WH'   :1.08, 'ZH'   :1.08, 'VH'   :1.08, 'qqH'   :1.04 ,
                                                  'ggH_SM':1.05, 'WH_SM':1.08, 'ZH_SM':1.08, 'VH_SM':1.08, 'qqH_SM':1.04}]
      else :
          nuisances['UEPS_hww_vbf'] = [ ['lnN'], {'ggH'   :1.20, 'qqH'   :1.10 ,
                                                  'ggH_SM':1.20, 'qqH_SM':1.10}]

    # UEPS for WW 0 jet
    if   jets == 0: nuisances['UEPS_WW'] = [ ['lnN'], {'WW':1.035}]   # see AN 2014/077

    # UEPS for offshell part
    #nuisances['UEPS_off'] = [ ['lnN'], {'ggH_sbi':1.10, 'ggH_b':1.10, 'ggH_s':1.10, 'qqH_sbi':1.10, 'qqH_b':1.10, 'qqH_s':1.10  }]
    if 'UEPS' in nuisances.keys() :
      nuisances['UEPS'][1].update( {'ggH_sbi':1.10, 'ggH_b':1.10, 'ggH_s':1.10, 'qqH_sbi':1.10, 'qqH_b':1.10, 'qqH_s':1.10  } )
    else :
      nuisances['UEPS'] = [ ['lnN'], {'ggH_sbi':1.10, 'ggH_b':1.10, 'ggH_s':1.10, 'qqH_sbi':1.10, 'qqH_b':1.10, 'qqH_s':1.10  }]


    #if ((not qqWWfromData) and (jets != 2)): nuisances['QCDscale_WW_EXTRAP'] = [ ['lnN'], {'WW':1.06}]
    # --- new ---
    # not needed with line-shape reweighting
    #nuisances['theoryUncXS_HighMH'] = [ ['lnN'], { 'ggH':1+1.5*pow(float(mass)/1000,3), 'vbfH':1+1.5*pow(float(mass)/1000,3), 'wzttH':1+1.5*pow(float(mass)/1000,3) } ]
    # -- Interference term ----------------
    # FIXME: Interference Error !!!!!
    if not newInterf :
      print '\033[0;31m FIXME: Interference Error !!!!!  \033[m'
      if mass>=400:
         nuisances['interf_ggH'] = [ ['lnN'], {'ggH':ggH_intf[mass]['intf']}] # --> high mass paper, new interference uncertainty (see below)
      else :
         nuisances['interf_ggH'] = [ ['lnN'], {'ggH':1.00}]
    else :
      if ewksinglet:
       if mass>=400:
         nuisances['interf_ggH125'] = [ ['lnN'], {'ggH':1.10, 'WW':1.15}]
       else :
         nuisances['interf_ggH125'] = [ ['lnN'], {'ggH':1.00}]

    # BR H > VV uncertainty
    BRunc    = 1.+GetYRVal(HWW_BR_vals,mass)/100.
    BRunc_SM = 1.+GetYRVal(HWW_BR_vals,mh_SM)/100.
    print "BR UNCERTAINIES: ",BRunc,BRunc_SM
    nuisances['BRhiggs_hvv'] = [ ['lnN'], {'ggH':BRunc, 'qqH':BRunc, 'ZH':BRunc , 'WH':BRunc , 'ggH_sbi':BRunc, 'ggH_b':BRunc, 'ggH_s':BRunc, 'qqH_sbi':BRunc, 'qqH_b':BRunc, 'qqH_s':BRunc , 'ggH_SM':BRunc_SM , 'qqH_SM':BRunc_SM , 'ZH_SM':BRunc_SM , 'WH_SM':BRunc_SM }]




    #if options.WJsub:
    #    nuisances['CMS_FakeRate_e'] = [ ['lnN'], { 'WJet': 1.0+options.WJsub } ]
    #    nuisances['CMS_FakeRate_m'] = [ ['lnN'], { 'WJet': 1.0+options.WJsub } ]
    #elif options.WJadd:
    #    nuisances['CMS_FakeRate_e'] = [ ['lnN'], { 'WJet': 1.0+options.WJadd } ]
    #    nuisances['CMS_FakeRate_m'] = [ ['lnN'], { 'WJet': 1.0+options.WJadd } ]
    addFakeRateSyst(nuisances, mass, channel, jets, shape, suffix)

    # FIXME : AM: "where are these numbers coming from?" 
    #if 'e' in channel:     nuisances['CMS'+suffix+'_eff_l'] = [ ['lnN'], dict([(p,pow(1.02,channel.count('e'))) for p in MCPROC])]
    ##elif channel == 'all': nuisances['CMS_eff_l'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC])]
    ##elif channel == 'sf':  nuisances['CMS_eff_l'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC])]
    ##elif channel == 'of':  nuisances['CMS_eff_l'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC])]
    #else :
        #nuisances['CMS'+suffix+'_eff_e'] = [ ['lnN'], dict([(p,1.04) for p in MCPROC])]
        #nuisances['CMS'+suffix+'_eff_m'] = [ ['lnN'], dict([(p,1.03) for p in MCPROC])]

    # just put a common one now
    if   channel == 'mumu': nuisances['CMS_p_scale_m'] = [ ['lnN'], dict([(p,1.015) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'elmu': nuisances['CMS_p_scale_m'] = [ ['lnN'], dict([(p,1.015) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'muel': nuisances['CMS_p_scale_e'] = [ ['lnN'], dict([(p,1.020) for p in MCPROC if p != 'DTT'] )]
    elif channel == 'elel': nuisances['CMS_p_scale_e'] = [ ['lnN'], dict([(p,1.020) for p in MCPROC if p != 'DTT'] )]
    elif channel in ['all', 'sf', 'of']: 
       nuisances['CMS'+suffix+'_p_scale_e'] = [ ['lnN'], dict([(p,1.020) for p in MCPROC if p != 'DTT'] )]
       nuisances['CMS'+suffix+'_p_scale_m'] = [ ['lnN'], dict([(p,1.015) for p in MCPROC if p != 'DTT'] )]
    nuisances['CMS'+suffix+'_met'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC if p!='DYTT'])]
    nuisances['CMS'+suffix+'_p_scale_j'] = [ ['lnN'], dict([(p,1.02) for p in MCPROC if p!='DYTT'])]
    if channel == 'of':
       nuisances['CMS'+suffix+'_norm_DYof'] = [ ['lnN'], { 'DYLL':2.0 } ]
    nuisances['CMS'+suffix+'_norm_DYTT'] = [ ['lnN'], { 'DYTT':1.3 } ]
    nuisances['CMS'+suffix+'_norm_Vg']   = [ ['lnN'], { 'Vg':1.3 } ]
    #nuisances['CMS'+suffix+'_norm_VVV']  = [ ['lnN'], { 'VVV':1.3 } ]

    return nuisances

def addFakeBackgroundSysts(nuisances, mass,channel,jets,errWW=0.2,errggWW=0.2,errDY=1.0,errTop0j=1.0,errTop1j=0.3,errWJ=0.5):
    if errWW:  nuisances['CMS_norm_WW'           ] = [ ['lnN'], { 'WW':(1+errWW)} ] 
    if errggWW:  nuisances['CMS_norm_ggWW'           ] = [ ['lnN'], { 'ggWW':(1+errggWW)} ] 
    if errTop0j and jets == 0: nuisances['CMS_norm_Top0j'] = [ ['lnN'], { 'Top':(1+errTop0j)}] 
    if errTop1j and jets == 1: nuisances['CMS_norm_Top1j'] = [ ['lnN'], { 'Top':(1+errTop1j)}] 
    if errWJ: nuisances['CMS_fake_%s'%channel[2]] = [ ['lnN'], { 'WJet':(1+errWJ) } ] 
    if errDY and channel in ['mumu', 'elel']:
        proc = "DY{0}{0}".format(channel[0].upper())
        nuisances['CMS_norm_DY_%s' % channel[0]]  = [ ['lnN'], { proc:(1+errDY) } ]

def addFakeRateSyst(nuisances, mass, channel, jets, shape, suffix=''):
    # some hard-coded numbers for now
    # define later in separate text + make it mH dependent for cut-based
    # do an asymmetric error, if possible...
    fake_m_up = 1.0
    fake_m_dn = 1.0
    fake_e_up = 1.0
    fake_e_dn = 1.0
    
    if (shape):
        if mass <= 300: 
            if jets == 0:
                fake_m_up = 0.75
                fake_m_dn = 1.18
                fake_e_up = 0.86
                fake_e_dn = 1.07
            elif jets == 1:
                fake_m_up = 0.78
                fake_m_dn = 1.13
                fake_e_up = 0.84
                fake_e_dn = 1.04
        else:
            if jets == 0:
                fake_m_up = 0.81
                fake_m_dn = 1.12
                fake_e_up = 0.73
                fake_e_dn = 1.02
            elif jets == 1:
                fake_m_up = 0.83
                fake_m_dn = 1.10
                fake_e_up = 0.76
                fake_e_dn = 1.03
    else:
        if jets == 0:
            fake_m_up = 0.69
            fake_m_dn = 1.23
            fake_e_up = 0.93
            fake_e_dn = 1.10
        elif jets == 1:
            fake_m_up = 0.71
            fake_m_dn = 1.18
            fake_e_up = 0.91
            fake_e_dn = 1.06

    nuisances['CMS'+suffix+'_hww_FakeRate']   = [ ['lnN'], { 'WJet': 1.15 } ] # from stat on FR (+) prompt contamination
    #nuisances['CMS'+suffix+'_hww_FakeRate']   = [ ['lnN'], { 'WJet': 1.20 } ] # from closure test
    # FR normalisation systematics are together in the shape variation
    #nuisances['CMS_hww_FakeRate_m'] = [ ['lnN'], { 'WJet': (fake_m_up,fake_m_dn) } ] # from jet ET variation
    #nuisances['CMS_hww_FakeRate_e'] = [ ['lnN'], { 'WJet': (fake_e_up,fake_e_dn) } ] # from jet ET variation


def floatNorm(process,jetcat):
    nuisances = {}
    if process in 'WW' and process != '':
        nuisances['CMS_norm_'+jetcat+'_'+process] = [ ['lnU'], { 'WW':2.00, 'ggWW':2.00 } ]
    elif process != '':
        nuisances['CMS_norm_'+jetcat+'_'+process] = [ ['lnU'], { process:2.00 } ]
    return nuisances

