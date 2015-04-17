import re
import HWWAnalysis.Misc.odict as odict

#  ___                         _
# | _ \__ _ _ _ __ _ _ __  ___| |_ ___ _ _ ___
# |  _/ _` | '_/ _` | '  \/ -_)  _/ -_) '_(_-<
# |_| \__,_|_| \__,_|_|_|_\___|\__\___|_| /__/
#
# flavors             = ['mm', 'ee', 'em', 'me']
# flavors				 = dict( sf=['ee','mm'],of=['em','me'] )
# masses               = [ 110 , 115 , 118 , 120 , 122 , 124 , 126 , 128 , 130 , 135 , 140 , 150 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600]
# jets 				 = [0,1]

# HCP
# ptllCut = 45. 
# Spin/Moriond
ptllCut = 30. 

class wwcutsB:
    wwcommon = odict.OrderedDict([
        ('trigger',                'trigger==1.'),
        ('min #slash{E}_{T}',      'pfmet>20.'),
        ('m_{ll}',                 'mll>12'),                       # ema7
        ('Z_{veto}',               '(zveto==1||!sameflav)'),
        ('min proj #slash{E}_{T}', 'mpmet>20.'),                    # ema9
        ('Soft #mu veto',          'bveto_mu==1'),
        ('Extra Lepton',           'nextra==0'),
        ('B veto',                 '(bveto_ip==1 && (nbjettche==0 || njet>3)  )'),
        ('p_{T}^{ll}',             'ptll>%f'%ptllCut),                     # ema 14
        #('Extra Jet',              'njet<4'),
    ])

    wwmin = odict.OrderedDict([
        ('trigger',                'trigger==1.'),
        ('min #slash{E}_{T}',      'pfmet>20.'),
        ('m_{ll}',                 'mll>12'),                       # ema7
        ('min proj #slash{E}_{T}', 'mpmet>20.'),                    # ema9
        ('Soft #mu veto',          'bveto_mu==1'),
        ('Extra Lepton',           'nextra==0'),
        #('Extra Jet',              'njet<4'),
    ])

    #dy cuts
    dphilljj   = '(njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  )'

    # met cuts 
    met  = '( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'
    mpmet  = '( !sameflav || ( (njet>1 || (mpmet>45 && dphiveto)) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'

    wwcommon.update([
        ('d#Phi_{lljj}', dphilljj),
        ('dymva',        met),
    ])
    
    wwmin.update([
        ('d#Phi_{lljj}', dphilljj),
    ])
        
    zerojet = 'njet == 0'
    onejet  = 'njet == 1'
    vbf     = '(njet >= 2 && njet <= 3 && (jetpt3 <= 30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) '
    atleastonejet = 'njet >= 1'

class wwcuts:
    wwcommon = [
        'pt1>20',
        'pt2>10',
        '(ch1*ch2)<0', # use 0.5 to get charge symmetric template  
        'trigger==1.',
        'pfmet>20.',
        'mll>12',                       # ema7
        '(zveto==1||!sameflav)',
        'mpmet>20.',                    # ema9
        'bveto_mu==1',
        'nextra==0',
        '(bveto_ip==1 && nbjettche==0)',
        'ptll>%f'%ptllCut,                     # ema 14
    ]

    wwcommon2011 = [
        'pt1>20',
        'pt2>10',
        '((ch1*ch2))<0', # use 0.5 to get charge symmetric template  
        'trigger==1.',
        'pfmet>20.',
        'mll>(12+8*sameflav)',
        '(zveto==1||!sameflav)',
        'mpmet>(20+(17+nvtx/2.)*sameflav)',
        '(dphiveto || ! sameflav)',
        'bveto_mu==1',
        'nextra==0',
        '(bveto_ip && nbjettche==0)',
        '(pt2 > 15||!sameflav)',
        'ptll>%f'%ptllCut,
    ]

    # minimum for skimming
    wwmin = [
        'pt1>20',
        'pt2>10',
#        '(ch1*ch2)<0.5',
        'trigger==1',
        'pfmet>20',
        'mll>12',
        'mpmet>20',
        'bveto_mu==1',
        'nextra==0',
        '(bveto_ip==1 && nbjettche==0)',
#        'njet<4',
    ]
    
    #dy cuts
    dphilljj   = '(njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  )'

    # met cuts
    met  = '( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'
    mpmet  = '( !sameflav || ( (njet>1 || (mpmet>45 && dphiveto)) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'

    # met or
    metor = '( !sameflav || ( (njet!=0 || dymva1>0.88 || mpmet>35) && (njet!=1 || dymva1>0.84 || mpmet>35) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'

    wwcommon = wwcommon+[dphilljj, met]
    wwmin    = wwmin   +[dphilljj, metor]
                
    zerojet = '(njet == 0)'
    onejet  = '(njet == 1)'
    twojet  = '(njet == 2)'
    zeroOrOnejet = '(njet == 0 || njet == 1)'
    loosevbf= '(njet >= 2 && njet <= 3) '
    vbf     = '(njet >= 2 && njet <= 3 && (jetpt3 <= 30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) '
    vbf2011 = '(njet >= 2 && njet <= 3 && njetvbf == 0) '
    vh      = '(njet >= 2)'
    atleastonejet = '(njet >= 1)'


# vbf
class vbfcuts:
    _massindep = [
        'abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5',
        'abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5',
    ]

    _cut   = [
        'detajj>3.5',
        'mjj>500',
    ]
    _shape = [
        'detajj>3.',
        'mjj>200',
    ]

    _jeteta2011 = [
        'abs(jeteta1)<4.5 && abs(jeteta2)<4.5',
    ]

    vbfcut     = wwcuts.wwcommon     + [wwcuts.vbf]     + _massindep + _cut
    vbf2011cut = wwcuts.wwcommon2011 + [wwcuts.vbf2011] + _massindep + _cut + _jeteta2011


    vbfshape = ['(ch1*ch2)<0 && pt1>20 && pt2>10 &&   trigger==1. && pfmet>20. && mll>12 && zveto==1 && mpmet>20. && (njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  ) && bveto_mu==1 && nextra==0 && (bveto_ip==1 &&  (nbjettche==0 || njet>3))  && ptll>45. &&   ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) ) && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>3.5     && mjj>500']

    specificvbfloshape = ['pt1>20   &&      pt2>0   &&      mth>30  &&  mth<280      &&    mll<200 ']
    specificvbfhishape = ['pt1>50   &&      pt2>0   &&      mth>30  &&  mth<680      &&    mll<600 ']

    vbfloshape = vbfshape + specificvbfloshape
    vbfhishape = vbfshape + specificvbfhishape



class vhcuts:
    _massindep = [
        'trigger==1.',
	    '((ch1*ch2)<0 && pt1>20 && pt2>10)',
        'pfmet>20.',
        'mpmet>20.',
        'mll>12 ',
        'zveto==1',
        '(dphilljetjet<pi/180.*165. || !sameflav )',
        'bveto_mu==1 ',
	    'bveto_ip==1',
	    'nbjettche==0',
        'nextra==0' ,
        'ptll>30.',
        '( !sameflav ||  (pfmet > 45.0))',
    	'njet>=2',
	    'abs(jeteta1)<2.5 && abs(jeteta2)<2.5',
	    '(mjj<105 && mjj>65)',
	    'detajj<1.5',
	    #'sqrt((jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2)))/sqrt((pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))*(pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))+(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))*(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))) > 0.7',
	    #'sqrt((jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2)))/sqrt((pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))*(pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))+(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))*(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))) < 1.4',
	    #'(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))>3.14159266)*(2*3.14159265-abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))))+(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))<3.14159266)*(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))) > 150./180*pi',
	    'mth>60',
        # OLD jet veto
        #'njet==2',
        # new btag cut
        'jetbjpb1<1.4 &&  jetbjpb2<1.4',
    ]
# achtung: njet<4 is applied here!!!! in wwcommon!!!
    #vhcut     = wwcuts.wwcommon + _massindep
    #vhcut      = _massindep + ['mll<50']
    vhcut      = _massindep

    #vhcutShape = _massindep + ['((sameflav && mll<50) || (!sameflav && mll<200))']
    vhcutShape = _massindep



# whsc: wh, h>ww, ww>lnjj same charge:  w(+/-)h>ww>l(+/-)njj
class whsccuts:
    # ---> ~30xSM sensitivity 
    #_massindep = [
        #'trigger==1.',
        #'((ch1*ch2)>0 && pt1>20 && pt2>10)',
        #'zveto==1',
        #'mtw1<100',
        #'mtw2<100',
        #'mjj>60',
        #'mjj<100',
        #'detajj<3.5',
        #'pfmet>20.',
        #'( !sameflav ||  (pfmet > 45.0) )',
        #'mpmet>20.',
        #'mll>12 ',
        #'bveto_mu==1 ',
        #'nextra==0' ,
        #'bveto_ip==1',
        #'nbjettche==0',
        #'njet>=2',
    #]

    _sqrtDependent2011 = [
        'abs(jeteta1)<4.5 && abs(jeteta2)<4.5',
        #'mll>(12+8*sameflav)',
        'mll>12',
        'mpmet>(20+(7+nvtx/2.)*sameflav)',
        #'(pt2 > 15||!sameflav)',
    ]
    _sqrtDependent2012 = [
        'mll>12 ',
        '( !sameflav ||  (pfmet > 35.0))',
    ]

    # sensitivity to 20x SM
    #_massindep = [
        #'trigger==1.',
        #'((ch1*ch2)>0 && pt1>20 && pt2>10)',
        #'zveto==1',
        ##'(mjj+(mtw1>mtw2)*mtw2+(mtw1<mtw2)*mtw1)<180',
        #'pt2<35',#'pt2<40',
        ##'(mjj+mtw1+mtw2)>180',
        #'abs(eta1-eta2) < 2.0',
        ##'(mtw2*(mtw1>mtw2)+mtw1*(mtw1<mtw2))<100',
        #'(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))>70',
        #'(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))<130',
        #'mjj>50',
        #'mjj<100',
        #'pfmet>20.',
        #'mpmet>20.',
        #'bveto_mu==1 ',
        #'nextra==0' ,
        #'bveto_ip==1',
        #'nbjettche==0',
        #'njet>=2',
    #]

    # x12 SM
    _massindep = [
        'trigger==1.',
        '((ch1*ch2)>0 && pt1>20 && pt2>10)',
        'zveto==1',
        'pt2>15 && ((channel != 1) || pt2 > 20)',
        'abs(eta1-eta2) < 2.0',

        '(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))>50', #   --> these selections seem to distort the mllijj peak!
        '(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))<160',
        '(mtw2*(mtw1>mtw2)+mtw1*(mtw1<mtw2))>10',

        #'(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))>70',
        #'(mtw1*(mtw1>mtw2)+mtw2*(mtw1<mtw2))<130',
        #'(mtw2*(mtw1>mtw2)+mtw1*(mtw1<mtw2))>20',
        #'(mjj*(njet>=2))<100',

        '(mjj*(njet>=2))<110',
        'pfmet>20.',
        'mpmet>20.',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'nbjettche==0',
        'njet>=1',
        #'minmlijj<400',
        'mindphimlijj<400',
    ]



    #whsccut      = _massindep
    whsccut      = _massindep + _sqrtDependent2012
    whsccut2011  = _massindep + _sqrtDependent2011


# ggH+2jet
class ggH2j:

    _sqrtDependent2011 = [
        'abs(jeteta1)<4.5 && abs(jeteta2)<4.5',
        'mll>(12+8*sameflav)',
        'mll>12',
        'mpmet>(20+(17+nvtx/2.)*sameflav)',
        '(pt2 > 15||!sameflav)',
    ]
    _sqrtDependent2012 = [
        'mll>12 ',
        '( !sameflav ||  (pfmet > 45.0))',
    ]

    _massindep = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>10)',
        'zveto==1',
        '(mjj<65 || mjj>105)',
        '(njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  )',
        'pfmet>20.',
        'mpmet>20.',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'ptll>30.',
        'detajj<3.5',
        'mth>50',
        'drll<2.0',
        'jetbjpb1<1.4 && jetbjpb2<1.4',
        '(nbjettche==0 || njet>3)',
        'njet>=2',
    ]

    ggH2jShape  = _massindep + _sqrtDependent2012



# WW 2j ewk
class wwewk:

    _sqrtDependent2012 = [
        'mll>12 ',
        '( !sameflav ||  (pfmet > 45.0))',
    ]

    _massindep = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>20)',  # pt2>20 !
        'zveto==1',
        'pfmet>20.',
        'mpmet>20.',
        'mll>50',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'nbjettche==0',
        'njet>=2',
        #'njetvbf==0',
        '(dphilljetjet<pi/180.*165. || !sameflav )',
# mild VBF
        'mjj>350', # >250 #>400
        'detajj>0.5',
    ]

    _cut = [
        'WWewkMVABDTG > 0.40'
    ]

    wwewkShape  = _massindep + _sqrtDependent2012
    wwewkCut    = _massindep + _sqrtDependent2012 + _cut
    wwewkShapeTopEnriched  = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>20)',  # pt2>20 !
        'zveto==1',
        'pfmet>20.',
        'mpmet>20.',
        'mll>50',
        #'bveto_mu==1 ',
        'nextra==0' ,
        #'bveto_ip==1',
        'nbjettche>=1',
        'njet>=2',
        '(dphilljetjet<pi/180.*165. || !sameflav )',
# mild VBF
        'mjj>250', # >350 # >250 #>400
        'detajj>0.5',
    ]




# WW 2j
class ww2j:

    _sqrtDependent2012 = [
        'mll>12 ',
        '( !sameflav ||  (pfmet > 60.0))',
    ]

    _massindep = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>10)',
        'zveto==1',
        'pfmet>20.',
        'mpmet>20.',
        'mll>70',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'nbjettche==0',
        'njet>=2',
        '(mth>40 || sameflav)',
        'abs(jeteta1)<2.5',
        'abs(jeteta2)<2.5',
        'jetbjpb1<0.5',
        'jetbjpb2<0.5',
        '(dphilljetjet<pi/180.*165. || !sameflav )',
    ]

    ww2jCut    = _massindep + _sqrtDependent2012




class ww:

    _massindep = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>20)',
        'nextra==0' ,
        'pfmet>20.',
        'mll>12',
        'mpmet>20.  && ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) ) ',
        'zveto==1',
        #'(njet >=2 || (dphiveto || ! sameflav))',  --> if dymva, no dphi ll-jet cut
        'bveto_mu==1 ',
        'bveto_ip==1',
        'nbjettche==0',
        'ptll>30 && (!sameflav || ptll>45) ',
    ]

    wwCut    = _massindep





# Higgs width
class Hwidth:

    _massindep = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>10)',
        'zveto==1',
        'pfmet>20.',
        'mpmet>20.',
        'mll>12',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'nbjettche==0',
        #'mth>30'
        #'ptll>10'
    ]

    HwidthCut    = _massindep


    _massindep7TeV = [
        'trigger==1.',
        '((ch1*ch2)<0 && pt1>20 && pt2>10)',
        'zveto==1',
        'pfmet>20.',
        'mpmet>20.',
        'mll>12',
        'bveto_mu==1 ',
        'nextra==0' ,
        'bveto_ip==1',
        'nbjettche==0',
        #'mth>30'
        #'ptll>10'
    ]

    HwidthCut7TeV    = _massindep7TeV




masses = [110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000]
#masses = [110, 115, 120, 125, 130, 135, 140, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]  

categoryCuts = {}
categoryCuts['0j'] = wwcuts.zerojet
categoryCuts['1j'] = wwcuts.onejet
categoryCuts['2jex'] = wwcuts.twojet
categoryCuts['01j'] = wwcuts.zeroOrOnejet
#categoryCuts['2j']   = wwcuts.vbf        # 2 or 3 jets, but the third not between the first two in \eta
#categoryCuts['2j']   = wwcuts.loosevbf   # 2 or 3 jets  ----> changed! Need to propagate to some analysis? VBF cuts already applied in specific VBF selections!
categoryCuts['2j']   = wwcuts.vh         # >=2 jets
categoryCuts['vh2j'] = wwcuts.vh         # >=2 jets
categoryCuts['whsc'] = wwcuts.atleastonejet     # >=1 jets
categoryCuts['2jtche05']   = wwcuts.vh         # >=2 jets
categoryCuts['2jtche05CJ'] = wwcuts.vh         # >=2 jets
categoryCuts['2jtche05FJ'] = wwcuts.vh         # >=2 jets


flavorCuts = {}
flavorCuts['all'] = '1'			       #'channel>-1'
flavorCuts['mm']  = 'channel == 0'     #'channel>-0.5 && channel<0.5'
flavorCuts['ee']  = 'channel == 1'     #'channel> 0.5 && channel<1.5'
flavorCuts['em']  = 'channel == 2'     #'channel> 1.5 && channel<2.5'
flavorCuts['me']  = 'channel == 3'     #'channel> 2.5 && channel<4.5'
flavorCuts['sf']  = 'channel < 1.5'
flavorCuts['of']  = 'channel > 1.5' 

flavors = {}
flavors['mm']=['mm']
flavors['ee']=['ee']
flavors['em']=['em']
flavors['me']=['me']
flavors['sf']=['mm','ee']
flavors['of']=['em','me']
flavors['ll']=['mm','ee','em','me']

channels = {}
channels['0j']     = ('0j','ll')
channels['1j']     = ('1j','ll')
channels['of_0j']  = ('0j','of')
channels['of_1j']  = ('1j','of')
channels['sf_0j']  = ('0j','sf')
channels['sf_1j']  = ('1j','sf')
channels['of_01j'] = ('01j','of')
channels['sf_01j'] = ('01j','sf')
channels['mm_0j']  = ('0j','mm')
channels['mm_1j']  = ('1j','mm')
channels['ee_0j']  = ('0j','ee')
channels['ee_1j']  = ('1j','ee')
channels['em_0j']  = ('0j','em')
channels['em_1j']  = ('1j','em')
channels['me_0j']  = ('0j','me')
channels['me_1j']  = ('1j','me')

channels['2j']    = ('2j','ll')
channels['of_2j'] = ('2j','of')
channels['sf_2j'] = ('2j','sf')
channels['mm_2j'] = ('2j','mm')
channels['ee_2j'] = ('2j','ee')
channels['em_2j'] = ('2j','em')
channels['me_2j'] = ('2j','me')

channels['of_2jtche05'] = ('2jtche05','of')
channels['sf_2jtche05'] = ('2jtche05','sf')

channels['of_2jtche05CJ'] = ('2jtche05CJ','of')
channels['sf_2jtche05CJ'] = ('2jtche05CJ','sf')
channels['of_2jtche05FJ'] = ('2jtche05FJ','of')
channels['sf_2jtche05FJ'] = ('2jtche05FJ','sf')

channels['2jex']    = ('2jex','ll')
channels['of_2jex'] = ('2jex','of')
channels['sf_2jex'] = ('2jex','sf')
channels['mm_2jex'] = ('2jex','mm')
channels['ee_2jex'] = ('2jex','ee')
channels['em_2jex'] = ('2jex','em')
channels['me_2jex'] = ('2jex','me')


channels['of_vh2j'] = ('vh2j','of')
channels['sf_vh2j'] = ('vh2j','sf')

channels['ll_whsc'] = ('whsc','ll')


#  __  __               ___     _      
# |  \/  |__ _ ______  / __|  _| |_ ___
# | |\/| / _` (_-<_-< | (_| || |  _(_-<
# |_|  |_\__,_/__/__/  \___\_,_|\__/__/
#
# smurfs index              1     1     3     6     8     9    10    10    11    11    12    13    14    15    16   20     21    22    23    24    25    26    27
# masses              = [ 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 , 700 , 800 , 900 ,1000]
#  normal index             0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26


#  normal index             1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23

cutmap = {}
cutmap['mtmin_vh']    = [ 60  , 60  , 60  , 60  , 60  , 60  , 60  ,  60 ,  60 ,  60 ,  60 ,  60 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ,  70 ]
cutmap['mllmax_vh']   = [ 60  , 60  , 60  , 60  , 60  , 60  , 60  ,  60 ,  60 ,  60 ,  60 ,  60 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ]
cutmap['mllmin_vh']   = [ 12  , 12  , 12  , 12  , 12  , 12  , 20  ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ,  20 ]
cutmap['drllmax_vh']  = [ 1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.6 , 1.6 , 1.7 , 1.8 , 1.9 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 , 2.0 ]

cutmap['mllmax_bdt']  = [ 70  , 70  , 70  , 80  , 80  , 90  , 90  , 100 , 100 , 100 , 100 , 100 , 110 , 120 , 130 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 , 700 , 800 , 900 ,1000 ]
cutmap['pt1min']      = [ 20  , 20  , 20  , 23  , 25  , 25  , 25  , 25  , 27  , 27  , 30  , 34  , 36  , 38  , 40  , 55  , 70  , 80  , 90  , 110 , 120 , 130 , 140 , 150 , 160 , 170 , 180 ]
cutmap['pt2min']      = [ 10  , 10  , 10  , 10  , 10  , 12  , 15  , 15  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  ]
cutmap['mllmax']      = [ 40  , 40  , 40  , 43  , 45  , 45  , 45  , 45  , 50  , 50  , 50  , 50  , 60  , 80  , 90  , 150 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 600 , 700 , 800 , 900 ]
cutmap['dphimax']     = [ 115 , 115 , 115 , 100 , 90  , 90  , 90  , 90  , 90  , 90  , 60  , 60  , 70  , 90  , 100 , 140 , 175 , 175 , 175 , 175 , 175 , 175 , 175 , 175 , 175 , 175 , 175 ]

cutmap['mtmin']       = [ 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 90  , 110 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 ]
cutmap['mtmax']       = [ 110 , 110 , 120 , 123 , 125 , 128 , 130 , 140 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 , 700 , 800 , 900 ,1000 ]

# for whsc                                                                                                                      [- random numbers from here on, I stop at 300 GeV
#cutmap['mlljjmin']    = [ 90  , 95  , 100  ,105  ,110  ,110 , 110 , 110 , 110 , 115 , 120 , 130 , 130 , 130 , 140 , 140 , 150 , 120 , 120 , 120 , 120 , 120 , 120 ]
cutmap['mlljjmin']    = [ 50  , 55  ,  60  , 65  , 70  , 75 ,  80 ,  80 ,  80 ,  85 ,  90 ,  95 , 100 , 120 , 140 , 140 , 150 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 ]
cutmap['mlljjmax']    = [ 140 , 145 , 150 , 165 , 170 , 175 , 175 , 175 , 180 , 185 , 190 , 200 , 210 , 220 , 220 , 270 , 320 , 190 , 190 , 190 , 190 , 190 , 190 , 190 , 190 , 190 , 190 ]
#cutmap['mlljjmax']    = [ 145 , 150 , 180 , 180 , 180 , 180 , 185 , 185 , 185 , 190 , 190 , 200 , 210 , 220 , 220 , 270 , 320 , 190 , 190 , 190 , 190 , 190 , 190 ]


# automatic conversion to radiants
cutmap['dphimax']     = [ str(phi)+'*pi/180' for phi in cutmap['dphimax'] ]

massDependantCutsbyVar = {}
for c,list in cutmap.iteritems():
    if len(list) != len(masses):
        raise RuntimeError('Wrong number of entries in {cut} ({list} vs. {masses})'.format( cut = c, list = len(list), masses = len(masses)) )
    massDependantCutsbyVar[c] = dict(zip(masses,list))
# cleanup
del cutmap
del c
del list
del phi



#  ___      _        _   _             
# / __| ___| |___ __| |_(_)___ _ _  ___
# \__ \/ -_) / -_) _|  _| / _ \ ' \(_-<
# |___/\___|_\___\__|\__|_\___/_||_/__/
#                                      

def massSelections(mass):

    # Spin
    #mthmin_2dlomass = 80.
    mthmin_2dlomass = 60. 

    mthmax_2dlomass = 280.
    mllmax_2dlomass = 200

    mthmin_bdt = 80.
    mthmin_vbf = 30.
    mthmin_vh  = 50.

    masscuts = dict([(cut,massDependantCutsbyVar[cut][mass]) for cut in massDependantCutsbyVar])

    sel = {}
    sel['ww-min']        = ' && '.join(wwcuts.wwmin)
    sel['ww-common']     = ' && '.join(wwcuts.wwcommon)
    sel['wwbtag-common'] = sel['ww-common'].replace('bveto_mu==1','bveto_mu>-1').replace('(bveto_ip==1 && nbjettche==0)','(bveto_ip>-1 && nbjettche>-1)')
    sel['ww2011-common'] = ' && '.join(wwcuts.wwcommon2011)
    sel['ww2011btag-common'] = sel['ww2011-common'].replace('bveto_mu==1','bveto_mu>-1').replace('(bveto_ip==1 && nbjettche==0)','(bveto_ip>-1 && nbjettche>-1)')

    sel['vbf-shape-2d-himass']    = ' && '.join(vbfcuts.vbfhishape)
    sel['vbf-shape-2d-lomass']    = ' && '.join(vbfcuts.vbfloshape)
    sel['vbf-shape-2d']           = sel['vbf-shape-2d-lomass'] if mass <= 250 else sel['vbf-shape-2d-himass'] 
    
    sel['shape-lomass'] = 'mth>%f && mth<%f && mll<%f '%(mthmin_2dlomass,mthmax_2dlomass,mllmax_2dlomass)
    sel['shape-himass'] = 'mth>80 && mth<380 && mll<450 && pt1>50'
    sel['shape-vhimass']= 'mth>80 && pt1>50'
    
    sel['vbf-level']     = ' && '.join(vbfcuts.vbfcut)
    sel['vbf2011-level'] = ' && '.join(vbfcuts.vbf2011cut)

    sel['ww-level']     = sel['ww-common']+'&& ptll>45'
    sel['wwbtag-level'] = sel['wwbtag-common']+'&& ptll>45'
    sel['ww2011-level'] = sel['ww2011-common']+'&& ptll>45'
    sel['ww2011btag-level'] = sel['ww2011btag-common']+'&& ptll>45'
    sel['bdt-specific'] = 'mll < {0} && (mth > {1:.0f} && mth < {2:.0f})'.format(masscuts['mllmax_bdt'], mthmin_bdt, int(mass))

    #sel['ww-xsec'] = sel['ww-common'].replace(wwcuts.met, wwcuts.mpmet)+' && pt2>20 && ptll>45'
    sel['ww-xsec'] = sel['ww-common']+' && pt2>20 && (!sameflav || ptll>45)'

    hwwlvl = {}
    hwwlvl['mll']    = 'mll < {0}'.format(masscuts['mllmax'])
    hwwlvl['pt1']    = 'pt1 > {0:.1f}'.format(masscuts['pt1min'])
    hwwlvl['pt2']    = 'pt2 > {0:.1f}'.format(masscuts['pt2min'])
    hwwlvl['dphill'] = 'dphill < {0}'.format(masscuts['dphimax'])
    hwwlvl['mth']    = '(mth > {0:.1f} && mth < {1:.1f})'.format(masscuts['mtmin'], masscuts['mtmax'])

    sel['wwxsec-selection']       = sel['ww-xsec']
    sel['ww-selection']           = sel['ww-level']
    sel['wwbtag-selection']       = sel['wwbtag-level']
    sel['ww2011-selection']       = sel['ww2011-level']
    sel['ww2011btag-selection']   = sel['ww2011btag-level']
    sel['wwr-selection']          = sel['ww-common']
    sel['wwrbtag-selection']      = sel['wwbtag-common']
    sel['wwr2011-selection']      = sel['ww2011-common']
    sel['wwr2011btag-selection']  = sel['ww2011btag-common']
    sel['wwtight-selection']      = sel['ww-level']+' && mth > %f'%mthmin_2dlomass
    sel['wwloose-selection']      = sel['ww-level'].replace('ptll>%f'%ptllCut,'ptll>20')+' && mth > %f'%mthmin_2dlomass
    sel['hww-selection']          = ' && '.join([sel['ww-level']]+[cut for cut     in hwwlvl.itervalues()])
    sel['hwwbtag-selection']      = ' && '.join([sel['wwbtag-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    sel['hww2011-selection']      = ' && '.join([sel['ww2011-level']]+[cut for cut     in hwwlvl.itervalues()])
    sel['hww2011btag-selection']  = ' && '.join([sel['ww2011btag-level']]+[cut for var,cut in hwwlvl.iteritems() if (var != 'mth' and var != 'mll' and var != 'pt1' and var != 'pt2')])
    sel['mll-selection']          = ' && '.join([sel['ww-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mll'])
    sel['mth-selection']          = ' && '.join([sel['ww-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    sel['dphi-selection']         = ' && '.join([sel['ww-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'dphill'])
    sel['gammaMRStar-selection']  = ' && '.join([sel['ww-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mll'])
    sel['bdt-selection']          = sel['ww-level']+' && '+sel['bdt-specific']
    sel['bdtl-selection']         = sel['bdt-selection']
    sel['vbf-shape-2d-selection'] = sel['vbf-shape-2d']


#    sel['ww-common-ss']           = sel['ww-common'].replace["(ch1*ch2)<0.5", "(ch1*ch2)>0.5"]+' && 1'
    sel['shape-withCR-selection']  = (sel['ww-common'].replace("zveto==1", "zveto>-1")).replace("(ch1*ch2)<0", "1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-himass']

    #### to look only at the same sign control region ####
    sel['shape-onlySSCR-selection']  = (sel['ww-common'].replace("zveto==1", "zveto>-1")).replace("(ch1*ch2)<0", "(ch1*ch2)>0")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-himass']


    sel['ww-common-ss']           = sel['ww-common'].replace("(ch1*ch2)<0.5", "(ch1*ch2)>0.5")
    sel['shape-ss-selection']     = sel['ww-common-ss'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common-ss'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-himass']

    sel['shape-selection']        = sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-himass']
    sel['shapebtag-selection']    = sel['wwbtag-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-lomass'] if mass <=250 else sel['wwbtag-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-himass']
    sel['shape2011-selection']     = sel['ww2011-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-himass']
    sel['shape2011btag-selection'] = sel['ww2011btag-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011btag-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-himass']

    # EWK Singlet 0/1-jet (2012 data)
    #if   mass <= 250 : sel['shape-ewks-selection'] = sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-lomass'] + ' && pt2>30.&&mth>80'
    if   mass <= 250 : sel['shape-ewks-selection'] = sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-lomass'] 
    elif mass <= 600 : sel['shape-ewks-selection'] = sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-himass']
    else             : sel['shape-ewks-selection'] = sel['ww-common'].replace("zveto==1", "zveto>-1").replace(wwcuts.met, wwcuts.mpmet)+' && '+sel['shape-vhimass']

    # EWK Singlet 0/1-jet (2011 data)
    #if   mass <= 250 : sel['shape2011-ewks-selection'] = sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-lomass'] + ' && pt2>30.&&mth>80'
    if   mass <= 250 : sel['shape2011-ewks-selection'] = sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-lomass'] 
    elif mass <= 600 : sel['shape2011-ewks-selection'] = sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-himass']
    else             : sel['shape2011-ewks-selection'] = sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-vhimass']


    sel['shapehcp-selection']      =  sel['shape-selection'].replace('ptll>%f'%ptllCut,'ptll>45').replace('mth>%f'%mthmin_2dlomass,'mth>80')
    sel['shapehcp2011-selection']  =  sel['shape2011-selection'].replace('ptll>%f'%ptllCut,'ptll>45').replace('mth>%f'%mthmin_2dlomass,'mth>80')
        
    #sel['shape-selection']        = sel['ww-common']+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common']+' && '+sel['shape-himass']
    #sel['shape2011-selection']     = sel['ww2011-common']+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011-common']+' && '+sel['shape-himass']
   
    
    # vbf #
    sel['vbf-selection-temp']     = ' && '.join([sel['vbf-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    sel['vbf-selection-temp']     = sel['vbf-selection-temp'].replace('ptll>30.', 'ptll>45.')
    sel['vbf-selection']          = sel['vbf-selection-temp'] + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(mthmin_vbf, masscuts['mtmax'])

    sel['vbf-selection-top-temp']  = ' && '.join([sel['vbf-level']])
    sel['vbf-selection-top']       = sel['vbf-selection-top-temp'].replace('ptll>30.', 'ptll>45.')


    # my AN 2013/097
    #sel['vbf-shape-temp']      = ' && '.join([sel['vbf-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth' and var != 'mll'])
    #sel['vbf-shape']           = sel['vbf-shape-temp'].replace('ptll>30.', 'ptll>45.')
    #sel['vbf-shape-selection'] = sel['vbf-shape'] + ' && ((sameflav && (mth > {0:.1f} && mth < {1:.1f}) && mll<{2:.1f}) || (!sameflav && (mth > {3:.1f} && mth < {4:.1f})))   '.format(mthmin_vbf, masscuts['mtmax'],masscuts['mllmax'],50,int(mass))


   # to synch with smurf # only emu
    sel['vbf-shape-temp']      = ' && '.join([sel['vbf-level']])
    sel['vbf-shape-selection'] = sel['vbf-shape-temp'] + ' && (!sameflav && (mth > {0:.1f} && mth < {1:.1f})) '.format(30,int(mass))
    sel['vbf-shape-selection'] = sel['vbf-shape-selection'].replace('ptll>30.', 'ptll>45.')

    if (mass > 250) :
        sel['vbf-shape-selection'] = sel['vbf-shape-selection'] + ' && pt1>50'



    # new cut mth > 30 everywhere to synch

    ##sel['vbf-shape']    = ' && '.join(vbfcuts.vbfshape)
    ##sel['vbf-shape']    =  sel['vbf-shape'].replace("ptll>30.", "ptll>45.")
    ##sel['vbf-shape-selection'] = sel['vbf-shape']+' && (mth > 50 && mth < {0:.0f})'.format(int(mass))

    ## vbf in 2011#
    #sel['vbf2011-selection-temp']     = ' && '.join([sel['vbf2011-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth' and var != 'mll' and var != 'pt1' and var != 'pt2' and var != 'dphill'])
    ##sel['vbf2011-selection-temp']     = ' && '.join([sel['vbf2011-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    #sel['vbf2011-selection-temp']     = sel['vbf2011-selection-temp'].replace('dphiveto', 'dphilljetjet<pi/180.*165.')
    #sel['vbf2011-selection-temp']     = sel['vbf2011-selection-temp'].replace('ptll>30.', 'ptll>45.')
    #sel['vbf2011-selection']          = sel['vbf2011-selection-temp'] + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(mthmin_vbf, masscuts['mtmax'])
    #if (mass <= 200) :
      #sel['vbf2011-selection']        = sel['vbf2011-selection']      + ' && (mll<100)'

    #sel['vbf2011-selection-top-temp']  = ' && '.join([sel['vbf2011-level']])
    #sel['vbf2011-selection-top']       = sel['vbf2011-selection-top-temp'].replace('ptll>30.', 'ptll>45.')
    #sel['vbf2011-selection-top']       = sel['vbf2011-selection-top'].replace('dphiveto', 'dphilljetjet<pi/180.*165.')

   # to synch with smurf
    # vbf in 2011#
    sel['vbf2011-selection-temp']     = ' && '.join([sel['vbf2011-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    sel['vbf2011-selection-temp']     = sel['vbf2011-selection-temp'].replace('dphiveto', 'dphilljetjet<pi/180.*165.')
    sel['vbf2011-selection-temp']     = sel['vbf2011-selection-temp'].replace('ptll>30.', 'ptll>45.')
    sel['vbf2011-selection']          = sel['vbf2011-selection-temp'] + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(mthmin_vbf, masscuts['mtmax'])

    sel['vbf2011-selection-top-temp']  = ' && '.join([sel['vbf2011-level']])
    sel['vbf2011-selection-top']       = sel['vbf2011-selection-top-temp'].replace('ptll>30.', 'ptll>45.')
    sel['vbf2011-selection-top']       = sel['vbf2011-selection-top'].replace('dphiveto', 'dphilljetjet<pi/180.*165.')


    # vbf in 2011#
    # to synch with smurf # only emu
    sel['vbf2011-shape-temp']      = ' && '.join([sel['vbf2011-level']])
    sel['vbf2011-shape-selection'] = sel['vbf2011-shape-temp'] + ' && (!sameflav && (mth > {0:.1f} && mth < {1:.1f}))   '.format(30,int(mass))
    sel['vbf2011-shape-selection'] = sel['vbf2011-shape-selection'].replace('ptll>30.', 'ptll>45.')
    if (mass > 250) :
        sel['vbf2011-shape-selection'] = sel['vbf2011-shape-selection'] + ' && pt1>50'



    # vbf fisher discriminant
    sel['vbf-shape-fish-temp']      = ' && '.join([sel['vbf-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth' and var != 'mll'])
    sel['vbf-shape-fish']           = sel['vbf-shape-fish-temp'].replace('ptll>30.', 'ptll>45.')
    #sel['vbf-shape-fish']           = sel['vbf-shape-fish-temp'].replace('ptll>45.', 'ptll>30.')
    sel['vbf-shape-fish']           = sel['vbf-shape-fish'].replace('mjj>500','mjj>150')
    sel['vbf-shape-fish']           = sel['vbf-shape-fish'].replace('detajj>3.5','(1.3491*log(mjj)-0.01163*(jeteta1*jeteta2)+1.433*detajj-12.647>3.0)')
    sel['vbf-shape-fish-selection'] = sel['vbf-shape-fish'] + ' && ((sameflav && (mth > {0:.1f} && mth < {1:.1f}) && mll<{2:.1f}) || (!sameflav && (mth > {3:.1f} && mth < {4:.1f})))   '.format(mthmin_vbf, masscuts['mtmax'],masscuts['mllmax'],50,int(mass))

    sel['vbf-selection-fish-top-temp']  = ' && '.join([sel['vbf-level']])
    sel['vbf-selection-fish-top']       = sel['vbf-selection-fish-top-temp'].replace('ptll>30.', 'ptll>45.')
    #sel['vbf-selection-fish-top']       = sel['vbf-selection-fish-top-temp'].replace('ptll>45.', 'ptll>30.')
    sel['vbf-selection-fish-top']       = sel['vbf-selection-fish-top'].replace('mjj>500','mjj>150')
    sel['vbf-selection-fish-top']       = sel['vbf-selection-fish-top'].replace('detajj>3.5','(1.3491*log(mjj)-0.01163*(jeteta1*jeteta2)+1.433*detajj-12.647>3.0)')


    # vbf banana #
    sel['vbf-banana-shape']    = ' && '.join(vbfcuts.vbfshape)
    sel['vbf-banana-shape']    = sel['vbf-banana-shape'].replace("(ch1*ch2)<0", "((ch1*ch2)<0 || !sameflav)").replace("ptll>30.", "ptll>45.")
    sel['vbf-banana-shape-selection'] = sel['vbf-banana-shape']+' && (mth > 50 && mth < {0:.0f})'.format(int(mass))






    # vh #
    sel['vh-level']    = ' && '.join(vhcuts.vhcut)
    sel['vh-selection'] = sel['vh-level']     + ' && (mth > {0:.1f} && mth < {1:.1f}) && drll<{2:.1f} '.format(masscuts['mtmin_vh'], int(mass), masscuts['drllmax_vh'])
    sel['vh-selection'] = sel['vh-selection'] + ' && (mll > {0:.1f} && mll < {1:.1f})'.format(masscuts['mllmin_vh'], masscuts['mllmax_vh'])

    sel['vh-shape-level']     = ' && '.join(vhcuts.vhcutShape)
    sel['vh-shape-selection'] = sel['vh-shape-level']
    sel['vh-shape-selection'] = sel['vh-shape-selection'] + ' && (mth > {0:.1f} && mth < {1:.1f}) && ((sameflav && drll<{2:.1f}) || (!sameflav && drll<2.5))'.format(masscuts['mtmin_vh'], int(mass),  masscuts['drllmax_vh'])
    sel['vh-shape-selection'] = sel['vh-shape-selection'] + ' && ((sameflav && mll>{0:.1f} && mll<{1:.1f}) || (!sameflav && mll<200))'.format(masscuts['mllmin_vh'], masscuts['mllmax_vh'])


    # vh banana #

    sel['vh-banana-shape-level']     = ' && '.join(vhcuts.vhcutShape)
    sel['vh-banana-shape-selection'] = sel['vh-banana-shape-level'].replace('(ch1*ch2)<0', '((ch1*ch2)<0 || !sameflav)')
    sel['vh-banana-shape-selection'] = sel['vh-banana-shape-selection'] + ' && (mth > {0:.1f} && mth < {1:.1f}) && ((sameflav && drll<1.5) || (!sameflav && drll<2.5))'.format(masscuts['mtmin_vh'], int(mass))
    sel['vh-banana-shape-selection'] = sel['vh-banana-shape-selection'] + ' && ((sameflav && mll>{0:.1f} && mll<{1:.1f}) || (!sameflav && mll<200)) && ((sameflav && drll<1.5) || (!sameflav && drll<2.5)) '.format(masscuts['mllmin_vh'], masscuts['mllmax_vh'])



    # whsc #

    sel['whsc-level']    = ' && '.join(whsccuts.whsccut) + ' && ( mindphimlijj > {0:.1f})  '.format(masscuts['mlljjmin'])   + ' && ( mindphimlijj < {0:.1f})  '.format(masscuts['mlljjmax'])
    #sel['whsc-level']    = ' && '.join(whsccuts.whsccut) + ' && ( minmlijj > {0:.1f})  '.format(masscuts['mlljjmin'])   + ' && ( minmlijj < {0:.1f})  '.format(masscuts['mlljjmax'])
    sel['whsc-selection'] = sel['whsc-level']


    sel['whsc2011-level']    = ' && '.join(whsccuts.whsccut2011) + ' && ((min(sqrt((2*pt1*cosh(eta1)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))*(2*pt1*cosh(eta1)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))-((2*pt1*sin(phi1)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(2*pt1*sin(phi1)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))+(2*pt1*cos(phi1)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(2*pt1*cos(phi1)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(2*pt1*(1-exp(-2*eta1))/(2.*exp(-eta1))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2)))*(2*pt1*(1-exp(-2*eta1))/(2.*exp(-eta1))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2))))),sqrt((2*pt2*cosh(eta2)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))*(2*pt2*cosh(eta2)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))-((2*pt2*sin(phi2)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(2*pt2*sin(phi2)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))+(2*pt2*cos(phi2)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(2*pt2*cos(phi2)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(2*pt2*(1-exp(-2*eta2))/(2.*exp(-eta2))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2)))*(2*pt2*(1-exp(-2*eta2))/(2.*exp(-eta2))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2))))))) > {0:.1f})  '.format(masscuts['mlljjmin'])   + ' && ((min(sqrt((2*pt1*cosh(eta1)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))*(2*pt1*cosh(eta1)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))-((2*pt1*sin(phi1)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(2*pt1*sin(phi1)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))+(2*pt1*cos(phi1)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(2*pt1*cos(phi1)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(2*pt1*(1-exp(-2*eta1))/(2.*exp(-eta1))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2)))*(2*pt1*(1-exp(-2*eta1))/(2.*exp(-eta1))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2))))),sqrt((2*pt2*cosh(eta2)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))*(2*pt2*cosh(eta2)+jetpt1*cosh(jeteta1)+jetpt2*cosh(jeteta2))-((2*pt2*sin(phi2)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(2*pt2*sin(phi2)+jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))+(2*pt2*cos(phi2)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(2*pt2*cos(phi2)+jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(2*pt2*(1-exp(-2*eta2))/(2.*exp(-eta2))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2)))*(2*pt2*(1-exp(-2*eta2))/(2.*exp(-eta2))+jetpt1*(1-exp(-2*jeteta1))/(2.*exp(-jeteta1))+jetpt2*(1-exp(-2*jeteta2))/(2.*exp(-jeteta2))))))) < {0:.1f})  '.format(masscuts['mlljjmax'])
    sel['whsc2011-selection'] = sel['whsc2011-level']

    sel['whscShape-level']     = ' && '.join(whsccuts.whsccut)
    sel['whscShape-selection'] = sel['whscShape-level']

    sel['whscShape2011-level']     = ' && '.join(whsccuts.whsccut2011)
    sel['whscShape2011-selection'] = sel['whscShape2011-level']


    sel['whscShape-ee-selection'] = sel['whscShape-level'] + ' && (channel == 1)'
    sel['whscShape-em-selection'] = sel['whscShape-level'] + ' && (channel == 2)'
    sel['whscShape-me-selection'] = sel['whscShape-level'] + ' && (channel == 3)'
    sel['whscShape-mm-selection'] = sel['whscShape-level'] + ' && (channel == 0)'

    # mumu #    channel == 0
    # mue #     channel == 3
    # emu #     channel == 2
    # ee #      channel == 1



    sel['whsc-ee-selection'] = sel['whsc-level'] + ' && (channel == 1)'
    sel['whsc-em-selection'] = sel['whsc-level'] + ' && (channel == 2)'
    sel['whsc-me-selection'] = sel['whsc-level'] + ' && (channel == 3)'
    sel['whsc-mm-selection'] = sel['whsc-level'] + ' && (channel == 0)'


    # ggH2j

    sel['ggH2jShape-level']     = ' && '.join(ggH2j.ggH2jShape) + ' && mth < {0:.1f}'.format(int(mass))
    sel['ggH2jShape-selection'] = sel['ggH2jShape-level']


    # WW 2j ewk

    sel['wwewkShape-level']     = ' && '.join(wwewk.wwewkShape)
    sel['wwewkShape-selection'] = sel['wwewkShape-level']

    sel['wwewkCut-level']     = ' && '.join(wwewk.wwewkCut)
    sel['wwewkCut-selection'] = sel['wwewkCut-level']

    sel['wwewkShape05-level']     = ' && '.join(wwewk.wwewkShape) + ' && jettche1 < 1.00   && jettche2 < 1.00 '
    sel['wwewkShape05-selection'] = sel['wwewkShape05-level']
    sel['wwewkShape21-level']     = ' && '.join(wwewk.wwewkShape) + ' && jettche1 >= 1.00  && jettche2 >= 1.00 '
    sel['wwewkShape21-selection'] = sel['wwewkShape21-level']

    sel['wwewkCut05-level']     = ' && '.join(wwewk.wwewkCut) + ' && jettche1 < 1.00   && jettche2 < 1.00 '
    sel['wwewkCut05-selection'] = sel['wwewkCut05-level']
    sel['wwewkCut21-level']     = ' && '.join(wwewk.wwewkCut) + ' && jettche1 >= 1.00  && jettche2 >= 1.00 '
    sel['wwewkCut21-selection'] = sel['wwewkCut21-level']

    sel['wwewkCut05CJ-level']     = ' && '.join(wwewk.wwewkCut) + ' && (((abs(jeteta1)<abs(jeteta2))  && (jettche1<1.00)) || ((abs(jeteta1)>=abs(jeteta2))  && (jettche2<1.00))) && (((abs(jeteta1)<abs(jeteta2))  && (jettche2>1.00)) || ((abs(jeteta1)>=abs(jeteta2))  && (jettche1>1.00))) '
    sel['wwewkCut05CJ-selection'] = sel['wwewkCut05CJ-level']
    sel['wwewkCut05FJ-level']     = ' && '.join(wwewk.wwewkCut) + ' && (((abs(jeteta1)<abs(jeteta2))  && (jettche2<1.00)) || ((abs(jeteta1)>=abs(jeteta2))  && (jettche1<1.00))) && (((abs(jeteta1)<abs(jeteta2))  && (jettche1>1.00)) || ((abs(jeteta1)>=abs(jeteta2))  && (jettche2>1.00))) '
    sel['wwewkCut05FJ-selection'] = sel['wwewkCut05FJ-level']



    sel['wwewkShape-top-enriched-level']     = ' && '.join(wwewk.wwewkShapeTopEnriched)
    sel['wwewkShape-top-enriched-selection'] = sel['wwewkShape-top-enriched-level']


    # WW 2j

    sel['ww2jShape-level']     = ' && '.join(wwewk.wwewkShape)
    sel['ww2jShape-selection'] = sel['ww2jShape-level']

    sel['ww2jCut-level']     = ' && '.join(ww2j.ww2jCut)
    sel['ww2jCut-selection'] = sel['ww2jCut-level']



    # WW 0/1 jet

    sel['CutWW-level']     = ' && '.join(ww.wwCut)
    sel['CutWW-selection'] = sel['CutWW-level']



    # Higgs width
    sel['Hwidth-level']     = ' && '.join(Hwidth.HwidthCut)
    sel['Hwidth7TeV-level']     = ' && '.join(Hwidth.HwidthCut7TeV)

    #sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * (1)) || '  + \
                          #' ((njet==1) * (1)) || ' +  \
                          #' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>3.5     && mjj>500    )) '  + \
                          #')'
    #sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * (1)) || '  + \
                          #' ((njet==1) * (1)) || ' +  \
                          #' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.0     && mjj>300    )) '  + \
                          #')'

    sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || '  + \
                          ' ((njet==1) * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || ' +  \
                          ' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          '  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          ')'


    sel['HwidthWWcontrolRegion-selection'] = sel['Hwidth-level']+ ' && (pt2>20) && ( (njet==0 * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || '  + \
                          ' ((njet==1) * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || ' +  \
                          ' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          '  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          ')'

    sel['HwidthWWcontrolRegionHighMllRemoved-selection'] = sel['Hwidth-level']+ ' && (pt2>20) && (mll<100)  && ( (njet==0 * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || '  + \
                          ' ((njet==1) * ( (mll<70 && pt2<40 && ptll>45) || (mll>=70 && pt2>20 && pt1>20)  )) || ' +  \
                          ' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          '  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          ')'




    #sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * ( (HwidthMVAggH>=0.0 && pt2>20) || (HwidthMVAggH<0.0 && mth>60 && mth<120 && ptll>45 && pt2<40)  )) || '  + \
                          #' ((njet==1) * ( (HwidthMVAggH>=0.0 && pt2>20) || (HwidthMVAggH<0.0 && mth>60 && mth<120 && ptll>45 && pt2<40)  )) || ' +  \
                          #' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          #'  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          #')'

    sel['Hwidthmthmll-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50)  )) || '  + \
                          ' ((njet==1) * ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50)  )) || ' +  \
                          ' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          '  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          ')'



    #sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * ( (HwidthMVAggH>=0.0) || (HwidthMVAggH<0.0 && mth<125 && ptll>45)  )) || '  + \
                          #' ((njet==1) * ( (HwidthMVAggH>=0.0) || (HwidthMVAggH<0.0 && mth<125 && ptll>45)  )) || ' +  \
                          #' ((njet>=2) * ( ptll>45 && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.0     && mjj>300 && ' + \
                          #'  (   (mll<70 && mth<125) || (mll>=70 && pt2>20 && pt1>20)  )  )) '  + \
                          #')'

    #sel['Hwidth-level']     = ' && '.join(Hwidth.HwidthCut)
    #sel['Hwidth-selection'] = sel['Hwidth-level']+ ' && ( (njet==0 * (1)) || '  + \
                          #' ((njet==1) * (1)) || ' +  \
                          #' ((njet>=2) * (  (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.0     && mjj>300    )) '  + \
                          #')'

    sel['Hwidthmthmll7TeV-selection'] = sel['Hwidth7TeV-level']+ ' && ( (njet==0 * ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50)  )) || '  + \
                          ' ((njet==1) * ( (mll<70 && pt2<35 && ptll>45) || (mll>=70 && pt2>20 && pt1>50)  )) || ' +  \
                          ' ((njet>=2) * ( ptll>45 && (njet >= 2 && njet <= 3 && njetvbf == 0)   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>2.5     && mjj>500 && ' + \
                          '  (   (mll<70 && pt2<40) || (mll>=70 && pt2>20 && pt1>20)  )  &&  '  + \
                          '  (  abs(jeteta1)<4.5 && abs(jeteta2)<4.5  )  )) '  + \
                          ')'


    return sel




