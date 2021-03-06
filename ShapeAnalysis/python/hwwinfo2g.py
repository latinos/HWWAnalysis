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
#ptllCut = 45. 
# Spin/Moriond
ptllCut = 30. 

class wwnamedcuts:
    wwcommon = odict.OrderedDict([
        ('trigger'  , 'trigger == 1'), 
        ('pt1'      , 'pt1>20'), 
        ('pt2'      , 'pt2>10'), 
        ('os'       , '(ch1*ch2)<0.5'), # don't use 0
        ('trigger'  , 'trigger==1.'), 
        ('pfmet'    , 'pfmet>20.'), 
        ('mllmin'   , 'mll>12'), # ema7
        ('zveto'    , 'zveto==1 || !sameflav'), 
        ('mpmet'    , 'mpmet>20.'), # ema9
        ('bveto_mu' , 'bveto_mu==1'), 
        ('nextra'   , 'nextra==0'), 
        ('bveto_ip' , 'bveto_ip==1 && nbjettche==0'), 
        ('ptll'     , 'ptll>%.2f'%ptllCut), # ema 14
        ('dphilljj' , 'njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )'),
        ('met'      , '!sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || pfmet > 45.0) )'),
    ])

    # commom cuts
    vbfcommon = wwcommon.copy()
    vbfcommon.update([
        ('lepcnt1','abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5'),
        ('lepcnt2','abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5'),
    ])

    vbfcutbased = vbfcommon.copy()
    vbfcutbased.update([
        ('detajj','detajj>3.5'),
        ('mjj','mjj>500'),
        ])

    vbfloshape = vbfcommon.copy()
    vbfloshape.update([
        ('minpt1','pt1>20'),
        ('minpt2','pt2>0'),
        ('mth'   ,'mth>30  &&  mth<280'),
        ('mllmax','mll<200'),
    ])

    vbfhishape = vbfcommon.copy()
    vbfhishape.update([
        ('minpt1','pt1>50'),
        ('minpt2','pt2>0'),
        ('mth'   ,'mth>30  &&  mth<680'),
        ('mllmax','mll<600'),
    
    ])
    
    zerojet = 'njet == 0'
    onejet  = 'njet == 1'
    vbf     = 'njet >= 2 && njet <= 3 && (jetpt3 <= 30 || !(jetpt3 > 30 && (  ( jeteta1 < jeteta3 && jeteta3 < jeteta2 ) || ( jeteta2 < jeteta3 && jeteta3 < jeteta1 ) ))  ) '

#   normal index              0     1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21    22
    masses               = [ 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]
    mcuts = {}
    mcuts['mtmin_vh']    = [  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  50 ,  60 ,  60 ,  60 ,  60 ,  60 , 100 , 100 , 100 , 100 , 100 , 100 ]
    mcuts['mllmax_vh']   = [  70 ,  70 ,  70 ,  80 ,  80 ,  90 ,  90 , 100 , 100 , 100 , 100 , 100 , 110 , 120 , 130 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

    mcuts['mllmax_bdt']  = [  70 ,  70 ,  70 ,  80 ,  80 ,  90 ,  90 , 100 , 100 , 100 , 100 , 100 , 110 , 120 , 130 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

    mcuts['pt1min']      = [  20 ,  20 ,  20 ,  23 ,  25 ,  25 ,  25 ,  25 ,  27 ,  27 ,  30 ,  34 ,  36 ,  38 ,  40 ,  55 ,  70 ,  80 ,  90 , 110 , 120 , 130 , 140 ]
    mcuts['pt2min']      = [  10 ,  10 ,  10 ,  10 ,  10 ,  12 ,  15 ,  15 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ,  25 ]
    mcuts['mllmax']      = [  40 ,  40 ,  40 ,  43 ,  45 ,  45 ,  45 ,  45 ,  50 ,  50 ,  50 ,  50 ,  60 ,  80 ,  90 , 150 , 200 , 250 , 300 , 350 , 400 , 450 , 500 ]
    mcuts['dphimax']     = [ 115 , 115 , 115 , 100 ,  90 ,  90 ,  90 ,  90 ,  90 ,  90 ,  60 ,  60 ,  70 ,  90 , 100 , 140 , 175 , 175 , 175 , 175 , 175 , 175 , 175 ]

    mcuts['mtmin']       = [  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  80 ,  90 , 110 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 ]
    mcuts['mtmax']       = [ 110 , 110 , 120 , 123 , 125 , 128 , 130 , 140 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

#                          [ 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

    # automatic conversion to radiants
    mcuts['dphimax']     = [ str(phi)+'*pi/180' for phi in mcuts['dphimax'] ]
    del phi

    @staticmethod
    def hww0j1jcb( mass ):
        
        mcuts = wwnamedcuts.mcuts
        i = wwnamedcuts.masses.index(mass)

        cuts = odict.OrderedDict()
        cuts['mll_max'] = 'mll < %.1f'                 % mcuts['mllmax'][i] 
        cuts['pt1_min'] = 'pt1 > %.1f'                 % mcuts['pt1min'][i]
        cuts['pt2_min'] = 'pt2 > %.1f'                 % mcuts['pt2min'][i]
        cuts['dphill']  = 'dphill < %s'                % mcuts['dphimax'][i]
        cuts['mth']     = '(mth > %.1f && mth < %.1f)' % (mcuts['mtmin'][i], mcuts['mtmax'][i])

        return cuts

    @staticmethod
    def vbfcb( mass ):
        
        mcuts = wwnamedcuts.mcuts
        mthmin = 30.
        i = wwnamedcuts.masses.index(mass)

        cuts = odict.OrderedDict()
        cuts['mll_max'] = 'mll < %s'                   % mcuts['mllmax'][i] 
        cuts['pt1_min'] = 'pt1 > %.1f'                 % mcuts['pt1min'][i]
        cuts['pt2_min'] = 'pt2 > %.1f'                 % mcuts['pt2min'][i]
        cuts['dphill']  = 'dphill < %s'                % mcuts['dphimax'][i]
        cuts['mth']     = '(mth > %.1f && mth < %.1f)' % (mthmin, mcuts['mtmax'][i])

        return cuts

    @staticmethod
    def hww0j1jcbfull( mass ):
        cuts  = wwnamedcuts.wwcommon.copy()
        cuts.update( wwnamedcuts.hww0j1jonly( mass ) )
        return cuts

    @staticmethod
    def vbfcbfull( mass ):
        cuts  = wwnamedcuts.vbfcutbased.copy()
        cuts.update( wwnamedcuts.vbfcb( mass ) )
        return cuts
        


class wwcuts:
    wwcommon = [
        'pt1>20',
        'pt2>10',
        '(ch1*ch2)<0.5', # don't use 0
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
        '(ch1*ch2)<0.5', # don't use 0
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
        'njet<4',
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
                
    zerojet = 'njet == 0'
    onejet  = 'njet == 1'
    vbf     = '(njet >= 2 && njet <= 3 && (jetpt3 <= 30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) '

# da rifare
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

#     vbflocut   = wwcuts.wwlo+_massindep+_cut
#     vbfhicut   = wwcuts.wwhi+_massindep+_cut
    vbfcut     = wwcuts.wwcommon+_massindep+_cut
    #vbfshape   = wwcuts.wwcommon+_massindep+_shape
    #vbfloshape = wwcuts.wwlo+_massindep+_shape
    #vbfhishape = wwcuts.wwhi+_massindep+_shape

    vbfshape = ['trigger==1. && pfmet>20. && mll>12 && zveto==1 && mpmet>20. && (njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  ) && bveto_mu==1 && nextra==0 && (bveto_ip==1 &&  (nbjettche==0 || njet>3))  && ptll>45. &&   ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) ) && (njet>=2 && njet<=3 && (jetpt3<=30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))))   && abs(eta1 - (jeteta1+jeteta2)/2)/detajj < 0.5 && abs(eta2 - (jeteta1+jeteta2)/2)/detajj < 0.5      && detajj>3.5     && mjj>500']

    specificvbfloshape = ['pt1>20   &&      pt2>0   &&      mth>30  &&  mth<280      &&    mll<200 ']
    specificvbfhishape = ['pt1>50   &&      pt2>0   &&      mth>30  &&  mth<680      &&    mll<600 ']

    vbfloshape = vbfshape + specificvbfloshape
    vbfhishape = vbfshape + specificvbfhishape



class vhcuts:
    _massindepincommonwithvbf = [
        'pfmet>20.',
        'mll>12 ',
        'zveto==1',
        'mpmet>20.',
        '(njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  )',
        'bveto_mu==1 ',
        'nextra==0' ,
        'ptll>45.',
        '( !sameflav ||  (pfmet > 45.0))'
    ]

    _massindep = [
        'abs(jeteta1)<2.5    && abs(jeteta2)<2.5',
        'drll < 1.3',
        'njet==2',
        '(mjj<105 && mjj>65)',
        'detajj<2.1',
        '(bveto_ip==1 &&   jettche1<1.6 && jettche2<1.6 )',
        'pt2>15',
        'sqrt((jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2)))/sqrt((pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))*(pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))+(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))*(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi)))>0.75',
        'sqrt((jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))*(jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))+(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2))*(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2)))/sqrt((pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))*(pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))+(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi))*(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi)))<1.5',
        '(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))>3.14159266)*(2*3.14159265-abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi))))+(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))<3.14159266)*(abs(atan2(jetpt1*sin(jetphi1)+jetpt2*sin(jetphi2),jetpt1*cos(jetphi1)+jetpt2*cos(jetphi2))-atan2(pt1*sin(phi1)+pt2*sin(phi2)+pfmet*sin(pfmetphi),pt1*cos(phi1)+pt2*cos(phi2)+pfmet*cos(pfmetphi)))) > 150*pi/180.'
    ]

    vhcut     = wwcuts.wwcommon+_massindepincommonwithvbf+_massindep




masses = [110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]
#masses = [110, 115, 120, 125, 130, 135, 140, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]  

categoryCuts = {}
categoryCuts['0j'] = wwcuts.zerojet
categoryCuts['1j'] = wwcuts.onejet
categoryCuts['2j']   = wwcuts.vbf        # 2 or 3 jets, but the third not between the first two in \eta
categoryCuts['vh2j'] = wwcuts.vbf        # 2 or 3 jets, but the third not between the first two in \eta


flavorCuts = {}
flavorCuts['all'] = '1'			       #'channel>-1'
flavorCuts['mm']  = 'channel == 0'     #'channel>-0.5 && channel<0.5'
flavorCuts['ee']  = 'channel == 1'     #'channel> 0.5 && channel<1.5'
flavorCuts['em']  = 'channel == 2'     #'channel> 1.5 && channel<2.5'
flavorCuts['me']  = 'channel == 3'     #'channel> 2.5 && channel<4.5'
flavorCuts['sf']  = 'channel < 1.5'
flavorCuts['of']  = 'channel > 1.5' 

flavors = {}
flavors['sf']=['mm','ee']
flavors['of']=['em','me']
flavors['ll']=['mm','ee','em','me']

channels = {}
channels['0j']    = ('0j','ll')
channels['1j']    = ('1j','ll')
channels['of_0j'] = ('0j','of')
channels['of_1j'] = ('1j','of')
channels['sf_0j'] = ('0j','sf')
channels['sf_1j'] = ('1j','sf')

channels['2j']    = ('2j','ll')
channels['of_2j'] = ('2j','of')
channels['sf_2j'] = ('2j','sf')

channels['of_vh2j'] = ('vh2j','of')
channels['sf_vh2j'] = ('vh2j','sf')


#  __  __               ___     _      
# |  \/  |__ _ ______  / __|  _| |_ ___
# | |\/| / _` (_-<_-< | (_| || |  _(_-<
# |_|  |_\__,_/__/__/  \___\_,_|\__/__/
#
# smurfs index              1     1     3     6     8     9    10    10    11    11    12    13    14    15    16   20     21    22    23    24    25    26    27
# masses              = [ 110 , 115 , 120 , 125 , 130 , 135 , 140 , 145 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600]
#  normal index             0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22


#  normal index             1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23

cutmap = {}
cutmap['mtmin_vh']    = [ 50  , 50  , 50  , 50  , 50  , 50  , 50  ,  50 ,  50 ,  50 ,  50 ,  50 ,  60 ,  60 ,  60 ,  60 ,  60 , 100 , 100 , 100 , 100 , 100 , 100 ]
cutmap['mllmax_vh']   = [ 70  , 70  , 70  , 80  , 80  , 90  , 90  , 100 , 100 , 100 , 100 , 100 , 110 , 120 , 130 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

cutmap['mllmax_bdt']  = [ 70  , 70  , 70  , 80  , 80  , 90  , 90  , 100 , 100 , 100 , 100 , 100 , 110 , 120 , 130 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]
cutmap['pt1min']      = [ 20  , 20  , 20  , 23  , 25  , 25  , 25  , 25  , 27  , 27  , 30  , 34  , 36  , 38  , 40  , 55  , 70  , 80  , 90  , 110 , 120 , 130 , 140 ]
cutmap['pt2min']      = [ 10  , 10  , 10  , 10  , 10  , 12  , 15  , 15  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  , 25  ]
cutmap['mllmax']      = [ 40  , 40  , 40  , 43  , 45  , 45  , 45  , 45  , 50  , 50  , 50  , 50  , 60  , 80  , 90  , 150 , 200 , 250 , 300 , 350 , 400 , 450 , 500 ]
cutmap['dphimax']     = [ 115 , 115 , 115 , 100 , 90  , 90  , 90  , 90  , 90  , 90  , 60  , 60  , 70  , 90  , 100 , 140 , 175 , 175 , 175 , 175 , 175 , 175 , 175 ]

cutmap['mtmin']       = [ 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 80  , 90  , 110 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 , 120 ]
cutmap['mtmax']       = [ 110 , 110 , 120 , 123 , 125 , 128 , 130 , 140 , 150 , 155 , 160 , 170 , 180 , 190 , 200 , 250 , 300 , 350 , 400 , 450 , 500 , 550 , 600 ]

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

    # HCP
    #mthmin_2dlomass = 80. 
    #mthmax_2dlomass = 280.
    #mllmax_2dlomass = 200

    # Spin
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
    
    sel['vbf-shape']    = ' && '.join(vbfcuts.vbfshape)
    sel['vbf-level']    = ' && '.join(vbfcuts.vbfcut)

    sel['ww-level']     = sel['ww-common']+'&& ptll>45'
    sel['wwbtag-level'] = sel['wwbtag-common']+'&& ptll>45'
    sel['ww2011-level'] = sel['ww2011-common']+'&& ptll>45'
    sel['ww2011btag-level'] = sel['ww2011btag-common']+'&& ptll>45'
    sel['bdt-specific'] = 'mll < {0} && (mth > {1:.0f} && mth < {2:.0f})'.format(masscuts['mllmax_bdt'], mthmin_bdt, int(mass))

    hwwlvl = {}
    hwwlvl['mll']    = 'mll < {0}'.format(masscuts['mllmax'])
    hwwlvl['pt1']    = 'pt1 > {0:.1f}'.format(masscuts['pt1min'])
    hwwlvl['pt2']    = 'pt2 > {0:.1f}'.format(masscuts['pt2min'])
    hwwlvl['dphill'] = 'dphill < {0}'.format(masscuts['dphimax'])
    hwwlvl['mth']    = '(mth > {0:.1f} && mth < {1:.1f})'.format(masscuts['mtmin'], masscuts['mtmax'])

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

    #sel['vbf-selection']          = sel['vbf-shape']+' && (mth > 50 && mth < {0:.0f})'.format(int(mass))
    #sel['vbf-selection']          = sel['vbf-shape']+' && (mth > 50 && mth < {0:.0f}) && mll > {0:.0f}    {0:.0f} {1:.0f} '.format(int(mass), float(12452))

    sel['vbf-selection-temp']     = ' && '.join([sel['vbf-level']]+[cut for var,cut in hwwlvl.iteritems() if var != 'mth'])
    sel['vbf-selection']          = sel['vbf-selection-temp'] + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(mthmin_vbf, masscuts['mtmax'])
    #sel['vbf-selection']          = sel['vbf-selection-temp'] + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(mthmin_vbf, int(mass))


    sel['shape-selection']        = sel['ww-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-himass']
    sel['shapebtag-selection']    = sel['wwbtag-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['wwbtag-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-himass']
    sel['shape2011-selection']     = sel['ww2011-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-himass']
    sel['shape2011btag-selection'] = sel['ww2011btag-common'].replace("zveto==1", "zveto>-1")+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011btag-common'].replace("zveto==1","zveto>-1")+' && '+sel['shape-himass']
        
    #sel['shape-selection']        = sel['ww-common']+' && '+sel['shape-lomass'] if mass <=250 else sel['ww-common']+' && '+sel['shape-himass']
    #sel['shape2011-selection']     = sel['ww2011-common']+' && '+sel['shape-lomass'] if mass <=250 else sel['ww2011-common']+' && '+sel['shape-himass']
   
    sel['vbf-shape-selection']    = sel['vbf-shape']+' && (mth > 50 && mth < {0:.0f})'.format(int(mass))

    # vh #
    sel['vh-level']    = ' && '.join(vhcuts.vhcut)
    sel['vh-selection'] = sel['vh-level']     + ' && (mth > {0:.1f} && mth < {1:.1f})'.format(masscuts['mtmin_vh'], int(mass))
    sel['vh-selection'] = sel['vh-selection'] + ' && (mll < {0:.1f})'.format(masscuts['mllmax_vh'])


    return sel




