#!/usr/bin/env python

import json
import sys
import ROOT
import optparse
import hwwinfo
import hwwsamples
import hwwtools
import os.path
import string
import logging
from HWWAnalysis.Misc.odict import OrderedDict
from HWWAnalysis.Misc.ROOTAndUtils import TH1AddDirSentry
import traceback
from array import array

pathWght = os.environ['CMSSW_BASE']+"/src/HWWAnalysis/ShapeAnalysis/ewksinglet/"
ROOT.gROOT.ProcessLineSync('.L '+pathWght+'getCPSWght.C+')
ROOT.gROOT.ProcessLineSync('.L '+pathWght+'getBWWght.C+')
ROOT.gROOT.ProcessLineSync('.L '+pathWght+'getIntWght.C+')

# ----------------------------------------------------- Read YR values from combination area --------------

def file2map(x):
        ret = {}; headers = []
        for x in open(x,"r"):
            cols = x.split()
            if len(cols) < 2: continue
            if "mH" in x:
                headers = [i.strip() for i in cols[1:]]
            else:
                fields = [ float(i) for i in cols ]
                ret[fields[0]] = dict(zip(headers,fields[1:]))
        return ret

# ----------------------------------------------------- Add point if needed to YR ------------------------

def GetYRVal(YRDic,iMass,Key):

    if iMass in YRDic :
       #print iMass,YRDic[iMass][Key]
       return YRDic[iMass][Key]
    else:
       n=len(YRDic.keys())
       x=[]
       y=[]
       for jMass in sorted(YRDic.keys()):
         x.append(jMass)
         y.append(YRDic[jMass][Key])
       gr = ROOT.TGraph(n,array('f',x),array('f',y));
       sp = ROOT.TSpline3("YR",gr);
       #print iMass,sp.Eval(iMass)
       return sp.Eval(iMass)

# ----------------------------------------------------- ShapeFactory --------------------------------------

class ShapeFactory:
    _logger = logging.getLogger('ShapeFactory')
 
    # _____________________________________________________________________________
    def __init__(self):
        self._stdWgt = 'baseW*puW*effW*triggW'
        self._systByWeight = {}

        ranges = {}
        ranges['counting']         = (1   , 0.  , 2.)
        ranges['bdtl']             = (400 , -1. , 1.)
        ranges['bdts']             = (400 , -1. , 1.)
        ranges['mth']              = (400 , 0.  , 200)
        ranges['dphill']           = (400 , 0.  , 3.15)
        ranges['detajj']           = (240 , 0.  , 6.)
        ranges['mll-vbf']          = (60  , 12  , 135)
        ranges['mll']              = self._getmllrange
        ranges['mllsplit']         = self._getmllsplitrange
        ranges['gammaMRStar']      = self._getGMstarrange
        ranges['vbf2D']            = self._getVBF2Drange
        ranges['mth-mll-hilomass'] = self._getMllMth2Drange
        ranges['mth-mll-hilospin'] = self._getMllMth2DSpinrange
        ranges['mth-mll-hilospin-withControlRegion']  = self._getMllMth2DSpinrangeWithControlRegion
        ranges['mth-mll-hilospin-withSSmirrorRegion'] = self._getMllMth2DSpinrangeWithSSmirrorRegion
        ranges['vbfMll-2011-range']  = self._getMll2011VBFrange
        ranges['vbfMll-range']       = self._getMllVBFrange
        ranges['vbfMll-fish-range']  = self._getMllVBFFishrange
        ranges['vhMll-range']        = self._getMllVHrange
        ranges['vhMllBanana-range']  = self._getMllVHrangeWithControlRegion
        ranges['vbfMllBanana-range'] = self._getMllVBFrangeWithControlRegion
        ranges['whsc-range']         = self._getWHSCrange
        ranges['wwewk-range']        = self._getWWewkrange
        ranges['wwewk-range-top']    = self._getWWewkrangeTop
        ranges['Hwidth-range']       = self._getHwidthrange
        ranges['HwidthSimple-range'] = self._getSimpleHwidthrange
        ranges['HwidthSimple7TeV-range'] = self._getSimpleHwidthrange7TeV

        self._ranges = ranges

        self._dataTag         = '2012A'
        self._sigTag          = 'SM'
        self._mcTag           = '0j1j'
        self._masses          = []
        self._channels        = {}
        # paths (to move out)
        self._outFileFmt      = ''
        self._paths           = {}
        self._range           = None
        self._splitmode       = None
        self._lumi            = 1
        self._muVal           = '1.'

        self._mh_SM           = 125.
        self._newcps          = False
        self._ewksinglet      = False
        self._cprimesq        = 1.
        self._brnew           = 0.
        self._approxewk       = False

        self._YR3rewght       = False
        self._YRValues        = {}

        variables = {}
        variables['2dWithCR']             = self._getMllMth2DSpinWithControlRegion
        variables['2dWithSSmirrorRegion'] = self._getMllMth2DSpinWithSSmirrorRegion
        variables['vhMllBanana']          = self._getMllVHWithControlRegion
        variables['vbfMllBanana']         = self._getMllVBFWithControlRegion
        variables['wwewk']                = self._getVarWWewk
        variables['Hwidth']               = self._getVarHwidth
        variables['HwidthSimple']         = self._getSimpleVarHwidth

        self._variables = variables


    # _____________________________________________________________________________
    def __del__(self):
        pass

# _____________________________________________________________________________
    def getvariable(self,tag,mass,cat):

        if tag in self._variables :
            try:
                theVariable = (self._variables[tag])(mass,cat)
            except KeyError as ke:
                self._logger.error('Variable '+tag+' not available. Possible values: '+', '.join(self._variables.iterkeys()) )
                raise ke
        else :
            theVariable = tag

        return theVariable

        #if isinstance(tag,tuple):
            #theVariable = tag
        #else:
            #try:
                #theVariable = self._variables[tag]
            #except KeyError as ke:
                #self._logger.error('Variable '+tag+' not available. Possible values: '+', '.join(self._variables.iterkeys()) )
                #raise ke

        #if isinstance(theVariable,tuple):
            #return theVariable
        #elif isinstance(theVariable,dict):
            #return theVariable[mass][cat]
        #elif callable(theVariable):
            #return theVariable(mass,cat)

    # _____________________________________________________________________________
    def _getHwidthrange(self,mass,cat):

        if cat not in ['0j','1j','2j','01j']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')

        if cat in ['1j'] :
          #   HwidthMVAggH in x,     HwidthMVAbkg in y
          #return ([12, 30, 45, 60, 70, 100,          300],[-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00])
          #return ([12, 30, 45,     70, 100,          300],[-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70,       1.00])
          return ([12, 30, 45,     70,               300],[-1.00,        -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70,       1.00])
        elif cat in ['0j'] :
          return ([12, 30, 45, 60, 70, 100,          300],[-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00])
          # test:
          #return ([12,          60, 70, 100,          300],[-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00])
        elif cat in ['01j'] :
          return ([12, 30, 45, 60, 70, 100,          300],[-1.00, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00])
        else :
          # mth:mll                    mll in x                    mth in y
          #return ([12,  45, 70,           100, 300],[0, 50,     100, 130, 170, 250, 450])
          return ([12,      70,                300],[0, 50,     100, 130, 170,      450])


    # _____________________________________________________________________________
    def _getSimpleHwidthrange(self,mass,cat):

        if cat not in ['0j','1j','2j','01j']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')


        if cat in ['0j','1j'] :
         # mth:mll                    mll in x                    mth in y
          #return ([12, 30, 45, 60, 70, 100, 150, 200, 300],[0, 50, 70, 90, 100, 110, 120, 130, 140, 160, 180, 200, 450])
          return ([12, 30, 45,     70, 100, 150, 200, 300],[0, 50, 70, 90, 100, 110, 120, 130, 140, 160,      200, 450])
        else :
          # mth:mll                    mll in x                    mth in y
          #return ([12,  45, 70,           100, 300],[0, 50,     100, 130, 170, 250,  450])
          return ([12,      70,                300],[0, 50,     100, 130, 170,      450])

    # _____________________________________________________________________________
    def _getSimpleHwidthrange7TeV(self,mass,cat):

        if cat not in ['0j','1j','2j','01j']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')


        if cat in ['0j','1j'] :
         # mth:mll                    mll in x                    mth in y
          #return ([12, 30, 45, 60, 70, 100, 150, 200, 300],[0, 50, 70, 90, 100, 110, 120, 130, 140, 160, 180, 200, 450])
          #return ([12, 30, 45,     70, 100, 150, 200, 300],[0, 50, 70, 90, 100, 110, 120, 130, 140, 160,      200, 450])
          return ([12, 30, 45,     70, 100, 150, 200, 300],[0, 50,      90,      110,      130,      160,      200, 450])
        else :
          # mth:mll                    mll in x                    mth in y
          #return ([12,  45, 70,           100, 300],[0, 50,     100, 130, 170, 250,  450])
          #return ([12,      70,                300],[0, 50,           130, 170,      450])
          #return ([12,      70,                300],[0, 80,           130, 170,      450])
          return ([12,      70,                300],[0, 80,           130,           450])



    # _____________________________________________________________________________
    def _getWWewkrange(self,mass,cat):

        if cat not in ['2j','2jtche05','2jtche05CJ','2jtche05FJ']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')

        #return ([-1.0, -0.5, 0.0, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],)
        #return ([-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],)
        #return ([-1.0, -0.75, -0.5, -0.25, 0.0, 0.30, 0.50, 0.70, 0.90, 1.0],)

        #return ([-1.0, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00],)
        #return ([-1.0, -0.50, -0.00, 0.30, 0.60, 1.00],)
        return ([-1.0, -0.50, -0.00, 0.50, 1.00],)

        #                 mva jets                             mll
        #return ([-1.0, -0.50, -0.00, 0.30, 0.60, 1.00],[50.00, 150.00, 250.00])



    # _____________________________________________________________________________
    def _getWWewkrangeTop(self,mass,cat):

        if cat not in ['2j','2jtche05','2jtche05CJ','2jtche05FJ']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')

        #return ([-1.0, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00],)
        #return ([-1.0, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],)
        #return ([-1.0, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00],)

        return ([-1.0, -0.50, -0.00, 0.50, 1.00],)


    # _____________________________________________________________________________
    def _getWHSCrange(self,mass,cat):

        if cat not in ['whsc']:
            print cat
            raise RuntimeError('range for '+str(cat)+' not defined. !?!?!?')

        #return ([0,20,40,60,80,100,120,140,160,180,250,400],)
        return ([0,40,80,120,160,200,240,300,400],)




    # _____________________________________________________________________________
    def _getVarWWewk(self,mass,cat):

        return 'WWewkMVABDTG'
        #return 'mll:WWewkMVABDTG'


    # _____________________________________________________________________________
    def _getMllVBFWithControlRegion(self,mass,cat):

        if cat not in ['2j']:
            print cat
            raise RuntimeError('mll range in VBF for '+str(cat)+' not defined. Must be 2')

        return 'mll*((ch1*ch2)<0)-5*((ch1*ch2)>0)'


    # _____________________________________________________________________________
    def _getMllVHWithControlRegion(self,mass,cat):

        if cat not in ['vh2j']:
            print cat
            raise RuntimeError('mll range in VH for '+str(cat)+' not defined. Must be 2')

        return 'mll*((ch1*ch2)<0)+220*((ch1*ch2)>0)'


    # _____________________________________________________________________________
    def _getMllMth2DSpinWithControlRegion(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return 'mll*((ch1*ch2)<0)+13*((ch1*ch2)>0):mth*((ch1*ch2)<0)+290*((ch1*ch2)>0)'
        else:
            return 'mll*((ch1*ch2)<0)+13*((ch1*ch2)>0):mth*((ch1*ch2)<0)+290*((ch1*ch2)>0)'


    # _____________________________________________________________________________
    def _getMllMth2DSpinWithSSmirrorRegion(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return 'mll:mth*((ch1*ch2)<0)-mth*((ch1*ch2)>0)'
        else:
            return 'mll:mth*((ch1*ch2)<0)-mth*((ch1*ch2)>0)'



    # _____________________________________________________________________________
    def _getVarHwidth(self,mass,cat):

        if cat not in ['0j','1j','2j','01j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if cat in ['0j','1j','01j'] :
           #return 'HwidthMVAbkg:HwidthMVAggH'
           #             (mth-65)/65 vs HwidthMVAggH in low HwidthMVAggH region                            HwidthMVAbkg vs HwidthMVAggH if high HwidthMVAggH
           #return '((HwidthMVAbkg*(HwidthMVAggH>=0.0))+((HwidthMVAggH<0.0)*(mll<50)*(mth-65.)/65.)):(HwidthMVAggH*((HwidthMVAggH>=0.0)+((HwidthMVAggH<0.0)&&(mll<50))))'
           #             (mth-65)/65 vs HwidthMVAggH in low HwidthMVAggH region                            HwidthMVAbkg vs HwidthMVAggH if high HwidthMVAggH
           #return '((HwidthMVAbkg*(HwidthMVAggH>=0.0))+((HwidthMVAggH<0.0)*(mth-65.)/65.)):HwidthMVAggH'
           #             (mll-50)/50 vs HwidthMVAggH in low HwidthMVAggH region                            HwidthMVAbkg vs HwidthMVAggH if high HwidthMVAggH
           #return '((HwidthMVAbkg*(HwidthMVAggH>=0.0))+((HwidthMVAggH<0.0)*(mll-40.)/40.)):HwidthMVAggH'

           #return '((HwidthMVAbkg*(HwidthMVAggH>=0.0))+((HwidthMVAggH<0.0)*(mth-120.)/120.)):mll'
           #return '((HwidthMVAbkg*(HwidthMVAggH>=0.0))+((HwidthMVAggH<0.0)*(mth-200.)/200.)):mll'
           #return '((HwidthMVAbkg*(mll>=70))+((mll<70)*(mth-200.)/200.)):mll'
           return '((HwidthMVAbkg*(mll>=70))+((mll<70)*(mth-160.)/160.)):mll'

        else :
           return 'mth:mll'
           #return '(mth*((mll>=&&pt2>20&&pt1>20)+(mth<125&&ptll>20))):mll'
        #return '(mth*((mth>=120&&pt2>20&&pt1>20)+(mth<120&&ptll>20))):mll'
        #return '(mth*((mth>=120&&pt2>20&&pt1>20)+(mth<120&&ptll>45&&pfmet>30))):mll'
        #return '(mth*((mth>=130&&pt2>20&&pt1>20)+(mth<130&&ptll>45&&pfmet>30))):mll'
        #return '(mth*((mth>=130&&pt2>20&&pt1>20)+(mth<130&&ptll>30&&pfmet>30))):mll'
        #return 'mth:(mll*(mth>=130&&pt2>20&&pt1>10)-mll*(mth<130&&ptll>30&&pfmet>30))'
        #return 'mll*(mth>=125)-mll*(mth<125)'

    # _____________________________________________________________________________
    def _getSimpleVarHwidth(self,mass,cat):

        if cat not in ['0j','1j','2j','01j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if cat in ['0j','1j'] :
           return 'mth:mll'
        else :
           return 'mth:mll'



    # _____________________________________________________________________________
    def getrange(self,tag,mass,cat):

        if isinstance(tag,tuple):
            theRange = tag
        else:
            try:
                theRange = self._ranges[tag]
            except KeyError as ke:
                self._logger.error('Range '+tag+' not available. Possible values: '+', '.join(self._ranges.iterkeys()) )
                raise ke

        if isinstance(theRange,tuple):
            return theRange
        elif isinstance(theRange,dict):
            return theRange[mass][cat]
        elif callable(theRange):
            return theRange(mass,cat)

    # _____________________________________________________________________________
    def _getMllMth2Drange(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return (10,80,280,8,0,200) 
        elif mass < 700:
            return (10,80,380,8,0,450) 
        else:
            return (13,80,600,10,0,600)

    # _____________________________________________________________________________
    def _getMllMth2DSpinrange(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return ([60,70,80,90,100,110,120,140,160,180,200,220,240,260,280],[12,30,45,60,75,100,125,150,175,200])
        elif mass < 700.:
            return (10,80,380,8,0,450)
        else:
            #return (12,80,580,10,0,600)
            return ([80,110,140,170,200,230,260,290,320,350,380,600],[0,45,90,145,180,225,270,315,360,405,450,600])


    # _____________________________________________________________________________
    def _getMllMth2DSpinrangeWithSSmirrorRegion(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return ([-280,-240,-200,-160,-120,-100,-80,-60,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280],[12,30,45,60,75,100,125,150,175,200])
            #return ([-280,-240,-200,-160,-120,-100,-80,-60,0,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280],[12,30,45,60,75,100,125,150,175,200])
            #return ([-280,-260,-240,-220,-200,-180,-160,-140,-120,-110,-100,-90,-80,-70,-60,0,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280],[12,30,45,60,75,100,125,150,175,200])
        else:
            return (10,80,380,8,0,450)

    # _____________________________________________________________________________
    def _getMllMth2DSpinrangeWithControlRegion(self,mass,cat):

        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        if mass < 300.:
            return ([60,70,80,90,100,110,120,140,160,180,200,220,240,260,280,300],[12,30,45,60,75,100,125,150,175,200])
        else:
            return (10,80,380,8,0,450)


    # _____________________________________________________________________________
    def _getMllVBFrangeWithControlRegion(self,mass,cat):

        if cat not in ['2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        #return ([12,32,52,72,92,112,132,152,172,192,212,232],)
        #return ([12,20,40,60,80,100,150,200,230],)

        #return ([-10,0,12,30,50,70,90,120,150,200,250,300,350,400,500,600],)
        return ([-10,0,12,35,60,90,120,160,200,250,300,350,400,500,600],)

    # _____________________________________________________________________________
    def _getMll2011VBFrange(self,mass,cat):

        if cat not in ['2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        return ([12,45,75,100,150,200,250,300,350,400,600],)

    # _____________________________________________________________________________
    def _getMllVBFrange(self,mass,cat):

        if cat not in ['2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        return ([12,30,45,60,75,100,125,150,175,200,250,300,350,400,600],)


    # _____________________________________________________________________________
    def _getMllVBFFishrange(self,mass,cat):

        if cat not in ['2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        #return ([12,45,75,100,125,150,175,200,250,300,350,400,600],)
        #return ([12,30,45,60,75,100,125,150,175,200,250,300,350,400,600],)
        return ([12,45,60,75,100,125,150,175,200,250,300,350,400,600],)


    # _____________________________________________________________________________
    def _getMllVHrange(self,mass,cat):

        if cat not in ['vh2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        return ([12,30,45,60,75,100,125,150,175,200],)
        #return ([12,30,50,75,100,150,200],)


    # _____________________________________________________________________________
    def _getMllVHrangeWithControlRegion(self,mass,cat):

        if cat not in ['vh2j']:
            print cat
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')

        #return ([12,32,52,72,92,112,132,152,172,192,212,232],)
        #return ([12,20,40,60,80,100,150,200,230],)
        return ([12,30,50,75,100,150,200,230],)



    # _____________________________________________________________________________
    def _getVBF2Drange(self,mass,cat):

        if cat not in ['2j']:
            raise RuntimeError('mth:mll range for '+str(cat)+' not defined. Can be 2')

        if mass<300 :
          return (4, 30, 280, 4, 0, 200)
        else :
          return (2, 30, 330, 3, 0, 450)

    # _____________________________________________________________________________
    def _getmllrange(self,mass,cat):
        
        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')
        # TODO: cleanup

        # xmin
        xmin   = 0. if mass<300 else 0.2*mass-20      #;--> changed "(mH<300)" and "0.2*float(mH) - 20"
        # xmax
        if mass < 200.:
            xmax = 200.
        elif mass < 400.:
            xmax = mass
        else:
            xmax = mass-50.
        # bins
        if mass < 300.:
            bins = 400 if cat == '0j' else 200
        else:
            bins = 300 if cat == '0j' else 150
        return (bins,xmin,xmax)

    # _____________________________________________________________________________
    def _getmllsplitrange(self,mass,cat):
        if cat not in ['0j','1j']:
            raise RuntimeError('mll range for '+str(cat)+' not defined. Can be 0 or 1')
        # xmax
        xmax = hwwinfo.massDependantCutsbyVar['mllmax_bdt'][mass]  
        # xmin
        xmin = -1.*xmax
        # bins
        if mass < 300.:
            bins = 400 if cat == 0 else 200
        else:
            bins = 300 if cat == 0 else 150
        return (bins,xmin,xmax)
    
    # _____________________________________________________________________________
    def _getGMstarrange(self,mass,cat):
        if cat not in ['0j','1j']:
            raise RuntimeError('GMstar range '+str(cat)+' not defined. Can be 0 or 1')
        # lower alwyas 50
        # upper 100+(mH-100)*0.5
        xmin=40
        xmax=90.+(mass-100.)*0.6

        if cat == '1j': xmax += 20

        if mass < 300.:
            bins = 200 if cat == '0j' else 200
        else:
            bins = 150 if cat == '0j' else 150
        return (bins,xmin,xmax)

    # _____________________________________________________________________________
    def makeNominals(self, var, sel, inputDir, outPath, **kwargs):

        print "======================"
        print "==== makeNominals ===="
        print "======================"
        
        ROOT.TH1.SetDefaultSumw2(True)
        shapeFiles = []

        # mass dependent sample list, can be in the mass loop
        for mass in self._masses:
            #samples = hwwsamples.samples(mass, self._dataTag, self._sigTag, self._mcTag)
            # mass and variable selection
            allCuts = hwwinfo.massSelections( mass )

            alias = var if not self._splitmode else var+'*(-1+2*('+allCuts[self._splitmode+'-selection']+') )'
            alias = alias if not 'btag' in alias else '(bveto_mu && bveto_ip && nbjettche==0)'

            try:
                varSelection = allCuts[sel+'-selection']
            except KeyError as ke:
                raise RuntimeError('Config error: '+str(ke))


            #inner  jet and flavor loops
            for chan,(category,flavor) in self._channels.iteritems():
#                 cat = hwwinfo.categories[category]
                flavors = hwwinfo.flavors[flavor]
                for flavor in flavors:
                    pars = dict([
                        ('mass',mass),
                        ('category',category),
                        ('flavor',flavor)
                    ])
                    print '-'*80
                    print ' Processing channel '+chan+': mass',mass,'category',category,'flavor',flavor
                    print '-'*80

                    # define samples here and remove DYLL from DF and DYTT from SF
                    samples = hwwsamples.samples(mass,self._energy,self._dataTag, self._sigTag, self._mcTag)
                    if (flavor=='em' or flavor=='me'):
                        if 'DYLL'              in samples: samples.pop('DYLL')
                        if 'DYLL-template'     in samples: samples.pop('DYLL-template')
                        if 'DYLL-templatesyst' in samples: samples.pop('DYLL-templatesyst')
                    if (flavor=='ee' or flavor=='mm'):
                        if 'DYTT'              in samples: samples.pop('DYTT')

                    #print '----> samples = ',samples
                    # - define the source paths 
                    activeInputPaths = ['base']
                    # - if the current var is listes among the known paths,
                    #   add it to the actives
                    if var in self._paths: activeInputPaths.append(var)

                    # - apply the pars of the current sample 
                    dirmap = {}
                    for path in activeInputPaths:
                        dirmap[path]=(self._paths[path]+'/'+inputDir).format( **pars )

                    print 'Input dir:',dirmap.values()

                    inputs = self._connectInputs(var,samples, dirmap)

                    # and the output path (might be par dependent as well)
                    output = outPath.format(**pars)
                    outdir = os.path.dirname(output)
                    if outdir:
                        self._ensuredir(outdir)

                    print '.'*80
                    print 'Output file:',output

                    # - now build the selection
                    # - make a separate function to contain the exceptions
                    catSel = hwwinfo.categoryCuts[category]
                    selection = varSelection+' && '+catSel+' && '+hwwinfo.flavorCuts[flavor]
                    selections = dict(zip(samples.keys(),[selection]*len(samples)))

                    # - check if there are "CHI" regions, where datadriven estimations have been performed
                    for vsample, vcut in selections.iteritems() :
                       if vsample.startswith('CHITOP-') :
                         self._logger.info('CHI-TOP changed')
                         selections[vsample] = hwwinfo.flavorCuts[flavor]


                    self._addweights(mass,var,'nominals',selections,category,sel,flavor)

                    print '.'*80
                    # - extract the histogram range
                    rng = self.getrange(opt.range,mass,category)

                    # - extract the histogram variable
                    doalias = self.getvariable(alias,mass,category)

                    # - to finally fill it
                    self._draw(doalias, rng, selections, output, inputs)
                    # - then disconnect the files
                    self._disconnectInputs(inputs)

                    shapeFiles.append(output)
        return shapeFiles

    # _____________________________________________________________________________
    def makeSystematics(self,var,sel,syst,mask,inputDir,outPath,**kwargs):

        print "========================="
        print "==== makeSystematics ===="
        print "========================="

        ROOT.TH1.SetDefaultSumw2(True)
        shapeFiles = []

        nicks = kwargs['nicks'] if 'nicks' in kwargs else None
        # mass dependent sample list, can be in the mass loop
        for mass in self._masses:
            # interference systematic only at high mass:
            if 'interferenceGGH' in syst and mass <  350 : continue
            if 'interferenceVBF' in syst and mass <  350 : continue
            if 'interference'    in syst and not self._newcps : continue
            if 'interference'    in syst and mass >= 350 : self._interfSyst = syst

            #samples = hwwsamples.samples(mass, self._dataTag, self._sigTag, self._mcTag)
            
            # mass and variable selection
            allCuts = hwwinfo.massSelections( mass )
            varSelection = allCuts[sel+'-selection']

            alias = var if not self._splitmode else var+'*(-1+2*('+allCuts[self._splitmode+'-selection']+') )'
            alias = alias if not 'btag' in alias else '(bveto_mu && bveto_ip && nbjettche==0)'

            #inner  jet and flavor loops
            for chan,(category,flavor) in self._channels.iteritems():
#                 cat = hwwinfo.categories[category]
                flavors = hwwinfo.flavors[flavor]
                for flavor in flavors:
                    print '-'*80
                    print ' Processing channel '+chan+': mass',mass,'category',category,'flavor',flavor
                    print '-'*80

                    pars = dict([
                        ('mass',mass),
                        ('category',category),
                        ('flavor',flavor),
                        ('syst',syst),
                    ])
                    pars['nick'] = nicks[syst] if nicks else syst

                    # define samples here and remove DYLL from DF and DYTT from SF
                    samples = hwwsamples.samples(mass, self._energy, self._dataTag, self._sigTag, self._mcTag)
                    # remove the dirname
                    for tag,files in samples.iteritems():
                        samples[tag] = map(os.path.basename,files)
                    if (flavor=='em' or flavor=='me'):
                        if 'DYLL'              in samples: samples.pop('DYLL')
                    if (flavor=='ee' or flavor=='mm'):
                        if 'DYTT'              in samples: samples.pop('DYTT')
                            
                    # - define the source paths 
                    activeInputPaths = ['base']
                    # - if the current var is listes among the known paths,
                    #   add it to the actives
                    if var in self._paths: activeInputPaths.append(var)

                    # - apply the pars of the current sample 
                    dirmap = {}
                    for path in activeInputPaths:
                        dirmap[path]=(self._paths[path]+'/'+inputDir).format( **pars )

                    #print 'Input dir:',dirmap.values()
                    inputs = self._connectInputs(var,samples, dirmap, mask)
                    # ---

                    # - and the output path (might be par dependent as well)
                    output = outPath.format(**pars)
                    outdir = os.path.dirname(output)
                    if output:
                        self._ensuredir(outdir)
                    print '.'*80
                    print 'Output file: ',output

                    # - now build the selection
                    catSel = hwwinfo.categoryCuts[category]
                    selection = varSelection+' && '+catSel+' && '+hwwinfo.flavorCuts[flavor]
                    selections = dict(zip(samples.keys(),[selection]*len(samples)))

                    # - check if there are "CHI" regions, where datadriven estimations have been performed
                    for vsample, vcut in selections.iteritems() :
                       if vsample.startswith('CHITOP-') :
                         self._logger.info('CHI-TOP changed')
                         selections[vsample] = hwwinfo.flavorCuts[flavor]


                    self._addweights(mass,var,syst,selections, category, sel,flavor)

                    print '.'*80
                    # - extract the histogram range
                    rng = self.getrange(opt.range,mass,category)

                    # - extract the histogram variable
                    doalias = self.getvariable(alias,mass,category) 

                    # - to finally fill it
                    self._draw(doalias, rng, selections ,output,inputs)
                    # - then disconnect the files
                    self._disconnectInputs(inputs)
                shapeFiles.append(output)
            self._interfSyst = 'NONE' 
        return shapeFiles
    
    # _____________________________________________________________________________
    def _draw(self, var, rng, selections, output, inputs):
        '''
        var :       the variable to plot
        selection : the selection to draw
        output :    the output file path
        inputs :    the process-input files map
        '''
        self._logger.info('Yields by process')
        print output
        outFile = ROOT.TFile.Open(output,'recreate')
        print outFile
        vdim = var.count(':')+1
#         hproto,hdim = ShapeFactory._projexpr(rng)
        # 3 items per dimention
        hdim = self._bins2dim( rng )

        if vdim != hdim:
            raise ValueError('The variable\'s and range number of dimensions are mismatching')

        #print 'selections = ', selections

        print 'var: '+var
        if 'WW' in selections : 
          print 'selection (for WW  as example): '+selections['WW']
        if 'WWlow' in selections : 
          print 'selection (for WWlow  as example): '+selections['WWlow']

        if 'WW1' in selections : 
          print 'selection (for WW1  as example): '+selections['WW1']

        if 'ggH' in selection :
          print 'selection (for ggH as example): '+selections['ggH']
        #print 'inputs = ', inputs

        for process,tree  in inputs.iteritems():
#             print ' '*3,process.ljust(20),':',tree.GetEntries(),
            print '    {0:<20} : {1:^9}'.format(process,tree.GetEntries()),
            # new histogram
            shapeName = 'histo_'+process
#             hstr = shapeName+hproto

            outFile.cd()

            # prepare a dummy to fill
            shape = self._makeshape(shapeName,rng)
            cut = selections[process]

            self._logger.debug('---'+process+'---')
            self._logger.debug('Formula: '+var+'>>'+shapeName)
            self._logger.debug('Cut:     '+cut)
            self._logger.debug('ROOTFiles:'+'\n'.join([f.GetTitle() for f in tree.GetListOfFiles()]))
            entries = tree.Draw( var+'>>'+shapeName, cut, 'goff')
#             print ' >> ',entries,':',shape.Integral()
            shape = outFile.Get(shapeName)
            shape.SetTitle(process+';'+var)


            if isinstance(shape,ROOT.TH2):
                shape2d = shape
                # puts the over/under flows in
                self._reshape( shape )
                # go 1d
                shape = self._h2toh1(shape2d)
                # rename the old
                shape2d.SetName(shape2d.GetName()+'_2d')
                shape2d.Write()
                shape.SetDirectory(outFile)

            print '>> {0:>9} : {1:>9.2f}'.format(entries,shape.Integral())
            shape.Write()
        outFile.Close()
        del outFile

    # _____________________________________________________________________________
    @staticmethod
    def _moveAddBin(h, fromBin, toBin ):
        if not isinstance(fromBin,tuple) or not isinstance(toBin,tuple):
            raise ValueError('Arguments must be tuples')

        dims = [h.GetDimension(), len(fromBin), len(toBin) ]

        if dims.count(dims[0]) != len(dims):
            raise ValueError('histogram and the 2 bins don\'t have the same dimension')
        
        # get bins
        b1 = h.GetBin( *fromBin )
        b2 = h.GetBin( *toBin )

        # move contents
        c1 = h.At( b1 )
        c2 = h.At( b2 )

        h.SetAt(0, b1)
        h.SetAt(c1+c2, b2)

        # move weights as well
        sumw2 = h.GetSumw2()

        w1 = sumw2.At( b1 )
        w2 = sumw2.At( b2 )

        sumw2.SetAt(0, b1)
        sumw2.SetAt(w1+w2, b2)


    def _reshape(self,h):
        if h.GetDimension() == 1:
            # nx = h.GetNbinsX()
            # ShapeFactory._moveAddBin(h, (0,),(1,) )
            # ShapeFactory._moveAddBin(h, (nx+1,),(nx,) )
            return
        elif h.GetDimension() == 2:
            nx = h.GetNbinsX()
            ny = h.GetNbinsY()

            for i in xrange(1,nx+1):
                ShapeFactory._moveAddBin(h,(i,0   ),(i, 1 ) )
                ShapeFactory._moveAddBin(h,(i,ny+1),(i, ny) )

            for j in xrange(1,ny+1):
                ShapeFactory._moveAddBin(h,(0,    j),(1, j) )
                ShapeFactory._moveAddBin(h,(nx+1, j),(nx,j) )

            # 0,0 -> 1,1
            # 0,ny+1 -> 1,ny
            # nx+1,0 -> nx,1
            # nx+1,ny+1 ->nx,ny

            ShapeFactory._moveAddBin(h, (0,0),(1,1) )
            ShapeFactory._moveAddBin(h, (0,ny+1),(1,ny) )
            ShapeFactory._moveAddBin(h, (nx+1,0),(nx,1) )
            ShapeFactory._moveAddBin(h, (nx+1,ny+1),(nx,ny) )

    @staticmethod
    def _h2toh1(h):
        import array
        
        if not isinstance(h,ROOT.TH2):
            raise ValueError('Can flatten only 2d hists')

        sentry = TH1AddDirSentry()

#         H1class = getattr(ROOT,h.__class__.__name__.replace('2','1'))

        nx = h.GetNbinsX()
        ny = h.GetNbinsY()

        h_flat = ROOT.TH1D(h.GetName(),h.GetTitle(),nx*ny,0,nx*ny)

        
        sumw2 = h.GetSumw2()
        sumw2_flat = h_flat.GetSumw2()

        for i in xrange(1,nx+1):
            for j in xrange(1,ny+1):
                # i,j must be mapped in 
                b2d = h.GetBin( i,j )
#                 b2d = h.GetBin( j,i )
#                 b1d = ((i-1)+(j-1)*nx)+1
                b1d = ((j-1)+(i-1)*ny)+1

                h_flat.SetAt( h.At(b2d), b1d )
                sumw2_flat.SetAt( sumw2.At(b2d), b1d ) 

        h_flat.SetEntries(h.GetEntries())
        
        stats2d = array.array('d',[0]*7)
        h.GetStats(stats2d)

        stats1d = array.array('d',[0]*4)
        stats1d[0] = stats2d[0]
        stats1d[1] = stats2d[1]
        stats1d[2] = stats2d[2]+stats2d[4]
        stats1d[3] = stats2d[3]+stats2d[5]

        h_flat.PutStats(stats1d)

        xtitle = h.GetXaxis().GetTitle()
        v1,v2 = xtitle.split(':') # we know it's a 2d filled by an expr like y:x
        xtitle = '%s #times %s bin' % (v1,v2)

        h_flat.GetXaxis().SetTitle(xtitle)

        return h_flat


    # _____________________________________________________________________________
    # add the weights to the selection
    def _addweights(self,mass,var,syst,selections,cat='',sel='',flavor='of'):
        sampleWgts =  self._sampleWeights(mass,var,cat,sel,flavor)
        print '--',selections.keys()
        for process,cut in selections.iteritems():
            wgt = self._stdWgt
            if process in sampleWgts:
                wgt = sampleWgts[process]
            
            if syst in self._systByWeight:
                wgt = wgt+'*'+self._systByWeight[syst]
   

            selections[process] = wgt+'*('+cut+')'

    def loadYR(self):
        if self._energy == '7TeV' or self._energy == '8TeV' :
          path = os.getenv('CMSSW_BASE')+'/src/HWWAnalysis/ShapeAnalysis/data/' 
          # load YR2
          self._YRValues['YR2'] = {}
          self._YRValues['YR2']['xs'] = {}
          self._YRValues['YR2']['br'] = {}
          self._YRValues['YR2']['xs']['ggH'] = file2map(path+'lhc-hxswg-YR2/sm/xs/'+self._energy+'/'+self._energy+'-ggH.txt')
          self._YRValues['YR2']['xs']['qqH'] = file2map(path+'lhc-hxswg-YR2/sm/xs/'+self._energy+'/'+self._energy+'-vbfH.txt')
          self._YRValues['YR2']['xs']['WH']  = file2map(path+'lhc-hxswg-YR2/sm/xs/'+self._energy+'/'+self._energy+'-WH.txt')
          self._YRValues['YR2']['xs']['ZH']  = file2map(path+'lhc-hxswg-YR2/sm/xs/'+self._energy+'/'+self._energy+'-ZH.txt')
          self._YRValues['YR2']['xs']['ttH'] = file2map(path+'lhc-hxswg-YR2/sm/xs/'+self._energy+'/'+self._energy+'-ttH.txt')
          self._YRValues['YR2']['br']['VV']  = file2map(path+'lhc-hxswg-YR2/sm/br/BR.txt')
          self._YRValues['YR2']['br']['ff']  = file2map(path+'lhc-hxswg-YR2/sm/br/BR1.txt')
          # load YR3
          self._YRValues['YR3'] = {}
          self._YRValues['YR3']['xs'] = {}
          self._YRValues['YR3']['br'] = {}
          self._YRValues['YR3']['xs']['ggH'] = file2map(path+'lhc-hxswg-YR3/sm/xs/'+self._energy+'/'+self._energy+'-ggH.txt')
          self._YRValues['YR3']['xs']['qqH'] = file2map(path+'lhc-hxswg-YR3/sm/xs/'+self._energy+'/'+self._energy+'-vbfH.txt')
          self._YRValues['YR3']['xs']['WH']  = file2map(path+'lhc-hxswg-YR3/sm/xs/'+self._energy+'/'+self._energy+'-WH.txt')
          self._YRValues['YR3']['xs']['ZH']  = file2map(path+'lhc-hxswg-YR3/sm/xs/'+self._energy+'/'+self._energy+'-ZH.txt')
          self._YRValues['YR3']['xs']['ttH'] = file2map(path+'lhc-hxswg-YR3/sm/xs/'+self._energy+'/'+self._energy+'-ttH.txt')
          self._YRValues['YR3']['br']['VV']  = file2map(path+'lhc-hxswg-YR3/sm/br/BR.txt')
          self._YRValues['YR3']['br']['ff']  = file2map(path+'lhc-hxswg-YR3/sm/br/BR1.txt')

    #
    # X. Janssen: Function to implement Higgs weigthing and replace old "kfw" from tree:
    #       * k-factor (2011 only)
    #       * New CPS   
    #       * EWK Singlet
    #       * Interference with WW (+systematics ?)
    #
    def _HiggsWgt(self,prodMode,mass,flavor):

       iSystGGH = 0
       if self._interfSyst == 'interferenceGGH_up'   : iSystGGH =  1
       if self._interfSyst == 'interferenceGGH_down' : iSystGGH = -1
       iSystVBF = 0
       if self._interfSyst == 'interferenceVBF_up'   : iSystVBF =  1
       if self._interfSyst == 'interferenceVBF_down' : iSystVBF = -1

       hWght = '1.'

       # Goto YR3 (from YR2)

       if self._YR3rewght :
         xs_Scale = 1.
         br_Scale = 1.  # YR2 stops at 600 for BR->WW !
         if self._energy == '7TeV' or self._energy == '8TeV' :
           if prodMode in ['ggH','qqH'] :
             print prodMode,mass 
             xs_Scale = GetYRVal(self._YRValues['YR3']['xs'][prodMode],mass,'XS_pb')/GetYRVal(self._YRValues['YR2']['xs'][prodMode],mass,'XS_pb')
             if mass <= 600 : br_Scale = GetYRVal(self._YRValues['YR3']['br']['VV'],mass,'H_WW')/GetYRVal(self._YRValues['YR2']['br']['VV'],mass,'H_WW')
           elif prodMode in ['WH','ZH','ttH'] :
             # Only up to 300 ofc
             if mass <= 300 :
               xs_Scale = GetYRVal(self._YRValues['YR3']['xs'][prodMode],mass, 'XS_pb')/GetYRVal(self._YRValues['YR2']['xs'][prodMode],mass, 'XS_pb')
               br_Scale = GetYRVal(self._YRValues['YR3']['br']['VV'],mass, 'H_WW')/GetYRVal(self._YRValues['YR2']['br']['VV'],mass ,'H_WW')
           elif prodMode == 'ggH_ALT' :
             xs_Scale = GetYRVal(self._YRValues['YR3']['xs']['ggH'],mass ,'XS_pb')/GetYRVal(self._YRValues['YR2']['xs']['ggH'],mass,'XS_pb')
             if mass <= 600 : br_Scale = GetYRVal(self._YRValues['YR3']['br']['VV'],mass,'H_WW')/GetYRVal(self._YRValues['YR2']['br']['VV'],mass,'H_WW')
           elif prodMode == 'qqH_ALT' :
             xs_Scale = GetYRVal(self._YRValues['YR3']['xs']['qqH'],mass, 'XS_pb')/GetYRVal(self._YRValues['YR2']['xs']['qqH'],mass, 'XS_pb')
             if mass <= 600 : br_Scale = GetYRVal(self._YRValues['YR3']['br']['VV'],mass, 'H_WW')/GetYRVal(self._YRValues['YR2']['br']['VV'],mass, 'H_WW')
           elif prodMode in ['ggH_SM','qqH_SM','WH_SM','ZH_SM','ttH_SM'] :
             prodModeStrip = prodMode.split('_')[0]
             massSM = 125
             xs_Scale = self._YRValues['YR3']['xs'][prodModeStrip][massSM]['XS_pb']/self._YRValues['YR2']['xs'][prodModeStrip][massSM]['XS_pb']
             br_Scale = self._YRValues['YR3']['br']['VV'][massSM]['H_WW']/self._YRValues['YR2']['br']['VV'][massSM]['H_WW']
           else:
             print 'YR3 reweight: UNKNW prodMode : '+prodMode
         hWght += '*'+str(xs_Scale)+'*'+str(br_Scale)

       # Scale _SM samples at exact mass  

       if prodMode in ['ggH_SM','qqH_SM','WH_SM','ZH_SM','ttH_SM'] :
         xs_Scale = 1.
         br_Scale = 1.
         if self._energy == '7TeV' or self._energy == '8TeV' :
           YR = 'YR2'
           if self._YR3rewght : YR = 'YR3'
           prodModeStrip = prodMode.split('_')[0]
           massSMSample  = 125.
           massMeasured  = self._mHSM
           xs_Scale = GetYRVal(self._YRValues[YR]['xs'][prodModeStrip],massMeasured, 'XS_pb')/GetYRVal(self._YRValues[YR]['xs'][prodModeStrip],massSMSample,'XS_pb')
           br_Scale = GetYRVal(self._YRValues[YR]['br']['VV'],massMeasured, 'H_WW')/GetYRVal(self._YRValues[YR]['br']['VV'],massSMSample,'H_WW')
         hWght += '*'+str(xs_Scale)+'*'+str(br_Scale)

       # New CPS (+ ptHiggs k-factor for 7 TeV MC)
       if self._newcps and mass >= 250 :
         if self._energy == '7TeV' : 
           hWght += '*kfW'
           if prodMode in ['ggH','ggH_ALT'] : 
             old_interf = {
               250 : 1.0,
               300 : 1.0,
               350 : 1.0,
               400 : 0.9979665200,
               450 : 1.0264436689,
               500 : 1.0932697009,
               550 : 1.1622655466,
               600 : 1.2697380560,
               700 : 1.5160605911,
               800 : 1.8482963461,
               900 : 2.2618724091,
               1000 : 2.3417787303
             }
             hWght += '/'+str(old_interf[mass])


         fileWght = os.environ['CMSSW_BASE']+"/src/HWWAnalysis/ShapeAnalysis/ewksinglet/data/cpsWght/cpsWght.root"
         ROOT.gROOT.ProcessLineSync('initCPSWght("'+fileWght+'","'+self._energy+'",'+str(mass)+','+str(int(self._mh_SM))+')')
 
         if prodMode in ['ggH','ggH_ALT'] : hWght += '*getCPSWght(0,MHiggs)' 
         if prodMode in ['qqH','qqH_ALT'] : hWght += '*getCPSWght(1,MHiggs)' 
       else:
         if prodMode in ['ggH','ggH_ALT','ggH_SM','qqH','qqH_ALT','qqH_SM'] : hWght += '*kfW'

       # EWK Singlet 
       if self._ewksinglet : 

         # ... Change high mass H width (only ggH and qqH for now)
         if mass >= 250 :
           fileBWParam = os.environ['CMSSW_BASE']+"/src/HWWAnalysis/ShapeAnalysis/ewksinglet/data/BWParam.json"
           jsf = open(fileBWParam,"r+")
           BWParam = (json.loads(jsf.read()))
           jsf.close()  
  
           GamSM = 0.
           if    prodMode == 'ggH' :
             Mass  = BWParam['ggH'][str(mass)]['Mass']
             GamSM = BWParam['ggH'][str(mass)]['Gamma']
  
           elif  prodMode == 'qqH' :
             Mass  = BWParam['qqH'][str(mass)]['Mass']
             GamSM = BWParam['qqH'][str(mass)]['Gamma']
   
           Gamma = GamSM * self._cprimesq / (1.-self._brnew)
           
           #if prodMode in ['ggH','qqH']  and not self._approxewk : hWght += '*getBWWght(MHiggs,%f,%f,%f)'%(Mass,GamSM,Gamma)
           #if prodMode in ['ggH','qqH']  : hWght += '*getBWWght(MHiggs,%f,%f,%f)'%(Mass,GamSM,Gamma)
           if prodMode in ['ggH'] : hWght += '*getBWWght(MHiggs,%f,%f,%f)'%(Mass,GamSM,Gamma)
           if prodMode in ['qqH'] and mass < 350 : hWght += '*getBWWght(MHiggs,%f,%f,%f)'%(Mass,GamSM,Gamma)

         # ... Change mu of both all H
         if prodMode in ['ggH','WH','ZH','ttH']                      :  hWght += '*'+str(self._cprimesq)+'*(1-'+str(self._brnew)+')'
         if prodMode in ['qqH'] and mass < 350                       :  hWght += '*'+str(self._cprimesq)+'*(1-'+str(self._brnew)+')'
         if prodMode in ['ggH_SM','qqH_SM','WH_SM','ZH_SM','ttH_SM'] :  hWght += '*(1-'+str(self._cprimesq)+')'


       # Inteference (Only meaningfull with new CPS as kfW already contains it otherwise)
       if self._newcps :
         EWKcase='0' 
         if self._ewksinglet : EWKcase='1' 
         if prodMode in ['ggH'] and mass >= 300 : 
            fileInt = os.environ['CMSSW_BASE']+'/src/HWWAnalysis/ShapeAnalysis/ewksinglet/data/Interference_ggH/1.0SMWidth/h_MWW_IonS_NNLO_'+str(mass)+'.root'
            ROOT.gROOT.ProcessLineSync('initIntWght("'+fileInt+'",0,'+str(iSystGGH)+','+str(mass)+','+str(self._cprimesq)+','+str(self._brnew)+','+EWKcase+')')
            hWght += '*getIntWght(0,MHiggs,'+str(self._cprimesq)+','+str(self._brnew)+')' 
         if prodMode in ['qqH'] and mass >= 350 : 
            if   flavor in ['of','em','me'] : iFlavor = '0'
            elif flavor in ['sf','mm','ee'] : iFlavor = '1'
            else :
              print '_HiggsWgt: Unknown flavor : ',flavor
              exit()
            EWKDir = os.environ['CMSSW_BASE']+'/src/HWWAnalysis/ShapeAnalysis/ewksinglet/'
            print EWKDir
            ROOT.gROOT.ProcessLineSync('initIntWght("'+EWKDir+'" ,1,'+str(iSystVBF)+','+str(mass)+','+str(self._cprimesq)+','+str(self._brnew)+','+EWKcase+')') 
            if not self._approxewk : hWght += '*getIntWght(1,MHiggs,'+str(self._cprimesq)+','+str(self._brnew)+','+iFlavor+')' 
            else                   : hWght += '*getIntWght(1,MHiggs,1.0,0.0,'+iFlavor+')'



       return hWght

    # _____________________________________________________________________________
    # this is too convoluted
    # define here the mass-dependent weights
    def _sampleWeights(self,mass,var,cat,sel,flavor):
        weights = {}
        #print ">>>> sel = ", sel
        if sel in ['CutWW','Hwidth','HwidthWWcontrolRegion','HwidthWWcontrolRegionHighMllRemoved'] : # only for WW xsec for the time being and for Higgs width analysis
        #if sel in ['CutWW'] : # only for WW xsec for the time being
            #print " WW xsec "
            #                                                            pow                      mc@nlo                   MG             NLO x-sec     nnll weight
            #weights['WW']              = self._stdWgt+'*(((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.8123*1000./baseW*nllW'
            #                                                            pow                      mc@nlo                   MG            NNLO x-sec     nnll weight
            weights['WW']              = self._stdWgt+'*(((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            #
            weights['WWlow']           = self._stdWgt+'*(mll<70)  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WWhigh']          = self._stdWgt+'*(mll>=70) * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'

            weights['WW1']      = self._stdWgt+'*(mll>=12 && mll<30)  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WW2']      = self._stdWgt+'*(mll>=30 && mll<45)  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WW3']      = self._stdWgt+'*(mll>=45 && mll<60)  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WW4']      = self._stdWgt+'*(mll>=60 && mll<70)  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WW5']      = self._stdWgt+'*(mll>=70 && mll<100) * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'
            weights['WW6']      = self._stdWgt+'*(mll>=100         )  * (((dataset==6)*1./999860.)+((dataset==2)*1./539594.)+((dataset==0)*1./1933232.))*5.984*1000./baseW*nllW'


        # tocheck
        weights['WJet']              = self._stdWgt+'*kfW*fakeW*(run!=201191)'
        weights['WJetFakeRate-nominal']  = self._stdWgt+'*kfW*fakeW*(run!=201191)'
        weights['WJetFakeRate-eUp']  = self._stdWgt+'*kfW*fakeWElUp*(run!=201191)'
        weights['WJetFakeRate-eDn']  = self._stdWgt+'*kfW*fakeWElDown*(run!=201191)'
        weights['WJetFakeRate-mUp']  = self._stdWgt+'*kfW*fakeWMuUp*(run!=201191)'
        weights['WJetFakeRate-mDn']  = self._stdWgt+'*kfW*fakeWMuDown*(run!=201191)'


        weights['WJetFakeRate-2j-template']              = self._stdWgt+'*kfW*fakeW'
        weights['WJetFakeRate-2j-eUp']  = self._stdWgt+'*kfW*fakeWElUp*(run!=201191)'
        weights['WJetFakeRate-2j-eDn']  = self._stdWgt+'*kfW*fakeWElDown*(run!=201191)'
        weights['WJetFakeRate-2j-mUp']  = self._stdWgt+'*kfW*fakeWMuUp*(run!=201191)'
        weights['WJetFakeRate-2j-mDn']  = self._stdWgt+'*kfW*fakeWMuDown*(run!=201191)'

        weights['WJetSS']            = self._stdWgt+'*fakeW*ssW*(run!=201191)'

        weights['WJet-template']              = self._stdWgt+'*kfW*fakeW'
        weights['WJetFakeRate-template']      = self._stdWgt+'*kfW*fakeWUp'
        weights['WJet-templatesyst']          = self._stdWgt+'*kfW*fakeWUp'

        weights['Data']              = '(run!=201191)'
        # problem with DYTT using embedded for em/me, for ee/mm it is inlcuded in DD DY estimate
        weights['DYTT']              = self._stdWgt
        #weights['DYLL']              = self._stdWgt+'*(1-(( dataset == 36 || dataset == 37 ) && mctruth == 2 ))*(channel<1.5)'
        weights['DYLL']              = self._stdWgt+'*(channel<1.5)'

        # beware:
        # mumu #    channel == 0
        # mue #     channel == 3
        # emu #     channel == 2
        # ee #      channel == 1
        weights['DYee']              = self._stdWgt+'*(channel<1.5 && channel>0.5)'
        weights['DYmm']              = self._stdWgt+'*(channel<0.5 && channel>-0.5)'
        weights['DYLL-template']     = self._stdWgt+'* dyW *(1-(( dataset == 36 || dataset == 37 ) && mctruth == 2 ))'
        weights['DYLL-templatesyst'] = self._stdWgt+'*dyWUp*(1-(( dataset == 36 || dataset == 37 ) && mctruth == 2 ))'
        #systematics
        weights['TopTW']             = self._stdWgt+'*(1+0.17*(dataset>=11 && dataset<=16))' # 17% on tW/tt ratio
        weights['TopCtrl']           = self._stdWgt+'*bvetoW'
        weights['Top-template']      = self._stdWgt+'*bvetoW'
        #filter and k-factor on Vg* done by kfW
        weights['VgS']               = self._stdWgt+'*kfW'
        weights['Vg']                = self._stdWgt+'*kfW'


        weights['ggH']               = self._stdWgt+'*'+self._HiggsWgt('ggH',mass,flavor)+'*'+self._muVal
        weights['qqH']               = self._stdWgt+'*'+self._HiggsWgt('qqH',mass,flavor)+'*'+self._muVal
        weights['ggHminlo']          = 'effW*triggW*kfW*puW*HEPMCweight/497500.*1000*0.108*0.108*9*0.216'+self._muVal

        weights['WH']                = self._stdWgt+'*'+self._HiggsWgt('WH',mass,flavor)+'*(mctruth == 26)*'+self._muVal
        weights['ZH']                = self._stdWgt+'*'+self._HiggsWgt('ZH',mass,flavor)+'*(mctruth == 24)*'+self._muVal
        weights['ttH']               = self._stdWgt+'*'+self._HiggsWgt('ttH',mass,flavor)+'*(mctruth == 121)*'+self._muVal

        weights['ggH_SM']            = self._stdWgt+'*'+self._HiggsWgt('ggH_SM',mass,flavor)+'*'+self._muVal
        weights['qqH_SM']            = self._stdWgt+'*'+self._HiggsWgt('qqH_SM',mass,flavor)+'*'+self._muVal

        weights['WH_SM']             = self._stdWgt+'*'+self._HiggsWgt('WH_SM',mass,flavor)+'*(mctruth == 26)*'+self._muVal
        weights['ZH_SM']             = self._stdWgt+'*'+self._HiggsWgt('ZH_SM',mass,flavor)+'*(mctruth == 24)*'+self._muVal
        weights['ttH_SM']            = self._stdWgt+'*'+self._HiggsWgt('ttH_SM',mass,flavor)+'*(mctruth == 121)*'+self._muVal


        weights['ggH_ALT']           = self._stdWgt+'*kfW*'+self._muVal
        weights['qqH_ALT']           = self._stdWgt+'*kfW*'+self._muVal
        weights['jhu_NORM']          = self._stdWgt+'*kfW*'+self._muVal
        weights['jhu_NLO']           = self._stdWgt+'*kfW*'+self._muVal


        # fix for 2011
        #print "sel = ",sel
        if (sel == "whscShape2011") or (sel == "whsc2011") or (sel == "vh2011") :
          if mass == 120:
            weights['VH_SM']             = 'puW*effW*triggW*0.0007238280'
          if mass == 125:
            weights['VH_SM']             = 'puW*effW*triggW*0.0011041573*0.941792557006409292'   ##-> 0011041573   @126 -> scaled to 125
          if mass == 130:
            weights['VH_SM']             = 'puW*effW*triggW*0.0011893260'
          if mass == 135:
            weights['VH_SM']             = 'puW*effW*triggW*0.0015132260'
          if mass == 140:
            weights['VH_SM']             = 'puW*effW*triggW*0.0015174284'
          if mass == 150:
            weights['VH_SM']             = 'puW*effW*triggW*0.0016502763'
          if mass == 160:
            weights['VH_SM']             = 'puW*effW*triggW*0.0016590295'
          if mass == 170:
            weights['VH_SM']             = 'puW*effW*triggW*0.0014509147'
          if mass == 180:
            weights['VH_SM']             = 'puW*effW*triggW*0.0011345486'
          if mass == 190:
            weights['VH_SM']             = 'puW*effW*triggW*0.0007897267'

          if mass == 120:
            weights['VH']             = 'puW*effW*triggW*0.0007238280'
          if mass == 125:
            weights['VH']             = 'puW*effW*triggW*0.0011041573*0.941792557006409292'
          if mass == 130:
            weights['VH']             = 'puW*effW*triggW*0.0011893260'
          if mass == 135:
            weights['VH']             = 'puW*effW*triggW*0.0015132260'
          if mass == 140:
            weights['VH']             = 'puW*effW*triggW*0.0015174284'
          if mass == 150:
            weights['VH']             = 'puW*effW*triggW*0.0016502763'
          if mass == 160:
            weights['VH']             = 'puW*effW*triggW*0.0016590295'
          if mass == 170:
            weights['VH']             = 'puW*effW*triggW*0.0014509147'
          if mass == 180:
            weights['VH']             = 'puW*effW*triggW*0.0011345486'
          if mass == 190:
            weights['VH']             = 'puW*effW*triggW*0.0007897267'


        # for Higgs width measurements
        #                                              2.1 from LO -> NNLO scaling
        #weights['ggH_sbi']            = self._stdWgt+'*2.1' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_s']              = self._stdWgt+'*2.1' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_b']              = self._stdWgt+'*2.1' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        # using powheg to normalize on-shell contribution
        #weights['ggH_sbi']            = self._stdWgt+'*(((njet==0) * (13.3258/5.85323)) + ((njet==1) * (5.78547/1.40855)) + ((njet>=2) * (1.79911/0.195922)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_s']              = self._stdWgt+'*(((njet==0) * (13.3258/5.85323)) + ((njet==1) * (5.78547/1.40855)) + ((njet>=2) * (1.79911/0.195922)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_b']              = self._stdWgt+'*(((njet==0) * (13.3258/5.85323)) + ((njet==1) * (5.78547/1.40855)) + ((njet>=2) * (1.79911/0.195922)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_sbi']            = self._stdWgt+'*2.739*(((njet==0) * (1.1986-4.9341/sqrt(mWW))) + ((njet==1) * (0.8415+68.0300/mWW)) + ((njet>=2) * (0.7655+165.6779/mWW)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_s']              = self._stdWgt+'*2.739*(((njet==0) * (1.1986-4.9341/sqrt(mWW))) + ((njet==1) * (0.8415+68.0300/mWW)) + ((njet>=2) * (0.7655+165.6779/mWW)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'
        #weights['ggH_b']              = self._stdWgt+'*2.739*(((njet==0) * (1.1986-4.9341/sqrt(mWW))) + ((njet==1) * (0.8415+68.0300/mWW)) + ((njet>=2) * (0.7655+165.6779/mWW)))' # +'*((dataset == 37) - (dataset == 37) + (dataset == 37))'

        # 1/sqrt(x) functional form
        #weights['ggH_sbi']            = self._stdWgt+'*2.739*(((njet==0) * (1.1983-4.9696/sqrt(mWW))) + ((njet==1) * (0.6711+7.0075/sqrt(mWW))) + ((njet>=2) * (0.3794+16.4944/sqrt(mWW))))'
        #weights['ggH_s']              = self._stdWgt+'*2.739*(((njet==0) * (1.1983-4.9696/sqrt(mWW))) + ((njet==1) * (0.6711+7.0075/sqrt(mWW))) + ((njet>=2) * (0.3794+16.4944/sqrt(mWW))))'
        #weights['ggH_b']              = self._stdWgt+'*2.739*(((njet==0) * (1.1983-4.9696/sqrt(mWW))) + ((njet==1) * (0.6711+7.0075/sqrt(mWW))) + ((njet>=2) * (0.3794+16.4944/sqrt(mWW))))'

        # 1/x functional form
        weights['ggH_sbi']            = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))'
        weights['ggH_s']              = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))'
        #weights['ggH_s']              = self._stdWgt
        weights['ggH_b']              = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))'

        #weights['ggH_sbi']            = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))       * (mll>70)'
        #weights['ggH_s']              = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))       * (mll>70)'
        ##weights['ggH_s']              = self._stdWgt
        #weights['ggH_b']              = self._stdWgt+'*2.739*(((njet==0) * (1.0753-48.9352/mWW)) + ((njet==1) * (0.8457+67.5294/mWW)) + ((njet>=2) * (0.7777+161.9353/mWW)))       * (mll>70)'

        #integral = 19.7036
        #integral = 7.19278
        # --> 2.739 as global scale factor



        # scale factors from h126 to h125 GeV
        # ggH ~ 1.035
        # qqH ~ 1.041
        #                                             only-offshell            1 sm                     9 sm                    25 sm
        #weights['qqH_sbi']            = self._stdWgt+'*(mWW>130)*( 1.000*(dataset == 150) - 0.000*(dataset == 151) + 0.000*(dataset == 152))'
        #weights['qqH_s']              = self._stdWgt+'*(mWW>130)*( 0.125*(dataset == 150) - 0.250*(dataset == 151) + 0.125*(dataset == 152))'
        #weights['qqH_b']              = self._stdWgt+'*(mWW>130)*( 1.875*(dataset == 150) - 1.250*(dataset == 151) + 0.375*(dataset == 152))'
        weights['qqH_sbi']            = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.578/1.568*2.15/2.31)*1.12*(mWW>130)*( 1.000*( (dataset == 160) ||  (dataset == 169) || (dataset == 176) || (dataset == 172)) - 0.000*( (dataset == 161) ||  (dataset == 175) || (dataset == 173) || (dataset == 170)) + 0.000*( (dataset == 162) ||  (dataset == 171) || (dataset == 174) || (dataset == 177)))'
        weights['qqH_s']              = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.578/1.568*2.15/2.31)*1.12*(mWW>130)*( 0.125*( (dataset == 160) ||  (dataset == 169) || (dataset == 176) || (dataset == 172)) - 0.250*( (dataset == 161) ||  (dataset == 175) || (dataset == 173) || (dataset == 170)) + 0.125*( (dataset == 162) ||  (dataset == 171) || (dataset == 174) || (dataset == 177)))'
        weights['qqH_b']              = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.578/1.568*2.15/2.31)*1.12*(mWW>130)*( 1.875*( (dataset == 160) ||  (dataset == 169) || (dataset == 176) || (dataset == 172)) - 1.250*( (dataset == 161) ||  (dataset == 175) || (dataset == 173) || (dataset == 170)) + 0.375*( (dataset == 162) ||  (dataset == 171) || (dataset == 174) || (dataset == 177)))'

        # fix for 2011
        #print "sel = ",sel
        # due to lack of mww variable, add and remove powheg sample dataset=8126
        if (sel == "Hwidthmthmll7TeV") :
          weights['qqH_sbi']            = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.222/1.211*2.15/2.31)*1.12*( 1.000*( (dataset == 272) ||  (dataset == 275) || (dataset == 278) || (dataset == 281)) - 0.000*( (dataset == 273) ||  (dataset == 276) || (dataset == 279) || (dataset == 282)) + 0.000*( (dataset == 274) ||  (dataset == 277) || (dataset == 280) || (dataset == 283)) -1.000 * (dataset == 8126))'
          weights['qqH_s']              = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.222/1.211*2.15/2.31)*1.12*( 0.125*( (dataset == 272) ||  (dataset == 275) || (dataset == 278) || (dataset == 281)) - 0.250*( (dataset == 273) ||  (dataset == 175) || (dataset == 279) || (dataset == 282)) + 0.125*( (dataset == 274) ||  (dataset == 277) || (dataset == 280) || (dataset == 283)) )'
          weights['qqH_b']              = self._stdWgt+'*((njet==0)*1.412499+(njet==1)*1.035659+(njet>=2)*0.9020033)*(1.222/1.211*2.15/2.31)*1.12*( 1.875*( (dataset == 272) ||  (dataset == 275) || (dataset == 278) || (dataset == 281)) - 1.250*( (dataset == 273) ||  (dataset == 175) || (dataset == 279) || (dataset == 282)) + 0.375*( (dataset == 274) ||  (dataset == 277) || (dataset == 280) || (dataset == 283)) -1.000 * (dataset == 8126))'

          # fix since we don't have gen-level information at 7 TeV
          #weights['ggH_sbi']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1986-4.9341/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1986-4.9341/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8415+68.030/    (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8415+68.030/(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7655+165.6779/  (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7655+165.6779/(400))) )'

          #weights['ggH_s']              = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1986-4.9341/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1986-4.9341/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8415+68.030/    (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8415+68.030/(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7655+165.6779/  (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7655+165.6779/(400))) )'

          #weights['ggH_b']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1986-4.9341/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1986-4.9341/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8415+68.030/    (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8415+68.030/(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7655+165.6779/  (200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7655+165.6779/(400))) )'

          #weights['ggH_sbi']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1983-4.9696/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1983-4.9696/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.6711+7.0075/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.6711+7.0075/sqrt(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.3794+16.4944/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.3794+16.4944/sqrt(400))) )'

          #weights['ggH_s']              = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1983-4.9696/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1983-4.9696/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.6711+7.0075/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.6711+7.0075/sqrt(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.3794+16.4944/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.3794+16.4944/sqrt(400))) )'

          #weights['ggH_b']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.1983-4.9696/sqrt(200)) \
                                                                           #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.1983-4.9696/sqrt(400))) + \
                                                                #(njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.6711+7.0075/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.6711+7.0075/sqrt(400))) + \
                                                                #(njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.3794+16.4944/sqrt(200)) \
                                                                          #+ ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.3794+16.4944/sqrt(400))) )'


          weights['ggH_sbi']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.0753-48.9352/200) \
                                                                           + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.0753-48.9352/400)) + \
                                                                (njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8457+67.5294/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8457+67.5294/400)) + \
                                                                (njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7777+161.9353/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7777+161.9353/400)) )'

          weights['ggH_s']              = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.0753-48.9352/200) \
                                                                           + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.0753-48.9352/400)) + \
                                                                (njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8457+67.5294/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8457+67.5294/400)) + \
                                                                (njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7777+161.9353/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7777+161.9353/400)) )'

          weights['ggH_b']            = self._stdWgt+'*2.739*((njet==0) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(1.0753-48.9352/200) \
                                                                           + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(1.0753-48.9352/400)) + \
                                                                (njet==1) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.8457+67.5294/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.8457+67.5294/400)) + \
                                                                (njet==2) * (((dataset == 253) || (dataset == 263) || (dataset == 270))*(0.7777+161.9353/200) \
                                                                          + ((dataset == 256) || (dataset == 266) || (dataset == 271))*(0.7777+161.9353/400)) )'



          # test to suppress cross-feed
          #weights['qqH_sbi'] = weights['qqH_sbi'] + '*(mll>70)'
          #weights['qqH_s']   = weights['qqH_s']   + '*(mll>70)'
          #weights['qqH_b']   = weights['qqH_b']   + '*(mll>70)'

          #weights['ggH_sbi'] = weights['ggH_sbi'] + '*((mll>70) || (mth>130))'
          #weights['ggH_s']   = weights['ggH_s']   + '*((mll>70) || (mth>130))'
          #weights['ggH_b']   = weights['ggH_b']   + '*((mll>70) || (mth>130))'

          #weights['qqH']   = weights['qqH'] + '*(mll<70)'
          #weights['ggH']   = weights['ggH'] + '*(mll<70)'



          # scale 125 GeV -> 125.6 GeV
          #weights['ggH'] = weights['ggH'] + '*1.035'
          #weights['qqH'] = weights['qqH'] + '*1.041'


   #Double_t S =  0.125 * P1 -0.250 *P9 + 0.125 * P25;
   #Double_t I = -1.000 * P1 +1.500 *P9 - 0.500 * P25;
   #Double_t B =  1.875 * P1 -1.250 *P9 + 0.375 * P25;


        if cat in ['2j','2jtche05','2jtche05CJ','2jtche05FJ']:
            #weights['WW']                = self._stdWgt+'*(1+(mjj>500)*(detajj>3.5))'
            if ("Hwidth" not in sel) :
              weights['WW']                = self._stdWgt
            #weights['WWewk']             = self._stdWgt+'*(numbLHE==0)'
            weights['WWewk']             = self._stdWgt+'*(abs(jetLHEPartonpid1)!=6 && abs(jetLHEPartonpid2)!=6 && abs(jetLHEPartonpid3)!=6)*(abs(jetLHEPartonpid1)!=5 && abs(jetLHEPartonpid2)!=5 && abs(jetLHEPartonpid3)!=5)'

            if (sel == 'vbf' or sel == 'vbf-shape') :
              weights['CHITOP-Top']        = self._stdWgt+'*('+hwwinfo.massSelections(mass)['vbf-selection-top']+')'

            if (sel == 'vbf-shape-fish') :
              weights['CHITOP-Top']        = self._stdWgt+'*('+hwwinfo.massSelections(mass)['vbf-selection-fish-top']+')'

            if (sel == 'vbf2011' or sel == 'vbf2011-shape') :
              weights['CHITOP-Top']        = self._stdWgt+'*('+hwwinfo.massSelections(mass)['vbf2011-selection-top']+')'


            weights['TopPt0']             = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<100)'
            weights['TopPt1']             = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=100)'

            if (var != 1) : # only if it's shape and not cut based
              #weights['TopPt2']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<50)'
              #weights['TopPt3']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=50)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<70)'
              #weights['TopPt1']           =  self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=70)'

              #weights['TopPt1']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<50)'
              #weights['TopPt2']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=50)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<70)'
              #weights['TopPt3']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=70)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<110)'
              #weights['TopPt4']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=110)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<150)'
              #weights['TopPt5']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=150)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<190)'
              #weights['TopPt6']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=190)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<230)'
              #weights['TopPt7']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=230)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<270)'
              #weights['TopPt8']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=270)'

              weights['TopPt1']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<50)'
              weights['TopPt2']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=50)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<70)'
              weights['TopPt3']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=70)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<110)'
              weights['TopPt4']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=110)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<150)'
              weights['TopPt5']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=150)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<200)'
              weights['TopPt6']           = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=200)'

            #weights['TopPt2']             =  self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=70)*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))<100)'
            #weights['TopPt3']             = self._stdWgt+'*(((abs(jeteta1)<abs(jeteta2))*(jetpt1)+((abs(jeteta1)>=abs(jeteta2))*(jetpt2)))>=100)'

        if var in ['bdts','bdtl']:
            weights['WW']       = self._stdWgt+'*2*(event%2 == 0)'
            weights['ggH']      = self._stdWgt+'*2*kfW*(event%2 == 0)*'+self._muVal
            weights['qqH']      = self._stdWgt+'*2*kfW*(event%2 == 0)*'+self._muVal
            weights['wzttH']    = self._stdWgt+'*2*(event%2 == 0)*'+self._muVal
            # TODO Signal injection weights, if available
            weights['ggH-SI']   = self._stdWgt+'*2*kfW*(event%2 == 0)*'+self._muVal
            weights['qqH-SI']   = self._stdWgt+'*2*kfW*(event%2 == 0)*'+self._muVal
            weights['wzttH-SI'] = self._stdWgt+'*2*(event%2 == 0)*'+self._muVal


        if ("Hwidth" in sel) :
          print " Hww width analysis "
          #weights['WW']   = self._stdWgt+'*((njet==0) * (1.10)  + (njet==1) * (1.20) + (njet>=2) * (1.0))'
          #weights['WW']   = self._stdWgt+'*((njet==0) * (0.80)  + (njet==1) * (0.80) + (njet>=2) * (0.80))'
          if ("7TeV" in sel) :
            print " 7 TeV"
            #weights['WW']   = self._stdWgt+'*((njet==0) * (0.80)  + (njet==1) * (0.80) + (njet>=2) * (0.80))'
            weights['WW']   = self._stdWgt+'*((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WWlow']    = self._stdWgt+'* (mll<70)  * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WWhigh']   = self._stdWgt+'* (mll>=70) * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'

            weights['WW1']      = self._stdWgt+'*(mll>=12 && mll<30)       * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WW2']      = self._stdWgt+'*(mll>=30 && mll<45)       * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WW3']      = self._stdWgt+'*(mll>=45 && mll<60)       * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WW4']      = self._stdWgt+'*(mll>=60 && mll<70)       * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WW5']      = self._stdWgt+'*(mll>=70 && mll<100)      * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'
            weights['WW6']      = self._stdWgt+'*(mll>=100         )       * ((njet==0) * (1.08) + (njet==1) * (0.88) + (njet>=2) * (1.0))'

        return weights

    # _____________________________________________________________________________
    def _ensuredir(self,directory):
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
            except OSError as e:
                if e.errno == 17:
                    pass
                else:
                    raise e

    # _____________________________________________________________________________
    def _connectInputs(self, var, samples, dirmap, mask=None):
        inputs = {}
        treeName = 'latino'
        for process,filenames in samples.iteritems():
            if mask and process not in mask:
                continue
            tree = self._buildchain(treeName,[ (dirmap['base']+'/'+f) for f in filenames])
            if 'bdt' in var:
                bdttreeName = 'latinobdt'
                bdtdir = self._paths[var]
                bdttree = self._buildchain(bdttreeName,[ (dirmap[var]+'/'+f) for f in filenames])
                
                if tree.GetEntries() != bdttree.GetEntries():
                    raise RuntimeError('Mismatching number of entries: '
                                       +tree.GetName()+'('+str(tree.GetEntries())+'), '
                                       +bdttree.GetName()+'('+str(bdttree.GetEntries())+')')
                logging.debug('{0:<20} - master: {1:<20} friend {2:<20}'.format(process,tree.GetEntries(), bdttree.GetEntries()))
                tree.AddFriend(bdttree)

            inputs[process] = tree

        return inputs

    # _____________________________________________________________________________
    def _disconnectInputs(self,inputs):
        for n in inputs.keys():
            friends = inputs[n].GetListOfFriends()
            if friends.__nonzero__():
                for fe in friends:
                    friend = fe.GetTree()
                    inputs[n].RemoveFriend(friend)
                    ROOT.SetOwnership(friend,True)
                    del friend
            del inputs[n]
    
    # _____________________________________________________________________________
    def _buildchain(self,treeName,files):
        tree = ROOT.TChain(treeName)
        for path in files:
            self._logger.debug('     '+str(os.path.exists(path))+' '+path)
            if not os.path.exists(path):
                raise RuntimeError('File '+path+' doesn\'t exists')
            tree.Add(path) 

        return tree


    # _____________________________________________________________________________
    @staticmethod
    def _bins2hclass( bins ):
        '''
        Fixed bin width
        bins = (nx,xmin,xmax)
        bins = (nx,xmin,xmax, ny,ymin,ymax)
        Variable bin width
        bins = ([x0,...,xn])
        bins = ([x0,...,xn],[y0,...,ym])
        
        '''

        from array import array
        if not bins:
            return name,0
        elif not ( isinstance(bins, tuple) or isinstance(bins,list)):
            raise RuntimeError('bin must be an ntuple or an arryas')

        l = len(bins)
        # 1D variable binning
        if l == 1 and isinstance(bins[0],list):
            ndim=1
            hclass = ROOT.TH1D
            xbins = bins[0]
            hargs = (len(xbins)-1, array('d',xbins))
        elif l == 2 and  isinstance(bins[0],list) and  isinstance(bins[1],list):
            ndim=2
            hclass = ROOT.TH2D
            xbins = bins[0]
            ybins = bins[1]
            hargs = (len(xbins)-1, array('d',xbins),
                    len(ybins)-1, array('d',ybins))
        elif l == 3:
            # nx,xmin,xmax
            ndim=1
            hclass = ROOT.TH1D
            hargs = bins
        elif l == 6:
            # nx,xmin,xmax,ny,ymin,ymax
            ndim=2
            hclass = ROOT.TH2D
            hargs = bins
        else:
            # only 1d or 2 d hist
            raise RuntimeError('What a mess!!! bin malformed!')
        
        return hclass,hargs,ndim

    @staticmethod
    def _bins2dim(bins):
        hclass,hargs,ndim = ShapeFactory._bins2hclass( bins )
        return ndim

    @staticmethod
    def _makeshape( name, bins ):
        hclass,hargs,ndim = ShapeFactory._bins2hclass( bins )
        return hclass(name, name, *hargs)


    
    # _____________________________________________________________________________
    @staticmethod
    def _projexpr( bins = None ):
        if not bins:
            return name,0
        elif not ( isinstance(bins, tuple) or isinstance(bins,list)):
            raise RuntimeError('bin must be an ntuple or an arrya')
            
        l = len(bins)
        if l in [1,3]:
            # nx,xmin,xmax
            ndim=1
        elif l in [4,6]:
            # nx,xmin,xmax,ny,ymin,ymax
            ndim=2
        else:
            # only 1d or 2 d hist
            raise RuntimeError('What a mess!!! bin malformed!')

        hdef = '('+','.join([ str(x) for x in bins])+')' if bins else ''
        return hdef,ndim

if __name__ == '__main__':
    print '''
--------------------------------------------------------------------------------------------------
  .-')    ('-. .-.   ('-.      _ (`-.    ('-.  _   .-')      ('-.    .-. .-')     ('-.  _  .-')   
 ( OO ). ( OO )  /  ( OO ).-. ( (OO  ) _(  OO)( '.( OO )_   ( OO ).-.\  ( OO )  _(  OO)( \( -O )  
(_)---\_),--. ,--.  / . --. /_.`     \(,------.,--.   ,--.) / . --. /,--. ,--. (,------.,------.  
/    _ | |  | |  |  | \-.  \(__...--'' |  .---'|   `.'   |  | \-.  \ |  .'   /  |  .---'|   /`. ' 
\  :` `. |   .|  |.-'-'  |  ||  /  | | |  |    |         |.-'-'  |  ||      /,  |  |    |  /  | | 
 '..`''.)|       | \| |_.'  ||  |_.' |(|  '--. |  |'.'|  | \| |_.'  ||     ' _)(|  '--. |  |_.' | 
.-._)   \|  .-.  |  |  .-.  ||  .___.' |  .--' |  |   |  |  |  .-.  ||  .   \   |  .--' |  .  '.' 
\       /|  | |  |  |  | |  ||  |      |  `---.|  |   |  |  |  | |  ||  |\   \  |  `---.|  |\  \  
 `-----' `--' `--'  `--' `--'`--'      `------'`--'   `--'  `--' `--'`--' '--'  `------'`--' '--' 
--------------------------------------------------------------------------------------------------
'''    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--tag'            , dest='tag'            , help='Tag used for the shape file name'           , default=None)
    parser.add_option('--selection'      , dest='selection'      , help='Selection cut'                              , default=None)
    parser.add_option('--sigset'         , dest='sigset'         , help='Signal samples [SM]'                        , default='SM')
    parser.add_option('--dataset'        , dest='dataset'        , help='Dataset to process'                         , default=None)
    parser.add_option('--mcset'          , dest='mcset'          , help='Mcset to process'                           , default=None)
    parser.add_option('--path_latino'    , dest='path_latino'    , help='Root of the master trees'                   , default=None)
    parser.add_option('--path_bdt'       , dest='path_bdt'       , help='Root of the friendly bdt trees'             , default=None)
    parser.add_option('--path_shape_raw' , dest='path_shape_raw' , help='Destination directory of nominals'          , default='raw')
    parser.add_option('--range'          , dest='range'          , help='Range (optional default is var)'            , default=None)
    parser.add_option('--splitmode'      , dest='splitmode'      , help='Split in channels using a second selection' , default=None)

    parser.add_option('--keep2d',        dest='keep2d',     help='Keep 2d histograms (no unrolling)',     action='store_true',    default=False)
    parser.add_option('--no-noms',       dest='makeNoms',   help='Do not produce the nominal',            action='store_false',   default=True)
    parser.add_option('--no-syst',       dest='makeSyst',   help='Do not produce the systematics',        action='store_false',   default=True)
    parser.add_option('--do-syst',       dest='doSyst',     help='Do only one systematic',                default=None)
#     parser.add_option('--skip-syst',     dest='skipSyst',   help='Skip set of systematics',               default='')
    parser.add_option('--skip-syst',     dest='skipSyst',   help='Skip set of systematics',               default=[] , type='string' , action='callback' , callback=hwwtools.list_maker('skipSyst'))
    parser.add_option('--mu'       ,     dest='muVal',   help='Initial signal strengh',               default='1.' , type='string' )
    parser.add_option('--newcps',        dest='newcps',  help='On/Off New CPS weights',               default=False , action='store_true')   
    # EWK Doublet Model
    parser.add_option('--ewksinglet',    dest='ewksinglet',  help='On/Off EWK singlet model',           default=False , action='store_true')   
    parser.add_option('--cprimesq'  ,    dest='cprimesq',    help='EWK singlet C\'**2 mixing value',    default=[1.]  , type='string'  , action='callback' , callback=hwwtools.list_maker('cprimesq',',',float))
    parser.add_option('--brnew'     ,    dest='brnew'   ,    help='EWK singlet BRNew values',           default=[0.]  , type='string'  , action='callback' , callback=hwwtools.list_maker('brnew',',',float))
    parser.add_option('--approxewk'  ,    dest='approxewk',    help='EWK not scaling width/interf',     default=False , action='store_true') 
    parser.add_option('--YR3'        ,    dest='YR3rewght',    help='Rewhgt from YR2 to YR3' ,          default=True , action='store_true')  
    parser.add_option('--mHSM'      ,    dest='mHSM',        help='Mass of the SM Higgs boson@125',     default=125.6 , type='float')

    hwwtools.addOptions(parser)
    hwwtools.loadOptDefaults(parser)
    (opt, args) = parser.parse_args()

    print 'EWK Singlet:' , opt.ewksinglet
    print 'CPrime**2  :' , opt.cprimesq
    print 'BRNew      :' , opt.brnew
    print 'Approx. EWK :' , opt.approxewk 

    sys.argv.append( '-b' )
    ROOT.gROOT.SetBatch()



    if not opt.debug:
        pass
    elif opt.debug == 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig(level=logging.DEBUG)
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig(level=logging.INFO)

    try:
#    if True:
        checks = [
            ('energy'         , 'Energy not defined')            , 
            ('sigset'         , 'Signal not defined')            , 
            ('mcset'          , 'MonteCarlo not defined')        , 
            ('dataset'        , 'Dataset not defined')           , 
            ('selection'      , 'Selection not defined')         , 
            ('path_latino'    , 'Master tree path not defined')  , 
            ('path_shape_raw' , 'Where shall I put the shapes?') , 
        ]
        
        for dest,msg in checks:
            if not getattr(opt,dest):
                parser.print_help()
                parser.error(msg)

        if not opt.range:
            opt.range = opt.variable

        tag       = opt.tag if opt.tag else opt.variable
        variable  = opt.variable
        selection = opt.selection

        latinoDir           = opt.path_latino
        bdtDir              = opt.path_bdt
        nomOutDir           = os.path.join(opt.path_shape_raw,'nominals/{mass}/')
        systOutDir          = os.path.join(opt.path_shape_raw,'systematics/{mass}/')

        nomInputDir         = ''
        systInputDir        = '{syst}/'

        nModel = 1
        if opt.ewksinglet : nModel = len(opt.cprimesq)*len(opt.brnew)

        for iModel in xrange(0,nModel):
          iCP2 = iModel%len(opt.cprimesq) 
          iBRn = (int(iModel/len(opt.cprimesq)))

          if opt.ewksinglet:
            nominalOutFile      = 'shape_Mh{mass}_{category}_'+tag+'_shapePreSel_{flavor}.EWKSinglet_CP2_'+str(opt.cprimesq[iCP2]).replace('.','d')+'_BRnew_'+str(opt.brnew[iBRn]).replace('.','d')+'.root'
            systematicsOutFile  = 'shape_Mh{mass}_{category}_'+tag+'_shapePreSel_{flavor}.EWKSinglet_CP2_'+str(opt.cprimesq[iCP2]).replace('.','d')+'_BRnew_'+str(opt.brnew[iBRn]).replace('.','d')+'_{nick}.root'
          else: 
            nominalOutFile      = 'shape_Mh{mass}_{category}_'+tag+'_shapePreSel_{flavor}.root'
            systematicsOutFile  = 'shape_Mh{mass}_{category}_'+tag+'_shapePreSel_{flavor}_{nick}.root'
          
          factory = ShapeFactory()
          factory._outFileFmt  = nominalOutFile
  
  #         masses = hwwinfo.masses[:] if opt.mass == 0 else [opt.mass]
  #         factory._masses   = masses
          factory._masses = opt.mass
  
          # go through final channels
          factory._channels = dict([ (k,v) for k,v in hwwinfo.channels.iteritems() if k in opt.chans])
          print factory._channels
  
          factory._paths['base']  = latinoDir
          factory._paths['bdtl']  = bdtDir
          factory._paths['bdts']  = bdtDir
  
          factory._energy    = opt.energy
          factory._dataTag   = opt.dataset
          factory._sigTag    = opt.sigset
          factory._mcTag     = opt.mcset
          factory._range     = opt.range
          factory._splitmode = opt.splitmode
          factory._lumi      = opt.lumi
          factory._muVal     = opt.muVal 
 
          factory._newcps    = opt.newcps 
          factory._interfSyst= "NONE" 
          factory._ewksinglet= opt.ewksinglet
          factory._approxewk = opt.approxewk
          if not opt.ewksinglet :
            factory._cprimesq  = 1.
            factory._brnew     = 0. 
          else:  
            factory._cprimesq  = opt.cprimesq[iCP2]
            factory._brnew     = opt.brnew[iBRn]
          factory._YR3rewght   = opt.YR3rewght
          factory._mHSM        = opt.mHSM

          # load YR if needed
          if opt.YR3rewght : factory.loadYR()

          if opt.makeNoms:
              # nominal shapes
              print factory.makeNominals(variable,selection,nomInputDir,nomOutDir+nominalOutFile)
  
          if opt.makeSyst:
              class Systematics:
                  def __init__(self,name,nick,indir,mask):
                      pass
              # systematic shapes
              systematics = OrderedDict([
                  ('electronResolution'      , 'p_res_e'),
                  ('electronScale_down'      , 'p_scale_eDown'),
                  ('electronScale_up'        , 'p_scale_eUp'),
                  ('jetEnergyScale_down'     , 'p_scale_jDown'),
                  ('jetEnergyScale_up'       , 'p_scale_jUp'),
                  ('leptonEfficiency_down'   , 'eff_lDown'),
                  ('leptonEfficiency_up'     , 'eff_lUp'),
                  ('muonEfficiency_down'     , 'eff_mDown'),
                  ('muonEfficiency_up'       , 'eff_mUp'),
                  ('electronEfficiency_down' , 'eff_eDown'),
                  ('electronEfficiency_up'   , 'eff_eUp'),
                  ('puW_up'                  , 'puModelUp'),
                  ('puW_down'                , 'puModelDown'),
                  ('metResolution'           , 'met'),
                  ('metScale_down'           , 'p_scale_metDown'),
                  ('metScale_up'             , 'p_scale_metUp'),
                  ('muonScale_down'          , 'p_scale_mDown'),
                  ('muonScale_up'            , 'p_scale_mUp'),
                  ('chargeResolution'        , 'ch_res'), 
                  ('interferenceGGH_up'      , 'interf_ggHUp'),
                  ('interferenceGGH_down'    , 'interf_ggHDown'),
                  ('interferenceVBF_up'      , 'interf_qqHUp'),
                  ('interferenceVBF_down'    , 'interf_qqHDown'),
                  ('JER_down'                , 'p_res_jDown'),
                  ('JER_up'                  , 'p_res_jUp'),
                  ('JER_down'                , 'p_res_jDown'),
                  ('NNLL_up'                 , 'nnllUp'),
                  ('NNLL_down'               , 'nnllDown'),
                  ('NNLLR_up'                , 'nnllRUp'),
                  ('NNLLR_down'              , 'nnllRDown'),
                  ('WWmodel_up'              , 'WWmodelUp'),
                  ('WWmodel_down'            , 'WWmodelDown'),
              ])


              # remove skip-syst list
  #             if opt.skipSyst!='':
  #                for s in opt.skipSyst.split(','):
              for s in opt.skipSyst:
                print 'skipping systematics: ',s
                if s in systematics.keys() :
                  systematics.pop(s)
                else :
                  print '>> Beware! you are trying to remove a nuisance that is not even there: ',s

              systByWeight = {}
              # use only leptonEfficiency or muonEfficiency+electronEfficiency
              # skipSyst = ['leptonEfficiency_down', 'leptonEfficiency_up']
              #AdditionalSkipSyst = ['JER_down', 'JER_up', 'chargeResolution', 'puW_up', 'puW_down']  # -> temporary fix
              AdditionalSkipSyst = ['chargeResolution', 'puW_up', 'puW_down']  # -> temporary fix
              for s in AdditionalSkipSyst:
                print 'skipping systematics: ',s
                if s in systematics.keys() :
                  systematics.pop(s)
                else :
                  print '>> Beware! you are trying to remove a nuisance that is not even there: ',s


              systByWeight['leptonEfficiency_down'] = 'effWDown/effW'
              systByWeight['leptonEfficiency_up']   = 'effWUp/effW'
              systByWeight['muonEfficiency_down'] = 'effWMuDown/effW'
              systByWeight['muonEfficiency_up']   = 'effWMuUp/effW'
              systByWeight['electronEfficiency_down'] = 'effWElDown/effW'
              systByWeight['electronEfficiency_up']   = 'effWElUp/effW'
              systByWeight['interferenceGGH_up']   = '1.0'
              systByWeight['interferenceGGH_down'] = '1.0'
              systByWeight['interferenceVBF_up']   = '1.0'
              systByWeight['interferenceVBF_down'] = '1.0'

              systByWeight['puW_down'] = 'puWup/puW'
              systByWeight['puW_up']   = 'puWdown/puW'

              #if selection in ['CutWW'] :
              if selection in ['CutWW','Hwidth','HwidthWWcontrolRegion','HwidthWWcontrolRegionHighMllRemoved'] :
                systByWeight['NNLL_down']  = 'nllW_Qdown/nllW'
                systByWeight['NNLL_up']    = 'nllW_Qup/nllW'
                systByWeight['NNLLR_down'] = 'nllW_Rdown/nllW'
                systByWeight['NNLLR_up']   = 'nllW_Rup/nllW'
              else :
                systematics.pop('NNLL_down')
                systematics.pop('NNLL_up')
                systematics.pop('NNLLR_down')
                systematics.pop('NNLLR_up')

              if selection in ['Hwidth','HwidthWWcontrolRegion','HwidthWWcontrolRegionHighMllRemoved'] :
                systByWeight['WWmodel_up']   = '(1.27 - 0.14 * dphill)'
                systByWeight['WWmodel_down'] = '(1)'
              else :
                systematics.pop('WWmodel_up')
                systematics.pop('WWmodel_down')



              factory._systByWeight = systByWeight

              processMask = ['ggH', 'ggH_ALT',  'qqH',  'qqH_ALT', 'wzttH', 'ZH', 'WH', 'ttH', 'ggWW', 'Top', 'TopPt0', 'TopPt1', 'TopPt2', 'TopPt3', 'TopPt4', 'TopPt5', 'TopPt6', 'TopPt7', 'TopPt8', 'WW', 'WWlow', 'WWhigh', 'WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6', 'VV', 'VgS', 'Vg', 'DYTT', 'Other', 'VVV', 'WWewk', 'CHITOP-Top' , 'ggH_SM', 'qqH_SM', 'wzttH_SM' , 'WH_SM','ZH_SM','ttH_SM','ggH_sbi','ggH_b','ggH_s','qqH_sbi','qqH_b','qqH_s']

              if '2011' in opt.dataset:
                  processMask = ['ggH', 'ggH_ALT', 'qqH', 'qqH_ALT', 'VH' , 'wzttH', 'ZH', 'WH', 'ttH', 'ggWW', 'Top', 'WW', 'WWlow', 'WWhigh', 'WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6', 'VV', 'CHITOP-Top', 'ggH_SM', 'qqH_SM','VH_SM', 'wzttH_SM', 'ZH_SM', 'WH_SM', 'ttH_SM']

              systMasks = dict([(s,processMask[:]) for s in systematics])
              # interference is only on signal samples:
              systMasks['interferenceGGH_up'  ] = ['ggH', 'ggH_ALT']
              systMasks['interferenceGGH_down'] = ['ggH', 'ggH_ALT']
              systMasks['interferenceVBF_up'  ] = ['qqH', 'qqH_ALT']
              systMasks['interferenceVBF_down'] = ['qqH', 'qqH_ALT']

              # NNLL reweight and unceratinty only if WW
              #if selection in ['CutWW'] :
              if selection in ['CutWW','Hwidth','HwidthWWcontrolRegion','HwidthWWcontrolRegionHighMllRemoved'] and '2011' not in opt.dataset:
                systMasks['NNLL_up']    = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']
                systMasks['NNLL_down']  = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']
                systMasks['NNLLR_up']   = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']
                systMasks['NNLLR_down'] = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']

              if selection in ['Hwidth','HwidthWWcontrolRegion','HwidthWWcontrolRegionHighMllRemoved'] and '2011' not in opt.dataset:
                systMasks['WWmodel_up']   = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']
                systMasks['WWmodel_down'] = ['WW','WWlow','WWhigh','WW1', 'WW2', 'WW3', 'WW4', 'WW5', 'WW6']



              # remove selected nuisances for some samples
              processMaskNoDYTT = processMask
              processMaskNoDYTT = filter(lambda a: a != 'DYTT', processMaskNoDYTT)
              for syst,mask in systMasks.iteritems():
                if syst == 'JER_down' or syst == 'JER_up' :
                  #print ' >>> Old ', mask
                  mask = processMaskNoDYTT
                  #print ' >>> New ', mask
                  #print 'processMaskNoDYTT  = ',processMaskNoDYTT
                  systMasks[syst] = mask



              systDirs  = dict([(s,systInputDir if s not in systByWeight else 'templates/' ) for s in systematics])
              #systDirs  = dict([(s,systInputDir if s not in systByWeight else 'nominals/' ) for s in systematics])
              print "systDirs = ",systDirs

              for syst,mask in systMasks.iteritems():
                  if opt.doSyst and opt.doSyst != syst:
                      continue
                  print '='*80
                  print ' Processing ',syst,' for samples ',' '.join(mask)
                  print '='*80
                  files = factory.makeSystematics(variable,selection,syst,mask,systDirs[syst],systOutDir+systematicsOutFile, nicks=systematics)

#          factory.Delete() 

    except Exception as e:
        print '*'*80
        print 'Fatal exception '+type(e).__name__+': '+str(e)
        print '*'*80
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback, file=sys.stdout)
#         traceback.print_tb(exc_traceback, limit=3, file=sys.stdout)
        print '*'*80
    finally:
        print 'Used options'
        print ', '.join([ '{0} = {1}'.format(a,b) for a,b in opt.__dict__.iteritems()])
