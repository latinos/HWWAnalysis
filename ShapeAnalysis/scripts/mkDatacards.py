#!/usr/bin/env python

import optparse
import ROOT
import sys
import re
import os.path
from HWWAnalysis.Misc.odict import OrderedDict
import hwwinfo
import datadriven
import fnmatch
from WWAnalysis.AnalysisStep.systematicUncertainties import getCommonSysts,addFakeBackgroundSysts

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class ShapeDatacardWriter:
    '''Dump a crappy datacard to file'''
    
    def __init__(self, mass, bin):
        self._mass = mass
        self._bin = bin

    def __del__(self):
        pass

    def write(self, yields, nuisances, path, fileFmt, signals = ['ggH', 'vbfH', 'wzttH']):

        cardPath = path.format(mass = self._mass, bin = self._bin)
        print cardPath
        card = open( cardPath ,"w")
        card.write('## Shape input card for H->WW analysis using 2.12/fb\n')
        
        card.write('imax 1 number of channels\n')
        card.write('jmax * number of background\n')
        card.write('kmax * number of nuisance parameters\n') 

        card.write('-'*100+'\n')
        card.write('bin         %s' % self._bin+'\n')
        if 'Data' not in yields:
            print yields.keys()
            raise RuntimeError('No Data found!')
        card.write('observation %.0f\n' % yields['Data']._N)
        card.write('shapes  *           * '+
                   fileFmt.format(mass=self._mass, bin=self._bin)+
                   '     histo_$PROCESS histo_$PROCESS_$SYSTEMATIC'+'\n')
        card.write('shapes  data_obs    * '+
                   fileFmt.format(mass=self._mass, bin=self._bin)+
                   '     histo_Data'+'\n')
        card.write('-'*100+'\n')
        
        bkgs = [ name for name in yields if name not in signals and name != 'Data']
        sigs = [ name for name in yields if name in signals ]
        keyline = []
        keyline.extend([ (-i,s,yields[s]._N) for i,s in enumerate(sigs) ])
        keyline.extend([ (i+1,b,yields[b]._N) for i,b in enumerate(bkgs) ])


        card.write('bin'.ljust(48)+''.join([self._bin.ljust(8)*len(keyline)])+'\n')
        card.write('process'.ljust(48)+''.join([n.ljust(8) for (i,n,N) in keyline])+'\n' )
        card.write('process'.ljust(48)+''.join([('%d' % i).ljust(8) for (i,n,N) in keyline])+'\n' )
        card.write('rate'.ljust(48)+''.join([('%-.2f' % N).ljust(8) for (i,n,N) in keyline])+'\n' )
        card.write('-'*100+'\n')

        for name in nuisances:
            (pdf,effect) = nuisances[name]
            if len(pdf) == 1: card.write('{0:<31} {1:<7}         '.format(name,pdf[0]))
            else:             card.write('{0:<31} {1:<7} {2:<7} '.format(name,pdf[0],pdf[1]))
            for i,p,y in keyline:
                if p in effect: 
                    if pdf[0]=='gmN':   card.write('%6.4f  ' % effect[p])
                    elif (pdf[0]=='shape' or pdf[0]=='shapeN2'): card.write('%-7d ' % effect[p])
                    else:               card.write('%5.3f   ' % effect[p] )
                else: card.write('-'.ljust(8))
            card.write('\n')

        card.close()

class Yield:
    def __init__(self,*args,**kwargs):
        if not args:
            raise RuntimeError('Specify number of entries')
        self._N = args[0]
        if 'name' in kwargs:
            self._name = kwargs['name']
        if 'entries' in kwargs:
            self._entries = kwargs['entries']


class ShapeLoader:
    '''Load the histogram data from the shape file
    + Yields
    + Nuisance shapes and parameters'''

    def __init__(self, path):
        self._systRegex = re.compile('^histo_([^_]+)_(.+)(Up|Down)$')
        self._nomRegex  = re.compile('^histo_([^_]+)$')
        self._src = ROOT.TFile.Open(path)
        self._yields = OrderedDict()

    def __del__(self):
        del self._src

    def yields(self):
        return self._yields.copy()

    def effects(self):
        return self._effects.copy()

    def load(self):
        # load the histograms and calculate the yields
        names = [ k.GetName() for k in self._src.GetListOfKeys()]
        self._nominals = sorted([ name for name in names if self._nomRegex.match(name) ]) 
#         print self._nominals
        self._systematics = sorted([ name for name in names if self._systRegex.match(name) ])
#         print '\n'.join(self._nominals)
#         print '\n'.join(self._systematics)
        for name in self._nominals:
            process = self._nomRegex.match(name).group(1)
            h = self._src.Get(name)
            N =  h.Integral(0,h.GetNbinsX())
            entries = h.GetEntries()

            # TODO: DYTT cleanup
#             if entries < 5: continue

            self._yields[process] = Yield( N, name=process, entries=entries ) 
#             self._yields[process] = Yield( N, name=process, entries=entries, shape=h ) 
#             print process, '%.3f' % h.Integral(0,h.GetNbins())
        
#         print self._systematics
        ups = {}
        downs = {}
        for name in self._systematics:
            # check for Up/Down
            (process,effect,var) = self._systRegex.match(name).group(1,2,3)
            if var == 'Up': 
                if effect not in ups: ups[effect]= []
                ups[effect].append(process)
            elif var == 'Down':
                if effect not in downs: downs[effect]= []
                downs[effect].append(process)

#         del ups['p_scale_j'][0]
#         del ups['p_scale_e'][1]
        # check 
        for effect in ups:
            if set(ups[effect]) != set(downs[effect]):
                sUp = set(ups[effect])
                sDown = set(downs[effect])
                raise RuntimeError('Some systematics shapes for '+effect+' not found in up and down variation: \n '+', '.join( (sUp | sDown) - ( sUp & sDown ) ))
        
        # all checks out, save only one
        self._effects = ups

class NuisanceMapBuilder:
    def __init__(self, ddPath):
        self._common    = OrderedDict()
        self._0jetOnly  = OrderedDict()
        self._1jetOnly  = OrderedDict()
        self._ddEstimates = OrderedDict()
        # to options
        self.ddPath     = ddPath
        self.statShapeVeto = []

        self.ddreader = datadriven.DDCardReader(ddPath)

        self._build()
 
    def _build(self):
        # common 0/1 jet systematics
        pureMC = [ 'Vg', 'VV', 'DYTT', 'ggH', 'vbfH', 'wzttH'] 
        dummy = {}
        dummy['CMS_fake_e']    = (1.50, ['WJet']) # take the average of ee/me 
#         dummy['CMS_fake_m']    = (1.42, ['WJet']) # take the average of mm/em
        dummy['CMS_eff_l']     = (1.04, pureMC)
        dummy['CMS_p_scale_e'] = (1.02, pureMC)
#         dummy['CMS_p_scale_m'] = (1.01, pureMC)
        dummy['CMS_p_scale_j'] = (1.02, pureMC)
        dummy['CMS_met']       = (1.02, pureMC)
        dummy['lumi']          = (1.04, pureMC)

        for k,v in dummy.iteritems():
            self._common[k] = (['lnN'], dict([( process, v[0]) for process in v[1] ]) )

        self._common['pdf_gg']    = (['lnN'],dict([('ggWW',1.0),('ggH',1.08)]) )
        self._common['pdf_qqbar'] = (['lnN'],dict([('WW',1.0),('VV',1.04),('vbfH',1.02)]) )
        self._common['pdf_assoc'] = (['lnN'],dict([('WW',1.04)]) )

        dummy = {} 
        # both 0/1 jets but different
        dummy['CMS_QCDscale_WW_EXTRAP'] = ([0.95, 1.21], ['WW'])
        dummy['QCDscale_VV']            = ([1.03, 1.03], ['VV'])
        dummy['QCDscale_ggH1in']        = ([0.89, 1.39], ['ggH'])
        dummy['QCDscale_ggH_ACEPT']     = ([1.02, 1.02], ['ggH'])
        dummy['QCDscale_ggVV']          = ([1.5,  1.5],  ['ggWW'])
        dummy['QCDscale_qqH']           = ([1.01, 1.01], ['vbfH'])
        dummy['QCDscale_qqH_ACEPT']     = ([1.02, 1.02], ['vbfH'])
        dummy['QCDscale_wzttH_ACEPT']   = ([1.02, 1.02], ['wzttH'])
        dummy['QCDscale_wzttH']         = ([1.01, 1.01], ['wzttH'])
        dummy['UEPS']                   = ([0.94, 1.11], ['ggH'])

        for k,v in dummy.iteritems():
            self._0jetOnly[k] = (['lnN'], dict([( process, v[0][0]) for process in v[1] ]) )
            self._1jetOnly[k] = (['lnN'], dict([( process, v[0][1]) for process in v[1] ]) )

        # 0 jets only
        dummy = {}
        
        dummy['CMS_fake_Vg']  = (2.00,['Vg']) # Vg, 0jet 
        dummy['QCDscale_Vg']  = (1.50,['Vg']) 
        dummy['QCDscale_ggH'] = (1.16,['ggH']) # 0 jets only
        for k,v in dummy.iteritems():
            self._0jetOnly[k] = (['lnN'], dict([( process, v[0]) for process in v[1] ]) )

        # 1 jet only
        dummy = {}
        dummy['QCDscale_ggH2in'] = (0.95,['ggH']) # 1 jey only
        for k,v in dummy.iteritems():
            self._1jetOnly[k] = (['lnN'], dict([( process, v[0]) for process in v[1] ]) )

    
    def _addDataDrivenNuisances(self, nuisances, yields, mass, jets, flavor):
        
        (estimates,dummy) = self.ddreader.get(mass, jets, flavor)

        pdf = 'lnN'

        mapping = {}
        mapping['WW']   = ['WW','ggWW']
        mapping['Top']  = ['Top']
        mapping['DYLL'] = ['DYLL']

        eff_bin1_tmpl = 'CMS_hww_{0}_{1}_{2}j_stat_bin1'
        for tag,processes in mapping.iteritems(): 
            extr_entries = {}
            stat_entries = {}
            eff_extr = 'CMS_hww_{0}_{1}j_extr'.format(tag,jets)
            eff_stat = 'CMS_hww_{0}_{1}j_stat'.format(tag,jets)
            
            available = [ p for p in processes if p in estimates ]
            if not available: continue

            # check the dd to have the same events in the ctr region before associating them
            listNctr = [ estimates[p].Nctr for p in available ]
            
            if len(available) != listNctr.count(estimates[available[0]].Nctr):
                raise RuntimeError('Mismatch between Nctr in the same systematic: '+', '.join([ '{0}{1}'.format(n[0],n[1]) for n in zip(available, listNctr) ]) )
            
            for process in available:
                e = estimates[process]
                extrUnc = 1+e.delta/e.alpha if pdf != 'gmM' else e.delta/e.alpha

                extr_entries[process] = extrUnc
                stat_entries[process] = e.alpha
                eff_bin1 = eff_bin1_tmpl.format(process,flavor,jets)
                if eff_bin1 in nuisances:
                    del nuisances[eff_bin1]


            nuisances[eff_extr] = ([pdf], extr_entries )
            nuisances[eff_stat] = (['gmN',e.Nctr], stat_entries)


    def _addStatisticalNuisances(self,nuisances, yields,jets,channel):
        for p,y in yields.iteritems():
            if p == 'Data':
                continue
            name  = 'CMS_hww_{0}_{1}_{2}j_stat_bin1'.format(p,channel,jets)
            if y._entries == 0.:
                continue
            value = 1+(1./ROOT.TMath.Sqrt(y._entries) if y._entries > 0 else 0.)
            nuisances[name] = ( ['lnN'], dict({p:value}) )

    def _addExperimentalMassDepNuisances(self, nuisances, mass, jets, flavor ):
        '''Mass dependent experimental nuisances [normalization only]'''
        channels = dict([('sf',['ee','mm']),('of',['em','me'])])

        cards = [ self._expcards[mass][jets][ch] for ch in channels[flavor]]
        # further x-check on the items we are going to average
        diffEffs = set(cards[0].keys()) ^ set(cards[1].keys())
        if diffEffs:
            raise RuntimeError('Different systematics found for '+'-'.join(channels[flavor])+': '+' '.join(diffEffs))

        effects = cards[0].keys()
        for eff in effects:
            # aslo check the same processes re available for the same eff
            diffProcs = set(cards[0][eff].keys()) ^ set(cards[1][eff].keys())
            if diffProcs:
                raise RuntimeError('Different systematics found for '+'-'.join(channels[flavor])+', '+eff+': '+' '.join(diffProcs))

            processes = cards[0][eff].keys()
    
            tag = 'CMS_'+eff
            if tag in nuisances: del nuisances[tag]
            values = [(cards[0][eff][p]+cards[1][eff][p])/2. for p in processes]
            nuisances[tag] = (['lnN'], dict(zip(processes,values)) )

    def _addWWShapeNuisances(self,nuisances, effects):
        '''Generator relates shape nuisances'''
        wwRegex  = re.compile('Gen_(.+)$')
        for eff,processes in effects.iteritems():
            # select the experimental effects only (starting with gen)
            if not wwRegex.match(eff):
                continue
            tag = eff
            if tag in nuisances: del nuisances[tag]
            nuisances[tag] = (['shapeN2'],dict([ (p,1) for p in processes]) )


    def _addExperimentalShapeNuisances(self, nuisances, effects):
        '''Experimental Shape-based nuisances'''
        # expr for CMS nuisances
        expRegex  = re.compile('CMS_(.+)')
        # expr for statistical nuisances
        statRegex = re.compile('CMS_hww_([^_]+)_.+_stat_shape')

#         mask = ['Vg','DYLL','DYTT']#'WJet',]#'Vg','DYLL','DYTT',]#]#,'Top',]#'WW','ggWW','VV',]#'ggH','vbfH']
        for eff in sorted(effects):
            # select the experimental effects only (starting with CMS)
            if not expRegex.match(eff):
                continue

            m = statRegex.match(eff)
            if m and m.group(1) in self.statShapeVeto:
                print 'Skipping',eff,' (vetoed, data driven)'
                continue
            tag = eff
            if tag in nuisances: del nuisances[tag]
            nuisances[tag] = (['shapeN2'],dict([ (p,1) for p in effects[eff] ]) )

    def _addShapeNuisances(self, nuisances, effects, opts):
        # local copy
        shapeNu = OrderedDict()

        self._addWWShapeNuisances(shapeNu, effects)
        self._addExperimentalShapeNuisances(shapeNu, effects)

        if 'shapeFlags' not in opts:
            sys.exit(-1)
        flags = opts['shapeFlags']

        nus = set(shapeNu.keys())
        dummy = nus.copy()
        for exp,flag in flags:
            subset = set(fnmatch.filter(nus,exp))
            if flag:
                dummy |= subset
            else:
                dummy -= subset
#         print dummy

        for eff in shapeNu:
            if eff not in dummy: continue
            if eff in nuisances: del nuisances[eff]
            nuisances[eff] = shapeNu[eff]

    #  _  _      _                          
    # | \| |_  _(_)___ __ _ _ _  __ ___ ___
    # | .` | || | (_-</ _` | ' \/ _/ -_|_-<
    # |_|\_|\_,_|_/__/\__,_|_||_\__\___/__/
    #                                      
    def nuisances(self, yields, effects, mass, jets, channel, opts):
        '''Add the nuisances according to the options'''
        allNus = OrderedDict()
#         allNus.update(self._common)
#         if jets == 0:
#             allNus.update(self._0jetOnly)
#         elif jets == 1:
#             allNus.update(self._1jetOnly)

        optMatt = mattOpts()
        optMatt.WJadd = 0.36
        optMatt.WJsub = 0.0

        qqWWfromData = int(mass) < 200

        CutBased = getCommonSysts(int(mass),channel,jets,qqWWfromData, optMatt)
        common = OrderedDict()
        for k in sorted(CutBased):
            common[k] = CutBased[k]
        allNus.update( common )
#         print channel
#         addFakeBackgroundSysts(allNus, mass,channel,jets)#,errWW=0.2,errggWW=0.2,errDY=1.0,errTop0j=1.0,errTop1j=0.3,errWJ=0.5)

        self._addStatisticalNuisances(allNus, yields, jets, channel)
        self._addDataDrivenNuisances(allNus, yields, mass, jets, channel)
        
        self._addShapeNuisances(allNus,effects, opts)

        if 'nuisFlags' not in opts:
            raise RuntimeError('nuisFlags not found among the allNus options')

        flags = opts['nuisFlags']

        finalNuisances = OrderedDict()
        nus = set(allNus.keys())
        dummy = nus.copy()
        for exp,flag in flags:
            subset = set(fnmatch.filter(nus,exp))
            if flag:
                dummy |= subset
            else:
                dummy -= subset

        nuisances = OrderedDict()
        for eff in allNus:
            if eff not in dummy: continue
            nuisances[eff] = allNus[eff]

        return nuisances

def incexc(option, opt_str, value, parser):
    print option
    print opt_str
    print value
    if not hasattr(parser.values,'shapeFlags'):
        setattr(parser.values,'shapeFlags',[])

    optarray = str(option).split('/')
    print optarray
    if '--Ish' in optarray:
        parser.values.shapeFlags.append((value,True))
    elif '--Xsh' in optarray:
        parser.values.shapeFlags.append((value,False))
    elif '-I' in optarray:
        parser.values.nuisFlags.append((value,True))
    elif '-X' in optarray:
        parser.values.nuisFlags.append((value,False))
    
class mattOpts: pass

if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--prefix','-p',dest='prefix',help='prefix',default=None)
    parser.set_defaults(shapeFlags=[])
    parser.set_defaults(nuisFlags=[])
    parser.add_option('--Xsh','--excludeShape',action='callback', dest='shapeFlags', type='string', callback=incexc)#)action='append',default=None)
    parser.add_option('--Ish','--includeShape',action='callback', dest='shapeFlags', type='string', callback=incexc)#)action='append',default=None)
    parser.add_option('-X','--exclude',action='callback', dest='nuisFlags', type='string', callback=incexc)#)action='append',default=None)
    parser.add_option('-I','--include',action='callback', dest='nuisFlags', type='string', callback=incexc)#)action='append',default=None)
#     parser.add_option('--ddpath', dest='ddpath', help='Data driven path', default=None)
    parser.add_option('--path_dd'           , dest='path_dd'           , help='Data driven path'                 , default=None)
    parser.add_option('--path_shape_merged' , dest='path_shape_merged' , help='Destination directory for merged' , default=None)
    hwwinfo.addOptions(parser)
    hwwinfo.loadOptDefaults(parser)

    (opt, args) = parser.parse_args()
    print 'ShapeFlags: ',opt.shapeFlags
    print 'NuisFlags:  ',opt.nuisFlags
#     sys.exit(0)

    # checks
    if not opt.var or not opt.lumi:
        parser.error('The variable and the luminosty must be defined')
#     lumi = float(opt.lumi)
    var = opt.var

    sys.argv.append('-b')
    ROOT.gROOT.SetBatch()

    flavors = ['sf','of']
    jetBins = hwwinfo.jets[:]
    masses = hwwinfo.masses[:] if opt.mass == 0 else [opt.mass]

#     histPath = 'merged/'
#     ddPath   = '/shome/thea/HWW/ShapeAnalysis/data/AnalFull2011_BDT/'
    mergedPath = opt.path_shape_merged
    outPath    = 'datacards/'
    if opt.prefix:
        if opt.prefix[0] == '/':
            raise NameError('prefix: Only subdirectories are supported')
        outPath = (opt.prefix if opt.prefix[-1] == '/' else opt.prefix+'/')+outPath
    shapeSubDir = 'shapes/'

    shapeDir = outPath+shapeSubDir[:-1]

    os.system('mkdir -p '+outPath)
    if os.path.exists(shapeDir):
        os.unlink(shapeDir)
    os.symlink(os.path.abspath(mergedPath), shapeDir)

    optsNuis = {}
    optsNuis['shapeFlags'] = opt.shapeFlags
    optsNuis['nuisFlags'] = opt.nuisFlags
    lumistr = '{0:.2f}'.format(opt.lumi)
    shapeTmpl = os.path.join(mergedPath,'hww-'+lumistr+'fb.mH{mass}.{flavor}_{jets}j_shape.root')
    mask = ['Vg','DYLL','DYTT']

    builder = NuisanceMapBuilder( opt.path_dd )
    builder.statShapeVeto = mask
    for mass in masses:
        for jets in jetBins:
            for flavor in flavors:
                print '- Processing',mass, jets, flavor
                loader = ShapeLoader(shapeTmpl.format(mass = mass, jets=jets, flavor=flavor) ) 
                loader.load()
    
                writer = ShapeDatacardWriter( mass,'{0}_{1}j'.format(flavor, jets) )
                print '   + loading yields'
                yields = loader.yields()

                # reshuffle the order
                order = [ 'vbfH', 'ggH', 'wzttH', 'ggWW', 'Vg', 'WJet', 'Top', 'WW', 'DYLL', 'VV', 'DYTT', 'Data']
                oldYields = yields.copy()
                yields = OrderedDict([ (k,oldYields[k]) for k in order if k in oldYields])
                
                effects = loader.effects()

                print '   + making nuisance map'
                nuisances = builder.nuisances( yields, effects , mass, jets, flavor, optsNuis)
    
                basename = 'hww-'+lumistr+'fb.mH{mass}.{bin}_shape'
                print '   + dumping all to file'
                writer.write(yields,nuisances,outPath+basename+'.txt',shapeSubDir+basename+'.root')
