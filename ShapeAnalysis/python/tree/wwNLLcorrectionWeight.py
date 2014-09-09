#
#
#
#   \ \        / \ \        /       \  |  |      |
#    \ \  \   /   \ \  \   /         \ |  |      |
#     \ \  \ /     \ \  \ /        |\  |  |      |
#      \_/\_/       \_/\_/        _| \_| _____| _____|
#
#
#
#


from tree.gardening import TreeCloner
import numpy
import ROOT
import sys
import optparse
import re
import warnings
import os.path
from array import array;

class wwNLLcorrectionWeightFiller(TreeCloner):
    def __init__(self):
       pass

    def help(self):
        return '''Add weight for WW NLL reweighting'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)
        group.add_option('-d', '--data'        , dest='datafile', help='Name of the input root file with reweight histograms',)
        group.add_option('-m', '--mcsample'    , dest='mcsample', help='Name of the mc sample to be considered. Possible options [powheg, mcatnlo, madgraph]',)
        parser.add_option_group(group)

        return group


    def checkOptions(self,opts):
        if (
             not hasattr(opts,'datafile') and
             not hasattr(opts,'mcsample')      ) :
            raise RuntimeError('Missing parameter')

        self.datafile = opts.datafile
        self.mcsample = opts.mcsample

    def process(self,**kwargs):

        # change this part into correct path structure... 
        cmssw_base = os.getenv('CMSSW_BASE')
        try:
            ROOT.gROOT.LoadMacro(cmssw_base+'/src/HWWAnalysis/ShapeAnalysis/python/tree/wwNLLcorrectionWeight.C+g')
        except RuntimeError:
            ROOT.gROOT.LoadMacro(cmssw_base+'/src/HWWAnalysis/ShapeAnalysis/python/tree/wwNLLcorrectionWeight.C++g')
        #----------------------------------------------------------------------------------------------------

        print " file = ",self.datafile
        wwNLL = ROOT.wwNLL(self.datafile,self.mcsample)



        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        # does that work so easily and give new variable itree and otree?
        self.connect(tree,input)
        newbranches = ['nllW', 'nllW_Rup', 'nllW_Qup', 'nllW_Rdown', 'nllW_Qdown']
        self.clone(output,newbranches)

        nllW    = numpy.ones(1, dtype=numpy.float32)
        nllW_Rup    = numpy.ones(1, dtype=numpy.float32)
        nllW_Qup    = numpy.ones(1, dtype=numpy.float32)
        nllW_Rdown    = numpy.ones(1, dtype=numpy.float32)
        nllW_Qdown    = numpy.ones(1, dtype=numpy.float32)

        self.otree.Branch('nllW'  , nllW  , 'nllW/F')
        self.otree.Branch('nllW_Rup'  , nllW_Rup  , 'nllW_Rup/F')
        self.otree.Branch('nllW_Qup'  , nllW_Qup  , 'nllW_Qup/F')
        self.otree.Branch('nllW_Rdown'  , nllW_Rdown  , 'nllW_Rdown/F')
        self.otree.Branch('nllW_Qdown'  , nllW_Qdown  , 'nllW_Qdown/F')

        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 

        #what is self.itree? what is self.otree?
        itree     = self.itree
        otree     = self.otree

        print '- Starting eventloop'
        step = 5000

        for i in xrange(nentries):
        #for i in xrange(100):

            itree.GetEntry(i)

            if i > 0 and i%step == 0.:
                print i,'events processed.'

            #ptl1 = itree.genVV_lepton1_LHE_pt
            #ptl2 = itree.genVV_lepton2_LHE_pt
            #phil1 = itree.genVV_lepton1_LHE_phi
            #phil2 = itree.genVV_lepton2_LHE_phi

            #ptv1 = itree.genVV_neutrino1_LHE_pt
            #ptv2 = itree.genVV_neutrino2_LHE_pt
            #phiv1 = itree.genVV_neutrino1_LHE_phi
            #phiv2 = itree.genVV_neutrino2_LHE_phi

            # after ISR but before QED FSR
            # because at calculation step these are not defined
            # and they don't know about photons
            # use lvlv status=3 (?? why not 1??) particles

            ptl1 = itree.genVV_lepton1_pt
            ptl2 = itree.genVV_lepton2_pt
            phil1 = itree.genVV_lepton1_phi
            phil2 = itree.genVV_lepton2_phi

            ptv1 = itree.genVV_neutrino1_pt
            ptv2 = itree.genVV_neutrino2_pt
            phiv1 = itree.genVV_neutrino1_phi
            phiv2 = itree.genVV_neutrino2_phi



            wwNLL.SetPTWW(ptl1, phil1, ptl2, phil2, ptv1, phiv1, ptv2, phiv2)

            nllW[0]   = wwNLL.nllWeight(0)
            nllW_Rup[0]   = wwNLL.nllWeight(1,1)
            nllW_Qup[0]   = wwNLL.nllWeight(1,0)
            nllW_Rdown[0] = wwNLL.nllWeight(-1,1)
            nllW_Qdown[0] = wwNLL.nllWeight(-1,0)

            otree.Fill()

        self.disconnect()
        print '- Eventloop completed'