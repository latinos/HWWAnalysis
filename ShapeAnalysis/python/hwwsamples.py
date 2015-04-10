import hwwtools
import re


signals = ['ggH','qqH','wzttH','ggH_ALT','qqH_ALT','WH','ZH','ttH','VH']

#--------------
# mcsets,
# list of samples and compact dictionary
#
# filter the list of samples
# and create association     label          -> label                  -> vector of root files
#                       used by mkShapes      just for association       blabla.root
#

mcsets = {
    '0j1j-JHU' : [
        #signals
        'ggH','ggH_ALT','jhu_NORM','jhu_NLO',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','VVV','DYTT','DYLL','WWnlo','WWnloUp','WWnloDown',
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template', 'Vg-template',
        # 0j1j specific
        ('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        ('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],
    '0j1j-JHUOthers' : [
        #signals
        'ggH','ggH_ALT','jhu_NORM','jhu_NLO','qqH','qqH_ALT','wzttH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','VVV','DYTT','DYLL','WWnlo','WWnloUp','WWnloDown',
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        ('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        ('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],
    '0j1j' : [
        #signals
        'ggH','qqH','WH','ZH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','VVV','DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],
    '0j1j-mH125' : [
        #signals
        'ggH','qqH','WH','ZH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','VVV','DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j') ,         #    mkmerged vuole "-template"
        # mH125 as background
        'ggH_SM', 'qqH_SM', 'WH_SM','ZH_SM',
    ],
     '0j1j-ss' : [
        #signals
        'ggH','qqH','WH','ZH',
        # bkgs
        'Other','VgS','Vg','WJet',
        #'WW',
        #'ggWW','VgS','Vg','WJet',
        #'Top',
        #'VV','DYTT','DYLL',
        #'WWnlo','WWnloUp','WWnloDown','TopTW','TopCtrl',
        #,'WJetSS',
        # systematics
        #'WJetFakeRate-eUp', 'WJetFakeRate-eDn','WJetFakeRate-mUp', 'WJetFakeRate-mDn'
        # templates
        #'VgS-template','Vg-template','Top-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],
    '0j1jWW' : [
        #signals
        'ggH','qqH','WH','ZH',
        # bkgs
        #('WW', 'WWmad'),
        ('WW', 'WWpowheg'),
        #('WW', 'WWmc'),
        #'ggWW',
        ('ggWW','ggWWnew'),
        'VgS','Vg','WJet','Top','VV',
        #'VVV',
        ('VVV','VVVall'),
        'DYTT','WWnlo','WWnloUp','WWnloDown',
        #
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],
    'cutbased' : [
        #signals
        'ggH','qqH','WH','ZH',
        # bkgs
        'WW', 'VgS','Vg','WJet','Top','VV','VVV','DYTT',
        # templates
        'VgS-template','Vg-template',
        # replce ggWW and DYLL
        ('ggWW',    'WW'),
        ('DYLL',    'WW'),
    ],
    'vbf_sf' : [
        #signals
        'ggH','qqH', #'wzttH',
        # bkgs
        #'WW','ggWW','VgS','Vg','WJet','Top','VV','DYTT','DYLL',
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'Top','VV',
        'WWewk',
        # dummy
        ('DYee', 'qqH'),
        ('DYmm', 'qqH'),
        # nuisance
        'WWpow',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # templates
        'VgS-template','Vg-template',
        # templates for Top estimation
        ('CHITOP-Top',     'Top'),
    ],
    'vbf_sf-mH125' : [
        #signals
        'ggH','qqH', #'wzttH',
        # bkgs
        #'WW','ggWW','VgS','Vg','WJet','Top','VV','DYTT','DYLL',
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'Top','VV',
        'WWewk',
        # dummy
        ('DYee', 'qqH'),
        ('DYmm', 'qqH'),
        # nuisance
        'WWpow',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # templates
        'VgS-template','Vg-template',
        # templates for Top estimation
        ('CHITOP-Top',     'Top'),
        # mH125 as background
        'ggH_SM', 'qqH_SM', # 'WH_SM','ZH_SM',
    ],
    'vbf_of' : [
        #signals
        'ggH','qqH', #'wzttH',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'Top','VV','DYTT',
        'WWewk',
        # nuisance
        'WWpow',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        ('CHITOP-Top',     'Top'),
    ],
   'vbf_of-mH125' : [
        #signals
        'ggH','qqH', #'wzttH',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'Top','VV','DYTT',
        'WWewk',
        # nuisance
        'WWpow',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        ('CHITOP-Top',     'Top'),
        # mH125 as background
        'ggH_SM', 'qqH_SM', # 'WH_SM','ZH_SM',
    ],
   'ggH2j_of' : [
        #signals
        'ggH','qqH', #'wzttH',
        #'ggHminlo',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'Top','VV','DYTT',
        'WWewk',
        # nuisance
        'WWpow',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        #('CHITOP-Top',     'Top'),
    ],
    'vh_sf' : [
        #signals
        'ggH','qqH','WH','ZH',
        #'ttH','VH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV',
        'VVV',
        # nuisance
        'WWpow',
        #'WWnloNorm',
        #'WWnlo','WWnloUp','WWnloDown',
        # templates
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal'),
        # templates
        'VgS-template','Vg-template',
         # dummy
        ('DYee', 'WW'),
        ('DYmm', 'WW'),
    ],
   'vh_of' : [
        #signals
        'ggH','qqH','WH','ZH',
        #'ttH','VH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','DYTT',
        'VVV',
        # nuisance
        'WWpow',
        #'WWnloNorm',
        #'WWnlo','WWnloUp','WWnloDown',
        # templates
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal'),
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        #('WJetFakeRate-nominal', 'WJetFakeRate-vh-template-nominal'),
        #('WJetFakeRate-eUp',     'WJetFakeRate-vh-template-nominal'),
        #('WJetFakeRate-eDn',     'WJetFakeRate-vh-template-nominal'),
        #('WJetFakeRate-mUp',     'WJetFakeRate-vh-template-nominal'),
        #('WJetFakeRate-mDn',     'WJetFakeRate-vh-template-nominal'),
    ],
    'zh4j_mm' : [
        #signals
        'ggH','qqH','WH','ZH','ttH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV',
        'VVV',
         # dummy
        ('DYee', 'WW'),
        ('DYmm', 'WW'),
        # templates
        'VgS-template','Vg-template',
    ],
   'zh4j_ee' : [
        #signals
        'ggH','qqH','WH','ZH','ttH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV',
        'VVV',
         # dummy
        ('DYee', 'WW'),
        ('DYmm', 'WW'),
        # templates
        'VgS-template','Vg-template',
    ],
   'whsc' : [
        #signals
        'ggH','qqH','WH','ZH','ttH','VH',
        #'WH','ZH','ttH','VH',
        #'VH',
        # bkgs
        'WW','ggWW','VgS','Vg','WJet','Top','VV','DYTT',
        'DYee',
        'DYmm',
        'VVV',
        # templates
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal'),
        'VgS-template','Vg-template',
    ],

   'wwewkshape_of' : [
        #signals
        'ggH','qqH', #'wzttH',
        #'ggHminlo',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'VV','DYTT',
        #('TopPt0','Top'),
        ('TopPt1','Top'),  # 30-50
        ('TopPt2','Top'),  # 50-70
        ('TopPt3','Top'),  # 70-110
        ('TopPt4','Top'),  # 110-150
        ('TopPt5','Top'),  # 150-200
        ('TopPt6','Top'),  # 200-
        #('TopPt7','Top'),
        #('TopPt8','Top'),
        'WWewk',
        'VVV',
        # nuisance
        'WWpow',
        'WWewkMG',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        #('CHITOP-Top',     'Top'),
    ],

   'wwewk_of' : [
        #signals
        'ggH','qqH', #'wzttH',
        #'ggHminlo',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'VV','DYTT',
        ('TopPt0','Top'),
        ('TopPt1','Top'),
        #('TopPt2','Top'),
        'WWewk',
        'VVV',
        # nuisance
        'WWpow',
        'WWewkMG',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        #('CHITOP-Top',     'Top'),
    ],
   'wwewk_sf' : [
        #signals
        'ggH','qqH', #'wzttH',
        #'ggHminlo',
        # bkgs
        'WW','ggWW','VgS','Vg',
        ('WJet', 'WJet-2j-fix'),
        'VV',
         # dummy
        ('TopPt0','WWewk'),
        ('TopPt1','WWewk'),
        #('TopPt2','WWewk'),
         # dummy
        ('DYee', 'WWewk'),
        ('DYmm', 'WWewk'),
        'WWewk',
        'VVV',
        # nuisance
        'WWpow',
        'WWewkMG',
        # systematics
        ('WJetFakeRate-2j-template','WJetFakeRate-nominal-2j-fix'), # here and in the following I put the "template" distributions (relaxed cuts)
        ('WJetFakeRate-2j-eUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-eDn',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mUp',     'WJetFakeRate-nominal-2j-fix'),
        ('WJetFakeRate-2j-mDn',     'WJetFakeRate-nominal-2j-fix'),
        # top shape ---> to be checked if needed
        #('TopTW',   'Top'),
        #('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 2j specific
        #('WJet-template',    'WJet-template-2j'),              #    A   <-   sorgente
        #('WJet-templatesyst','WJet-templatesyst-2j')           #    mkmerged vuole "-template"
        # templates for Top estimation
        #('CHITOP-Top',     'Top'),
    ],

   'topestimate' : [
        # signals
        'ggH','qqH','wzttH',
        # backgrounds
        'WW','ggWW','VgS','Vg','WJet','Top','VV','DYTT','DYLL',
   ],
   'topplots' : [
        # signals
        'ggH','qqH','wzttH',
        # backgrounds
        'WW','ggWW','VgS','Vg','WJet','ttbar','tW','VV','DYTT','DYLL',
   ],

   'Hwidth_01j' : [
        #signals
        'ggH','qqH',
        'ggH_sbi','qqH_sbi',
        'ggH_b','qqH_b',
        'ggH_s','qqH_s',
        #'ggHminlo',
        # bkgs
        #('WW','WWmad'),
        ('WWlow', 'WWmad'),
        ('WWhigh','WWmad'),
        #'WW',
        #'ggWW',
        'VgS','Vg','WJet','Top','VV',
        #'VVV',
        ('VVV','VVVnoEWK'),
        'DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],

   'Hwidth_2j' : [
        #signals
        'ggH','qqH',
        'ggH_sbi','qqH_sbi',
        'ggH_b','qqH_b',
        'ggH_s','qqH_s',
        #'ggHminlo',
        # bkgs
        #('WW','WWmad'),
        ('WWlow', 'WWmad'),
        ('WWhigh','WWmad'),
        #'WW',
        #'ggWW',
        'VgS','Vg','WJet','Top','VV',
          #'VVV',
        ('VVV','VVVnoEWK'),
        'DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        #
        'Top'
        #('CHITOP-Top',     'Top'), -> to be fixed for data-driven estimation
        # templates
        'VgS-template','Vg-template',
    ],
   'Hwidth7TeV_01j' : [
        #signals
        'ggH','qqH',
        'ggH_sbi','qqH_sbi',
        'ggH_b','qqH_b',
        'ggH_s','qqH_s',
        #'ggHminlo',
        # bkgs
        'WW',
        #'ggWW',
        'VgS','Vg','WJet','Top','VV',
        #'VVV',
        #('VVV','VVVnoEWK'),
        'DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        ('TopTW',   'Top'),
        ('TopCtrl', 'Top'),
        # templates
        'VgS-template','Vg-template',
        # 0j1j specific
        #('DYLL-template',    'DYLL-template-0j1j'),              #    A   <-   sorgente
        #('DYLL-templatesyst','DYLL-templatesyst-0j1j')           #    mkmerged vuole "-template"
    ],

   'Hwidth7TeV_2j' : [
        #signals
        'ggH','qqH',
        'ggH_sbi','qqH_sbi',
        'ggH_b','qqH_b',
        'ggH_s','qqH_s',
        #'ggHminlo',
        # bkgs
        'WW',
        #'ggWW',
        'VgS','Vg','WJet','Top','VV',
          #'VVV',
        #('VVV','VVVnoEWK'),
        'DYTT','WWnlo','WWnloUp','WWnloDown',
        ('DYLL','ggH'),
        # systematics
        'WJetFakeRate-nominal',
        ('WJetFakeRate-eUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-eDn', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mUp', 'WJetFakeRate-nominal'),
        ('WJetFakeRate-mDn', 'WJetFakeRate-nominal'),
        #
        'Top'
        #('CHITOP-Top',     'Top'), -> to be fixed for data-driven estimation
        # templates
        'VgS-template','Vg-template',
    ],
}

#--------------
def samples(mass, energytag, datatag='Data2012', sigtag='SM', mctag='all'):
    '''
    mass: mass for the higgs samples'
    datatag: tag for the dataset to be included
    sigtag: kind of signal (SM or Graviton)
    mctag: tag for the set of mc to be included
    '''


    if energytag == '8TeV':
        import samples.samples8TeV as sampledb
    elif energytag == '7TeV':
        import samples.samples7TeV as sampledb
    else:
        raise ValueError('Energy misdefined: '+energytag)

    signals = sampledb.signalSamples(sigtag,mass)
    print "*** signals = ",signals

    mcsamples = {}
    mcsamples.update(signals)
    if ('vbf' in mctag or 'wwewk' in mctag or 'vh' in mctag) and  datatag == 'Data2012' :
      sampledb.backgrounds['WJet'] = ['wjets/latino_RunA_892pbinv_LooseLoose.root', 
                                      'wjets/latino_RunB_4404pbinv_LooseLoose.root', 
                                      'wjets/latino_RunC_7032pbinv_LooseLoose.root', 
                                      'wjets/latino_RunD_7274pbinv_LooseLoose.root']
      sampledb.backgrounds['WJetFakeRate-nominal'] = ['wjets/latino_RunA_892pbinv_LooseLoose.root',
                                                      'wjets/latino_RunB_4404pbinv_LooseLoose.root',
                                                      'wjets/latino_RunC_7032pbinv_LooseLoose.root',
                                                      'wjets/latino_RunD_7274pbinv_LooseLoose.root']
      print sampledb.backgrounds['WJet']
      print sampledb.backgrounds['WJetFakeRate-nominal'] 

    mcsamples.update(sampledb.backgrounds)


    if 'mH' in mctag:
        print mctag
        # get the background-higgs mass from the tag
        if   '0j1j'   in mctag : mHbkg = int(re.match('0j1j-mH(\d+)', mctag).group(1))
        elif 'vbf_of' in mctag : mHbkg = int(re.match('vbf_of-mH(\d+)', mctag).group(1))
        elif 'vbf_sf' in mctag : mHbkg = int(re.match('vbf_sf-mH(\d+)', mctag).group(1))
        print 'signal as background', mHbkg
        #signalbkg = sampledb.signalSamples(sigtag, mHbkg, str(mHbkg))
        signalbkg = sampledb.signalSamples(sigtag, mHbkg, "_SM")
        print signalbkg
        mcsamples.update(signalbkg)

    if isinstance(mctag,list):
        mclabels = mctag
    else:
        if mctag not in mcsets:
            raise ValueError('MCtag '+mctag+' not supported')
        mclabels = mcsets[mctag]

    selectedMc = hwwtools.filterSamples( mcsamples, mclabels )

    # add data
    selectedData = {}
    if datatag == 'NoData':
        pass
    elif 'Data' in datatag:
        selectedData['Data'] = sampledb.data[datatag]
    elif 'SI' in datatag:
        # add the signal samples of the given mass
        m = re.match('SI(\d+)',datatag)
        print 'SI with mH =', int(m.group(1))
        if not m:
            raise ValueError('Signal injection must have the format SImmm where mmm is the mass')
        simass = int(m.group(1))
        siSamples = sampledb.signalSamples(sigtag, simass,suffix='-SI')
        selectedData.update(siSamples)
    else:
        raise ValueError('Data tag '+datatag+' not supported')

    samples = {}
    samples.update(selectedMc)
    samples.update(selectedData)

    return samples



