lumi=5.064
rebin=5
chans=['of_0j','of_1j','sf_0j','sf_1j']

# 'mll' selection also exist, but use bdt for this one
variable='mll'
selection='bdt'
range='mllsplit'
split='mll'

tag='mllsplit_bdtsel'
xlabel='Split in hi/lo S/B, m(ll) [GeV]'
dataset='Data2012'
mcset='0j1j'

# directories
path_latino='/shome/mtakahashi/HWW/Tree/ShapeAna/tree_skim_wwcommon'
path_bdt='/shome/mtakahashi/HWW/Tree/ShapeAna/bdt_skim_wwcommon/mva_MH{mass}_{category}'
path_dd='/shome/mtakahashi/HWW/Data/dd/bdt_2012_51fb'
#path_latino_dd='/shome/mtakahashi/HWW/Tree/ShapeAna/tree_skim_ddskim'

# other directories
path_shape_raw='raw'
path_shape_merged='merged'

