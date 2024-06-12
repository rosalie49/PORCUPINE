import sys

sys.path.append('./')

from init import *

edges = load_edges_file('rand_genes.RData')
pathways = load_gmt('c2.cp.reactome.v7.1.symbols.gmt')

def test_filter_pathways():
    pathways_filt = filter_pathways(pathways, edges)
    assert len(pathways_filt) == 1042  # R result

    

def test_filter_pathways_size():
    pathways_filt = filter_pathways(pathways, edges)
    pathways_filt_size = filter_pathways_size(pathways_filt,minSize=5,maxSize=150)
    assert len(pathways_filt_size) == 210 # R result 