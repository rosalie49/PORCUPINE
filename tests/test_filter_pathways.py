import sys

sys.path.append('./')

from init import *

def test_filter_pathways():
    edges = load_edges_file('rand_genes.RData')
    pathways = load_gmt('c2.cp.reactome.v7.1.symbols.gmt')
    pathways_filt = filter_pathways(pathways, edges)
    assert len(pathways_filt) == 1042  # R result