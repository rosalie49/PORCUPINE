import pandas as pd
import numpy as np
from run_pca import *

def pca_pathway(pathways_list, reg_net, edges, scale_data=True, center_data=True):
    results = []
    for index, pathway in pathways_list.iterrows():
        result = process_pathway(pathway, reg_net, edges, scale_data=scale_data, center_data=center_data)
        results.append(result)
    res_pca = pd.DataFrame(results, columns=['pathway', 'pc1', 'n_edges', 'pathway_size'])
    return res_pca

def process_pathway(pathway, reg_net, edges, scale_data=True, center_data=True):
    idx = edges[edges['tar'].isin(pathway['genes'])].index
    subnet = reg_net.loc[idx]
    pca_result = run_pca(subnet, scale_data=scale_data, center_data=center_data)
    return [pathway, pca_result['pc1'].iloc[0], pca_result['n_edges'].iloc[0], len(pathway['genes'])]