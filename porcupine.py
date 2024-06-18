import pandas as pd
import numpy as np
from scipy.stats import stats

def calculate_statistics(res_pca_pathway, res_pca_rndm_set):
    pc1_rndm = res_pca_rndm_set['pc1'].dropna().values
    pathway = res_pca_pathway["pathway"]
    path_size = res_pca_pathway["pathway_size"]
    pc1_pathway = float(res_pca_pathway["pc1"])
    
    _, pvalue = stats.ttest_1samp(pc1_rndm, popmean=pc1_pathway, alternative='less')

    if len(pc1_rndm) > 1:
        effect_size = (np.mean(pc1_rndm) - pc1_pathway) / np.std(pc1_rndm, ddof=1)
    else:
        effect_size = np.nan
    
    stat_res = pd.DataFrame({
        "pathway": [pathway],
        "pathway_size": [path_size],
        "pc1": [pc1_pathway],
        "pval": [pvalue],
        "es": [effect_size]
    })
    
    return stat_res


def porcupine(res_pca_pathways, res_pca_rndm):
    pathways_size = res_pca_pathways['pathway_size'].unique()
    res_all = pd.DataFrame()
    for path_size in pathways_size:
        print(path_size)
        res_pathway_set = res_pca_pathways[res_pca_pathways['pathway_size'] == path_size]
        print(res_pathway_set)
        res_rndm_set = res_pca_rndm[res_pca_rndm['pathway_size'] == path_size]
        res = res_pathway_set.apply(lambda x: calculate_statistics(x, res_rndm_set), axis=1)
        res_all = pd.concat([res_all, pd.concat(res.tolist(), ignore_index=True)], ignore_index=True)    
        

    return res_all
