import sys

sys.path.append('./')

import pandas as pd
from run_pca import * 
from pca_pathway import process_pathway
import pytest 

data_df = pd.read_csv("random_dataset_run_pca.csv")

def test_run_pca():

    result = run_pca(data_df)
    
    assert result['pc1'].iloc[0] == pytest.approx(40.076, abs=0.001) #R result 
    assert result['n_edges'].iloc[0] == 10


def test_process_pathway():
    net_data = data_df
    pathway_data = {'pathway':'PATHWAY1', 'genes' : ['gene1','gene2','gene3','gene4','gene5']}
    edges_data = {
    'reg': ['geneA', 'geneB', 'geneC', 'geneD', 'geneE'], 
    'tar': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']}
    edges_df = pd.DataFrame(edges_data)
    
    res = process_pathway(pathway_data,net_data, edges_df)

    assert res[1] == pytest.approx(56.856, abs=0.001) # 
    assert res[2] == 5
    assert res[3] == 5

