import sys

sys.path.append('./')

import pandas as pd
from run_pca import * 
import pytest 

data_df = pd.read_csv("random_dataset_run_pca.csv")

def test_run_pca():
    result = run_pca(data_df)
    
    assert result['pc1'].iloc[0] == pytest.approx(40.076, abs=0.001)
    assert result['n_edges'].iloc[0] == 10