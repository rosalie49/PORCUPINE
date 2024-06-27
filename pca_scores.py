import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def get_ind_scores(data, scale_data=True, center_data=True):
    """This function extracts patient heterogeneity scores on the first two principal components.

    Args:
        data (table): numeric matrix with samples in columns and features in rows
        scale_data (bool, optional): Whether to scale the data (TRUE) or not (FALSE). Defaults to True.
        center_data (bool, optional): Whether to center the data (TRUE) or not (FALSE). Defaults to True.

    Returns:
        pd.DataFrame: patient heterogeneity scores for PC1 and PC2
    """
    data_t = np.transpose(data)

    # Scale and/or center the data if specified
    if scale_data or center_data:
        # Initialize the StandardScaler object with mean centering and scaling
        scaler = StandardScaler(with_mean=center_data, with_std=scale_data)
        # Scale and center the data using the scaler object
        data_t = scaler.fit_transform(data_t)

    pca = PCA(n_components=2) # Initialize PCA object
    individual_scores = pca.fit_transform(data_t)  # Fit PCA to data

    # Return a dataframe with the first two principal components
    return pd.DataFrame(individual_scores, columns=['pc1', 'pc2'], index = data.columns)

def get_pathway_ind_scores(pathways_list, reg_net, edges, scale_data=True, center_data=True):
    """This function extracts patient heterogeneity scores on the first two principal components for a list of pathways.

    Args:
        pathways_list (pd.DataFrame): list of pathways
        reg_net (pd.DataFrame): Numeric matrix with samples in columns, features in rows
        edges (pd.DataFrame): containing information on "reg" and "tar"
        scale_data (bool, optional): Whether to scale the data (TRUE) or not (FALSE). Defaults to True.
        center_data (bool, optional): Whether to center the data (TRUE) or not (FALSE). Defaults to True.

    Returns:
        pd.DataFrame: patient heterogeneity scores for each pathway
    """
    # Initialize an empty list to store results for each pathway
    res = []
    # Iterate over each pathway in the pathways list
    for _, pathway in pathways_list.iterrows():
        # Identify the indices of edges where the target genes are in the pathway's gene list
        all_genes = pathway['genes']
        idx = edges[edges['tar'].isin(all_genes)].index
        # Extract the subset of the regulatory network corresponding to these indices
        subnet = reg_net.loc[idx]
        # Perform PCA analysis on the subset
        individual_scores = get_ind_scores(subnet, scale_data=scale_data, center_data=center_data)
        # Store the results in a DataFrame with patient as index, a pathway column and scores of the first two principal components
        res_pca = pd.DataFrame({'pathway':pathway.iloc[0], 'pc1':individual_scores['pc1'], 'pc2':individual_scores['pc2']})
        res_pca.index = individual_scores.index
        res.append(res_pca)

    res_all = pd.concat(res)   

    return res_all