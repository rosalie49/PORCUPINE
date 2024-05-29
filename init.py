def load_data(net_file,edges_file):
    '''
    Loads network data and edge information from files specified by the user

    Parameters:
    - net_file (str): path to the input networks data file
    - edges_file (str) : path to the input edges data file 

    Returns:
        - net (DataFrame): The network data
        - edges (DataFrame): The edge data
    '''
    import pandas as pd
    import pyreadr
    
    data = pyreadr.read_r(net_file)
    net = data['net']  # Network data is stored in net
    
    data_edges = pyreadr.read_r(edges_file)
    edges = data_edges['edges']  # Edge data is stored in edges
    
    return net, edges


def load_gmt(gmt_file):
    '''
    Load a .gmt file containing information about biological pathways

    Parameters:
    gmt_file (str): Path to the input .gmt file. Gmt file is downloaded
        from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

    Returns:
    A list of pathways from the .gmt file
    '''
    
    pathways_list = []
    
    # Reading the .gmt file and extracting pathway and gene information
    with open(gmt_file, 'r') as file:
        for line in file:
            fields = line.split('\t')  # Fields are tab-separated 
            pathway = fields[0]  # Pathway name is the first field
            genes = fields[2:]  # Genes start from the third field
            pathways_list.append((pathway, genes))  # Adding pathway name and gene list to pathways_list
    
    return pathways_list


def filter_pathways(pathways_list, edges):
    """
    Filter a list of pathways to include only genes present in networks

    Parameters:
    pathways_list : list of pathways
    edges : dataframe containing information on "reg" and "tar"

    Returns:
    A list of filtered pathways
    """
    import pandas as pd
    
    # Creating a dataframe from the pathways list with columns 'pathway' and 'genes'
    pathways_filt = pd.DataFrame(pathways_list, columns=['pathway', 'genes'])
    
    # Extracting unique genes from edge data
    unique_edges = edges['tar'].unique()
    
    # Filtering pathways by keeping only those whose genes are present in the unique genes from edges
    pathways_filt = pathways_filt[pathways_filt['genes'].apply(lambda x: any(gene in unique_edges for gene in x))]
    
    return pathways_filt


def filter_pathways_size(pathways_filt, minSize=5, maxSize=150):
    '''
    Filter a list of pathways based on specified minimum and maximum size for number of genes in a pathway

    Parameters:
    pathways_list : list of pathways
    minSize (int): minimum size for number of genes in a pathway (default = 5)
    maxSize (int): maximum size for number of genes in a pathway (default = 150)

    Returns:
    A list of filtered pathways
    '''
    # Filtering pathways based on the length of the gene list
    pathways_to_use = pathways_filt[
        (pathways_filt['genes'].apply(len) >= minSize) &  # Filtering pathways with specified minimum size
        (pathways_filt['genes'].apply(len) <= maxSize)    # Filtering pathways with specified maximum size
    ]
    return pathways_to_use








    