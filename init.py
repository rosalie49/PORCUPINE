#Function init permits to load the data
def load_data():
    import pandas as pd
    import pyreadr
    #Loading of the network data
    net_file_path = input("Network file path => ")
    data = pyreadr.read_r(net_file_path)
    net = data['net']

    #Loading of the edges information
    edges_file_path = input("Edges file path => ")
    data_edges = pyreadr.read_r(edges_file_path)
    edges = data_edges['edges']

    return net, edges

#This function permits to load the gmt file with pathways
def load_gmt():
    pathways_file_path = input("Pathways file path => ")
    pathways_list = []
    with open(pathways_file_path, 'r') as file:
        for line in file:
            fields = line.split('\t')  
            pathway = fields[0]  
            genes = fields[2:]  
            pathways_list.append((pathway, genes))
    return pathways_list


#This function filters the list of pathways to include only genes in pathways present in networks
def filter_pathways(pathways_list,edges):
    import pandas as pd
    
    pathways_filt = pd.DataFrame(pathways_list, columns=['pathway', 'genes'])
    unique_edges = edges['tar'].unique()
    pathways_filt = pathways_filt[pathways_filt['genes'].apply(lambda x: any(gene in unique_edges for gene in x))]

    return pathways_filt


#This function filters a list of pathways based on specified minimum and maximum size for number of genes in a pathway
def filter_pathways_size(pathways_filt,minSize=5,maxSize=150):
    pathways_to_use = pathways_filt[
    (pathways_filt['genes'].apply(len) >= minSize) &
    (pathways_filt['genes'].apply(len) <= maxSize)
]
    return pathways_to_use