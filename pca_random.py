from random import *
from pca_pathway import *

def create_gene_set(universe, psize, n_perm=1000):
    universe_list = list(universe) 
    gene_set = pd.DataFrame({
        'pathway': [i + 1 for i in range(n_perm)],
        'genes': [random.sample(universe_list, psize) for _ in range(n_perm)]})
    return gene_set


def pca_random(reg_net, edges, res_pca_pathways, pathways_list, n_perm=1000, ncores=1, scale_data = True, center_data = True):
    pathways_size = res_pca_pathways['pathway_size'].unique()
    universe = set(gene for pathway, genes in pathways_list for gene in genes)
    res_pca_random = []
    for psize in pathways_size:
        print("Pathways with size", psize)
        random_genes = create_gene_set(universe, psize, n_perm=n_perm)
        res_pca = pca_pathway(random_genes, reg_net, edges,scale_data=scale_data, center_data=center_data)
        res_pca_random.append(res_pca)

    res_pca_random_all = pd.concat(res_pca_random, ignore_index=True)
    return res_pca_random_all


    

