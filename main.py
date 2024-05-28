from init import *


def main():
    net, edges = load_data()
    pathways=load_gmt()
    print(len(pathways))
    pathways_filt=filter_pathways(pathways,edges)
    print(len(pathways_filt))
    pathways_to_use=filter_pathways_size(pathways_filt,5,150)
    print(len(pathways_to_use))

main()