# PORCUPINE
**P**rincipal Components Analysis to **O**btain **R**egulatory **C**ontributions **U**sing **P**athway-based Interpretation of **N**etwork **E**stimates is an R package to identify biological pathway which drive inter-tumour heterogeneity in a population of gene regulatory networks. It is a Principal Components Analysis (PCA)-based approach that can be used to determine whether a specific set of variables—for example a set of genes in a specific pathway—have coordinated variability in their regulation.

## Method
PORCUPINE uses as input individual patient networks, for example networks modeled using PANDA and LIONESS, as well as a .gmt file that includes biological pathways and the genes belonging to them. For each pathway, it extracts all edges connected to the genes belonging to that pathway and scales each edge across individuals. It then performs a PCA analysis on these edge weights, as well as on a null background that is based on random pathways. To identify significant pathways, PORCUPINE applies a one-tailed t-test and calculates the effect size (ES). Here we provide an example of how we analyzed heterogeneity among single gene regulatory asmple networks in Leiomyosarcomas (LMS).  Networks were obtained with PANDA and LIONESS algorithms. Our  dataset contains data for 80 TCGA LMS patient-specific gene regulatory networks. 

## Setup
The requirements are provided in a requirements.txt file.
Install package 
```{r}
pip install git+https://github.com/yourusername/PORCUPINE
```

## Usage
```{r}
import porcupine as pcp
```
First, we load the network data and edges information. 
In this example, we have patient-specific gene regulatory networks for 80 TCGA leiomyosarcoma patients.
The first three columns in the data provide information on the regulators (TFs) and target genes. 
Edges information includes three columns: reg (the transcription factor's gene symbol),tar (Ensembl ID), prior (whether an edge is prior (1) or not (0)).

```{r}
net, edges = pcp.load_data(net_file_path, edges_file_path)
print(net.head())

   0E244FE2-7C17-4642-A51F-2CCA796D9C70  75435ED8-93E8-45FB-8480-98D8EB2EF8CB  \
0                                  0.76                                  0.10   
1                                  0.94                                  1.43   
2                                  1.09                                  2.78   
3                                  1.13                                  2.60   
4                                 -0.71                                 -1.42   

   B6D11678-15A9-4F43-A0A2-225067DCAF1C  B7F5A41E-9559-4329-81F5-1B88A74730B7  \
0                                 -1.27                                  0.01   
1                                  0.30                                  0.91   
2                                  1.01                                  2.13   
3                                  1.66                                  1.71   
4                                  0.02                                  0.27   

   04823F53-A12D-4852-8F34-77B9DCBB7DF0  49684C2B-D31C-4B45-A400-3497C3CCEC01  \
0                                 -7.18                                  2.74   
1                                 -5.69                                  2.43   
2                                 -6.09                                  3.23   
3                                 -6.04                                  2.88   
4                                  5.31                                 -1.02   

   FFDD7A12-DDEF-4974-8D60-64B7EEAAC994  830DFA6F-A85A-4317-82B2-791FAB998A01  \
0                                 -0.80                                 -0.31   
1                                 -1.55                                 -2.56   
2                                  1.04                                 -2.07   
...
3                                  0.13  
4                                 -0.06 


print(edges.head())

    reg    tar  prior
0  AIRE  AAGAB    1.0
1  ALX1  AAGAB    0.0
2  ALX3  AAGAB    0.0
3  ALX4  AAGAB    0.0
4    AR  AAGAB    0.0


```
