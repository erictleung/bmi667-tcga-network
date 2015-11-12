# BMI667 TCGA Gene Network Analysis

Network analyses on 127 significantly mutated genes from well-known and emerging
cellular processes in cancer using Python and Shell scripting. The paper
[Mutational landscape and significance across 12 major cancer
types](http://www.ncbi.nlm.nih.gov/pubmed/24132290) by Kandoth et. al. published
about the 127 significantly mutated genes.

## Execute Analysis

```shell
$ bash extract_and_save.sh
$ python network_analysis.py
$ python distance_and_paths.py
```

## Analysis Instructions

### Part 1

- Download Reactome interaction files in PSI-MITAB format from
  http://www.reactome.org/download with the link: `Human protein-protein
  interaction pairs in PSI-MITAB format`. 
- Extract interactions into pair-wise gene symbols: one interaction per line
  (e.g. NFS1\tITGA7). 
- Map the 127 TCGA PanCancer genes to create a sub-network. 
- Calculate four centralities (degree, eigenvector, PageRank, Katz) each gene in
  the sub-network.

### Part 2

- Implement breadth first search (BFS) algorithm and apply it to the sub-network
  for 127 cancer genes you have generated.
- Find how many components.
- Calculate pair-wise distances.
- Calculate average shortest path and the diameter of the sub-network.

## Contents

- `bin`: scripts directory
- `data`: data repository
- `report`: written report for the analysis

## Scripts

TBD

## Data

TBD
