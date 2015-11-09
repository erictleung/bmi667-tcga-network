#!/usr/bin/env python

# import packages
import networkx as nx
import operator
import matplotlib.pyplot as plt

# import interactions as edge list
def get_prot_inter():
    ppList = []
    with open("prot_interaction.tsv", "r") as fh:
        for line in fh:
            splitLine = line.rstrip().split("\t")
            ppList.append(splitLine)
    return ppList

# make larger graph
def make_net(ppList):
    ppNet = nx.Graph() # create empty network
    for pp in ppList:
        ppNet.add_edge(pp[0], pp[1]) # add interaction to network
    return ppNet

# import TCGA nodes
def get_TCGA():
    tcga = []
    with open("tcga_genes.txt", "r") as fh:
        for line in fh:
            tcga.append(line.rstrip())
    return tcga

# make subgraph with TCGA genes
def tcga_subgraph(ppNet, tcga):
    return ppNet.subgraph(tcga)

# calculate centrality measures
def calculate_center(tcgaSubgraph):
    centers = {}
    centers["eigen"] = nx.eigenvector_centrality(tcgaSubgraph)
    centers["degree"] = nx.degree_centrality(tcgaSubgraph)
    centers["katz"] = nx.katz_centrality_numpy(tcgaSubgraph)
    centers["pagerank"] = nx.pagerank(tcgaSubgraph)
    return centers

# write file with centrality measures for nodes
def write_central(centrality, centralName):
    fileName = centralName + "Rounded.tsv"
    sortedCenter = sorted(centrality.items(), key = operator.itemgetter(1),
            reverse = True)
    with open(fileName, "w") as fh:
        for pair in sortedCenter:
            fh.write(pair[0] + "\t" + str(round(pair[1], 4)) + "\n")

# write file for each of the centralities
def write_all(centers):
    for center in centers.keys():
        write_central(centers[center], center)

# plot and save network
def plot_net(graph, graphName):
    nx.draw(graph)
    saveName = graphName + ".png"
    plt.savefig(saveName)

# main function
def main():
    ppList = get_prot_inter()
    ppNet = make_net(ppList)
    tcga = get_TCGA()
    tcgaSubgraph = tcga_subgraph(ppNet, tcga)
    centers = calculate_center(tcgaSubgraph)
    write_all(centers)
    plot_net(tcgaSubgraph, "tcga_subgraph")

if __name__ == '__main__':
    main()
