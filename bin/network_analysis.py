#!/usr/bin/env python

# Python: 2.7.8

# import packages
import networkx as nx
import operator
import matplotlib.pyplot as plt
import numpy.linalg

def get_prot_inter():
    """
    DESCRIPTION: Import interactions as an edge list
    INPUT: None
    OUTPUT: List of edges in separate lists
    """
    ppList = []
    print "Importing protein interaction pairs into a list."
    with open("../data/prot_interaction.tsv", "r") as fh:
        for line in fh:
            splitLine = line.rstrip().split("\t")
            ppList.append(splitLine)
    print "Successfully imported protein interaction pairs into a list."
    print "There are %d edges in this network" % (len(ppList))
    return ppList

def make_net(ppList):
    """
    DESCRIPTION: Make a graph object with adjancency list
    INPUT: List of edges in separate lists
    OUTPUT: Simple NetworkX graph object
    """
    ppNet = nx.Graph() # create empty network
    for pp in ppList:
        ppNet.add_edge(pp[0], pp[1]) # add interaction to network
    return ppNet

def get_TCGA():
    """
    DESCRIPTION: Import TCGA genes
    INPUT: None
    OUTPUT: List of imported TCGA genes
    """
    tcga = []
    with open("../data/tcga_genes.txt", "r") as fh:
        for line in fh:
            tcga.append(line.rstrip())
    return tcga

def tcga_subgraph(ppNet, tcga):
    """
    DESCRIPTION: Make a subgraph with the TCGA genes
    INPUT: Graph object and list of TCGA genes
    OUTPUT: Subgraph (graph) object made with only TCGA genes
    """
    return ppNet.subgraph(tcga)

def largest_component(ppNet):
    """
    DESCRIPTION: Find largest connected component of subgraph
    INPUT: Graph object
    OUTPUT: Largest component subgraph
    """
    subnetGraphs = nx.connected_component_subgraphs(ppNet)
    largest = nx.Graph() # empty graph
    greatN = 0
    for sub in subnetGraphs:
        if len(sub.nodes()) > greatN:
            largest = sub
            greatN = len(sub.nodes())
    return largest

def calculate_center(tcgaSubgraph):
    """
    DESCRIPTION: Calculate centrality measures
    INPUT: Graph object
    OUTPUT: Dictionary of dictionaries, each being a different centrality
    measure
    """
    # calculate maximum eigenvalue of graph
    denseMat = nx.adjacency_matrix(tcgaSubgraph).todense() # convert to adj mat
    eigs = numpy.linalg.eig(denseMat)[0] # calculate eigenvalues
    maxEig = max(eigs)
    alpha = 1 / maxEig.real

    # calculate centrality measures
    centers = {}
    centers["eigen"] = nx.eigenvector_centrality(tcgaSubgraph)
    centers["degree"] = nx.degree_centrality(tcgaSubgraph)
    centers["katz"] = nx.katz_centrality(tcgaSubgraph, alpha=alpha-0.01, beta=1.0)
    centers["pagerank"] = nx.pagerank(tcgaSubgraph)
    return centers

# write file with centrality measures for nodes
def write_central(centrality, centralName):
    """
    DESCRIPTION: Write a file with the centrality measures for nodes as a .tsv
    file
    INPUT: Dictionary of centrality and name of centrality
    OUTPUT: A tab-delimited file of centrality measure
    """
    fileName = "../data/" + centralName + "Rounded.tsv"
    sortedCenter = sorted(centrality.items(), key = operator.itemgetter(1),
            reverse = True)
    with open(fileName, "w") as fh:
        for pair in sortedCenter:
            fh.write(pair[0] + "\t" + str(round(pair[1], 4)) + "\n")

def write_all(centers):
    """
    DESCRIPTION: Write file for each of the centralities
    INTPUT: Dictionary with each of the centralities
    OUTPUT: Separate file for each centrality measure
    """
    for center in centers.keys():
        write_central(centers[center], center)

def plot_net(graph, graphName):
    """
    DESCRIPTION: Plot and save network as image
    INPUT: Graph object and desired name of image file
    OUTPUT: Plotted network visualization as .png file
    """
    nx.draw(graph)
    saveName = "../images/" + graphName + ".png"
    plt.savefig(saveName)

def main():
    """
    DESCRIPTION: Main function to run entire analysis
    INPUT: None
    OUTPUT: The network analysis and plot of TCGA subgraph
    """
    ppList = get_prot_inter()
    ppNet = make_net(ppList)
    tcga = get_TCGA()
    tcgaSubgraphRaw = tcga_subgraph(ppNet, tcga)
    tcgaSubgraph = largest_component(tcgaSubgraphRaw)
    centers = calculate_center(tcgaSubgraph)
    write_all(centers)
    plot_net(tcgaSubgraph, "tcga_subgraph")

if __name__ == '__main__':
    main()
