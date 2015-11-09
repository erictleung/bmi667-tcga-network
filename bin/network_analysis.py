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

# make larger graph
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

# import TCGA nodes
def get_TCGA():
    tcga = []
    with open("../data/tcga_genes.txt", "r") as fh:
        for line in fh:
            tcga.append(line.rstrip())
    return tcga

# make subgraph with TCGA genes
def tcga_subgraph(ppNet, tcga):
    return ppNet.subgraph(tcga)

# find largest connected component of subgraph
def largest_component(ppNet):
    subnetGraphs = nx.connected_component_subgraphs(ppNet)
    largest = nx.Graph() # empty graph
    greatN = 0
    for sub in subnetGraphs:
        if len(sub.nodes()) > greatN:
            largest = sub
            greatN = len(sub.nodes())
    return largest

# calculate centrality measures
def calculate_center(tcgaSubgraph):
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
    fileName = "../data/" + centralName + "Rounded.tsv"
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
    saveName = "../images/" + graphName + ".png"
    plt.savefig(saveName)

def bfs_edges(subnet, source):
    """
    INPUT: NetworkX undirected graph and string of node
    OUTPUT: generator of edges using the breadth first search algorithm
    Templated from: http://networkx.readthedocs.org/en/stable/_modules/networkx/algorithms/traversal/breadth_first_search.html#bfs_edges
    """
    neighbors = subnet.neighbors_inter
    visited = set([source])
    queue = [[source, neighbors(source)]]
    while queue:
        # children is an iterator
        parent, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent, child # add to iterator generator
                visited.add(child)
                queue.append([child, neighbors(child)])
        except StopIteration:
            queue.reverse()
            queue.pop()
            queue.reverse()

# main function
def main():
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
