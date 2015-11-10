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

def bfs_edges(subnet, source):
    """
    INPUT: NetworkX undirected graph and string of node
    OUTPUT: generator of edges using the breadth first search algorithm
    Templated from: http://networkx.readthedocs.org/en/stable/_modules/networkx/algorithms/traversal/breadth_first_search.html#bfs_edges
    """
    neighbors = subnet.neighbors_iter
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

def num_components(subgraph, centrality):
    """
    DESCRIPTION: Find number of components in graph
    INPUT: Graph object and dictionary of centrality measures
    OUTPUT: Number of components
    """
    # get gene with largest centrality measure
    start = max(centrality.iteritems(), key=operator.itemgetter(1))[0]

    comp = 1 # count number of components
    nodes = set(subgraph.nodes())
    unvisited = set(nodes)

    paths = bfs_edges(subgraph, start) # get paths from source

    mapPaths = map(set, paths) # convert inner lists to sets
    reducePath = reduce(lambda x, y: x.union(y), mapPaths) # reduce to one set

    unvisited = unvisited.difference(reducePath) # find left over nodes

    # continue doing breadth first searches on left over nodes to find
    # components
    while unvisited:
        node = unvisited.pop() # take random node to examine
        paths = bfs_edges(subgraph, node)
        mapPaths = map(set, paths)
        reducePath = reduce(lambda x, y: x.union(y), mapPaths)
        unvisited = unvisited.difference(reducePath) # find left over nodes
        comp += 1

    print "There are %d components in the subgraph." % (comp)
    return comp

def pairwise_dist(graph):
    """
    DESCRIPTION: Find all pairwise distances between genes
    INPUT: Graph object
    OUTPUT: dictionary with keys being genes and values being dictionaries of
    distances
    """
    print "Calculating all pairwise distances."
    allDist = {} # dictionary to hold pairwise distances for all nodes
    for node in graph.nodes():
        paths = bfs_edges(graph, node) # breadth first search algorithm
        dist = {} # dictionary to carry total distance from source
        for edge in paths:
            # need to start distance calculation with source node
            if node in edge:
                neighbor = set(edge).difference(set([node])).pop()
                dist[neighbor] = 1
            else:
                if edge[0] in dist.keys():
                    dist[edge[1]] = dist[edge[0]] + 1
                else:
                    dist[edge[0]] = dist[edge[1]] + 1
        allDist[node] = dist
    print "Finishing calculating all pairwise distances"
    return allDist

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
    comp = num_components(tcgaSubgraph, centers["degree"])
    allDist = pairwise_dist(tcgaSubgraph)

if __name__ == '__main__':
    main()
