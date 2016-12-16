__author__ = 'orsegal'


import os
# from GenAssignmentModified import *
from tabulate import tabulate
from hackingtrees.sequences import seq_reader
# from hackingtrees.trees.node import Node
from hackingtrees.trees.tree_printer import *
import time
import sys
import pickle
import operator
import copy


class Node:
    def __init__(self):
        self.label = ""
        self.length = 0.0
        self.time_length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False
        self.height = 0
        self.id = None

    def add_child(self,child):
        # make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self

    def remove_child(self,child):
        #make sure that the child is in there
        assert child in self.children
        self.children.remove(child)
        child.parent = None

    def leaves(self,v=None):
        if v == None:
            v = []
        if len(self.children) == 0:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self):
        return [n for n in self.iternodes() if n.istip ]

    def lvsnms(self):
        return [n.label for n in self.iternodes() if n.istip ]

    def iternodes(self,order="preorder"):
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order.lower() == "postorder":
            yield self

    def prune(self):
        p = self.parent
        if p != None:
            p.remove_child(self)
        return p

    def get_newick_repr(self,showbl=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
        if showbl == True:
            ret += ":" + str(self.length)
        return ret

    def set_height(self):
        if len(self.children) == 0:
            self.height = 0
        else:
            tnode = self
            h = 0
            while len(tnode.children) > 0:
                if tnode.children[1].length < tnode.children[0].length:
                    tnode = tnode.children[1]
                else:
                    tnode = tnode.children[0]
                h += tnode.length
            self.height = h



def readStrings(filename):
    # Using the hackingtrees function for reading fasta files, import the strings.
    # Shorten the string names.
    seqs = seq_reader.read_fasta_file(filename)
    for i in seqs:
        strtInd = i.name.index('(')
        stringName = i.name[strtInd+3:]
        endInd = stringName.index('(')
        stringName = stringName[:endInd]
        i.name = stringName
    return seqs


def importNewick(filename):
    with open(filename, 'r') as f:
        for line in f:
            tree = line
    return tree


def seqsToNodes(seqs):
    allNodes = []
    tempID = 0
    for seq in seqs:
        tempNode = Node()
        tempNode.label = seq.name
        tempNode.data['seq'] = seq.seq
        tempNode.data['seqNoAlign'] = seq.seq.translate(None, '-')
        tempNode.data['NN'] = None
        tempNode.data['distToNN'] = None
        tempNode.id = tempID
        tempNode.istip = True
        tempID += 1
        allNodes.append(tempNode)
    return allNodes


def EditDistance(A, B):
    m = len(A)
    n = len(B)
    C = [[0] * (n+1) for _ in range(m+1)]

    # Distance of string to the other when it's empty
    for i in xrange(m+1):
        C[i][0] = i
    for j in xrange(n+1):
        C[0][j] = j

    for j in xrange(1, n+1):
        for i in xrange(1, m+1):
            if A[i-1] == B[j-1]:
                C[i][j] = C[i-1][j-1]
            else:
                deletion = C[i-1][j] + 1
                insertion = C[i][j-1] + 1
                substitution = C[i-1][j-1] + 1
                C[i][j] = min(deletion, insertion, substitution)
    return C[m][n]

def distanceTableToPrint(strings):
    Nstrings = len(strings)
    DTable = [[0] * (Nstrings) for _ in range(Nstrings)]
    for i in xrange(Nstrings):
        for j in xrange(i+1, Nstrings):
            DTable[i][j] = DTable[j][i] = EditDistance(strings[i], strings[j])
    return DTable

def StringsTableAnalysis(STable, strings, NHeader = None):
    Nstrings = len(strings)
    if NHeader == None:
        NHeader = range(Nstrings)
    print tabulate(STable, showindex='always', headers=NHeader, floatfmt=".3f")


def storeTable(nodes, filename, tableFunc):
    # Store the distance table to file. Building the table is the most costly function so saving it saves time.

    STable = tableFunc(nodes)
    with open(filename, 'wb') as fp:
        pickle.dump(STable, fp)


def openTable(filename):
    with open(filename, 'rb') as fp:
        Table = pickle.load(fp)
    return Table



def distanceTable(nodes):
    # Build a distance table based on the edit distances. The table is built as a dictionary of dictionaries
    # (instead of a list of lists) in order to be able to edit it later and change the names of the nodes

    DTable = {}
    for node1 in nodes:
        DTable[node1.id] = {}
        for node2 in nodes:
            if node1 != node2:
                DTable[node1.id][node2.id] = 0
    keepTrack = []
    for node1 in nodes:
        keepTrack.append(node1.id)
        for node2 in nodes:
            if node2.id not in keepTrack:
                DTable[node1.id][node2.id] = DTable[node2.id][node1.id] = EditDistance(node1.data['seqNoAlign'],
                                                                                       node2.data['seqNoAlign'])
    return DTable



def updateNN(clusters, DTable):
    # Goes over all nodes and searches for the closest node or cluster of nodes (which is represented
    # as a node with children).

    for node in clusters:
        node.data['NN'] = None
        node.data['distToNN'] = None
        tempDist = float('inf')
        NNID = 0
        for key in DTable[node.id]:
            if DTable[node.id][key] < tempDist:
                tempDist = DTable[node.id][key]
                NNID = key
        node.data['NN'] = [x for x in clusters if x.id == NNID][0]
        node.data['distToNN'] = tempDist


def getClosestPair(clusters):
    minDist = float('inf')
    for node in clusters:
        if node.data['distToNN'] < minDist:
            minDist = node.data['distToNN']
            c1 = node
            c2 = node.data['NN']
            assert c2 is not None
            assert c2 is not c1
    return (c1, c2)


def clusterDistance(newC, c2, originalDTable):
    newClen = len(newC.data['nodesInBranch'])
    if len(c2.children) == 0:
        c2len = 1
    else:
        c2len = len(c2.data['nodesInBranch'])

    dist = 0
    for i in [n for n in newC.iternodes() if n.istip]:
        if c2len == 1:
            dist += originalDTable[i.id][c2.id]
        else:
            for j in [n for n in c2.iternodes() if n.istip]:
                dist += originalDTable[i.id][j.id]
    dist = dist/(newClen * c2len)
    return dist


def createOTU(c1, c2, DTable, originalDTable, clusters):
    # remove the two clusters that we will combine
    clusters.remove(c1)
    clusters.remove(c2)

    # create a new middle node
    newCluster = Node()
    newCluster.label = ("(" + str(c1.id) + ", " + str(c2.id) + ")")
    newCluster.data['NN'] = None
    newCluster.data['distToNN'] = None
    if len(clusters) != 0:
        newCluster.id = max(node.id for node in clusters) + 1
    newCluster.add_child(c1)
    newCluster.add_child(c2)
    newCluster.data['nodesInBranch'] = [n.id for n in newCluster.iternodes() if n.istip]

    clusters.append(newCluster)

    # Update the distance table
    ncDist = {}                     # dictionary of distances from the new cluster to append later
    for key in DTable:
        if (key != c1.id) and (key != c2.id):
            Ncomp = [x for x in clusters if x.id == key][0]
            ncDist[key] = clusterDistance(newCluster, Ncomp, originalDTable)

    DTable.pop(c1.id)
    DTable.pop(c2.id)
    for key in DTable:
        DTable[key].pop(c1.id)
        DTable[key].pop(c2.id)
        DTable[key][newCluster.id] = ncDist[key]
    DTable[newCluster.id] = ncDist


def UPGMA(clusters, originalDTable):
    DTable = copy.deepcopy(originalDTable)   # The distance table copy that we actually manipulate as the algorithm runs

    for i in xrange(len(clusters) - 1):
        updateNN(clusters, DTable)
        c1, c2 = getClosestPair(clusters)
        createOTU(c1, c2, DTable, originalDTable, clusters)

    return clusters[0]


def main():
    # strings = importStrings("seqdump.txt")
    # strLst = []
    # for strng in strings:
    #     strLst.append(strng[1])
    #
    # # storeTable(strLst, 'STable.txt', stringsTable)
    # STable = openSTable('STable.txt')
    # StringsTableAnalysis(STable, strLst)
    # ratioSTable = stringsRatioTable(STable, strLst)
    # ratioSTableAnalysis(ratioSTable, strLst, NHeader=None)
    #
    # storeTable(strLst, 'DTable.txt', distanceTable)
    # DTable = openSTable('DTable.txt')
    # StringsTableAnalysis(DTable, strLst)


    # Read strings from file
    seqs = readStrings("FluGeneSeq-ClustalWAlligned10.fas")

    # Convert strings to nodes
    nodes = seqsToNodes(seqs)

    # Print a distance table
    strLst = []
    for node in nodes:
        strLst.append(node.data['seqNoAlign'])
    storeTable(strLst, 'DTableToPrint.txt', distanceTable)
    DTableToPrint = openTable('DTable.txt')
    StringsTableAnalysis(DTableToPrint, strLst)

    #
    #
    #
    # storeTable(nodes, 'DTable.txt', distanceTable)
    # DTable = openTable('DTable.txt')
    #
    # tree = UPGMA(nodes, DTable)
    # print tree2ascii(tree, scaled=0, show_internal_labels=0)



if __name__ == "__main__":
    main()