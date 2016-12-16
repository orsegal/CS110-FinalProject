import sys,os
import tree_reader
import node
from copy import *

"""
This is going to only assume nucleotide sequences. It could easily be extended
for amino acids but for the purposes of demonstration, it currently only 
considers nucleotide data
"""


position = {"A":['A'],"C":['C'],"G":['G'],"T":['T'],"-":['A','C','G','T'],"N":['A','C','G','T'],"Y":['C','T'],"R":['A','G'],"W":['A','T'],"M":['A','C'],"B":['G','C','T'],"V":['G','A','C'],"S":['G','C'],"K":['G','T'],"H":['A','T','C']}


def match_tips_and_seqs(tree,seqs):
    lvs = tree.leaves()
    for i in lvs:
        test = False
        for j in seqs:
            if j.name == i.label:
                i.data['seq'] = j.seq.upper()
                test = True
                break
        if test == False:
            print "can't find "+i.label+" in seqs"
            return False
    return True

"""
this calculates the fitch cost for a tree
given a list of sequences that are contained in the tree

"""
def calc_fitch_cost(tree,seqs):
    #first check to make sure that the seqs match the tips
    ret = match_tips_and_seqs(tree,seqs)
    if ret == False:
        return
    sites = len(seqs[0].seq)
    #now calculate the cost
    totalcost = 0
    #print "site cost"
    for i in range(sites):
        cost = 0
        for j in tree.iternodes(order="postorder"):
            if len(j.children) > 0:
                child1 = j.children[0].data['positions']
                child2 = j.children[1].data['positions']
                tset = set(child1).intersection(set(child2))
                #print tset,child1,child2,j.children[0].label,j.children[1].label
                if len(tset) == 0:
                    cost += 1
                    tset = set(child1).union(set(child2))
                j.data['positions'] = tset
            else:
                tset = position[j.data['seq'][i]]
                j.data['positions'] = tset
        #print " ",i," ",cost
        totalcost += cost
    return totalcost


"""
This calculates the Sankoff parsimony with the costs 
so at each node we record the score for each position
"""
def calc_sankoff_dna_cost(tree,seqs,cost_matrix_dict):
    ret = match_tips_and_seqs(tree,seqs)
    if ret == False:
        return
    sites = len(seqs[0].seq)
    totalcost = 0
    for i in range(sites):
        for j in tree.iternodes(order="postorder"):
            if len(j.children) > 0:
                child1 = j.children[0].data['scores']
                child2 = j.children[1].data['scores']
                nddict = {'A':float("inf"),'C':float("inf"),'G':float("inf"),'T':float("inf")}
                for s in nddict:
                    minh = float("inf")
                    for h in child1:
                        tscore = cost_matrix_dict[s][h] + child1[h]
                        if tscore < minh:
                            minh = tscore
                    mink = float("inf")
                    for k in child2:
                        tscore = cost_matrix_dict[s][k] + child2[k]
                        if tscore < mink:
                            mink = tscore
                    nddict[s] = minh+mink
                j.data['scores'] = nddict
            else:
                nddict = {'A':float("inf"),'C':float("inf"),'G':float("inf"),'T':float("inf")}
                tset = position[j.data['seq'][i]]
                for k in tset:
                    nddict[k] = 0
                j.data['scores'] = nddict
        mint = float("inf")
        for j in tree.data['scores']:
            if tree.data['scores'][j] < mint:
                mint = tree.data['scores'][j]
        totalcost += mint
    return totalcost

"""
the general case for discrete
"""
def calc_sankoff_discrete_cost(tree,seqs,cost_matrix_dict):	
    totalcost = 0
    maindict = {}
    for i in cost_matrix_dict:
        maindict[i] = float("inf")
    sites = len(tree.leaves()[0].data['seq'])
    for i in range(sites):
        for j in tree.iternodes(order="postorder"):
            if len(j.children) > 0:
                child1 = j.children[0].data['scores']
                child2 = j.children[1].data['scores']
                nddict = deepcopy(maindict)
                for s in nddict:
                    minh = float("inf")
                    for h in child1:
                        tscore = cost_matrix_dict[s][h] + child1[h]
                        if tscore < minh:
                            minh = tscore
                    mink = float("inf")
                    for k in child2:
                        tscore = cost_matrix_dict[s][k] + child2[k]
                        if tscore < mink:
                            mink = tscore
                    nddict[s] = minh+mink
                j.data['scores'] = nddict
            else:
                nddict = deepcopy(maindict)
                nddict[j.data['seq'][i]] = 0
                j.data['scores'] = nddict
        mint = float("inf")
        for j in tree.data['scores']:
            if tree.data['scores'][j] < mint:
                mint = tree.data['scores'][j]
        totalcost += mint
    return totalcost
