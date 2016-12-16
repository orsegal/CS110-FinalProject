import sys,os
from sequences import seq_reader,sequence
from trees import tree_reader,tree_likelihood_calculator
from trees.tree_likelihood_calculator import *
import copy 
from alignments import aln_utils
from alignments import aln_reader
from models import model_calc

def calc_sel(params,tree,codonfreqs,sites):
    K = params[0]
    w = params[1]
    if K < 0.0001 or w < 0.0001 or K > 999.0 or w > 999.0:
        return LARGE
    like = calc_tree_likelihood_M0_selection(tree,codonfreqs,K,w,sites,verbose=False)
    #print tree.get_newick_repr(True)
    print like
    if like < 0 or isnan(like):
        return LARGE
    return like


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" tree fastaaln"
        sys.exit(0)
    infile = open(sys.argv[1],"r")
    trees = infile.readline()
    print trees
    tree = tree_reader.read_tree_string(trees)
    infile.close()
    seqs = aln_reader.read_phylip_file(sys.argv[2])
    sites = len(seqs[0].seq)/3
    match_tips_and_codons(tree,seqs)
    
    Q = model_calc.calc_simple_goldman_yang_q(1,1,1)
    print Q
    print model_calc.calc_p_matrix(Q,0.1)
    codonfreqs = [1./61] * 61
    w = 1
    K = 1
    print calc_tree_likelihood_M0_selection(tree,codonfreqs,K,w,sites,verbose=True)

    a = [K,w]
    res =  fmin_powell(calc_sel,a,args=(tree,codonfreqs,sites),full_output=True)
    print res
