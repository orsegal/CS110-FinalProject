import math
from models.model_calc import *
import numpy as np
from numpy import *
from scipy.linalg import *
from scipy.optimize import fmin_powell

#move to datatype
DNAposition = {"A":[0],"C":[1],"G":[2],"T":[3],"-":[0,1,2,3],"-":[0,1,2,3],"N":[0,1,2,3],"Y":[1,3],"R":[0,2],"W":[0,3],"M":[0,1],"B":[2,1,3],"V":[2,0,1],"S":[2,1],"K":[2,3],"H":[0,3,1]}
BINARYposition = {"0":[0],"1":[1],"?":[0,1]}
CODONposition = {'CTT': [13], 'ATG': [[32]], 'ACA': [35], 'ACG': [36], 'ATC': [30], 'ATA': [31], 'AGG': [44], 'CCT': [17], 'AGC': [42], 'AGA': [43], 'ATT': [29], 'CTG': [16], 'CTA': [15], 'ACT': [33], 'CCG': [20], 'AGT': [41], 'CCA': [19], 'CCC': [18], 'TAT': [8], 'GGT': [57], 'CGA': [27], 'CGC': [26], 'CGG': [28], 'GGG': [60], 'GGA': [59], 'GGC': [58], 'TAC': [9], 'CGT': [25], 'GTA': [47], 'GTC': [46], 'GTG': [48], 'GAG': [56], 'GTT': [45], 'GAC': [54], 'GAA': [55], 'AAG': [40], 'AAA': [39], 'AAC': [38], 'CTC': [14], 'CAT': [21], 'AAT': [37], 'CAC': [22], 'CAA': [23], 'CAG': [24], 'TGT': [10], 'TCT': [4], 'GAT': [53], 'TTT': [0], 'TGC': [11], 'TGG': [12], 'TTC': [1], 'TCG': [7], 'TTA': [2], 'TTG': [3], 'TCC': [5], 'ACC': [34], 'TCA': [6], 'GCA': [51], 'GCC': [50], 'GCG': [52], 'GCT': [49],'---':range(0,61)}
LARGE = 10000000

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

def match_tips_and_codons(tree,seqs):
	lvs = tree.leaves()
	for i in lvs:
		test = False
		for j in seqs:
			if j.name == i.label:
                                a = j.seq.upper()
				i.data['codon'] = [a[k:k+3] for k in range(0, len(a), 3)]
				test = True
				break
		if test == False:
			print "can't find "+i.label+" in seqs"
			return False

"""
for binary characters, simplified model
"""
def calc_bin_tree_likelihood(tree,a,b,sites,datatypeposition,equil_root = False,verbose=False):
	q = calc_asym_q_matrix(a,b)
	ps = {}#key is length, value is pmatrix
	loglike = 0
	for s in range(sites):
		tempretlike = 0
		for c in tree.iternodes(order = "postorder"):
			if len(c.children) > 0: #this is an internal node
				c.data['probs'] = ones(2)
				c.data['margex'] = {}
				for j in c.children:
					c.data['margex'][j] = ones(2)
				for j in range(2):
					for i in c.children:
						if i.length not in ps:
							p = calc_p_matrix(q,i.length)
							ps[i.length] = p
						p = ps[i.length]
						templike = 0
						for k in range(2):
							templike += (p[j][k] * i.data['probs'][k])
						c.data['probs'][j] *= templike
						c.data['margex'][i][j] = templike
				#print c.get_newick_repr(),c.data['probs']
			else: #this is a tip
				c.data['probs'] = zeros(2)
				positions = datatypeposition[c.data['seq'][s]]
				for j in positions:
					c.data['probs'][j] = 1
		#the tempretlike will be the sum(root prob * basefreq) for all bases
		tempretlike = 0
		root0 = 0.5
		root1 = 0.5
		if equil_root == False:
			root0 = (b/(a+b))
			root1 = 1-root0
		tempretlike  += tree.data['probs'][0]*root0#0.5#add root*basefreq[i]
		tempretlike  += tree.data['probs'][1]*root1#0.5#add root*basefreq[i]

		if verbose:
			print "site",s,"log(L):",log(tempretlike),"like:",tempretlike
		loglike -= log(tempretlike)
	return loglike

def calc_tree_likelihood(tree,rmatrix,basefreq,sites,datatypeposition,verbose=False):
	q = calc_q_matrix(rmatrix,basefreq)
	ps = {}#key is length, value is pmatrix
	loglike = 0
	for s in range(sites):
		tempretlike = 0
		for c in tree.iternodes(order = "postorder"):
			if len(c.children) > 0: #this is an internal node
				c.data['probs'] = ones(len(basefreq))
				for j in range(len(basefreq)):
					for i in c.children:
						if i.length not in ps:
							p = calc_p_matrix(q,i.length)
							ps[i.length] = p
						p = ps[i.length]
						templike = 0
						for k in range(len(basefreq)):
							templike += (p[j][k] * i.data['probs'][k])
						c.data['probs'][j] *= templike
				#print c.get_newick_repr(),c.data['probs']
			else: #this is a tip
				c.data['probs'] = zeros(len(basefreq))
				positions = datatypeposition[c.data['seq'][s]]
				for j in positions:
					c.data['probs'][j] = 1
		#the tempretlike will be the sum(root prob * basefreq) for all bases
		tempretlike = 0
		for i in range(len(basefreq)):
			tempretlike  += tree.data['probs'][i]*basefreq[i]
		if verbose:
			print "site",s,"log(L):",log(tempretlike),"like:",tempretlike
		loglike -= log(tempretlike)
	return loglike

def calc_tree_likelihood_M0_selection(tree,codonfreq,K,w,sites,verbose=False):
	q = calc_simple_goldman_yang_q(codonfreq,K,w)
	ps = {}#key is length, value is pmatrix
	loglike = 0
	for s in range(sites):
		tempretlike = 0
		for c in tree.iternodes(order = "postorder"):
			if len(c.children) > 0: #this is an internal node
				c.data['probs'] = ones(len(codonfreq))
				for j in range(len(codonfreq)):
					for i in c.children:
						if i.length not in ps:
							p = calc_p_matrix(q,i.length)
							ps[i.length] = p
						p = ps[i.length]
						templike = 0
						for k in range(len(codonfreq)):
							templike += (p[j][k] * i.data['probs'][k])
						c.data['probs'][j] *= templike
				#print c.get_newick_repr(),c.data['probs']
			else: #this is a tip
				c.data['probs'] = zeros(len(codonfreq))
				positions = CODONposition[c.data['codon'][s]]
				for j in positions:
					c.data['probs'][j] = 1
		#the tempretlike will be the sum(root prob * basefreq) for all bases
		tempretlike = 0
		for i in range(len(codonfreq)):
			tempretlike  += tree.data['probs'][i]*codonfreq[i]
		if verbose:
			print "site",s,"log(L):",log(tempretlike),"like:",tempretlike
		loglike -= log(tempretlike)
	return loglike

def calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats,verbose=False):
	q = calc_q_matrix(rmatrix,basefreq)
	ps = {}#key is length, value is pmatrix
	cs = get_gamma_cats(alpha,cats,False)
	loglike = 0
	for s in range(sites):
		tempretlikes = 0
		for x in cs:
			tempretlike = 0
			for c in tree.iternodes(order = "postorder"):
				if len(c.children) > 0: #this is an internal node
					c.data['probs'] = ones(4)
					for j in range(4):
						for i in c.children:
							if i.length not in ps:
								p = calc_p_matrix(q,i.length*x)
								ps[i.length*x] = p
							p = ps[i.length*x]
							templike = 0
							for k in range(4):
								templike += p[j][k] * i.data['probs'][k]
							c.data['probs'][j] *= templike
				else: #this is a tip
					c.data['probs'] = zeros(4)
					#could add here for ambiguity
					positions = DNAposition[c.data['seq'][s]]
					for j in positions:
						c.data['probs'][j] = 1
			#at this point the tempretlike will be the sum(root prob * basefreq) for all bases
			tempretlike = 0
			for i in range(4):
				tempretlike  += tree.data['probs'][i]*basefreq[i]
			tempretlikes += (tempretlike*(1/float(len(cs))))
			if verbose: 
				print "site",s,"log(L):",log(tempretlike),log(tempretlikes)
		loglike -= log(tempretlikes)
		if verbose:
			print log(tempretlikes)
	return loglike


def set_binary_rmatrix_basefreq(params):
	rmatrix = ones((2,2))
	rmatrix[0][1] = params[0]
	rmatrix[1][0] = params[1]
	fill_diagonal(rmatrix,0)
	basecomp = ones(2)	
	basecomp[0] = params[2]
	basecomp[1] = 1-params[2]
	#sc = sum(basecomp)
	#basecomp /= sc
	return rmatrix,basecomp

"""
this scales so that the sixth rate (G->T) is 1 
the matrix will be
  A C G T
A   1 2 3
C	 4 5
G 
T 

the three base composition rates will be parameters
6, 7, 8. The last base composition (T) will be 1- sum(6,7,8)
and the total will be scaled to 1
"""
def set_gtr_rmatrix_basefreq(params):
	rmatrix = ones((4,4))
	rmatrix[0][1] = params[0]
	rmatrix[1][0] = params[0]
	rmatrix[0][2] = params[1]
	rmatrix[2][0] = params[1]
	rmatrix[0][3] = params[2]
	rmatrix[3][0] = params[2]
	rmatrix[1][2] = params[3]
	rmatrix[2][1] = params[3]
	rmatrix[1][3] = params[4]
	rmatrix[3][1] = params[4]
	rmatrix[2][3] = 1.
	rmatrix[3][2] = 1.
	fill_diagonal(rmatrix,0)
	basecomp = ones(4)	
	basecomp[0] = params[5]
	basecomp[1] = params[6]
	basecomp[2] = params[7]
	basecomp[3] = 1-sum(params[5:8])
	sc = sum(basecomp)
	basecomp /= sc
	return rmatrix,basecomp


"""
params contain the rmatrix and basefreq params in the form
  A C G T
A   1 2 3
C	 4 5
G 
T 
the three base composition rates will be parameters
6, 7, 8. The last base composition (T) will be 1- sum(6,7,8)
and the total will be scaled to 1
"""
def calc_params_rateparams(params,tree,sites,alpha=None,cats=None):
	for i in params: 
		if i < 0:
			return LARGE
	if sum(params[5:8]) > 1:
		return LARGE
	rmatrix,basefreq = set_gtr_rmatrix_basefreq(params)
	like = -1
	if alpha != None:
		like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
	else:
		like = calc_nuc_tree_likelihood(tree,rmatrix,basefreq,sites)
	print like
	if like < 0 or isnan(like):
		return LARGE
	return like

"""
params contains parameter values for everything but the tree topology and gamma
"""
def calc_params_allparams(params,tree,sites,alpha,cats):
	for i in params: 
		if i < 0:
			return LARGE
	if sum(params[5:8]) > 1:
		return LARGE
	rateparams = params[0:8]
	count = 8
	for i in tree.iternodes():
		if i != tree:
			if params[count] < 0:
				return LARGE
			i.length = params[count]
			count += 1
	rmatrix,basefreq = set_gtr_rmatrix_basefreq(rateparams)
	like = -1
	if alpha != None:
		like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
	else:
		like = calc_nuc_tree_likelihood(tree,rmatrix,basefreq,sites)
	print like
	if like < 0 or isnan(like):
		return LARGE
	return like

"""
params contains just the alpha shape parameter
"""
def calc_params_gamma(params,tree,sites,rmatrix,basefreq,cats):
	alpha = params[0]
	if alpha < 0:
		return LARGE
	like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
	if like < 0 or isnan(like):
		return LARGE
	print like
	return like

"""
params contains the branch lengths
"""
def calc_params_treebl(params,tree,sites,rmatrix,basefreq,alpha=None,cats=None):
	count = 0
	for i in tree.iternodes():
		if i != tree:
			if params[count] < 0:
				return LARGE
			i.length = params[count]
			count += 1
	like = -1
	if alpha != None:
		like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
	else:
		like = calc_tree_likelihood(tree,rmatrix,basefreq,sites,DNAposition)
	#print tree.get_newick_repr(True)
	print "like:",like
	if like < 0 or isnan(like):
		return LARGE
	return like
		

"""
params contains the internal node dates
"""
def calc_params_treebl_dates(params,tree,sites,rmatrix,basefreq,alpha=None,cats=None):
	count = 0
	for i in tree.iternodes(order="postorder"):
		if len(i.children) > 0: #not internal
			if params[count] < 0:
				return LARGE
			i.height = params[count]
			for j in i.children:
				if len(j.children) == 0:
					j.height = 0
				j.length = i.height - j.height
				if j.length < 0:
					return LARGE
			count += 1
	like = -1
	if alpha != None:
		like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
	else:
		like = calc_tree_likelihood(tree,rmatrix,basefreq,sites,DNAposition,False)
	#print tree.get_newick_repr(True)
	print "clock:",like
	if like < 0 or isnan(like):
		return LARGE
	return like

"""
params contains just the two rates (alpha and beta) for the asym binary likelihood
"""
def calc_params_bin_like(params,tree,sites,datatypeposition,equal_root):
	a = params[0]
	b = params[1]
	if a < 0 or b < 0:
		return LARGE
	like = calc_bin_tree_likelihood(tree,a,b,sites,datatypeposition,equal_root, False)
	#print tree.get_newick_repr(True)
	print like
	if like < 0 or isnan(like):
		return LARGE
	return like

def optimize_bin_like(a,tree,sites,datatypeposition,equal_root):
	assert len(a) == 2
	res =  fmin_powell(calc_params_bin_like,a,args=(tree,sites,datatypeposition,equal_root),full_output=True)
	return res

"""
optimizing params
"""

"""
this should really be done with derivatives, but 
for demonstration purposes we will just use powell 
derivative free things
"""
def optimize_brlen(tree,sites,rmatrix,basefreq,alpha=None,cats=None):
	blstart = []
	for i in tree.iternodes():
		if i != tree:
			blstart.append(i.length)
	res = fmin_powell(calc_params_treebl,blstart,args=(tree,sites,rmatrix,basefreq,alpha,cats),full_output=True)
	return res
	
def optimize_brlen_clock(tree,sites,rmatrix,basefreq,alpha=None,cats=None):
	hts = []
	for i in tree.iternodes(order="postorder"):
		i.set_height()
		if len(i.children) > 0:
			hts.append(i.height)
	res = fmin_powell(calc_params_treebl_dates,hts,args=(tree,sites,rmatrix,basefreq,alpha,cats),full_output=True)
	return res	

def optimize_rateparams(a,tree,sites,alpha=None,cats=None):
	assert len(a) == 8
	res =  fmin_powell(calc_params_rateparams,a,args=(tree,sites,alpha,cats),full_output=True)
	return res

"""
#maximize rate params
res = fmin_powell(optimize_rateparams,a,args=(tree,sites),full_output=True)
print res[0][0:6]
print res[0][5:]
print sum(res[0][5:])
rmatrix,basefreq = set_gtr_rmatrix_basefreq(res[0])
print rmatrix
print basefreq
print calc_nuc_q_matrix(rmatrix,basefreq)
"""


def optimize_allparams(a,tree,sites,alpha=None,cats=None):
	assert len(a) == 8
	blstart = []
	for i in tree.iternodes():
		if i != tree:
			blstart.append(i.length)
	a = a + blstart
	res =  fmin_powell(calc_params_allparams,a,args=(tree,sites,alpha,cats),full_output=True)
	return res

def optimize_gammaparam(alphastart,tree,sites,rmatrix,basefreq,cats):
	res =  fmin_powell(calc_params_gamma,[alphastart],args=(tree,sites,rmatrix,basefreq,cats),full_output=True)
	return res
