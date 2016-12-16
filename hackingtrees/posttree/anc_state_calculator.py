import math
from models.model_calc import *
import numpy as np
from numpy import *
from scipy.linalg import *
from trees import tree_parsimony_calculator
from trees.tree_likelihood_calculator import *
from scipy.optimize import fmin_powell
from copy import *
import random


"""
calculatates marginal anc states
"""
def calc_marginal_bin_anc_states(tree,a,b,equal_root):
	q = calc_asym_q_matrix(a,b)
	root0 = 0.5
	root1 = 0.5
	if equal_root == False:
		root0 = (b/(a+b))
		root1 = 1 - root0
	for i in tree.iternodes():
		if i.istip == True:
			continue	
		print i.get_newick_repr()
		if i == tree:
			if equal_root == False:
				st = i.data['probs'][0] * root0
				st += i.data['probs'][1] * root1
				print "0:",root0 * i.data['probs'][0] / st,"1:",root1 * i.data['probs'][1] / st
			else:
				print i.data['probs'] / sum(i.data['probs']) 
				print i.data['probs']*0.5
			i.data['revb'] = [1,1]
			i.data['revb'][0] = root0
			i.data['revb'][1] = root1
		else: 
			parent = i.parent
			for j in parent.children:
				if j != i:
					p = calc_p_matrix(q,i.length)
					i.data['revb'] = [1,1]
					tk0 = ((parent.data['margex'][j][0] * p[0][0] *parent.data['revb'][0] ) + (parent.data['margex'][j][1] * p[1][0] *parent.data['revb'][1]))
					k0 =  (i.data['probs'][0] * tk0)
					tk1 = ((parent.data['margex'][j][0] * p[0][1] *parent.data['revb'][0]) + (parent.data['margex'][j][1] * p[1][1] *parent.data['revb'][1]))
					k1 =  (i.data['probs'][1] * tk1)
					i.data['revb'][0] = tk0
					i.data['revb'][1] = tk1
					print "0:",k0 / (k0+k1) ,"1:",k1 / (k0+k1) 

"""
calculates square change parsimony for continuous characters
"""
def calc_square_change_anc_states(tree,rate):
	df = 0
	nodenum = {}
	count = 0
	for i in tree.iternodes(order="postorder"):
		if i.istip:
			i.data['val'] = float(i.data['cont_values'][0])
			i.data['valse'] = float(i.data['cont_values'][0])
		else:
			nodenum[i] = count
			count += 1
			df += 1
			i.data['val'] = 0.
			i.data['valse'] = 0.
	df -= 1
	#compute the mlest of the root
	fullMcp = zeros((df+1,df+1))
	fullVcp = zeros(df+1)
	count = 0
	for i in tree.iternodes(order="postorder"):
		if i.istip == False:
			nni = nodenum[i]
			for j in i.children:
				tbl = 2./j.length
				fullMcp[nni][nni] += tbl;
				if j.istip:
					fullVcp[nni] += (j.data['val'] * tbl)
				else:
					nnj = nodenum[j]
					fullMcp[nni][nnj] -= tbl;
					fullMcp[nnj][nni] -= tbl;
					fullMcp[nnj][nnj] += tbl;
			count += 1
	b = cho_factor(fullMcp)
	#these are the ML estimates for the ancestral states
	mle = cho_solve(b,fullVcp)
	sos = 0
	for i in tree.iternodes(order="postorder"):
		if i.istip == False:
			i.data['val'] = mle[nodenum[i]]
			print i.data['val']
			i.label = str(mle[nodenum[i]])
			for j in i.children:
				temp = (i.data['val'] - j.data['val'])
				sos += temp*temp / j.length
	print "Square Length: ",sos
	#calcSE
	for i in tree.iternodes(order="postorder"):
		if i.istip == False:
			qpq = fullMcp[nodenum[i]][nodenum[i]]
			tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
			tm = np.delete(tm1,(nodenum[i]),axis=1)
			b = cho_factor(tm)
			sol = cho_solve(b,tm1[:,nodenum[i]])
			tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
			i.data['valse'] = math.sqrt(2*sos/(df*tempse))
	return 0


"""
this calculates ancestral states for discrete (multistate) characters using parsimony
the trees are assumed to have been match seq and tips
"""
def calc_discrete_parsimony_anc_states(tree,seqs,cost_matrix_dict):
	totalcost = tree_parsimony_calculator.calc_sankoff_discrete_cost(tree,seqs,cost_matrix_dict)
	tree.data["parstate"] = None
	minv = min(tree.data["scores"].values())
	for i in tree.data["scores"]:
		if tree.data["scores"][i] == minv:
			tree.data["parstate"] = i
			break
	for i in tree.iternodes(order="preorder"):
		ancstate = i.data["parstate"]
		for j in i.children:
			bestscore = float("inf")
			bestval = None
			for k in j.data["scores"]:
				tscore = cost_matrix_dict[ancstate][k] + j.data["scores"][k]
				if tscore < bestscore:
					bestscore = tscore
					bestval = k
			j.data["parstate"] = bestval
	
	for i in tree.iternodes():
		if i.istip == True:
			continue	
		print i.get_newick_repr()
		print i.data["parstate"]
	
	return totalcost

def sliding_window(x,w):
	return random.uniform(x-w/2.,x+w/2.)

"""
reps : number of replicates
sampleiter : sample every iter
priors : array of prior distributions in the order of the params
params : array of the params with rate matrix, basefreq, branch lengths
nsites : number of sites
tree : a starting tree
freeparams : 0=rateparams,1=branch lengths,2=both
"""

def mcmc_marginal_bin_anc_states(reps,sampleiter,priors,params,nsites,tree,BINARYposition,equal_rates,outfilename):
	#calc the first prior
	curprior = 1.
	ppdfs = []
	for i in range(len(priors)):
		curprior *= priors[i].pdf(params[i])
		ppdfs.append(priors[i].pdf(params[i]))
	#calc the first likelihood
	a = params[0]
	b = params[1]
	curcost = calc_bin_tree_likelihood(tree,a,b,nsites,BINARYposition,equal_rates,False)
	curpost = log(curprior)+(-curcost)
	print "initial prior:",curprior
	print "initial likelihood:",curcost
	print "initial posterior:",curpost
	paramindex = range(len(params))
	sw = 0.1 # sliding window
	outfile = open(outfilename,"w") 
	for j in range(reps):
		#choose which parameter to update
		tp = random.sample(paramindex,1)[0]
		#calc the new proposal
		oldp = params[tp]
		params[tp] = sliding_window(oldp,sw)
		if params[tp] < 0:
			params[tp] *= -1
		
		#calc the new prior
		oldpr = ppdfs[tp]
		ppdfs[tp] = priors[tp].pdf(params[tp])
		prior = 1.
#		for i in range(len(ppdfs)):
#			prior *= ppdfs[i]
		#calc the new likelihood
		a = params[0]
		b = params[1]
		newcost = calc_bin_tree_likelihood(tree,a,b,nsites,BINARYposition,equal_rates,False)
		if isinf(newcost) or isnan(newcost) or prior == 0:
			params[tp] = oldp
			ppdfs[tp] = oldpr
			continue
		newpost = (-newcost) + log(prior)
		if exp(newpost-curpost) >= 1:
			#keep
			curpost = newpost
			curcost = newcost
			curprior = prior
		else:
			r = random.random()
			if r <= exp(newpost-curpost):
				#keep
				curpost = newpost
				curcost = newcost
				curprior = prior
			else:
				params[tp] = oldp
				ppdfs[tp] = oldpr
		#calc whether to keep the value
		if j % sampleiter == 0:
			print j,curcost
			outfile.write(str(j)+"\t"+str(curcost)+"\t"+str(params[0])+"\t"+str(params[1])+"\n")
	outfile.close()

