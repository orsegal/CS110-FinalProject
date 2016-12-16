import math
from sequences.sequence import Sequence
import numpy as np
from numpy import *
from scipy.linalg import *
from scipy.stats import *
from scipy.special import *

codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
codons_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
#codons.remove('TAA')
#codons.remove('TAG')
#codons.remove('TGA')


"""
this takes seqstrings for two ALIGNED sequences
if you have sequence objects then it would be seqobject.seq
"""
def prop_diff(seqstring1,seqstring2):
	totallen = len(seqstring1)
	assert totallen == len(seqstring2)
	totaldiff = 0
	for i in range(totallen):
		if seqstring1[i] != seqstring2[i]:
			totaldiff += 1
	return totaldiff/float(totallen)	
		
"""
given the distance calculated by the prop_diff
this will calculate the corrected difference
"""
def jc_diff(propdiff):
	return -0.75*math.log(1.-((4./3.)*propdiff))

def jc_loglike(dist,diff,numsites):
	#can't have below 0
	if dist < 0:
		return 10000000
	first = diff * math.log((1/16.)-((1/16.)*math.exp(-4*(dist/3.))))
	second = (numsites-diff)*math.log((1/16.)+((3/16.)*math.exp(-4*(dist/3.))))
	return -(first + second)

def transition_transversion_diff(seqstring1,seqstring):
	totallen = len(seqstring1)
	assert totallen == len(seqstring2)
	S = 0
	V = 0
	for i in range(totallen):
		if seqstring1[i] == 'T' and seqstring2[i] == 'C':
			S += 1
		elif seqstring1[i] == 'C' and seqstring2[i] == 'T':
			S += 1
		elif seqstring1[i] == 'A' and seqstring2[i] == 'G':
			S += 1
		elif seqstring1[i] == 'G' and seqstring2[i] == 'A':
			S += 1
		elif seqstring1[i] != seqstring2[i]:
			V += 1
	return S/float(totallen),V/float(totallen)

def k80_diff(S,V):
	d = -(1/2.) * math.log(1-(2*S)-V) - (1/4.)* math.log(1-(2*V))
	k = ((2. * math.log(1-(2*S)-V))/(math.log(1-(2*V))))-1
	return d,k

def k80_loglike(x,S,V,numsites):
	dist = x[0]
	kappa = x[1]
	if dist < 0 or kappa < 0:
		return 1000000000
	p0 = (1/4.)+((1/4.)*math.exp(-4*dist/(kappa+2.)))+((1/2.)*math.exp(-2*dist*(kappa+1.)/(kappa+2.)))
	p1 = (1/4.)+((1/4.)*math.exp(-4*dist/(kappa+2.)))-((1/2.)*math.exp(-2*dist*(kappa+1.)/(kappa+2.)))
	p2 = (1/4.)-((1/4.)*math.exp(-4*dist/(kappa+2.)))
	return -((numsites-S-V)*math.log(p0/4.)+S*math.log(p1/4.)+V*math.log(p2/4.))

"""
for codons
"""
"""
this would be the Q as in 


     | 0, more than one change
     | pj, if i and j differ in syn transver 
qij =| Kpj, if i and j differ in syn transit
     | wpj, if i and j differ in nonsyn transv
     | wkpj, if i and j differ in nonsyn transit
ps are the pjs 
K is kapp 
w is omega
"""
def test_transition(a,b):
        ret = False
        if (a == 'A' and b == 'G' ) or (a == 'C' and b == 'T') or (a == 'G' and b == 'A' ) or (a == 'T' and b == 'C'):
                ret = True
        return ret

def test_nonsynonymous(a,b):
        ret = False
        if codons_dict[a] != codons_dict[b]:
                ret = True
        return ret

def calc_simple_goldman_yang_q(ps, K, w):
        bigpi = ones((len(codons),len(codons)))
        basefreq = [1./len(codons)] * len(codons)
        count = 0
        for i in codons:
                count2 = 0
                for j in codons:
                        diff = 0
                        transit = False
                        nonsyn = False
                        for m in range(3):
                                if i[m] != j[m]:
                                        diff += 1
                                        transit = test_transition(i[m],j[m])
                                        nonsyn = test_nonsynonymous(i,j)
     
                        if diff > 1:
                                bigpi[count][count2] = 0
                        else:
                                bigpi[count][count2] = basefreq[count2]
                                if transit:
                                        bigpi[count][count2] *= K
                                if nonsyn:
                                        bigpi[count][count2] *= w
                        count2 += 1
                count += 1
        fill_diagonal(bigpi,0)
        tscale = sum(bigpi)
        t = bigpi/tscale
        count = 0
	#make it so the diags make the rows sum to 0
	for i in t:
		t[count][count] = 0-sum(i)
		count += 1
	t = t/basefreq
	t = transpose(t)
	return t

"""
for binary characters
"""
def calc_asym_q_matrix(a,b):
	bigpi = ones((2,2))
	bigpi[0][0] = -a
	bigpi[0][1] = a
	bigpi[1][0] = b
	bigpi[1][1] = -b
	return bigpi

def calc_q_matrix(rmatrix,basefreq):
	bigpi = ones((len(basefreq),len(basefreq)))
	count = 0
	for i in bigpi:
		count2 = 0
		for j in i:
			if count != count2:
				bigpi[count][count2] = basefreq[count] * basefreq[count2]
			else:
				bigpi[count][count2] = basefreq[count]
			count2 += 1
		count += 1
	t = rmatrix * bigpi
	#print t
	count = 0
	fill_diagonal(t,0)
	tscale = sum(t)
	t= t/tscale
	#make it so the diags make the rows sum to 0
	for i in t:
		t[count][count] = 0-sum(i)
		count += 1
	t = t/basefreq
	t = transpose(t)
	return t

def calc_p_matrix(q,bl):
	return expm2(q*bl)

def calc_jc_p_matrix(time):
	rmatrix = ones((4,4))
	#one way to calculate it
	#rmatrix = rmatrix * 0.33333
	#fill_diagonal(rmatrix,-1)
	#return expm3(rmatrix*time)
	#second way to calculate it
	#i -> i =  (1/4.)+(3/4.)*exp(-4*float(time)/3.)
	#i -> j =  (1/4.)+(1/4.)*exp(-4*float(time)/3.)
	rmatrix = rmatrix * (1/4.)-((1/4.)*exp(-4*float(time)/3.))
	fill_diagonal(rmatrix,(1/4.)+(3/4.)*exp(-4*float(time)/3.))
	return rmatrix
	#third way to calculate it
	#basefreq = ones(4)
	#basefreq = basefreq * 0.25
	#q = calc_nuc_q_matrix(rmatrix,basefreq)
	#return expm3(q*time)

def dist_like_jc2seq(dist,seq1,seq2):
	mat = calc_jc_p_matrix(dist)
	position = {"A":0,"C":1,"G":2,"T":3}
	p = 1
	for i in range(len(seq1)):
		p *= mat[position[seq1[i]]][position[seq2[i]]] 
		p *= 0.25
	return p

def dist_like_jc2seq(dist,seq1,seq2):
	mat = calc_jc_p_matrix(dist)
	position = {"A":0,"C":1,"G":2,"T":3}
	p = 0
	for i in range(len(seq1)):
		p += math.log(mat[position[seq1[i]]][position[seq2[i]]]) 
		p += math.log(0.25)
	return p

def point_gamma(p,a, b):
	return chi2.ppf(p,2*(a))/(2.0*(b))

def get_gamma_cats(alpha,cats,median):
	K=cats
	a,b = alpha,alpha
	factor=a/b*K
	rK = [0]*K
	if median:
		gap05=1.0/(2.0*K)
		for i in range(K):
			rK[i] = point_gamma((i*2.0+1)*gap05,a,b)
		t=0
		for i in range(K):
			t+=rK[i]
		for i in range(K):
			rK[i] *= factor/t
		return rK
	else:
		freqK = [0]*K
		for i in range(K-1):
			freqK[i]=point_gamma((i+1.0)/K, a, b);
		for i in range(K-1):
			freqK[i]=gammainc(a+1,freqK[i]*b)
		rK[0] = freqK[0]*factor;
		rK[K-1] = (1-freqK[K-2])*factor;
		for i in range(1,K-1):
			rK[i] = (freqK[i]-freqK[i-1])*factor;
		return rK
