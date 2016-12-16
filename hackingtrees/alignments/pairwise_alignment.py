import sys,os
from numpy import *
from sequences.sequence import Sequence
from time import clock,time
import random

"""
this reads matrices as they are formatted here
ftp://ftp.ncbi.nih.gov/blast/matrices/
"""
def read_scoring_matrix(filename):
	infile = open(filename,"r")
	sc_mat = {}
	order = []
	first = True
	for i in infile:
		if i[0] == "#":
			continue
		else:
			spls = i.strip().split()
			if first == True:
				first = False
				order = spls
				for j in order:
					sc_mat[j] = {}
				continue
			for j in range(len(order)):
				sc_mat[spls[0]][order[j]] = float(spls[j+1]) #changed from int to float
	infile.close()
	return sc_mat

def get_simple_score(aln1,aln2):
	#just some simple calculation of match and mismatch
	match = 0
	mismatch = 0
	print aln1
	print aln2
	for i in range(len(aln1)):
		if aln1[i] == aln2[i]:
			match += 1
		else:
			mismatch += 1
	return mismatch,match

def nwg_pairwise(seq1,seq2,gap_penalty,scoringmatrix,verbose = False):
	seq1 = seq1.upper()
	seq2 = seq2.upper()
	F = []
	for i in range(len(seq2)+1):
		b = []
		for j in range(len(seq1)+1):
			b.append(0)
		F.append(b)
	#print_latex_matrix(seq1,seq2,F)
	d = gap_penalty
	for i in range(len(seq1)+1):
		F[0][i] = d*i
	for j in range(len(seq2)+1):
		F[j][0] = d*j
	#print_latex_matrix(seq1,seq2,F)
	for j in range(1,len(seq2)+1):
		for i in range(1,len(seq1)+1):
			match = F[j-1][i-1] + scoringmatrix[seq1[i-1]][seq2[j-1]]
			delete = F[j-1][i] + d
			insert = F[j][i-1] + d
			F[j][i] = max(match,delete,insert)
 #   print_latex_matrix(seq1,seq2,F)
	aln1 = ""
	aln2 = ""
	i = len(seq1)
	j = len(seq2)
	if verbose:
		print "unalign: "+seq1
		print "unalign: "+seq2
	while( i > 0 or j > 0):
		score = F[j][i]
		scorediag = F[j-1][i-1]
		scoreup = F[j-1][i]
		scoreleft = F[j][i-1]
		if (i > 0 and j > 0 and score == scorediag + scoringmatrix[seq1[i-1]][seq2[j-1]]):
			aln1 += seq1[i-1]
			aln2 += seq2[j-1]
			i = i - 1
			j = j - 1
		elif ( i > 0 and score == scoreleft + d):
			aln1 += seq1[i-1]
			aln2 += "-"
			i = i - 1
		elif (j > 0 and score == scoreup + d):
			aln1 += "-"
			aln2 += seq2[j-1]
			j = j - 1
	#print F
	aln1 = aln1[::-1]
	aln2 = aln2[::-1]
	if verbose:
		print "aligned: "+aln1
		print "aligned: "+aln2
	#just some simple calculation of match and mismatch
	match = 0
	mismatch = 0
	score=0;
	for i in xrange(0,len(aln1)):
		if aln1[i]!='-' and aln2[i]!='-':
			score=score+scoringmatrix[aln1[i]][aln2[i]]
		else:
			score=score+gap_penalty

	for i in range(len(aln1)):
		if aln1[i] == aln2[i]:
			match += 1
		else:
			mismatch += 1
	if verbose:
		print aln1,aln2,score
	return aln1,aln2,score

def sw_pairwise(seq1,seq2,gap_penalty,scoringmatrix,verbose=False):
	seq1 = seq1.upper()
	seq2 = seq2.upper()
	F = zeros((len(seq2)+1,len(seq1)+1))
	d = gap_penalty
	for i in range(len(seq1)+1):
		F[0][i] = 0
	for j in range(len(seq2)+1):
		F[j][0] = 0
#	print F
#	print_latex_matrix(seq1,seq2,F)
	for j in range(1,len(seq2)+1):
		for i in range(1,len(seq1)+1):
			match = F[j-1][i-1] + scoringmatrix[seq1[i-1]][seq2[j-1]]
			delete = F[j-1][i] + d
			insert = F[j][i-1] + d
			F[j][i] = max(match,delete,insert,0)
	if verbose:
		print F
	# fancy way to get the j and i with the highest score
	bestj,besti = unravel_index(F.argmax(), F.shape)
	if verbose:
		print "best i:",besti,"best j:",bestj
 #   print_latex_matrix(seq1,seq2,F)
	aln1 = ""
	aln2 = ""
	i = besti
	j = bestj
	if verbose:
		print "unalign: "+seq1
		print "unalign: "+seq2
	while( i > 0 and j > 0):
		score = F[j][i]
		if score == 0:
			break
		scorediag = F[j-1][i-1]
		scoreup = F[j-1][i]
		scoreleft = F[j][i-1]
		if (score == scorediag + scoringmatrix[seq1[i-1]][seq2[j-1]]):
			aln1 += seq1[i-1]
			aln2 += seq2[j-1]
			i = i - 1
			j = j - 1
		elif (score == scoreleft + d):
			aln1 += seq1[i-1]
			aln2 += "-"
			i = i - 1
		elif (score == scoreup + d):
			aln1 += "-"
			aln2 += seq2[j-1]
			j = j - 1
	if verbose:
		print "aligned: "+aln1[::-1]
		print "aligned: "+aln2[::-1]
	#just some simple calculation of match and mismatch
	match = 0
	mismatch = 0
	for i in range(len(aln1)):
		if aln1[i] == aln2[i]:
			match += 1
		else:
			mismatch += 1
	return mismatch,match


def print_latex_matrix(seq1,seq2,mat):
	print "\\begin{matrix}"
	print " & & "+" & ".join(seq1) + " \\\\"
	count = 0
	for i in mat:
		ar = [str(st) for st in i]
		if count == 0:
			print " & " +" & ".join(ar)+" \\\\"
		else:
			print seq2[count-1]+" & " +" & ".join(ar)+" \\\\"
		count += 1
	print "\\end{matrix}"

