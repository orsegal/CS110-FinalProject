import pairwise_alignment
from sequences.sequence import Sequence
from numpy import *
from random import *

#calculates the score of random seqs
def srand(s1,s2,penalty,matrix):
	s_rand_1,s_rand_2="",""
	#amino acid
	if len(matrix) == 24:
		for i in range(len(s1)):
			s_rand_1 += choice('ARNDCQEGHILKMFPSTWYV')
		for i in range(len(s2)):
			s_rand_2 += choice('ARNDCQEGHILKMFPSTWYV')
	else:
		for i in range(len(s1)):
			s_rand_1 += choice('ACGT')
		for i in range(len(s2)):
			s_rand_2 += choice('ACGT')	
	al1, al2, sc = pairwise_alignment.nwg_pairwise(s_rand_1,s_rand_2,penalty,matrix)
	return sc
	
def matrixmin(A):#return the coordinate of min in a 2-D matrix
	minnum=10000;
	for i in range(len(A)):
		for j in range(i+1,len(A)):
			if A[i][j]<=minnum and A[i][j]!=0:
				min_x=i
				min_y=j
				minnum=A[i][j]
	return min_x,min_y

#string search in list
def str_find(list_name,substr):
	found=-1
	for i in range(len(list_name)):
		if list_name[i].find(substr)==0:
			found=i
			break
		else: 
			found=-1
	return found

#function to build guide tree
def guide_tree(distm,guide,pair):
	min_x,min_y=matrixmin(distm)
	marker=[]
	marker.append(str(min_x))
	marker.append(str(min_y))
	newline=[marker,guide[min_x][1]+guide[min_x][2],guide[min_y][1]+guide[min_y][2],str(len(guide))+':'+str(str(min_x))+' '+(str(min_y))]
	guide.append(newline)

	dim=len(distm)
	distm=vstack((distm,zeros((1,dim))))
	distm=hstack((distm,zeros((dim+1,1))))

	for k in range(len(distm)-1):#original index is i
		tmpdis=0
		for t in range(len(guide[k][1])):
			i=int(guide[k][1][t])
			for j in range(len(newline[1])):
				index=int(newline[1][j])
				tmpdis=tmpdis+pair[index][i]
				tmpdis=tmpdis+pair[i][index]
			for j in range(len(newline[2])):
				index=int(newline[2][j])
				tmpdis=tmpdis+pair[index][i]
				tmpdis=tmpdis+pair[i][index]
		for t in range(len(guide[k][2])):
			i=int(guide[k][2][t])
			for j in range(len(newline[1])):
				index=int(newline[1][j])
				tmpdis=tmpdis+pair[index][i]
				tmpdis=tmpdis+pair[i][index]
			for j in range(len(newline[2])):
				index=int(newline[2][j])
				tmpdis=tmpdis+pair[index][i]
				tmpdis=tmpdis+pair[i][index]
		tmpdis=float(tmpdis)/(len(newline[1])*len(newline[2]))
		distm[k][len(distm)-1]=tmpdis

	for k1 in range(len(guide)):#original index is i
		for k2 in range(len(pair),len(guide)):#original index is i
			if str_find(guide[k2][0],str(k1))!=-1:		
				distm[:,k1]=0
				distm[k1,:]=0
	return distm,guide

def str_insert(original, new, pos):#insert a character into a string
	return original[:pos] + new + original[pos:]

"""
This is an implementation of the Feng Doolittle algorithm
s = seq strings in a list
"""
def feng_doolittle(seqs,penalty,matrix,verbose = False):
	s  = []
	for i in seqs:
		s.append(i.seq)
	#generate score matrix
	print "calculating distance matrix"
	score_matrix=zeros((len(s),len(s)));
	for i in range(0,len(s)):
		for j in range(i,len(s)):
			al1, al2, tscore =pairwise_alignment.nwg_pairwise(s[i],s[j],penalty,matrix)
			score_matrix[i][j]= tscore
	#Convert pairwise alignment score into distance
	distm=zeros((len(s),len(s)));
	for i in xrange(0,len(s)):
		for j in xrange(i+1,len(s)):
			s_obs=score_matrix[i][j]
			s_rand=float(srand(s[i],s[j],penalty,matrix))
			s_max=(float(score_matrix[i][i]+score_matrix[j][j]))/2.
			#changed this to abs so as to not be negative
			distm[i][j]=-log(abs((s_obs-s_rand)/(s_max-s_rand)))
			whereAreNaNs = isnan(distm)
			if sum(whereAreNaNs) > 0:
				print "found nans, could be a weird result"
			distm[whereAreNaNs] = 0
			whereAreInfs = isinf(distm)
			if sum(whereAreInfs) > 0:
				print "found infs, could be a weird result"
			distm[whereAreInfs] = 0
	#Initialize guide tree
	guide=[]
	for i in range(0,len(s)):
		guide.append([[str(i)],str(i),'',str(i)+':'+str(i)])
	#build guide tree
	pair=distm
	while sum(distm)!=0:
		distm,guide=guide_tree(distm,guide,pair)
	#for i in range(len(pair),len(guide)):
	#	print guide[i][3], " MSA1=",guide[i][1]," MSA2=",guide[i][2];
	
	print "aligning"
	#align MSAs
	for i in range(len(pair),len(guide)):
		if verbose:
			print "####################"
		#find highest alignemtn pair
		ali1,ali2,AlignMax=pairwise_alignment.nwg_pairwise(s[int(guide[i][1][0])],s[int(guide[i][2][0])],penalty,matrix);
		for msa1 in range(len(guide[i][1])):
			seq1_index=int(guide[i][1][msa1])
			seq1=s[seq1_index]
			for msa2 in range(len(guide[i][2])):
				seq2_index=int(guide[i][2][msa2])
				seq2=s[seq2_index]
				align1,align2,alignscore=pairwise_alignment.nwg_pairwise(seq1,seq2,penalty,matrix)
				if alignscore>=AlignMax:
					template_num_1=seq1_index #seq index template 1
					template_num_2=seq2_index #seq index template 2
					template1=align1
					template2=align2
					AlignMax=alignscore
		s[template_num_1]=template1
		s[template_num_2]=template2
		
		replaceletter = "N"
		#for amino acid
		if len(matrix) == 24:
			replaceletter = "X"
		if verbose:
			print template1.replace(replaceletter,'-')
		#insert gap into all other sequences in same MSA
		for msa1 in xrange(0,len(guide[i][1])):#all sequences in same MSA group
			if int(guide[i][1][msa1])!=template_num_1:#but not itself
				replacing_seq=s[int(guide[i][1][msa1])]
				for j in xrange(0,len(template1)):
					if template1[j]=='-':
						#insert a gap in replacing_seq
						replacing_seq=str_insert(replacing_seq,'-',j)
				replacing_seq=replacing_seq.replace('-',replaceletter)			
				s[int(guide[i][1][msa1])]=replacing_seq
				if verbose:
					print(replacing_seq.replace(replaceletter,'-'))
		if verbose:
			print template2.replace(replaceletter,'-')
		for msa2 in xrange(0,len(guide[i][2])):#all sequences in same MSA group
			if int(guide[i][2][msa2])!=template_num_2:#but not itself
				replacing_seq=s[int(guide[i][2][msa2])]
				for j in xrange(0,len(template2)):
					if template2[j]=='-':
						#insert a gap in replacing_seq
						replacing_seq=str_insert(replacing_seq,'-',j)
				replacing_seq=replacing_seq.replace('-',replaceletter)			
				s[int(guide[i][2][msa2])]=replacing_seq
				if verbose:
					print(replacing_seq.replace(replaceletter,'-'))
		#replace - by X
		#only need to do this for templates 
		s[template_num_1]=s[template_num_1].replace('-',replaceletter)
		s[template_num_2]=s[template_num_2].replace('-',replaceletter)										
	
	print "finished calculation"
	for i in range(len(s)):
		seqs[i].seq = s[i].replace(replaceletter,"-")
	#calculate SP score
	spscore=0;
	for i in range(len(s[0])):#i-th bit
		for j in range(len(s)-1):#sequence j
			for k in range(j+1,len(s)):#all sequences after it
				spscore=spscore+matrix[s[j][i]][s[k][i]]
	if verbose:
		print "####################\n",spscore
	return seqs
