import sys,os
from sequences.sequence import Sequence
import copy
import random
import png

def sample_with_repl(pop, nitems):
    n = len(pop)
    result = [None] * nitems
    for i in xrange(nitems):
        j = int(random.random() * n)
        result[i] = pop[j]
    return result

def sample_without_repl(pop, nitems):
    return random.sample(pop,nitems)

def bootstrap_aln(seqs):
	bsseqs = copy.deepcopy(seqs)
	res = sample_with_repl(range(0,len(seqs[0].seq)),len(seqs[0].seq))
	#print res
	for i in range(len(bsseqs)):
		s = ""
		for j in res:
			s += seqs[i].seq[j]
		bsseqs[i].seq = s
	return bsseqs

"""
the nsites here should be the number of sites we want in the jackknife
alignment
"""
def jackknife_aln(seqs,nsites):
	bsseqs = copy.deepcopy(seqs)
	res = sample_without_repl(range(0,len(seqs[0].seq)),nsites)
	for i in range(len(bsseqs)):
		s = ""
		for j in res:
			s += seqs[i].seq[j]
		bsseqs[i].seq = s
	return bsseqs


def write_aln_to_fig(seqs,outfilename,sizeincrease=10,sizexincrease=5,sizeyincrease=20,breakup=1):
    #all amb are grey, a red, t blue, c yellow, g green
    position = {"A":1,"C":2,"G":3,"T":4,"-":0,"N":0,"Y":5,"R":5,"W":5,"M":5,"B":5,"V":5,"S":5,"K":5,"H":5}
    order = {}
    count = 0
    for i in seqs:
        order[i.name] = count
        count += 1
    array = ['']*len(order)
    palette=[(0xff,0xff,0xff), (0xff,0x00,0x00), (0xff,0xff,0x00), (0x00,0xff,0x00), (0x00,0xff,0xff), (0x80,0x80,0x80)]
    for i in seqs:
        #duplicateit
        std = ""
        for j in i.seq:
            for k in range(sizexincrease):
                std += str(position[j.upper()])
        array[order[i.name]] = "".join(std)
    start = 0
    end = 0
    for i in range(breakup):
        nend = end
        start = nend
        end = min(start + (len(array)/breakup),len(array))+1
    #size increase in blocks
        order2 = []
        for j in array[start:end]:
            for k in range(sizeyincrease):
                order2.append(j)
        tarray = map(lambda x: map(int, x), order2)
        outfile = None
        if breakup == 1:
            outfile = open(outfilename+".png","wb")
        else:
            outfile = open(outfilename+"."+str(i)+".png","wb")
        w = png.Writer(len(tarray[0]),len(tarray),palette=palette,bitdepth=4)
        w.write(outfile,tarray)
        outfile.close()

#seqsarray should be an array of seqs
def write_concat_genes_to_phylip(seqsarray,outfilename):
    specieslist = []
    genespecies = []
    genespeciesd = []
    count = 0
    nsites = 0
    genelens = []
    for i in seqsarray:
        nsites += len(i[0].seq)
        genelens.append(len(i[0].seq))
        gsp = []
        gspd = {}
        for j in i:
            specieslist.append(j.name)
            if j.name in gsp:
                print "error",j.name,"duplicated in alignment",count
            gsp.append(j.name)
            gspd[j.name]= j
        count += 1
        genespecies.append(gsp)
        genespeciesd.append(gspd)
    specieslist = set(specieslist)
    
    outfile = open(outfilename,"w")
    seqstring = {}
    for j in specieslist:
        seqstring[j] = j+"\t"

    outfile.write(str(len(specieslist))+"\t"+str(nsites)+"\n")
    for i in range(len(seqsarray)):
        for j in specieslist:
            if j in genespecies[i]:
                seqstring[j] += genespeciesd[i][j].seq
            else:
                seqstring[j] += ("-"*genelens[i])
    
    for i in seqstring:
        outfile.write(seqstring[i]+"\n")

    outfile.close()

    outfile = open(outfilename+".models","w")
    start = 0
    end = 1
    for i in range(len(seqsarray)):
        start += 1
        end = start + genelens[i] - 1
        outfile.write("DNA, gene"+str(i)+" = "+str(start)+"-"+str(end)+"\n")
        start = end
    outfile.close()
