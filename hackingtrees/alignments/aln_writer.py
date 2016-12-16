from sequence import Sequence

"""
Phylip writer
You could write a phylip output file pretty easily
but this will make it a little more convenient

The input is expected to be a list of seqs
"""
def phylip_writer(seqs,outfilename):
	numtax = str(len(seqs))
	numsites = str(len(seqs[0].seq))
	outfile = open(outfilename,"w")
	outfile.write(numtax+"\t"+numsites+"\n")
	for i in seqs:
		outfile.write(i.name+" "+i.seq+"\n")
	outfile.close()


"""
NEXUS writer
"""

