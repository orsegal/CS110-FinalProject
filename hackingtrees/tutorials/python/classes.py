
class Sequence:
	def __init__(self, inlabel,inseq):
		self.label = inlabel
		self.seq = inseq
	
	def rev_seq(self):
		self.seq = self.seq[::-1]
	
	def comp_seq(self):
		comp = {'A':'T','C':'G','G':'C','T':'A'}
		newseq = ""
		for i in self.seq: #there are better ways
			newseq += comp[i]
		self.seq = newseq

oneseq = Sequence("seq1","ACGAATTTCT")
print ">"+oneseq.label
print oneseq.seq
oneseq.rev_seq()
print "reverse: "+oneseq.seq
oneseq.rev_seq() #turn it back around
print "and back: "+oneseq.seq
oneseq.comp_seq()
print "comp    : "+oneseq.seq
