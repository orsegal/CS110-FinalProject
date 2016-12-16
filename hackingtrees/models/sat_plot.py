import math
from sequences.sequence import Sequence
from model_calc import *

"""
this will return two lists one with distances
the other with corrected distances as calculated by 
the JC correction
raw distances are calculated as the proporiton of 
differences between A->C and A->G and A->T and so on
"""
def get_sat_plot_data(seqs):
	lenseqs = len(seqs)
	for i in range(lenseqs):
		for j in range(lenseqs):
			if j > i:
				df = prop_diff(seqs[i].seq,seqs[j].seq)
				jdf = jc_diff(df)
				yield df,jdf
