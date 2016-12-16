import string, sys
from hackingtrees.trees.node import Node


"""
this takes a newick string as instr
and reads the string and makes the 
nodes and returns the root node
"""
def read_tree_string(instr):
	root = None
	index = 0
	nextchar = instr[index]
	start = True
	keepgoing = True
	curnode = None
	while keepgoing == True:
		if nextchar == "(":
			if start == True:
				root = Node()
				curnode = root
				start = False
			else:
				newnode = Node()
				curnode.add_child(newnode)
				curnode = newnode
		elif nextchar == ',':
			curnode = curnode.parent
		elif nextchar == ")":
			curnode = curnode.parent
			index += 1
			nextchar = instr[index]
			name = ""
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			curnode.label = name
			index -= 1
		elif nextchar == ';':
			keepgoing = False
			break
		elif nextchar == ":":
			index += 1
			nextchar = instr[index]
			brlen = ""
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				brlen += nextchar
				index += 1
				nextchar = instr[index]
			curnode.length = float(brlen)
			index -= 1
		elif nextchar == ' ':
			index += 1
			nextchar = instr[index]
		else: # this is an external named node
			newnode = Node()
			curnode.add_child(newnode)
			curnode = newnode
			curnode.istip = True
			name = ""
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			curnode.label = name
			index -= 1
		if index < len(instr) - 1:
			index += 1
		nextchar = instr[index]
	return root

if __name__ == "__main__":
	# s = "(a:3,(b:1e-05,c:1.3)int_|_and_33.5:5)root;"
	# s = "(A,B,(C,D)E)F;"
	s = "((('Influenza A virus (A/Argentina/3509/2006(H3N2))','Influenza A virus (A/Argentina/AG-3509/2006(H3N2))')0.9600,NewFluStrainGeneSequence)0.5700,'Influenza A virus (A/MEXICO/2519/2006(H3N2))',('Influenza A Virus (A/Iquitos/NFI310/2007(H3N2))',('Influenza A virus (A/Peru/3355/2006(H3N2))','Influenza A virus (A/Peru/PE-FLU-3355/2006(H3N2))')0.7000)0.8000,('Influenza A virus (A/BRISBANE/35/2006(H3N2))','Influenza A virus (A/NORWAY/566-1/2006(H3N2))','Influenza A virus (A/AUCKLAND/78/2006(H3N2))','Influenza A virus (A/WAIKATO/49/2006(H3N2))','Influenza A virus (A/AUCKLAND/47/2006(H3N2))','Influenza A virus (A/NORWAY/622/2006(H3N2))')0.9500,('Influenza A Virus (A/Piura/FSP1230/2006(H3))','Influenza A virus (A/Piura/FLU3482/2006(H3))')0.6000,'Influenza A Virus (A/Piura/FLU3646/2006(H3N2))','Influenza A virus (A/Tumbes/AST0136/2006(H3))''Influenza A virus (A/Tumbes/AST0129/2006(H3))','Influenza A Virus (A/Piura/FSP1248/2006(H3))','Influenza A Virus (A/Piura/FLU3638/2006(H3N2))','Influenza A virus (A/Tumbes/AST0074/2006(H3))');"
	# n2 = read_tree_string(s)
	# print n2.get_newick_repr(True)
	# print n2.children[1].label
	print s
	print s[50:100]
