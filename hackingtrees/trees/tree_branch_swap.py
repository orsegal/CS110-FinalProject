from node import *
import random

"""
this is going to perform a nearest neighbor interchange
on a rooted tree
this works by 
1) picking an internal branch of the tree,
2) removing the four nodes that attach to this
3) replacing them in one of the four spots randomly
"""
def nni_rooted(tree,nnodes):
	#pick a node
	#while a node is a tip continue picking
	cnode = None
	count = 0
	stop = random.randint(0,nnodes)
	for i in tree.iternodes():
		if count == stop:
			cnode = i
		count += 1
	while cnode.istip or cnode.parent == None:
		count = 0
		stop = random.randint(0,nnodes)
		for i in tree.iternodes():
			if count == stop:
				cnode = i
			count += 1
	ndli = []
	#get the parent
	pnode = cnode.parent
	# if the parent is not the root get the parent of the parent
	# if the parent is the root, there should be three children, two other than the child
	#	then just change which goes to which , list = [nodes],use random.shuffle(list)
	if pnode.parent == None:
		ndli.append(random.sample(cnode.children,1)[0])
		for i in pnode.children:
			if i != cnode:
				pnode = i
				break
		ndli.append(random.sample(pnode.children,1)[0])
		for i in ndli:
			i.parent.remove_child(i)
		random.shuffle(ndli)
		cnode.add_child(ndli[0])
		pnode.add_child(ndli[1])
	else:
		ndli.append(random.sample(cnode.children,1)[0])
		ndli[0].parent.remove_child(ndli[0])
		tnode = None
		p = False
		if random.sample([0,1],1)[0] == 0:
			tnode = random.sample(pnode.children,1)[0]
			while tnode == cnode:
				tnode = random.sample(pnode.children,1)[0]
			tnode.parent.remove_child(tnode)
		else:
			tnode = pnode.parent
			tnode.remove_child(pnode)
			p = True
		ndli.append(tnode)
		random.shuffle(ndli)
		cnode.add_child(ndli[0])
		if p == False:
			pnode.add_child(ndli[1])
		if p == True:
			ndli[1].add_child(pnode)
	return tree

def nni_rooted2(tree,nnodes):
	#pick a node
	#while a node is a tip continue picking
	cnode = None
	count = 0
	stop = random.randint(0,nnodes)
	for i in tree.iternodes():
		if count == stop:
			cnode = i
		count += 1
	while cnode.istip:
		count = 0
		stop = random.randint(0,nnodes)
		for i in tree.iternodes():
			if count == stop:
				cnode = i
			count += 1
	ndli = []
	# if the parent is not the root get the parent of the parent
	# if the parent is the root, there should be three children, two other than the child
	#	then just change which goes to which , list = [nodes],use random.shuffle(list)
	children = cnode.children
	swap1 = None
	swap2 = None
	count = 0
	ndlip = []
	for i in children:
		if count >= 2:
			break
		if len(i.children) > 0:
			ndli.append(random.sample(i.children,1)[0])
			ndlip.append(i)
		else:
			ndli.append(i)
			ndlip.append(cnode)
	for i in range(len(ndlip)):
		ndlip[i].remove_child(ndli[i])
	for i in range(len(ndlip)):
		ndlip[-len(ndlip)+1-i].add_child(ndli[i])
	return tree

