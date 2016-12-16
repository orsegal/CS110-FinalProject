import sys

infile = open(sys.argv[1],"r")
for i in infile:
	print i.strip()
infile.close()
