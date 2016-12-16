import sys

infile = open(sys.argv[1],"r")
countries = []
lats = []
longs = []
for i in infile:
	i = i.strip() #take off whitespace
	spls = i.split(",") #split with comma
	countries.append(spls[0])
	lats.append(spls[1])
	longs.append(spls[2])
infile.close()
print "number of records: "+str(len(countries))
