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
uniq_countries = set(countries) # a set only holds unique records
print "number of records: "+str(len(countries))
print "number of unique records: "+str(len(uniq_countries))
print "here are the first ten:"
count = 0
for i in uniq_countries:
	if count < 10:
		print i
	else:
		break
	count += 1
