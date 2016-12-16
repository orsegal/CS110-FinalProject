def sum(numlist):
	num = 0
	for i in numlist:
		num = num + i
	return num

def maxnum(numlist):
	curmax = numlist[0]
	for i in numlist:
		if i > curmax:
			curmax = i
	return curmax

a = [1,2,3,4]
b = sum(a)
print "sum: "+str(b)
c = maxnum(a)
print "max: "+str(c)
