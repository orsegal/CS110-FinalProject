
a = [1,2,3]
print len(a)

if len(a) == 3:
	print "it is three"

if len(a) != 2:
	print "it is not equal to two"

if len(a) == 3 and a[0] == 2:
	print "length is 3 and first is 2"
elif len(a) == 3 and a[0] == 1:
	print "length is 3 and first is 1"
else:
	print "none of these are correct"
		
