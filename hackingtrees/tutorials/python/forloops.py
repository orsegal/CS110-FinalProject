a = "ACCGTT"
b = [1,2,3,4,5]

print "print the a string: "+a
for i in a:
	print i

print "print the b list: "+str(b)
for i in b:
	print i

print "print the range of a number: "+str(range(4))
print "print the length of b: "+str(len(b))
print "print the range of len(b): "+str(range(len(b)))
#remember we add str() just to print
print range(len(b))

for i in range(len(b)):
	print i

for i in range(len(b)):
	print b[i]
