#numbers

a = 5 
b = 7
c = a+b #a and b are ints, so is c
print c

a = 5.5 #a is now holding a floating point
c = a+b #b is still an int, but c is a floating point
print c

a = 1
b = 5
c = a/b
print c

a = 1
b = 5
c = float(a)/b
print c

a = 1.
b = 5
c = a/b
print c

#strings

a = "programming"
print a

prea = "i am doing a bit of "
a = "programming"
print prea+a

prea = "i like "
two = 2
a = " program"
print prea+str(two)+a

#lists
a = [1,2,3]
print a
a.append(4)
print a
print a[0]
print a[0]+a[1]
print a[-1]
print a[:2]
print a[-2:]
a.remove(2)
print a

#dictionary
a = {}
a['candy1'] = 10
a['candy2'] = 15
print a
a['anotherproduct'] = 20
print a
print a['anotherproduct']
print a['candy2']
