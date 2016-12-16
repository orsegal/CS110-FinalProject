import sys

print sys.argv
word = sys.argv[1]
count = int(sys.argv[2]) #we want it to be a number

for i in range(count):
	print word
