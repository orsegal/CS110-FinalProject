import sys

def run_count(inword,incount):
	for i in range(incount):
		print inword

if __name__ == "__main__":
	print sys.argv
	word = sys.argv[1]
	count = int(sys.argv[2]) #we want it to be a number
	run_count(word,count)
