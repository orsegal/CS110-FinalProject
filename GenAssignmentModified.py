__author__ = 'orsegal'

from tabulate import tabulate
import numpy as np
from GenPhylogenyBuild import *
import sys
import pickle

def LCSTable(A, B):
    m = len(A)
    n = len(B)
    C = [[0] * (n + 1) for _ in range(m + 1)]
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            if A[i-1] == B[j-1]:
                C[i][j] = C[i-1][j-1] + 1
            else:
                C[i][j] = max(C[i-1][j], C[i][j-1])
    return C

def nicePrintLCSMatrix(C, A, B):
    NHeader = [i for i in A]
    print tabulate(C, showindex='always', headers=NHeader, floatfmt=".3f")

def findLCS(C, A, B, i = None, j = None):
    if i == None or j == None:
        i = len(A)
        j = len(B)
    if i == 0 or j == 0:
        return ""
    elif A[i-1] == B[j-1]:
        return findLCS(C, A, B, i-1, j-1) + A[i-1]
    else:
        if C[i][j-1] > C[i-1][j]:
            return findLCS(C, A, B, i, j-1)
        else:
            return findLCS(C, A, B, i-1, j)

def LCSDetailed(A, B):
    C = LCSTable(A, B)
    nicePrintLCSMatrix(C, A, B)
    LCS = findLCS(C, A, B)
    print "LCS: ", LCS
    print "LCS Length: ", len(LCS)

def LCSLen(A, B):
    C = LCSTable(A, B)
    LCS = findLCS(C, A, B)
    return len(LCS)

def EditDistance(A, B):
    m = len(A)
    n = len(B)
    C = [[0] * (n+1) for _ in range(m+1)]

    # Distance of string to the other when it's empty
    for i in xrange(m+1):
        C[i][0] = i
    for j in xrange(n+1):
        C[0][j] = j

    for j in xrange(1, n+1):
        for i in xrange(1, m+1):
            if A[i-1] == B[j-1]:
                C[i][j] = C[i-1][j-1]
            else:
                deletion = C[i-1][j] + 1
                insertion = C[i][j-1] + 1
                substitution = C[i-1][j-1] + 1
                C[i][j] = min(deletion, insertion, substitution)
    return C[m][n]

def distanceTable(strings):
    Nstrings = len(strings)
    DTable = [[0] * (Nstrings) for _ in range(Nstrings)]
    for i in xrange(Nstrings):
        for j in xrange(i+1, Nstrings):
            DTable[i][j] = DTable[j][i] = EditDistance(strings[i], strings[j])
    return DTable

def stringsTable(strings):
    Nstrings = len(strings)
    STable = [[0] * (Nstrings) for _ in range(Nstrings)]
    for i in xrange(Nstrings):
        for j in xrange(i+1, Nstrings):
            STable[i][j] = STable[j][i] = LCSLen(strings[i], strings[j])
    return STable

def StringsTableAnalysis(STable, strings, NHeader = None):
    Nstrings = len(strings)
    if NHeader == None:
        NHeader = range(Nstrings)
    print tabulate(STable, showindex='always', headers=NHeader, floatfmt=".3f")

def stringsRatioTable(STable, strings):
    Nstrings = len(strings)
    tempSTable = STable[:]
    for i in xrange(Nstrings):
        for j in xrange(Nstrings):
            # tempSTable[i][j] /= max(float(len(strings[i])), float(len(strings[j])))
            tempSTable[i][j] /= float(len(strings[i]))
    return tempSTable

def ratioSTableAnalysis(ratioSTable, strings, NHeader = None):
    Nstrings = len(strings)
    if NHeader == None:
        NHeader = range(Nstrings)
    print tabulate(ratioSTable, showindex='always', headers=NHeader, floatfmt=".6f")
    print "\nThe most and least similar strings:"
    for i in range(Nstrings):
        ratioAve = sum(ratioSTable[i])/len(ratioSTable[i])
        print NHeader[i], "- Ratio Average: %.4f" %ratioAve,
        print "    Minimal Ratio: %.4f" %min(ratioSTable[i][:i] + ratioSTable[i][i+1:])

def treeInference(strings, NChildren = 2):
    Nstrings = len(strings)
    familyTree = []
    fam = {}
    for i in xrange(Nstrings):
        fam[i] = strings[i]
    membersLeft = range(Nstrings)

    # First find the root by looking at the maxmin (maximal minimal ratio of closeness)
    STable = stringsTable(strings)
    ratioSTable = stringsRatioTable(STable, strings)
    minRatios = []
    for i in range(Nstrings):
        minRatios.append(min(ratioSTable[i][:i] + ratioSTable[i][i+1:]))
    tempMaxmin = max(minRatios)
    root = [a for a, b in enumerate(minRatios) if b == tempMaxmin] #check if there is more than one max
    if len(root) > 1:
        raise ValueError('There is no clear Root. One was chosen at random. Manually analysis is required!')
        root = root[:1]
        familyTree.append(root)
    else:
        familyTree.append(root)
    membersLeft.remove(root[0])

    # Now let's find the children, grandchildren, etc. until the last generation (in case that this one is not full)
    for i in xrange(1, 2**int(np.log2(len(strings)))-2):
        generation = []
        for currentMother in familyTree[-1]:
            branch = []
            tempRatio = ratioSTable[currentMother][:]
            tempRatio.sort()
            for _ in xrange(NChildren):
                currentChildRatio = tempRatio.pop()
                currentChild = [a for a, b in enumerate(ratioSTable[currentMother]) if b == currentChildRatio] #check if there is more than one max
                if len(currentChild) > 1:
                    raise ValueError('There is no clear child in branch height = %d. One was chosen at random. Manually analysis is required!')%i
                currentChild = currentChild[0]
                branch.append(currentChild)
                membersLeft.remove(currentChild)
            generation.append(branch)
        familyTree.append(generation)

    # Now for the last generation
    generation = []
    lastParents = familyTree[-1][:]
    while type(lastParents[0]) == list: #Make one list out of the list of lists in the last generation
        tempLst = []
        for innerLst in lastParents:
            for itm in innerLst:
                tempLst.append(itm)
        lastParents = tempLst[:]
    for _ in xrange(len(lastParents)):
        generation.append([])

    for currentChild in membersLeft:
        tempConnection = []
        for potentialParent in lastParents:
            tempConnection.append(ratioSTable[potentialParent][currentChild])
        trueParentRatio = max(tempConnection)
        childRatios = [] #Create a vector from a column in the matrix of ratios
        for p in xrange(Nstrings):
            childRatios.append(ratioSTable[p][currentChild])
        trueParent = childRatios.index(trueParentRatio)
        trueParentLoc = lastParents.index(trueParent)
        generation[trueParentLoc].append(currentChild)
    familyTree.append(generation)

    return familyTree

def drawNiceTree(tree):
    treeLst = tree[:]
    while type(treeLst[-1]) == list: #Make one list out of the list of lists in the last generation
        tempLst = []
        for innerLst in treeLst:
            if type(innerLst) == list:
                for itm in innerLst:
                    tempLst.append(itm)
            else:
                tempLst.append(innerLst)
        treeLst = tempLst[:]
    Tlen = len(treeLst)
    levels = int(np.log2(Tlen))
    space = 3
    width = 2**levels + (2**levels - 1)*space
    sideSpace = 0
    lines = []
    for lev in xrange(levels, -1, -1):
        line = ' ' * sideSpace
        for leafInd in xrange(2**lev, 0, -1):
            leaf = treeLst.pop(-leafInd)
            line += (str(leaf) + ' '*space)
        sideSpace = space - 1
        space = space*2 + 1
        lines.append(line)
    for i in xrange(len(lines)):
        print lines.pop()

def storeTable(strings, filename, tableFunc):
    STable = tableFunc(strings)
    with open(filename, 'wb') as fp:
        pickle.dump(STable, fp)

def openSTable(filename):
    with open(filename, 'rb') as fp:
        STable = pickle.load(fp)
    return STable

def importStrings(filename):
    strings = []
    i = 0
    currentString = ''
    stringName = ''
    with open(filename, 'r') as f:
        for line in f:
            thisLine = line.decode('utf-8-sig')
            if thisLine[0] == ">":
                if stringName != '':
                    strings.append((stringName, currentString))
                    currentString = ''
                lastIndx = thisLine.index(')')
                stringName = thisLine[12:(lastIndx+2)]
            else:
                currentString += thisLine[0:-1]
    return strings

def main():
    sys.setrecursionlimit(2000)
    
    #### Full Analysis
    # STable = stringsTable(strings)
    # StringsTableAnalysis(STable, strings, NHeaderlst=None)
    # ratioSTable = stringsRatioTable(STable, strings)
    # ratioSTableAnalysis(ratioSTable, strings, NHeader=None)

    #### Specific Analysis
    # LCSDetailed(strings[0], strings[1])
    # strings = strings[0:1] + strings[2:5] + strings[6:7]
    # NHeader = range(0, 1) + range(2, 5) + range(6, 7)
    # STable = stringsTable(strings)
    # StringsTableAnalysis(STable, strings, NHeaderlst=NHeader)
    # ratioSTable = stringsRatioTable(STable, strings)
    # ratioSTableAnalysis(ratioSTable, strings, NHeader=NHeader)

    strings = importStrings("seqdump.txt")
    strLst = []
    for strng in strings:
        strLst.append(strng[1])
    # STable = stringsTable(strLst)
    #
    A = strLst[3]
    B = strLst[2]
    print LCSLen(A, B)
    # # C = LCSTable(A, B)
    # # nicePrintLCSMatrix(C, A, B)
    print len(A), A
    print len(B), B
    # StringsTableAnalysis(STable, strLst)

    # storeTable(strLst, 'STable.txt', stringsTable)
    STable = openSTable('STable.txt')
    StringsTableAnalysis(STable, strLst)
    ratioSTable = stringsRatioTable(STable, strLst)
    ratioSTableAnalysis(ratioSTable, strLst, NHeader=None)
    # C = EditDistanceTable(A, B)
    # print C
    # DTable = distanceTable(strLst)
    # StringsTableAnalysis(DTable, strLst)

    # storeTable(strLst, 'DTable.txt', distanceTable)
    DTable = openSTable('DTable.txt')
    StringsTableAnalysis(DTable, strLst)


    # familyTree = treeInference(strings)
    # print familyTree
    # drawNiceTree(familyTree)


if __name__ == "__main__":
    main()





