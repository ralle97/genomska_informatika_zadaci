import sys
import psutil
import os
import time
from Bio import SeqIO
from pydivsufsort import divsufsort

def SuffixArrayImprovedSort(seq):
    return list(divsufsort(seq))

def BWTViaSAImprovedSort(seq, suffixArray):
    bwt = []
    for si in suffixArray:
        if si == 0:
            bwt.append('$')
        else:
            bwt.append(seq[si - 1])
    return ''.join(bwt)

def RankBWT(bw):
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots:
            tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots

def FirstColumn(tots):
    first = {}
    totc = 0
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count)
        totc += count
    return first

def SetRank(ranks, lColumn, char, lowerIndex, higherIndex):
    indexesOfChar = []
    
    for i in range(lowerIndex, higherIndex):#lower index is inclusive, higher is exclusive
        if lColumn[i] == char:
            indexesOfChar.append(i)
        
    if indexesOfChar:
        lowerIndex = min(indexesOfChar)
        higherIndex = max(indexesOfChar)
    else:
        return (-1, -1)
    
    return (ranks[lowerIndex], ranks[higherIndex])

def SetIndex(fColumn, char, lowerRank, higherRank):
    lowerIndex = fColumn[char][0] + lowerRank
    higherIndex = fColumn[char][0] + higherRank + 1 #+1 is to make higherIndex exclusive
    return (lowerIndex, higherIndex)

def SearchViaImprovedSort(sequence, pattern):
    positions = SuffixArrayImprovedSort(sequence)
    lColumn = BWTViaSAImprovedSort(sequence, positions)
    ranks, tots = RankBWT(lColumn)
    fColumn = FirstColumn(tots)
    
    lowerIndex = 0
    higherIndex = 0
    lowerRank = 0
    higherRank = 0
    firstIteration = True
    
    for char in reversed(pattern):
        if firstIteration:
            firstIteration = False
            if fColumn.get(char) is None:
                return [-1]
            else:
                (lowerIndex, higherIndex) = fColumn[char]
                continue
        (lowerRank, higherRank) = SetRank(ranks, lColumn, char, lowerIndex, higherIndex)
        if lowerRank == -1 or higherRank == -1:
            return [-1]
        else:
            (lowerIndex, higherIndex) = SetIndex(fColumn, char, lowerRank, higherRank)
    
    return [positions[i] for i in range(lowerIndex, higherIndex)]

def GetWholeSequenceFromFile(file):
    # all genome records of given FASTA file
    records = list(SeqIO.parse(file, "fasta"))
    sequence = ""
    # iterate over each record element
    for i in range(0, len(records)):
        sequence += records[i].seq
    # sequence is now whole
    return "".join(str(sequence).split())

def SearchOverGenomeWithImprovedSort(genome, pattern, stepSize):
    indexes = []
    totalTime = 0
    patternLength = len(pattern)
    
    subString = genome[:stepSize] + "$"
    
    startTime = time.time()
    tempIndexes = SearchViaImprovedSort(subString, pattern)
    endTime = time.time()
    totalTime += endTime - startTime
    
    indexes.append(list(filter(lambda x:x>-1, tempIndexes)))
    
    for i in range(1, (len(genome)//stepSize)+1):
        subString = genome[(i*stepSize)-patternLength+1:(i*stepSize)+patternLength-1] + "$"
        
        startTime = time.time()
        tempIndexes = SearchViaImprovedSort(subString, pattern)
        endTime = time.time()
        totalTime += endTime - startTime
        
        indexes.append(list(map(lambda x:x+(i*stepSize)-patternLength+1,filter(lambda x:x>-1, tempIndexes))))
        
        subString = genome[i*stepSize:i*stepSize+stepSize] + "$"
        
        startTime = time.time()
        tempIndexes = SearchViaImprovedSort(subString, pattern)
        endTime = time.time()
        totalTime += endTime - startTime
        
        indexes.append(list(map(lambda x:x+i*stepSize,filter(lambda x:x>-1, tempIndexes))))
        
    finalIndexes = [indexes[i] for i in range(0,len(indexes)) if indexes[i] != []]
    if finalIndexes:
        return ([item for sublist in finalIndexes for item in sublist], totalTime)
    else:
        return ([-1], totalTime)

seq = GetWholeSequenceFromFile("./data/13443_ref_Cara_1.0_chr1c.fa")
pattern = sys.argv[1]
stepSize = int(sys.argv[2])
SearchOverGenomeWithImprovedSort(seq, pattern, stepSize)

process = psutil.Process(os.getpid()) 
mem = process.memory_info().rss / 1024 / 1024 
print("Used this much memory: " + str(mem) + ' Mb')
