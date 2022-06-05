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

class FMIndexImproved():
    @staticmethod
    def SampleSuffixArray(suffixArray, step = 32):
        sampledSA = {}
        for index, suffix in enumerate(suffixArray):
            if suffix % step == 0:
                sampledSA[index] = suffix
        return sampledSA
    
    def __init__(self, seq, suffixArray = None, checkpointStep = 128, sampledSAStep = 32):
        if seq[-1] != '$':
            seq += '$'
        if suffixArray == None:
            suffixArray = SuffixArrayImprovedSort(seq)
        self.bwt = BWTViaSAImprovedSort(seq, suffixArray)
        self.sampledSA = self.SampleSuffixArray(suffixArray, sampledSAStep)
        self.length = len(self.bwt)
        
        self.CreateCheckpoints(checkpointStep)
        
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count
    
    def CreateCheckpoints(self, checkpointStep = 128):
        self.checkpoints = {}
        self.checkpointStep = checkpointStep
        tally = {}
        
        for char in self.bwt:
            if char not in tally:
                tally[char] = 0
                self.checkpoints[char] = []
        
        for index, char in enumerate(self.bwt):
            tally[char] += 1
            if index % checkpointStep == 0:
                for c in tally.keys():
                    self.checkpoints[c].append(tally[c])
    
    def Rank(self, bwt, char, row):
        if row < 0 or char not in self.checkpoints:
            return 0
        index, numOccurences = row, 0
        
        while index % self.checkpointStep != 0:
            if bwt[index] == char:
                numOccurences += 1
            index -= 1
        return self.checkpoints[char][index // self.checkpointStep] + numOccurences
    
    def Range(self, pattern):
        left, right = 0, self.length - 1
        for i in range(len(pattern) - 1, -1, -1):
            left = self.Rank(self.bwt, pattern[i], left - 1) + self.Count(pattern[i])
            right = self.Rank(self.bwt, pattern[i], right) + self.Count(pattern[i]) - 1
            if right < left:
                break
        return left, right + 1
    
    def Resolve(self, row):
        def StepLeft(row):
            char = self.bwt[row]
            return self.Rank(self.bwt, char, row - 1) + self.Count(char)
        
        numSteps = 0
        while row not in self.sampledSA:
            row = StepLeft(row)
            numSteps += 1
        return self.sampledSA[row] + numSteps
    
    def Count(self, char):
        if char not in self.first:
            for cc in sorted(self.first.keys()):
                if char < cc:
                    return self.first[cc]
            return self.first[cc]
        else:
            return self.first[char]
    
    def HasSubstring(self, pattern):
        left, right = self.Range(pattern)
        return right > left
    
    def HasSuffix(self, pattern):
        left, right = self.Range(pattern)
        if left >= self.length:
            return False
        offset = self.Resolve(left)
        return right > left and offset + len(pattern) == self.length - 1
    
    def Search(self, pattern):
        left, right = self.Range(pattern)
        return [self.Resolve(x) for x in range(left, right)]

def GetWholeSequenceFromFile(file):
    # all genome records of given FASTA file
    records = list(SeqIO.parse(file, "fasta"))
    sequence = ""
    # iterate over each record element
    for i in range(0, len(records)):
        sequence += records[i].seq
    # sequence is now whole
    return "".join(str(sequence).split())

seq = GetWholeSequenceFromFile("./data/13443_ref_Cara_1.0_chr1c.fa")
checkpointStep = int(sys.argv[1])
sampledSAStep = int(sys.argv[2])
pattern = sys.argv[3]
fm = FMIndexImproved(seq, None, checkpointStep, sampledSAStep)
fm.Search(pattern)

process = psutil.Process(os.getpid())
mem = process.memory_info().rss / 1024 / 1024
print("Used this much memory: " + str(mem) + ' Mb')
