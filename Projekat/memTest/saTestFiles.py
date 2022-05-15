import sys
import psutil
import os
from Bio import SeqIO

def SuffixArray(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # Extract and return just the offsets
    return map(lambda x: x[1], satups)

def BWTViaSA(t):
    bw = []
    for si in SuffixArray(t):
        if si == 0:
            bw.append('$')
        else:
            bw.append(t[si - 1])
    return ''.join(bw) # returns string version of list bw

def GetWholeSequenceFromFile(file):
    # all genome records of given FASTA file
    records = list(SeqIO.parse(file, "fasta"))
    sequence = ""
    # iterate over each record element
    for i in range(0, len(records)):
        sequence += records[i].seq
    # sequence is now whole
    return "".join(str(sequence).split())

seq = ""
if int(sys.argv[1]) == 1:
	seq = GetWholeSequenceFromFile("./data/13443_ref_Cara_1.0_chr1c.fa")
elif int(sys.argv[1]) == 2:
	seq = GetWholeSequenceFromFile("./data/10093_ref_PAHARI_EIJ_v1.1_chrX.fa")
elif int(sys.argv[1]) == 3:
	seq = GetWholeSequenceFromFile("./data/144034_ref_Pbar_UMD_V03_chrUn.fa")
BWTViaSA(seq[:50000] + '$')

process = psutil.Process(os.getpid()) 
mem = process.memory_info().rss / 1024 / 1024 
print("Used this much memory: " + str(mem) + ' Mb')
