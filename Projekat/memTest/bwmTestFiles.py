import sys
import psutil
import os
from Bio import SeqIO

def Rotations(t):
    tt = t * 2
    return [tt[i : i + len(t)] for i in range(0, len(t))]

def BWM(t):
    return sorted(Rotations(t))

def BWTViaBWM(t):
    return ''.join(map(lambda x: x[-1], BWM(t)))

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
BWTViaBWM(seq[:50000] + '$')

process = psutil.Process(os.getpid())
mem = process.memory_info().rss / 1024 / 1024
print("Used this much memory: " + str(mem) + ' Mb')
