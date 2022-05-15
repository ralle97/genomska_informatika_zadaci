import sys
import psutil
import os

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

BWTViaSA(sys.argv[1])
process = psutil.Process(os.getpid()) 
mem = process.memory_info().rss / 1024 / 1024 
print("Used this much memory: " + str(mem) + ' Mb')
