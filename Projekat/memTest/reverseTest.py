import sys
import psutil
import os

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

def ReverseBWT(bw):
    ranks, tots = RankBWT(bw)
    first = FirstColumn(tots)
    rowi = 0   # first row
    t = '$'    # rightmost character
    while bw[rowi] != '$':
        c = bw[rowi]
        t = c + t    # prepend to answer
        # jump to row that starts with c of same rank
        rowi = first[c][0] + ranks[rowi]
    return t

ReverseBWT(sys.argv[1])
process = psutil.Process(os.getpid()) 
mem = process.memory_info().rss / 1024 / 1024 
print("Used this much memory: " + str(mem) + ' Mb')
