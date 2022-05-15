import sys
import psutil
import os

def Rotations(t):
    tt = t * 2
    return [tt[i : i + len(t)] for i in range(0, len(t))]

def BWM(t):
    return sorted(Rotations(t))

def BWTViaBWM(t):
    return ''.join(map(lambda x: x[-1], BWM(t)))

BWTViaBWM(sys.argv[1])
process = psutil.Process(os.getpid()) 
mem = process.memory_info().rss / 1024 / 1024 
print("Used this much memory: " + str(mem) + ' Mb')
