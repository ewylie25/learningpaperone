#!/usr/bin/env python
import sys
from random import shuffle
import itertools
from math import factorial
import datetime 

def time(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

def binomial(n,m): 
    return factorial(n) / (factorial(m) * factorial(n-m))

def LongestCommonSubstring(S1, S2):
   M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))]
   longest, x_longest = 0, 0
   for x in xrange(1,1+len(S1)):
       for y in xrange(1,1+len(S2)):
           if S1[x-1] == S2[y-1]:
               M[x][y] = M[x-1][y-1] + 1
               if M[x][y]>longest:
                   longest = M[x][y]
                   x_longest  = x
           else:
               M[x][y] = 0
   return S1[x_longest-longest: x_longest]

def DirectionlessLongestCommonSubstring(S1,S2):
    a = LongestCommonSubstring(S1, S2)
    b = LongestCommonSubstring(S1, S2[::-1])
    if len(a) >= len(b): return a
    else : return b

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        lines = f.readlines()
        shuffle(lines)
        select = 400
        percent = -1
        i = 0
        count = float(binomial(select,2))
        for s1,s2 in itertools.combinations(lines[:select],2):
            p = int(i/count * 100)
            if p > percent:
                percent = p 
                print >>sys.stderr, "{percent}%".format(percent=percent), time()
            i += 1
            print DirectionlessLongestCommonSubstring(s1.strip(),s2.strip())
