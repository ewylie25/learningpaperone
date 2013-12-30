#!/usr/bin/env python

import sys
from random import shuffle

def remove_spaces(s):
    i = 0
    out = []
    spaces = []
    for c in s:
        if c == ' ': 
            spaces.append(i)
        else:
            i += 1
            out.append(c)
    return (''.join(out), spaces)

def add_spaces(w,spaces):
    out = []
    for i in range(len(w)):
        if i in spaces: out.append(' ')
        out.append(w[i])
    return ''.join(out)

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        lines = f.readlines()
        shuffle(lines)
        for l in lines[:100000]:
            blob, segments = remove_spaces(l.strip())
            print blob, '\t', segments

    #test = "Hello I am a sentence"
    #print test
    #test2, spaces = remove_spaces(test)
    #print test2
    #test3 = add_spaces(test2, spaces)
    #print test3
