#!/usr/bin/env python

import sys
from random import shuffle
from deflate_sentences import remove_spaces
from collections import defaultdict
from longest_common_substring import time

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        frags = [frag.strip() for frag in f]
        frags = [f for f in frags if f]
    counts = defaultdict(int)
    with open(sys.argv[2]) as f:
        lines = f.readlines()
        shuffle(lines)
        i = 0
        percent = -1
        count = 10000
        for l in lines[:count]:
            p = int(i/float(count) * 100)
            if p > percent:
                percent = p 
                print >>sys.stderr, "{percent}%".format(percent=percent), time()
            i += 1
            blob, _ = remove_spaces(l.strip())
            for frag in frags:
                if frag in blob or frag[::-1] in blob:
                    counts[frag] += 1
    for k,v in counts.items():
        print v, '\t', k
