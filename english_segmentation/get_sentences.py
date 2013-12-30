#!/usr/bin/env python

#import os
#_curdir=os.path.dirname(os.path.realpath(__file__))
#_activate=os.path.join(_curdir,"venv","bin","activate_this.py")
#execfile(_activate, dict(__file__=_activate))
import bz2
import nltk.data
from glob import glob
from random import shuffle
import sys
import re

def striphtml(data):
    p = re.compile(r'<.*?>')
    return p.sub('', data)

def remove_whitespace(w):
    return re.sub(' +',' ',w.strip())

if __name__ == "__main__":
    sent_detector = nltk.data.load('nltk:tokenizers/punkt/english.pickle')
    files = glob(sys.argv[1])
    shuffle(files)
    for fn in files[:100]:
        with bz2.BZ2File(fn) as f:
            for sentence in sent_detector.tokenize(striphtml(f.read())):
                print remove_whitespace(sentence.replace('\n',' '))
