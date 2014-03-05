#!/usr/bin/env python
import sys
import os, errno
import shutil

from rdkit import Chem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.Crippen as Crippen

from pprint import pprint
from collections import defaultdict, OrderedDict
from operator import itemgetter

import struct
from math import log, floor

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def maximal_words(words,selection):
    for p in selection:
        if not any(words[q].HasSubstructMatch(words[p]) and
                   not words[p].HasSubstructMatch(words[q])
                   for q in selection): yield p

def weight_bonds(words,idfs,mol):
    # Compute the word matches for bonds in the molecule
    bond_words=defaultdict(set)
    # Keep track of counts for use in TF-IDF later
    doc_word_counts = defaultdict(float)
    for i,w in enumerate(words):
        #mol_matches = mol.GetSubstructMatches(w,uniquify=False)
        mol_matches = mol.GetSubstructMatches(w,uniquify=False)
        if mol_matches: 
            doc_word_counts[i] += len(mol_matches)
            for m in mol_matches:
                cmap = dict(enumerate(m))
                for b in w.GetBonds():
                    start = b.GetBeginAtomIdx()
                    end = b.GetEndAtomIdx()
                    bond_words[frozenset([cmap[start],cmap[end]])].add(i)

    # Compute the maximal words
    #words_to_use = set(maximal_words(words,doc_word_counts.keys()))
    words_to_use = doc_word_counts.keys()

    # Compute the TF-IDF scores for each word
    maxtf = float(max(doc_word_counts[t] for t in words_to_use))
    score = defaultdict(float,((t,doc_word_counts[t]/maxtf*idfs[t])
                               for t in words_to_use))
    # Alternate, modified TF-IDF score
    #score = defaultdict(float,((t,log(1+doc_word_counts[t]/maxtf)*idfs[t])
    #                          for t in maxwords))

    # Get the combined TF-IDF scores for each bond
    bond_weights = dict((k,sum(score[t] for t in v)) for k,v in bond_words.items())
    # Return the bond values
    return bond_weights

def bond_weights_to_hmap(bond_weights,mol):
    """Converts bond weights to atom weights"""
    d=defaultdict(float)
    for k,v in bond_weights.items():
        for atom in k:
           if not d.has_key(atom): d[atom] = v
           else: d[atom] = min(d[atom],v) #why min??
           #else: d[atom] = max(d[atom],v) #why min??
    maxweight = float(max(d.values()))
    #return dict((i,oppenheimer_color(d[i]/maxweight)) for i in range(mol.GetNumAtoms()))
    return OrderedDict((i, color_tuple_clamp(d[i],maxweight)) for i in range(mol.GetNumAtoms()))

def color_tuple_clamp(val,maxweight):
    return (0,(val/maxweight)**4,0)

def label_mol(mol):
    for a in mol.GetAtoms():
        print >>sys.stderr, "IDX", a.GetIdx() + 1
        a.SetProp('molAtomMapNumber',str(a.GetIdx() + 1))

def get_idfs(counts,docs):
    """Computes the IDF scores given a table of counts"""
    return [log(float(docs)/float(c)) for c in counts]

def tup_to_hex(r,g,b):
    return struct.pack('BBB',
                int(r*255),
                int(g*255),
                int(b*255)).encode('hex')

def print_hmap(hmap):
    #print >>sys.stderr, "HMap size", len(hmap)
    return "000000#"+ "#".join( tup_to_hex(r,g,b)
                                for r,g,b in hmap.values())

def clamp_hex(v,maxw):
    return tup_to_hex(*color_tuple_clamp(v,maxw))

def get_color_pairs(mol,bw):
    cp=[]
    print "ahmap: ", bw
    maxweight = max( bw.values())
    for x in mol.GetBonds():
        start=x.GetBeginAtomIdx()
        end=x.GetEndAtomIdx()
        #1)colortupleclamp,2)tohex
        col=clamp_hex(bw.get(frozenset([start,end]),(0,0,0)), maxweight)
        cp.append({'a':start,'b':end,'color':col})
    return cp

if __name__ == "__main__":
    # Parse the arguments
    smart_stats_fn=sys.argv[1]
    numdocs=float(sys.argv[2])
    mol=Chem.MolFromSmiles(sys.argv[3])
    #dest_file=sys.argv[4]

    # Get all the mols with counts from smarts file
    words = []
    counts = []
    smarts = []
    with open(smart_stats_fn) as smart_stats:
        for line in smart_stats:
            c,s = line.strip().split('\t')
            words.append(Chem.MolFromSmarts(s,mergeHs=True))
            counts.append(float(c))
            smarts.append(s)

    # Compute the idf scores
    idfs = get_idfs(counts,numdocs)
    bweights = weight_bonds(words,idfs,mol)
    print bweights
    ahmap = bond_weights_to_hmap(bweights,mol)

    # Label the atoms
    label_mol(mol)
    smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
    maxweights = sorted(ahmap.values())[-2:]
    maxatoms=[idx+1 for idx,weight in ahmap.items() 
                    if weight in maxweights]
    print >>sys.stderr, maxatoms, smiles

    #sorted_bweights=sorted(bweights.iteritems(), key=itemgetter(1))
    #pprint(sorted_bweights)

    # Label the mol
    print 500, 500, '"'+ smiles + '"',  '"' + print_hmap(ahmap) + '"'
    import json
    colorstring = print_hmap(ahmap)
    print json.dumps({
        'width':500,
        'height':500,
        'smiles': smiles,
        'pairs': get_color_pairs(mol,bweights),
        'colors': colorstring,
    })
    # Print the moleculue
    #Draw.MolToFile(mol,fileName=dest_file,imageType="svg", highlightMap=bond_weights_to_hmap(bweights,mol))
