#!/usr/bin/env python
from math import log
from rdkit import Chem
from collections import defaultdict
import sys

def load_lib(smarts_stats_text, count):
    """Loads a library of chemical data"""
    mol_stats=[]
    with open(smarts_stats_text) as smart_stats:
        for line in smart_stats:
            c,s = line.strip().split('\t')
            mol_stats.append((s, Chem.MolFromSmarts(s,mergeHs=True), -log(float(c)/count)))
    return mol_stats

def entropy(mol,lib):
    """Computes the TF-IDF entropy of a molecule"""
    raw_tf_idf=[(len(mol.GetSubstructMatches(frag)),idf) 
                for _,frag,idf in lib]
    maxtf = float(max(tf for tf,_ in raw_tf_idf))
    return sum(tf/maxtf*idf for tf,idf in raw_tf_idf)

def symmetry(mol,lib):
    """Uses a library of chemical fragments to estimate the symmetry"""
    # Compute the word matches for bonds in the molecule
    atom_words=defaultdict(set)
    for i,entry in enumerate(lib):
        _,w,_ = entry
        #mol_matches = mol.GetSubstructMatches(w,uniquify=False)
        mol_matches = mol.GetSubstructMatches(w,uniquify=False)
        if mol_matches: 
            for m in mol_matches:
                cmap = dict(enumerate(m))
                for b in w.GetBonds():
                    atom_words[cmap[b.GetBeginAtomIdx()]].add(i)
                    atom_words[cmap[b.GetEndAtomIdx()]].add(i)
    atom_classes = defaultdict(set)
    for a,c in atom_words.items():
        atom_classes[frozenset(c)].add(a)
    return max(len(v) for v in atom_classes.values())

if __name__ == "__main__":
    lib = load_lib(sys.argv[1], int(sys.argv[2]))
    with open(sys.argv[3]) as data:
        for line in data:
            a,b,smiles = line.strip().split('\t')
            try:
                m = Chem.MolFromSmiles(smiles)
                print '\t'.join(map(str,[symmetry(m,lib),a,b,smiles,entropy(m,lib)]))
            except: pass # nothing to see here...
