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
from subprocess import call

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

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
           else: d[atom] = min(d[atom],v)
    maxweight = float(max(d.values()))
    #return dict((i,oppenheimer_color(d[i]/maxweight)) for i in range(mol.GetNumAtoms()))
    return OrderedDict((i,(0,(d[i]/maxweight),0)) for i in range(mol.GetNumAtoms()))
        
def label_mol(mol):
    for a in mol.GetAtoms():
        print >>sys.stderr, "IDX", a.GetIdx() + 1
        a.SetProp('molAtomMapNumber',str(a.GetIdx() + 1))

def get_idfs(counts,docs):
    """Computes the IDF scores given a table of counts"""
    return [log(float(docs)/float(c)) for c in counts]

def print_hmap(hmap):
    print >>sys.stderr, "HMap size", len(hmap)
    return "000000#"+ "#".join(struct.pack('BBB',
                                int(r*255), 
                                int(g*255), 
                                int(b*255)).encode('hex') 
                                for r,g,b in hmap.values())

if __name__ == "__main__":
    # Parse the arguments
    smart_stats_fn=sys.argv[1]
    numdocs=float(sys.argv[2])
    molecules_file=sys.argv[3]
    svgs_dir=sys.argv[4]
    svg_data_file=sys.argv[5]

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

    svg_data=open(svg_data_file,'w')
    with open(molecules_file) as f:
        for line in f:
            _,_,id_,H,_,smiles,brn =  line.split('\t')
            print H, id_, smiles 
            try:
                mol = Chem.MolFromSmiles(smiles)
            except: continue
            # Compute the idf scores
            idfs = get_idfs(counts,numdocs)
            bweights = weight_bonds(words,idfs,mol)
            ahmap = bond_weights_to_hmap(bweights,mol)

            # Label the atoms
            label_mol(mol)
            new_smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
            maxweights = max(ahmap.values())
            maxatoms=[idx+1 for idx,weight in ahmap.items() 
                            if weight == maxweights]
            bestahmap={}
            for k,v in ahmap.items():
                if k+1 in maxatoms : bestahmap[k]=v
                else : bestahmap[k] = (0,0,0)
            print >>sys.stderr, maxatoms, new_smiles

            #sorted_bweights=sorted(bweights.iteritems(), key=itemgetter(1))
            #pprint(sorted_bweights)

            # Make the svg
            # If xvfb-run exists, use that
            if which("xvfb-run"):
                print ["xvfb-run", "java",  "-cp", "..", "-jar", "../colorer.jar",
                      os.path.join(svgs_dir,"{id}.svg".format(entropy=H,id=id_)), 
                      str(500), str(500), new_smiles , print_hmap(bestahmap)]
                call(["xvfb-run", "java",  "-cp", "..", "-jar", "../colorer.jar",
                      os.path.join(svgs_dir,"{id}.svg".format(entropy=H,id=id_)), 
                      str(500), str(500), new_smiles , print_hmap(bestahmap)])
            else:
                call(["java",  "-cp", "..", "-jar", "../colorer.jar", 
                      os.path.join(svgs_dir,"{id}.svg".format(entropy=H,id=id_)), 
                      str(500), str(500), new_smiles , print_hmap(bestahmap)])

            # Write out where the algorithm thought the atom seperated
            print >>svg_data,"\t".join(map(str,[id_,smiles,new_smiles,maxatoms,H,brn]))

            # Print the molecule
            #Draw.MolToFile(mol,fileName='test.svg',imageType="svg", highlightMap=bestahmap)
