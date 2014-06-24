#!/usr/bin/env python
import sys
import os
import errno
from rdkit import Chem
from collections import defaultdict, OrderedDict
import struct
from math import log


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def maximal_words(words, selection):
    for p in selection:
        if not any(words[q].HasSubstructMatch(words[p]) and
                   not words[p].HasSubstructMatch(words[q])
                   for q in selection):
            yield p


def weight_bonds(words, idfs, mol):
    # Compute the word matches for bonds in the molecule
    bond_words = defaultdict(set)
    # Keep track of counts for use in TF-IDF later
    doc_word_counts = defaultdict(float)
    for i, w in enumerate(words):
        mol_matches = mol.GetSubstructMatches(w, uniquify=False)
        if mol_matches:
            doc_word_counts[i] += len(mol_matches)
            for m in mol_matches:
                cmap = dict(enumerate(m))
                for b in w.GetBonds():
                    start = b.GetBeginAtomIdx()
                    end = b.GetEndAtomIdx()
                    bond_words[frozenset([cmap[start], cmap[end]])].add(i)

    # Compute the maximal words
    words_to_use = doc_word_counts.keys()

    # Compute the TF-IDF scores for each word
    maxtf = float(max(doc_word_counts[t] for t in words_to_use))
    score = defaultdict(float, ((t, doc_word_counts[t] / maxtf * idfs[t])
                                for t in words_to_use))
    # Alternate, modified TF-IDF score
    # score = defaultdict(float,((t,log(1+doc_word_counts[t]/maxtf)*idfs[t])
    #                          for t in maxwords))

    # Get the combined TF-IDF scores for each bond
    bond_weights = dict((k, sum(score[t] for t in v))
                        for k, v in bond_words.items())
    # Return the bond values
    return bond_weights


def bond_weights_to_hmap(bond_weights, mol):
    """Converts bond weights to atom weights"""
    d = defaultdict(float)
    for k, v in bond_weights.items():
        for atom in k:
            if atom not in d:
                d[atom] = v
            else:
                d[atom] = min(d[atom], v)
    maxweight = float(max(d.values()))
    return OrderedDict((i, color_tuple_clamp(d[i], maxweight)) for i in range(mol.GetNumAtoms()))


def color_tuple_clamp(val, maxweight):
    return (0, (val / maxweight) ** 4, 0)


def label_mol(mol):
    for a in mol.GetAtoms():
        a.SetProp('molAtomMapNumber', str(a.GetIdx() + 1))


def get_idfs(counts, docs):
    """Computes the IDF scores given a table of counts"""
    return [log(float(docs) / float(c)) for c in counts]


def tup_to_hex(r, g, b):
    return struct.pack('BBB',
                       int(r * 255),
                       int(g * 255),
                       int(b * 255)).encode('hex')


def print_hmap(hmap):
    return "000000#" + "#".join(tup_to_hex(r, g, b)
                                for r, g, b in hmap.values())


def clamp_hex(v, maxw):
    return tup_to_hex(*color_tuple_clamp(v, maxw))


def get_color_pairs(mol, bw):
    cp = []
    maxweight = max(bw.values())
    for x in mol.GetBonds():
        start = x.GetBeginAtomIdx()
        end = x.GetEndAtomIdx()
        col = clamp_hex(bw.get(frozenset([start, end]), (0, 0, 0)), maxweight)
        cp.append({'a': start + 1, 'b': end + 1, 'color': col})
    return cp

if __name__ == "__main__":
    # Parse the arguments
    smart_stats_fn = sys.argv[1]
    numdocs = float(sys.argv[2])
    mol = Chem.MolFromSmiles(sys.argv[3])

    # Get all the mols with counts from smarts file
    words = []
    counts = []
    smarts = []
    with open(smart_stats_fn) as smart_stats:
        for line in smart_stats:
            c, s = line.strip().split('\t')
            words.append(Chem.MolFromSmarts(s, mergeHs=True))
            counts.append(float(c))
            smarts.append(s)

    # Compute the idf scores
    idfs = get_idfs(counts, numdocs)
    bweights = weight_bonds(words, idfs, mol)
    # print bweights
    ahmap = bond_weights_to_hmap(bweights, mol)

    # Label the atoms
    label_mol(mol)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    maxweights = sorted(ahmap.values())[-2:]
    maxatoms = [idx + 1 for idx, weight in ahmap.items()
                if weight in maxweights]
    import json

    print 500, 500, json.dumps({
        'smiles': smiles,
        'pairs': get_color_pairs(mol, bweights),
    })
