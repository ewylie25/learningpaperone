#!/usr/bin/env python
import sys
from rdkit import Chem
from collections import defaultdict
from math import log
import json


def weight_bonds(words_list, idf_list, molecule):
    # Compute the word matches for bonds in the molecule
    bond_words = defaultdict(set)
    # Keep track of counts for use in TF-IDF later
    doc_word_counts = defaultdict(float)
    for i, w in enumerate(words_list):
        mol_matches = molecule.GetSubstructMatches(w, uniquify=False)
        if mol_matches:
            doc_word_counts[i] += len(mol_matches)
            for m in mol_matches:
                c_map = dict(enumerate(m))
                for b in w.GetBonds():
                    start = b.GetBeginAtomIdx()
                    end = b.GetEndAtomIdx()
                    bond_words[frozenset([c_map[start], c_map[end]])].add(i)

    # Compute the maximal words
    words_to_use = doc_word_counts.keys()

    # Compute the TF-IDF scores for each word
    max_tf = float(max(doc_word_counts[t] for t in words_to_use))
    score = defaultdict(float, ((t, doc_word_counts[t] / max_tf * idf_list[t])
                                for t in words_to_use))

    # Get the combined TF-IDF scores for each bond
    bond_weights_dict = dict((k, sum(score[t] for t in v))
                             for k, v in bond_words.items())
    # Return the bond values
    return bond_weights_dict


def label_mol(molecule):
    for a in molecule.GetAtoms():
        a.SetProp('molAtomMapNumber', str(a.GetIdx() + 1))


def get_idfs(counts_list, docs):
    """Computes the IDF scores given a table of counts"""
    return [log(float(docs) / float(c)) for c in counts_list]


def process_scores(molecule, bw):
    scores = []
    for x in molecule.GetBonds():
        start = x.GetBeginAtomIdx()
        end = x.GetEndAtomIdx()
        score = bw.get(frozenset([start, end]))
        scores.append({'a': start + 1, 'b': end + 1, 'score': score})
    return scores

if __name__ == "__main__":
    # Parse the arguments
    smart_stats_fn = sys.argv[1]
    num_docs = float(sys.argv[2])
    mol = Chem.MolFromSmiles(sys.argv[3])
    dest_file = sys.argv[4]

    # Get all the mols with counts from smarts file
    words = []
    counts = []
    smarts = []
    with open(smart_stats_fn) as smart_stats:
        for line in smart_stats:
            count, smart = line.strip().split('\t')
            words.append(Chem.MolFromSmarts(smart, mergeHs=True))
            counts.append(float(count))
            smarts.append(smart)

    # Compute the idf scores
    idf_scores = get_idfs(counts, num_docs)
    bond_weights = weight_bonds(words, idf_scores, mol)

    # Label the atoms
    label_mol(mol)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    with open(dest_file, 'w') as f:
        json.dump({
            'smiles': smiles,
            'pairs': process_scores(mol, bond_weights),
        }, f, indent=4)