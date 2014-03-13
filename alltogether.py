
from collections import defaultdict, OrderedDict
from math import log, floor

def label_mol(mol):
    for a in mol.GetAtoms():
        a.SetProp('molAtomMapNumber',str(a.GetIdx() + 1))

def weight_bonds(words,idfs,mol):
    # Compute the word matches for bonds in the molecule
    bond_words=defaultdict(set)
    # Keep track of counts for use in TF-IDF later
    doc_word_counts = defaultdict(float)
    for i,w in enumerate(words):
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
    words_to_use = doc_word_counts.keys()

    # Compute the TF-IDF scores for each word
    maxtf = float(max(doc_word_counts[t] for t in words_to_use))
    score = defaultdict(float,((t,doc_word_counts[t]/maxtf*idfs[t])
                               for t in words_to_use))

    # Get the combined TF-IDF scores for each bond
    bond_weights = dict((k,sum(score[t] for t in v)) for k,v in bond_words.items())
    # Return the bond values
    return bond_weights


def get_idfs(counts,docs):
    """Computes the IDF scores given a table of counts"""
    return [log(float(docs)/float(c)) for c in counts]

def get_color_pairs(mol,bw, maxes, mins):
    cp=[]
    palette={
                mins[0]:"FF0000",
                mins[1]:"C76939",
                mins[2]:"FFDA45",
                maxes[0]:"00FF00",
                maxes[1]:"00AA00",
                maxes[2]:"005500"
            }
    for x in mol.GetBonds():
        start=x.GetBeginAtomIdx()
        end=x.GetEndAtomIdx()
        score=bw.get(frozenset([start,end]), 0.0)
        cp.append({'a':start+1,'b':end+1,'color': palette.get(score, "000000")})

    return cp


if __name__ == "__main__":
    # arguments:
    # 1) fragment dictionary
    # 2) number of documents
    # 3) input file (list of smiles)
    # 4) output file.json

    import json
    import sys
    from rdkit import Chem

    smart_stats_fn=sys.argv[1]
    numdocs=float(sys.argv[2])

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
    data = []

    with open(sys.argv[3],'r') as inp:
        i=0
        for line in inp.readlines():
            i+=1
            try:
                mol=Chem.MolFromSmiles(line)
                bweights = weight_bonds(words,idfs,mol)

                allweights = sorted(list(set(bweights.values()))) #needs at least 6 diff kind of bonds!!
                maxweights = allweights[-3:] #top3
                minweights = allweights[:3] #bottom 3

                # Label the atoms
                label_mol(mol)
                smiles=Chem.MolToSmiles(mol, isomericSmiles=True)

                data.append({
                "format":{
                    'width':500,
                    'height':500,
                    'fname' : str(i) + ".svg"
                    },
                "data":{
                    'smiles': smiles,
                    'pairs': get_color_pairs(mol,bweights,maxweights,minweights),
                    }
                })
            except Exception as e:
                print "line", i , line, "failed."

    with open(sys.argv[4],'w') as outfile:
        json.dump({'entries':data}, outfile, sort_keys=True,indent=4, separators=(',', ': '))