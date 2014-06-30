def get_color_pairs(mol,bw, maxes, mins):
    cp=[]
    maxweight = max( bw.values()) # this is for normalization.. might as well throw it away:

    for x in mol.GetBonds():
        start=x.GetBeginAtomIdx() #this is rdkit only!
        end=x.GetEndAtomIdx()

        #1)colortupleclamp,2)tohex
        score=bw.get(frozenset([start,end]), 0.0) #returns the scoretuple.. right?
        #print score
        palette={
                mins[0]:"FF0000",
                mins[1]:"C76939",
                mins[2]:"FFDA45",
                maxes[0]:"00FF00",
                maxes[1]:"00AA00",
                maxes[2]:"005500"
            }
        cp.append({'a':start+1,'b':end+1,'color':palette.get(score, "000000")})

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
    import weight_bonds_to_json as wb

    smart_stats_fn=sys.argv[1]
    #numdocs=float(sys.argv[2])

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

    numdocs = float(len(words))

    # Compute the idf scores
    idfs = wb.get_idfs(counts,numdocs)

    data = []

    with open(sys.argv[2],'r') as inp:
        i=0
        for line in inp.readlines():
            i+=1
            try:
                mol=Chem.MolFromSmiles(line)
                bweights = wb.weight_bonds(words,idfs,mol)

                allweights = sorted(list(set(bweights.values()))) #needs at least 6 diff kind of bonds!!
                maxweights = allweights[-3:] #top3
                minweights = allweights[:3] #bottom 3



                # Label the atoms
                wb.label_mol(mol)
                smiles=Chem.MolToSmiles(mol, isomericSmiles=True)

                # maxatoms=[idx+1 for idx,weight in ahmap.items()
                #             if weight in maxweights]
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

    with open(sys.argv[3],'w') as outfile:
        json.dump({'entries':data}, outfile, sort_keys=True,indent=4, separators=(',', ': '))
