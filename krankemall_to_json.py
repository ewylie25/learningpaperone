

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
    numdocs=float(sys.argv[2])

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
    idfs = wb.get_idfs(counts,numdocs)

    data = []

    with open(sys.argv[3],'r') as inp:
        i=0
        for line in inp.readlines():
            i+=1
            try:
                mol=Chem.MolFromSmiles(line)
                bweights = wb.weight_bonds(words,idfs,mol)
                ahmap = wb.bond_weights_to_hmap(bweights,mol)

                # Label the atoms
                wb.label_mol(mol)
                smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
                maxweights = sorted(ahmap.values())[-2:]
                maxatoms=[idx+1 for idx,weight in ahmap.items()
                            if weight in maxweights]
                data.append({
                "format":{
                    'width':500,
                    'height':500,
                    'fname' : str(i) + ".svg"
                    },
                "data":{
                    'smiles': smiles,
                    'pairs': wb.get_color_pairs(mol,bweights),
                    }
                })
            except Exception as e:
                print "line", i , line, "failed."

    with open(sys.argv[4],'w') as outfile:
        json.dump({'entries':data}, outfile, sort_keys=True,indent=4, separators=(',', ': '))