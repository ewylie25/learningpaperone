#!/usr/bin/env python
import rdkit
from rdkit import Chem
import pymongo
from entropy import load_lib, entropy, symmetry
import sys

if __name__ == "__main__":
    db=pymongo.MongoClient().new_data
    cursor=db.chemical.find({"$and":[{"mol_weight":{"$gt":100}}
                                     , {"mol_weight":{"$lt":250}}
                                     #, {"in_degree":{"$gt":5}}
                                     , {"out_degree":{"$gt":5}}
                                    ]})
    lib = load_lib(sys.argv[1], int(sys.argv[2]))
    for entry in cursor:
        try:
            m = Chem.MolFromSmiles(str(entry["smiles"]))
            out = '\t'.join(map(str,[entry["in_degree"], 
                                     entry["out_degree"],
                                     entry["_id"],
                                     entropy(m,lib),
                                     symmetry(m,lib), 
                                     entry["smiles"],
                                     entry["brn"]
                                    ]))
            print out
        except: pass
