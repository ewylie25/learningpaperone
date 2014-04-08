#!/usr/bin/env python
import rdkit
from rdkit import Chem
import pymongo
import sys

if __name__ == "__main__":
    db=pymongo.MongoClient('129.105.205.35').new_data
    cursor=db.chemical.find({"$and":[{"mol_weight":{"$gt":100}}
                                     , {"mol_weight":{"$lt":850}}
                                     #, {"in_degree":{"$gt":5}}
                                     , {"out_degree":{"$gt":5}}
                                    ]})
    for entry in cursor[:int(sys.argv[1])]:
        print entry["smiles"]
           
          
       
