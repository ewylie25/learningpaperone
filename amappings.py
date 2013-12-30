import rdkit
from rdkit import Chem

class Mapper():
    def __init__(self,mol):
        self.mol=mol

    def getMappedAtoms(self):
        '''returns a list of the atoms which have a mapping assigned'''
        self.mappedatoms  = [atom for atom in self.mol.GetAtoms() if atom.HasProp("molAtomMapNumber") ]
        self.maps = [int(atom.GetProp("molAtomMapNumber"))for atom in self.mol.GetAtoms() if atom.HasProp("molAtomMapNumber")]
        return [self.mappedatoms, self.maps]
    
    def usedIds(self):
        '''returns a list of the mappings already used'''
        return self.maps

    def getUnmappedAtoms(self):
        '''returns the a list of unmapped atoms present in the molecule'''
        self.unmappedatoms=[atom for atom in self.mol.GetAtoms() if not atom.HasProp("molAtomMapNumber")]
        return self.unmappedatoms
    def getMapList(self):
        '''returns a list of indexes, without the ones already used in the mappings'''
        self.maplist=[i for i in range(self.mol.GetNumAtoms()) if i not in self.usedIds() ]
        return self.maplist
    
    def setMapList(self, atomlist, idlist):
        '''assign to each atom in atom list the corresponding id in idlist'''
        i=0
        for each in atomlist:
            if not each.HasProp("molAtomMapNumber"):
                each.SetProp("molAtomMapNumber",str(idlist[i]))
                i=i+1
        