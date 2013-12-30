import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def idToIsotope(mol):
	for atom in mol.GetAtoms():
		atom.SetIsotope(atom.GetIdx())

def protectBondsMatchingTemplate(mol,template):
	for match in mol.GetSubstructMatches(template):
            base.GetAtomWithIdx(match[0]).SetProp('_protected','1')

def protectBondsNotMatchingTemplate(mol,template):
	matches = [match for match in mol.GetSubstructMatches(template)]
	rest = [i for i in m.GetAtoms() if i not in matches]
        base.GetAtomWithIdx(rest[0]).SetProp('_protected','1')

def mol_from_bond_somehow(b):
	a1 = b.GetBeginAtomIdx()
	a2 = b.GetEndAtomIdx()
	f=Chem.EditableMol(Chem.MolFromSmiles("[C:1]"))
	f.ReplaceAtom(0,b.GetBeginAtom())
	f.AddAtom(b.GetEndAtom())
	f.AddBond(a1,a2,b.GetBondType())
	return f.GetMol()

def mol_equals(m1,m2):
	return m1.HasSubstructureMatch(m2) and  m2.HasSubstructureMatch(m1) 

def whatChanged(rxn,mol):
	"""Runs a reaction (as in rxn.RunReactants), and with each result associates the bonds in the original mol that were modified"""
	#annotate mol with IsotopeID numbers
	idToIsotope(mol)
	d = dict((frozenset([b.GetBeginAtomIdx(),b.GetEndAtomIdx()]),mol_from_bond_somehow(b))
                 for b in mol.GetBonds())
	results = rxn.RunReactants([mol])

	for res in results:
		missing_bonds=set(d.keys())
		for b in res.GetBonds():
			bond_name=frozenset([b.GetStartAtom().GetIsotope(),
                                             b.GetEndAtom().GetIsotope()])
			if not d.has_key(bond_name): continue
			b_mol = mol_from_bond_somehow(b)
			old_b_mol= d[bond_name]
			if mol_equals(b_mol,old_b_mol): missing_bonds.remove(bond_name)
		yield (missing_bonds,res)
	

