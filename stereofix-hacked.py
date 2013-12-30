import copy
import re
import itertools
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondDir, BondType, ChiralType, BondStereo
from environments import compact


def run_reactants(rxnSmarts, reactantSmiles,
                  needsExplicitHs=False, needsEnvironments=True,
                  needsStereoTunneling=True,
                  error=None):
    """Return a list of unique reaction products.

    Function is a wrapper around the RDKit built-in RunReaction() which attempts to
    address its limitations related to handling stereochemistry in reaction cores.
    """
    if error == None : error=logging.getLogger("stereotest.stereofix").error

    # By default, rewrite the transform using SMARTS environments.
    if needsEnvironments:
        rxnSmarts = compact(rxnSmarts)

    # Generate reaction from input SMARTS.
    try:
        rxn = AllChem.ReactionFromSmarts(rxnSmarts)
        rxn.Initialize()
    except:
        error('Invalid SMARTS: {0}.'.format(rxnSmarts))
        return None

    # Convert reactants SMILES to molecules and add explicit H atoms need if
    # requested.
    reactants = [Chem.MolFromSmiles(s) for s in reactantSmiles.split('.')]
    if needsExplicitHs:
        reactants = [add_Hs(m) for m in reactants]

    if None not in reactants:
        [m.UpdatePropertyCache() for m in reactants]
    else:
        error('Invalid SMILES among reactants.')
        return None

    # Number atoms in reactants.
    number_atoms(reactants)

    # Create reactant and product atom maps.
    reactantAtomMap = create_atom_map(reactants)

    # Temporarily remove chiral centers from the reactants.
    if needsStereoTunneling:
        backup = copy.deepcopy(reactants)
        for reactant in reactants:
            for rAtom in reactant.GetAtoms():
                if rAtom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                    rAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)

    # Run reaction transform.
    try:
        productList = rxn.RunReactants(tuple(reactants))
    except:
        error('Applying transform failed.')
        return None

    # Restore chiral centers to reactants and introduce them in relevant atoms
    # of products.
    if needsStereoTunneling:
        reactants = backup

        for products in productList:
            for product in products:
                for pAtom in product.GetAtoms():
                    if pAtom.GetIsotope() not in reactantAtomMap.keys():
                        continue

                    rIdx, rAtomIdx = reactantAtomMap[pAtom.GetIsotope()]
                    rAtom = reactants[rIdx].GetAtomWithIdx(rAtomIdx)
                    if rAtom.GetChiralTag() != pAtom.GetChiralTag():
                        copy_chirality(rAtom, pAtom)

    # Create template atom maps.
    [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps, IsMatchChiral] = create_template_map(reactants, rxn)
    if len(IsMatchChiral) != len(productList):
        error('Template maps does not match products.')
        return None

    # Remove non-chiral matches AND check for ring-opening reactions.
    newProductList = []
    rTemplateListIdx = 0
    for idx, products in enumerate(productList):
        if IsMatchChiral[idx]:
            newProducts = list(products)

            while True:
                [mergeIdx, newProducts] = merge(newProducts)
                if mergeIdx < 0:
                    break

            # Fix bonds
            #
            # Note:
            # It is probably an overkill, but I am not sure at the moment
            # how to pass relevant reactant and reactant template in case
            # of multi reactants reaction. Thus, using a brute force
            # approach, I iterate over all of them.
            for newProductIdx in range(len(newProducts)):
                for reactant, rTemplate in zip(reactants, rTemplateList[rTemplateListIdx]):
                    newProducts[newProductIdx] = fix_bonds(reactant, rTemplate, newProducts[newProductIdx])

            newProductList.append(tuple(newProducts))
            rTemplateListIdx += 1
    productList = tuple(newProductList)

    # Parse products to correct chirality.
    for matchIdx, products in enumerate(productList):
        for product in products:
            for productAtom in product.GetAtoms():
                if not productAtom.HasProp('Core'):
                    # Unique atom ID
                    atomIdx = productAtom.GetIsotope()

                    # Map Atoms
                    tmp = reactantAtomMap[atomIdx] if atomIdx in reactantAtomMap.keys() else None
                    reactantAtom = reactants[tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None
                    tmp = rTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in rTemplateAtomMaps[matchIdx].keys() else None
                    rTemplateAtom = rTemplateList[matchIdx][tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None
                    tmp = pTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in pTemplateAtomMaps[matchIdx].keys() else None
                    pTemplateAtom = pTemplateList[matchIdx][tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None

                    # Correct Chirality and charge of Product Atom
                    fix_chirality(reactantAtom, productAtom, rTemplateAtom, pTemplateAtom)
                    fix_charge(reactantAtom, productAtom, rTemplateAtom, pTemplateAtom)

    # Remove Hs if necessary
    if needsExplicitHs:
        tempProductList = []
        for products in productList:
            mols = [remove_Hs(p) for p in products]
            if None in mols:
                return None
            tempProductList.append(mols)
        productList = tempProductList

    # Remove atom mappings.
    for products in productList:
        for product in products:
            for atom in product.GetAtoms():
                atom.SetIsotope(0)
                atom.ClearProp('Core')

    # Get unique product sets.
    uniqueProductSets = []
    for products in productList:
        uniqueProductSets.append('.'.join([Chem.MolToSmiles(p, isomericSmiles=True) for p in products]))
    uniqueProductsList = uniquify(uniqueProductSets)

    # Convert to molecules.
    productList = []
    for s in uniqueProductsList:
        products = [Chem.MolFromSmiles(smi) for smi in s.split('.')]
        productList.append(products)

    # Get transpose matrix of the productList.
    newProductList = map(list, zip(*productList))

    # Merge product list if any isomers are present.
    finalProducts = []
    needsMerging = False
    for products in newProductList:
        if are_isomers(products):
            needsMerging = True
            break

    needsMerging = False
    if needsMerging:
        for products in newProductList:
            mi = products[0]
            for mj in products:
                mi = combine_stereoisomers(mi, mj)
                if mi is None:
                    return None
            finalProducts.append(Chem.MolToSmiles(mi, isomericSmiles=True))
        finalProductSets = ['.'.join(finalProducts)]

    else:
        finalProductSets = uniqueProductsList

    return finalProductSets


def uniquify(seq):
    """ remove duplicates from list """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]


def are_perms_equal_parity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length
    and have the same elements. No need to check for these
    conditions.
    """
    perm1 = perm1[:]  # copy this list so we don't mutate the original

    transCount = 0
    for loc in range(len(perm0) - 1):           # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0) + loc  # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1    # Swap in perm1
            transCount += 1

    # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False


def number_atoms(reactants):
    """Number atoms using the Isotope field."""
    atomIdx = 1
    for reactant in reactants:
        for atom in reactant.GetAtoms():
            atom.SetIsotope(atomIdx)
            atom.SetProp('Core', '')
            atomIdx += 1


def create_atom_map(mols):
    """Create atom map.

    Creates a dictionary where the keys are the unique atom number stored in the
    Isotope field and values are 2-element lists containing the molecule number
    and the atom index.
    """
    molIdx = 0
    atomMap = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            atomMap[atom.GetIsotope()] = [molIdx, atom.GetIdx()]
        molIdx += 1
    return atomMap


def create_template_map(reactants, rxn):
    """Generate atom  maps for reactant and product templates.

    Atom maps are stored in a list of dictionaries. Each element of the list
    corresponds to single match between reactants and reactant templates. The
    dictionary keys are the unique atom indices stored in the Isotope field.
    The dictionary values contain lists of two elements: 1) the template index
    (0...number of templates) and 2) the atom index (0...number of atoms in the
    template).  Given a unique atom index, the atom map data structure allows
    one to rapidly locate that atom within the reactant and product templates.
    NOTE: In the case of product templates only "mapped atoms" (atoms with atom
    map numbers in the templates) are included.
    """

    # Extract templates.
    rTemplates = [rxn.GetReactantTemplate(x)
                  for x in range(rxn.GetNumReactantTemplates())]
    pTemplates = [rxn.GetProductTemplate(x)
                  for x in range(rxn.GetNumProductTemplates())]

    # Generate substructure matches between reactants and reactant templates
    rMatches = []
    rMatchesChiral = []
    for rTemplate, reactant in itertools.izip(rTemplates, reactants):
        matches = reactant.GetSubstructMatches(rTemplate, uniquify=0)
        rMatches.append(matches)
        rMatchesChiral.append([(is_chiral_match(reactant, rTemplate, match)
                                and is_charge_match(reactant, rTemplate, match)) for match in matches])

    # Generate all permutations of matches.
    IsMatchChiral = []
    for matches in itertools.product(*rMatchesChiral):
        IsMatchChiral.append(False if False in matches else True)

    # Generate atom map for reactant templates.
    rTemplateAtomMaps = []
    rTemplateList = []
    matchIdx = 0
    for idx, matches in enumerate(itertools.product(*rMatches)):
        if IsMatchChiral[idx]:
            rTemplateAtomMaps.append({})
            rTemplateList.append(copy.deepcopy(rTemplates))  # make a copy (NOTE: DOESN'T COPY CHIRALITY)
            for rTemplateIdx, match in enumerate(matches):
                for rTemplateAtomIdx, rAtomIdx in enumerate(match):
                    rAtom = reactants[rTemplateIdx].GetAtomWithIdx(rAtomIdx)
                    rTemplateAtom = rTemplateList[matchIdx][rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx)
                    rTemplateAtom.SetIsotope(rAtom.GetIsotope())
                    rTemplateAtom.SetChiralTag(rTemplates[rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx).GetChiralTag())
                    rTemplateAtomMaps[matchIdx][rAtom.GetIsotope()] = [rTemplateIdx, rTemplateAtomIdx]
            matchIdx += 1

    # Generate atom map for product templates.
    pTemplateAtomMaps = []
    pTemplateList = []
    for matchIdx, rTemplates in enumerate(rTemplateList):
        pTemplateAtomMaps.append({})
        pTemplateList.append(list(pTemplates))  # make a copy
        for rTemplate in rTemplates:
            for rTemplateAtom in rTemplate.GetAtoms():
                if rTemplateAtom.HasProp('molAtomMapNumber'):
                    for pTemplateIdx, pTemplate in enumerate(pTemplateList[matchIdx]):
                        for pTemplateAtomIdx, pTemplateAtom in enumerate(pTemplate.GetAtoms()):
                            if (pTemplateAtom.HasProp('molAtomMapNumber') and
                                (pTemplateAtom.GetProp('molAtomMapNumber') ==
                                 rTemplateAtom.GetProp('molAtomMapNumber'))):
                                pTemplateAtom.SetIsotope(rTemplateAtom.GetIsotope())
                                pTemplateAtomMaps[matchIdx][rTemplateAtom.GetIsotope()] = [pTemplateIdx, pTemplateAtomIdx]

    return [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps, IsMatchChiral]


def is_charge_match(molecule, template, match):
    """Return True if a molecule and a template has the same charge.

    Given a molecule, a template, and a match between them (as generated by
    GetSubstructureMatch()), this function checks whether or not the template
    has the same charge as the molecule.
    """

    #Set charge of template and molecule atoms correctly.
    for atom in molecule.GetAtoms():
        smarts = re.sub('\$\(.*\)','',atom.GetSmarts())
        matched = re.search('([-+])(\d*)', smarts)
        if matched:
            #print matched.group(0)
            number = matched.group(2)
            if not number:
                number = '1'
            atom.SetFormalCharge(int(matched.group(1) + number))

    for rTAtom in template.GetAtoms():
        if rTAtom is not None:
            smarts = re.sub('\$\(.*\)','',rTAtom.GetSmarts())
            matched = re.search('([-+])(\d*)', smarts)
            if matched:
                number = matched.group(2)
                if not number:
                    number = '1'
                rTAtom.SetFormalCharge(int(matched.group(1) + number))

    # Initialize list for all atom.
    matchlist = []

    # Check atoms
    for atomIdx, templateAtom in enumerate(template.GetAtoms()):
        fragmentAtom = molecule.GetAtomWithIdx(match[atomIdx])
        if templateAtom.GetFormalCharge() == fragmentAtom.GetFormalCharge():
            is_match = 1
        else:
            is_match = 0
        matchlist.append(is_match)

    if matchlist == []:
        return False
    elif all([v == 1 for v in matchlist]):
        return True
    else:
        return False


def is_chiral_match(molecule, template, match):
    """Return True if a molecule and a template has the same chirality.

    Given a molecule, a template, and a match between them (as generated by
    GetSubstructureMatch()), this function checks whether or not the template
    has the same chirality as the molecule.

    Note
    This version works only with "development" (i.e. from github or svn repo) version
    of RDKit. Current stable release (2012_12_1) does not support stereochemistry in
    substructure matching.
    """

    # Check for stereochemistry
    fragmentHasStereochem = False
    templateHasStereochem = False
    # Check atoms
    for atomIdx, templateAtom in enumerate(template.GetAtoms()):
        fragmentAtom = molecule.GetAtomWithIdx(match[atomIdx])
        if templateAtom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
            templateHasStereochem = True
        if fragmentAtom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
            fragmentHasStereochem = True
    # Check bonds
    for templateBond in template.GetBonds():
        fragmentBond = molecule.GetBondBetweenAtoms(match[templateBond.GetBeginAtomIdx()],match[templateBond.GetEndAtomIdx()])
        if templateBond.GetBondDir() != BondDir.NONE:
            templateHasStereochem = True
        if fragmentBond.GetBondDir() != BondDir.NONE:
            fragmentHasStereochem = True
    if (not fragmentHasStereochem) or (not templateHasStereochem):
        return True

    #### Extract the match as a separate molecule
    # label the match atoms
    for atomIdx in match:
        molecule.GetAtomWithIdx(atomIdx).SetProp('match', '')
    # create an editable molecule
    fragmentMol = Chem.EditableMol(molecule)
    # delete all the unlabeled atoms
    atomIdxList = range(molecule.GetNumAtoms())
    for a in molecule.GetAtoms():
        if not a.HasProp('match'):
            Idx2 = next(i2 for i2, i1 in enumerate(atomIdxList) if i1 == a.GetIdx())
            fragmentMol.RemoveAtom(Idx2)
            del atomIdxList[Idx2]
    # convert to normal molecule
    fragmentMol = fragmentMol.GetMol()
    # Assign stereochemistry (DO NOT SANITIZE)
    Chem.AssignStereochemistry(fragmentMol)
    # remove labels from molecule
    molecule.ClearProp('match')

    #### Generate a molecule from the template
    templateMol = Chem.EditableMol(Chem.Mol())
    # Copy atoms and bonds
    for atomIdx, atom in enumerate(template.GetAtoms()):
        templateMol.AddAtom(molecule.GetAtomWithIdx(match[atomIdx]))
    for bondIdx, bond in enumerate(template.GetBonds()):
        templateMol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    # Convert to molecule
    templateMol = templateMol.GetMol()
    # Copy atom info
    for atomIdx, atom in enumerate(templateMol.GetAtoms()):
        oldAtom = template.GetAtomWithIdx(atomIdx)
        atom.SetChiralTag(oldAtom.GetChiralTag())
    # Copy bond info
    for bondIdx, bond in enumerate(templateMol.GetBonds()):
        oldBond = template.GetBondWithIdx(bondIdx)
        bond.SetBondType(oldBond.GetBondType())
        bond.SetBondDir(oldBond.GetBondDir())
    # Assign Stereochemistry (DO NOT SANITIZE!)
    Chem.AssignStereochemistry(templateMol)

    # Check for chiral match
    return templateMol.HasSubstructMatch(fragmentMol, useChirality=True)


def contains_any(str, set):
    for c in set:
        if c in str:
            return True
    return False


def merge(products):
    """Merge products sharing the same atoms.

    Given a list of products, this function looks through all pairs of products to
    see if they share common atoms (i.e., atoms with the same unique atom labels
    stored in the Isotope field).  If so, it merges these products into a single
    product.
    """
    for i in range(len(products)):
        for j in range(i + 1, len(products)):
            idxi = [atom.GetIsotope() for atom in products[i].GetAtoms()]
            idxj = [atom.GetIsotope() for atom in products[j].GetAtoms()]
            if len(set(idxi).intersection(set(idxj))) > 0:
                products[i] = combine_fragments(products[i], products[j])
                del products[j]
                return [i, products]
    return [-1, products]


def combine_fragments(m1, m2):
    """Stich two molecules together.

    Given two molecules that share common atoms (i.e., atoms with the same unique
    atom labels stored in the Isotope field), this function stitches them back
    together using common atoms and common bonds.
    """

    # build dictionary of atomidx and atommaps
    atomidx = {}
    for atom in m1.GetAtoms():
        atomidx[atom.GetIsotope()] = atom.GetIdx()

    # make new editable molecule from m1
    newmol = Chem.EditableMol(m1)

    # add atoms from m2
    for atom in m2.GetAtoms():
        if atom.GetIsotope() not in atomidx:
            # add the atom
            idx = newmol.AddAtom(atom)
            atomidx[atom.GetIsotope()] = idx

    # add bonds from m2
    m = newmol.GetMol()
    for bond in m2.GetBonds():
        atom1 = atomidx[bond.GetBeginAtom().GetIsotope()]
        atom2 = atomidx[bond.GetEndAtom().GetIsotope()]
        if m.GetBondBetweenAtoms(atom1, atom2) is None:
            newmol.AddBond(atom1, atom2, bond.GetBondType())

    m = newmol.GetMol()
    # Chem.SanitizeMol(m)
    return m


def fix_bonds(reactant, rTemplate, product):
    """Add missing bonds.

    Bonds which atoms belongs to different product templates are lost when
    RDKit's RunReactants() is used. It causes problems in ring opening
    reactions.

    Function iterates through the reactants bonds making sure that any bond
    which:
      1) is NOT already in a product,
      2) is NOT in a reactant template,
      3) bond atoms are mapped atoms,
    will be added to the product.
    """

    # Create an editible copy of the product.
    copy = Chem.EditableMol(product)

    # Make a map between reactant and product atoms utilizing unique numbers
    # stored in the isotope field of reactant atoms.
    pAtomMap = {}
    for pAtom in product.GetAtoms():
        for rAtom in reactant.GetAtoms():
            if pAtom.GetIsotope() == rAtom.GetIsotope():
                pAtomMap[rAtom.GetIsotope()] = pAtom.GetIdx()

    # Make a map between reactant and reactant template atoms for each match
    # of reactant and reactant template.
    tAtomMap = {}
    for rAtom in reactant.GetAtoms():
        for tAtom in rTemplate.GetAtoms():
            if rAtom.GetIsotope() == tAtom.GetIsotope():
                tAtomMap[rAtom.GetIsotope()] = tAtom.GetIdx()

    # Initialize a map between bonds and a pair of corresponding atoms of the
    # product.
    bondMap = {}

    # Iterate through the bonds of the reactant adding bonds which may be
    # missing in the product due to fact that their atoms belong to
    # different product templates.
    for rBond in reactant.GetBonds():
        # Find reactant atoms connected by the bond and their unique indices.
        rAtomi, rAtomj = rBond.GetBeginAtom(), rBond.GetEndAtom()
        i, j = rAtomi.GetIsotope(), rAtomj.GetIsotope()

        # Find corresponding atoms in the product. Proceed to next bond, if
        # not found.
        if i in pAtomMap.keys() and j in pAtomMap.keys():
            pAtomi = product.GetAtomWithIdx(pAtomMap[i])
            pAtomj = product.GetAtomWithIdx(pAtomMap[j])
        else:
            continue

        # If the bond already present in the product. continue with a next bond.
        if product.GetBondBetweenAtoms(pAtomi.GetIdx(), pAtomj.GetIdx()):
            #print "Bond already exists in the product."
            continue

        # Find corresponding atoms in the reactant template and check if they
        # are connected by the bond.
        isBond = False
        # Check if atom pair is present in reactant template match. If no,
        # proceed to a next one.
        isInMatch = True
        if i not in tAtomMap.keys() or j not in tAtomMap.keys():
            isInMatch = False
            continue

        # If a corresponding pair of atoms can be found in a reactant
        # template, check if it is connected by a bond.
        rTAtomi = rTemplate.GetAtomWithIdx(tAtomMap[i])
        rTAtomj = rTemplate.GetAtomWithIdx(tAtomMap[j])
        if rTemplate.GetBondBetweenAtoms(rTAtomi.GetIdx(), rTAtomj.GetIdx()):
            isBond = True
            continue

        # If we are here, we've found a missing bond. Add it to the product
        # and save corresponding product atom indecies.
        copy.AddBond(pAtomi.GetIdx(), pAtomj.GetIdx())
        bondMap[rBond] = (pAtomi.GetIdx(), pAtomj.GetIdx())

    # Copy properties of added bonds from the reactant.
    product = copy.GetMol()
    for rBond, pAtomIndexes in bondMap.items():
        addedBond = product.GetBondBetweenAtoms(*pAtomIndexes)
        addedBond.SetBondType(rBond.GetBondType())
        addedBond.SetBondDir(rBond.GetBondDir())
        addedBond.SetIsAromatic(rBond.GetIsAromatic())

    return product


def fix_charge(rAtom, pAtom, rTAtom, pTAtom):
    """Correct charge.

    Assign the appropriate charge to the product atom, given the
    information present in the reaction atom, the reactant template atom, and the
    product template atom.
    """
    #Set charge of template atoms correctly.
    if rTAtom is not None:
        smarts = re.sub('\$\(.*\)','',rTAtom.GetSmarts())
        match = re.search('([-+])(\d*)', smarts)
        if match:
            number = match.group(2)
            if not number:
                number = '1'
            rTAtom.SetFormalCharge(int(match.group(1) + number))

    if pTAtom is not None:
        smarts = re.sub('\$\(.*\)','',pTAtom.GetSmarts())
        match = re.search('([-+])(\d*)', smarts)
        if match:
            number = match.group(2)
            if not number:
                number = '1'
            pTAtom.SetFormalCharge(int(match.group(1) + number))

    if rAtom is None:
        if rTAtom is None:
            if pTAtom is None:
                # This case corresponds to an atom in the product that was introduced
                # by the template.  RDkit handles this appropriately, so we do nothing.
                return
            else:
                error('Should not be here (1) in fix_charge().')
                return None
        else:
            error('Should not be here (2) in fix_charge().')
            return None
    elif rAtom.GetFormalCharge() == 0:
        if rTAtom is None:
            # rAtom-Uncharged, rTAtom-None
            error('Should not be here (3) in fix_charge().')
            return None
        elif rTAtom.GetFormalCharge() == 0:
            if pTAtom is None:
                # rAtom-Uncharged, rTAtom-Uncharged, pTAtom-None
                error('Should not be here (4) in fix_charge().')
                return None
            elif pTAtom.GetFormalCharge() == 0:
                # rAtom-Uncharged, rTAtom-Uncharged, pTAtom-Uncharged --> UNSPECIFIED
                pAtom.SetFormalCharge(0)
                return
            else:
                # rAtom-Uncharged, rTAtom-Uncharged, pTAtom-Charged --> FORCE pTAtom charge
                pAtom.SetFormalCharge(pTAtom.GetFormalCharge())
                return
        else:  # rTAtom is charged
            # rAtom-Uncharged, rTAtom-Charged --> UNSPECIFIED
            pAtom.SetFormalCharge(0)
            return
    else:  # rAtom is charged
        if rTAtom is None:
            if pTAtom is None:
                # Preserve charge.
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                return
            elif pTAtom.GetFormalCharge() == 0:
                # Preserve charge.
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                return
            else:
                # Add charge.
                pAtom.SetFormalCharge(pTAtom.GetFormalCharge())
                pTAtom.SetFormalCharge(0)
                return
        elif rTAtom.GetFormalCharge() == 0:
            if pTAtom is None:
                # Preserve charge.
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                return
            elif(pTAtom.GetFormalCharge() == 0):
                # rAtom-Charged, rTAtom-Uncharged, pTAtom-Uncharged --> PRESERVE rAtom charge
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                return
            else:
                # rAtom-Charged, rTAtom-Uncharged, pTAtom-Charged --> FORCE pTAtom charge
                pAtom.SetFormalCharge(pTAtom.GetFormalCharge())
                pTAtom.SetFormalCharge(0)
                return
        else:  # rTAtom is charged
            if pTAtom is None:
                # Preserve charge.
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                return
            elif pTAtom.GetFormalCharge() == 0:
                # rAtom-Charged, rTAtom-Charged, pTAtom-Uncharged --> REMOVE charge
                pAtom.SetFormalCharge(0)
                return
            elif pTAtom.GetFormalCharge() == rTAtom.GetFormalCharge():  # same charge
                # Preserve charge.
                pAtom.SetFormalCharge(rAtom.GetFormalCharge())
                pTAtom.SetFormalCharge(0)
                rTAtom.SetFormalCharge(0)
                return
            else:
                #Don't know yet what should be possible behavior. Different charge.
                pAtom.SetFormalCharge(pTAtom.GetFormalCharge())
                pTAtom.SetFormalCharge(0)
                rTAtom.SetFormalCharge(0)
                return


def fix_chirality(rAtom, pAtom, rTAtom, pTAtom):
    """Correct chirality.

    Assign the appropriate tetrahedral chirality to the product atom, given the
    information present in the reaction atom, the reactant template atom, and the
    product template atom.
    """
    if rAtom is None:
        if rTAtom is None:
            if pTAtom is None:
                # This case corresponds to an atom in the product that was introduced
                # by the template.  RDkit handles this appropriately, so we do nothing.
                return
            else:
                error('Error: Should not be here (1) in fix_chirality().')
                return None
        else:
            error('Error: Should not be here (2) in fix_chirality().')
            return None
    elif rAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
        if rTAtom is None:
            # rAtom-Achiral, rTAtom-None
            error('Error: Should not be here (3) in fix_chirality().')
            return None
        elif rTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            if pTAtom is None:
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-None
                error('Error: Should not be here (4) in fix_chirality().')
                return None
            elif(pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED):
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-Achiral --> UNSPECIFIED
                pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                return
            else:
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom CHIRALITY
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        else:  # rTAtom is chiral
            # rAtom-Achiral, rTAtom-Chiral --> UNSPECIFIED
            pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
            return
    else:  # rAtom is chiral
        if rTAtom is None:
            if pTAtom is None:  # 19
                # Preserve chirality.
                copy_chirality(rAtom, pAtom)
                return
            elif pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                # Preserve chirality.
                copy_chirality(rAtom, pAtom)
                return
            else:  # 21
                # Add chirality.
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        elif rTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            if pTAtom is None:
                # Preserve chirality.
                copy_chirality(rAtom, pAtom)
                return
            elif(pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED):
                # rAtom-Chiral, rTAtom-Achiral, pTAtom-Achiral --> PRESERVE rAtom CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            else:
                # rAtom-Chiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom CHIRALITY
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        else:  # rTAtom is chiral
            if pTAtom is None:
                # Preserve chirality.
                copy_chirality(rAtom, pAtom)
                return
            elif pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                # rAtom-Chiral, rTAtom-Chiral, pTAtom-Achiral --> REMOVE CHIRALITY
                pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                return
            elif pTAtom.GetChiralTag() == rTAtom.GetChiralTag():  # same chirality
                # Preserve chirality.
                copy_chirality(rAtom, pAtom)
                return
            else:
                # Different chiralitiy, invert it.
                copy_chirality(rAtom, pAtom)
                pAtom.InvertChirality()
                return


def copy_chirality(a1, a2):
    a2.SetChiralTag(a1.GetChiralTag())

    neighbors1 = [atom.GetIsotope() for atom in a1.GetNeighbors()]
    neighbors2 = [atom.GetIsotope() for atom in a2.GetNeighbors()]

    # Check if chirality was lost
    if len(neighbors2) < 3:
        return

    # Add Hydrogens
    while len(neighbors1) < 4:
        neighbors1.append(-1)
    while len(neighbors2) < 4:
        neighbors2.append(-1)

    # Check for ambiguous chirality.
    d1 = list(set(neighbors1) - set(neighbors2))
    if len(d1) > 1:
        return
    elif len(d1) == 1:
        d2 = list(set(neighbors2) - set(neighbors1))
        neighbors2 = [d1[0] if x == d2[0] else x for x in neighbors2]

    # Check Parity of Permutation
    if not are_perms_equal_parity(neighbors1, neighbors2):
        a2.InvertChirality()


def copy_chirality_Hs(a1, a2):
    a2.SetChiralTag(a1.GetChiralTag())

    neighbors1 = [(atom.GetIsotope() if ((atom.GetAtomicNum() != 1) and (atom.GetIsotope() != 0)) else atom.GetSymbol()) for atom in a1.GetNeighbors()]
    neighbors2 = [(atom.GetIsotope() if ((atom.GetAtomicNum() != 1) and (atom.GetIsotope() != 0)) else atom.GetSymbol()) for atom in a2.GetNeighbors()]

    # Add Hydrogens
    while len(neighbors1) < 4:
        neighbors1.append('H')
    while len(neighbors2) < 4:
        neighbors2.append('H')

    # Check for ambiguous chirality.
    d1 = list(set(neighbors1) - set(neighbors2))
    if len(d1) > 1:
        return
    elif len(d1) == 1:
        d2 = list(set(neighbors2) - set(neighbors1))
        neighbors2 = [d1[0] if x == d2[0] else x for x in neighbors2]

    # Check Parity of Permutation
    if not are_perms_equal_parity(neighbors1, neighbors2):
        a2.InvertChirality()


def add_Hs(mol):
    """Add explicit hydrogens.

    Returns a new molecule with H atoms added with atoms with implicit valence.
    """
    # Get a list of implicit H atoms to be added for each atom in the molecule.
    implicitH = []
    for atom in mol.GetAtoms():
        implicitValence = atom.GetImplicitValence()
        numExplicitHs = atom.GetNumExplicitHs()
        if numExplicitHs > 0:
            implicitValence += numExplicitHs
            atom.SetNumExplicitHs(0)
        implicitH.append(implicitValence)

    # Check to see if we need to do anything
    if implicitH == []:
        return mol
    elif all([v == 0 for v in implicitH]):
        return mol

    #Define H atom.
    mol_H = Chem.MolFromSmiles('[H]')
    newatom = mol_H.GetAtoms()[0]

    # Define an editable molecule similar to given molecule (Don't just copy
    # because then we are unable to change stereochemistry information)
    rdemol = Chem.EditableMol(Chem.Mol())
    for atom in mol.GetAtoms():
        rdemol.AddAtom(atom)
    for bond in mol.GetBonds():
        rdemol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Add H atom in place of implicit H atoms in the editable molecule.
    for idx, i in enumerate(implicitH):
        if i != 0:
            for j in range(i):
                newatomidx = rdemol.AddAtom(newatom)
                rdemol.AddBond(idx, newatomidx)

    # Convert editable molecule to normal molecule.
    rdmol = rdemol.GetMol()
    for atom in rdmol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atom.GetBonds()[0].SetBondType(BondType.SINGLE)
    # Add Chirality and aromaticity information
    for atomIdx, atom in enumerate(mol.GetAtoms()):
        rdmol.GetAtomWithIdx(atomIdx).SetChiralTag(atom.GetChiralTag())
    # Add Bond information
    for bond in mol.GetBonds():
        newBond = rdmol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        newBond.SetBondType(bond.GetBondType())
        newBond.SetBondDir(bond.GetBondDir())
        newBond.SetIsAromatic(bond.GetIsAromatic())
    try:
        Chem.SanitizeMol(rdmol)
        return rdmol
    except:
        error('Invalid molecule after adding explicit Hs.\n')
        return None


def remove_Hs(mol):
    """Remove explicit hydrogens from a molecule.

    Returns a new molecule after removing all the explicit H atoms present in
    the molecule and keeping stereo-chemistry of molecule preserved.
    """
    # Add stereobonds, but not to hydrogen atoms.
    for bondIdx, bond in enumerate(mol.GetBonds()):
        if bond.GetBondDir() != Chem.BondDir.NONE:
            # check for hydrogen
            if bond.GetBeginAtom().GetAtomicNum() == 1:
                otherAtom = bond.GetEndAtom()
            elif bond.GetEndAtom().GetAtomicNum() == 1:
                otherAtom = bond.GetBeginAtom()
            else:
                mol.GetBondWithIdx(bondIdx).SetBondDir(bond.GetBondDir())
                continue

            # Move bond directionality to the non-hydrogen atom
            for otherBond in otherAtom.GetBonds():
                if ((otherBond.GetBondType() == BondType.SINGLE)
                    and (otherBond.GetOtherAtom(otherAtom).GetAtomicNum() != 1)):
                        if bond.GetBondDir() == BondDir.ENDUPRIGHT:
                            mol.GetBondWithIdx(otherBond.GetIdx()).SetBondDir(BondDir.ENDDOWNRIGHT)
                        elif bond.GetBondDir() == BondDir.ENDDOWNRIGHT:
                            mol.GetBondWithIdx(otherBond.GetIdx()).SetBondDir(BondDir.ENDUPRIGHT)
                        break

    # Create a new molecule from scrach
    rdemol = Chem.EditableMol(Chem.Mol())

    # add all non-hydrogen atoms
    atomMap = {}
    for atomIdx, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() != 1:
            atomMap[atomIdx] = rdemol.AddAtom(atom)
    for bond in mol.GetBonds():
        if ((bond.GetBeginAtom().GetAtomicNum() != 1) and
            (bond.GetEndAtom().GetAtomicNum() != 1)):
            rdemol.AddBond(atomMap[bond.GetBeginAtomIdx()], atomMap[bond.GetEndAtomIdx()])

    # Convert editable molecule to normal molecule.
    rdmol = rdemol.GetMol()
    # Add Chirality information
    for atomIdx, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() != 1:
            if atom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                copy_chirality_Hs(atom, rdmol.GetAtomWithIdx(atomMap[atomIdx]))
                if len(atom.GetNeighbors()) == 3:
                    atom.SetNumExplicitHs(1)
    # Add Bond information
    for bond in mol.GetBonds():
        if ((bond.GetBeginAtom().GetAtomicNum() != 1) and
            (bond.GetEndAtom().GetAtomicNum() != 1)):
            newBond = rdmol.GetBondBetweenAtoms(atomMap[bond.GetBeginAtomIdx()], atomMap[bond.GetEndAtomIdx()])
            newBond.SetBondType(bond.GetBondType())
            newBond.SetBondDir(bond.GetBondDir())
            newBond.SetIsAromatic(bond.GetIsAromatic())

    try:
        Chem.SanitizeMol(rdmol)
        return rdmol
    except:
        error('Invalid molecule after removing explicit Hs.')
        return None


def are_isomers(mols):
    """Returns True if a given list of molecules are isomers and vice versa."""

    # We need at least two moleacules to talk about isomers.
    if len(mols) == 1:
        return False

    mi = mols[0]
    if all(mi.HasSubstructMatch(mj) and mj.HasSubstructMatch(mi) for mj in mols):
        smis = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
        return smis[0] != smis[1]
    else:
        return False


def combine_stereoisomers(mol1, mol2):
    """
    Given two stereoisomers,returns a new molecule with ambiguous chirality by merging their stereoisomerism"""
    Chem.AssignStereochemistry(mol1)
    Chem.AssignStereochemistry(mol2)

    #Label atoms with different chirality.
    for a in mol1.GetAtoms():
        b = mol2.GetAtomWithIdx(a.GetIdx())
        if a.GetChiralTag() != b.GetChiralTag():
            a.SetProp('removeChiralTag', '')

    # Label bonds with same chirality.
    for a in mol1.GetBonds():
        b = mol2.GetBondWithIdx(a.GetIdx())
        if (((a.GetStereo() == BondStereo.STEREOE) and (b.GetStereo() == BondStereo.STEREOE))
            or ((a.GetStereo() == BondStereo.STEREOZ) and (b.GetStereo() == BondStereo.STEREOZ))):
            a.SetProp('sameStereo', '')

    # Label bond directions for removal.
    for a in mol1.GetBonds():
        if a.GetBondDir() != BondDir.NONE:
            removeBond = True
            # check it's neighbor bonds
            for b in a.GetBeginAtom().GetBonds():
                if b.HasProp('sameStereo'):
                    removeBond = False
                    break
            for b in a.GetEndAtom().GetBonds():
                if b.HasProp('sameStereo'):
                    removeBond = False
                    break
            if removeBond: a.SetProp('removeBondDir','')

    # Create an editable molecule from scratch.
    rdemol = Chem.EditableMol(Chem.Mol())

    # Copy atoms.
    for atom in mol1.GetAtoms():
        rdemol.AddAtom(atom)

    # Copy bonds.
    for bond in mol1.GetBonds():
        rdemol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Convert into normal molecule.
    rdmol = rdemol.GetMol()

    # Copy chiral information of atoms and remove chiral information of labeled atoms.
    for atomIdx, atom in enumerate(mol1.GetAtoms()):
        if atom.HasProp('removeChiralTag'):
            rdmol.GetAtomWithIdx(atomIdx).SetChiralTag(ChiralType.CHI_UNSPECIFIED)
        else:
            rdmol.GetAtomWithIdx(atomIdx).SetChiralTag(atom.GetChiralTag())

    # Copy bond type and directionality and remove bond directionality of labeled bonds.
    for bond in mol1.GetBonds():
        newbond = rdmol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        newbond.SetBondType(bond.GetBondType())
        newbond.SetIsAromatic(bond.GetIsAromatic())
        if bond.HasProp('removeBondDir'):
            newbond.SetBondDir(Chem.BondDir.NONE)
        else:
            newbond.SetBondDir(bond.GetBondDir())

    # Remove atom and bond label.
    for atom in rdmol.GetAtoms():
        atom.ClearProp('removeChiralTag')
    for bond in rdmol.GetBonds():
        bond.ClearProp('removeBondDir')
        bond.ClearProp('sameStereo')

    try:
        Chem.SanitizeMol(rdmol)
        return rdmol
    except:
        error('Output molecule of combine_stereoisomers() is invalid\n')
        return None


if __name__ == '__main__':
    print run_reactants.__doc__
