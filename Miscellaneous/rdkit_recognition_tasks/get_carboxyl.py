from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_carboxyl(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_carboxyl = Fragments.fr_COO(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("[#6]C(=O)[O;H]") # we filter out the O-
    matches = mol.GetSubstructMatches(SMARTS_pattern)

    SMARTS_pattern_1 = Chem.MolFromSmarts("C(=O)[O;H]") # we filter out the O-
    matches_1 = mol.GetSubstructMatches(SMARTS_pattern_1)

    if len(matches) != len(matches_1): # filter out the carboxylic acids not attached to a carbon
        print(f"Number of carboxylic acids from the SMARTS pattern: {len(matches)}, from the SMARTS pattern_1: {len(matches_1)}, SMILES: {smiles}")
        return None, None

    
    num_carboxyl_from_pattern = len(matches)
    #assert num_carboxyl == num_carboxyl_from_pattern, f"Number of carboxylic acids from the function: {num_carboxyl}, from the SMARTS pattern: {num_carboxyl_from_pattern}, SMILES: {smiles}"
    if num_carboxyl != num_carboxyl_from_pattern:
        print(f"Number of carboxylic acids from the function: {num_carboxyl}, from the SMARTS pattern: {num_carboxyl_from_pattern}, SMILES: {smiles}")
    
        return None, None
    
    #print(matches)
    new_matches = []
    for match in matches:
        # for each match, there are two carbon and two oxygen atoms
        # we keep two oxygen atoms and one carbon atom that is connected to all other atoms
        new_match = []
        #print(f"match: {match}")
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "O":
                new_match.append(atom_idx)
            elif atom.GetSymbol() == "C":
                # check the connected atom indices should be within the match
                connected_atom_indices = [bond.GetOtherAtomIdx(atom_idx) for bond in atom.GetBonds()]
                if all([idx in match for idx in connected_atom_indices]):
                    new_match.append(atom_idx)
                else:
                    pass
            else:
                raise ValueError(f"Unexpected atom symbol: {atom.GetSymbol()}")
        
        #print(f"new_match: {new_match}")
        new_matches.append(new_match)
        
    # convert the list of lists to a tuple of tuples
    matches = tuple([tuple(match) for match in new_matches])
    #print(f"matches: {matches}")
        
        
    
    return num_carboxyl, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "C1C(CCC(C(O)=O)C(O)=O)CCC(C)C1CCC(C(O)=O)C(=O)O"
    print(get_carboxyl(SMILES))
