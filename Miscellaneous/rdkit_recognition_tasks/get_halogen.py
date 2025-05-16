from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_halogen(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_halogen = Fragments.fr_halogen(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("[#9,#17,#35,#53]")
    matches = mol.GetSubstructMatches(SMARTS_pattern)
    num_halogen_from_pattern = len(matches)
    assert num_halogen == num_halogen_from_pattern, f"Number of halogen rings from SMARTS pattern: {num_halogen_from_pattern}, Number of halogen rings from rdkit: {num_halogen}"
    
    return num_halogen, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "C(Br)(CC[C@@H](C(F)(F)F)O[Si](c1ccccc1)(c1ccccc1)C(C)(C)C)Br"
    print(get_halogen(SMILES))
