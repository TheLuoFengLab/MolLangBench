from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_thiophene(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_thiophene = Fragments.fr_thiophene(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("s1cccc1")
    matches = mol.GetSubstructMatches(SMARTS_pattern)
    num_thiophene_from_pattern = len(matches)
    assert num_thiophene == num_thiophene_from_pattern, f"Number of thiophene rings from SMARTS pattern: {num_thiophene_from_pattern}, Number of thiophene rings from rdkit: {num_thiophene}"
    
    return num_thiophene, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "c1c(CN=C(COc2ccccc2C=O)O)ccc(OCC(C)C)c1"
    print(get_thiophene(SMILES))
