from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_ester(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_ester = Fragments.fr_ester(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2H0][#6]")
    matches = mol.GetSubstructMatches(SMARTS_pattern)
    num_ester_from_pattern = len(matches)
    assert num_ester == num_ester_from_pattern, f"Number of ester rings from SMARTS pattern: {num_ester_from_pattern}, Number of ester rings from rdkit: {num_ester}"
    
    return num_ester, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "c1c(CN=C(COc2ccccc2C=O)O)ccc(OCC(C)C)c1"
    print(get_ester(SMILES))
