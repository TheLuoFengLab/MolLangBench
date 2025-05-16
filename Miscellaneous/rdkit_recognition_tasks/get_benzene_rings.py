from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_benzene_rings(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_benzene_rings = Fragments.fr_benzene(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("c1ccccc1")
    matches = mol.GetSubstructMatches(SMARTS_pattern)
    num_benzene_rings_from_pattern = len(matches)
    assert num_benzene_rings == num_benzene_rings_from_pattern, f"Number of benzene rings from SMARTS pattern: {num_benzene_rings_from_pattern}, Number of benzene rings from rdkit: {num_benzene_rings}"
    
    return num_benzene_rings, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "c1c(CN=C(COc2ccccc2C=O)O)ccc(OCC(C)C)c1"
    print(get_benzene_rings(SMILES))
