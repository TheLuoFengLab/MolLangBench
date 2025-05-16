from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random
from rdkit.Chem import Fragments

def get_amide(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # get the benzene rings
    num_amide = Fragments.fr_amide(mol)

    # get the benzene rings from the SMARTS pattern
    SMARTS_pattern = Chem.MolFromSmarts("C(=O)-N")
    matches = mol.GetSubstructMatches(SMARTS_pattern)
    num_amide_from_pattern = len(matches)
    assert num_amide == num_amide_from_pattern, f"Number of amide rings from SMARTS pattern: {num_amide_from_pattern}, Number of amide rings from rdkit: {num_amide}"
    
    return num_amide, matches

if __name__ == "__main__":
    #SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"
    SMILES = "c1c(CN=C(COc2ccccc2C=O)O)ccc(OCC(C)C)c1"
    print(get_amide(SMILES))
