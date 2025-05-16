from rdkit import Chem
from rdkit.Chem import Draw

def count_quaternary_carbons(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, []
    
    quaternary_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            num_carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if num_carbon_neighbors == 4:
                quaternary_carbons.append(atom.GetIdx())
    
    return len(quaternary_carbons), quaternary_carbons

if __name__ == "__main__":
    SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"

    num_quaternary_carbons, quaternary_carbons = count_quaternary_carbons(SMILES)
    print(f"Number of quaternary carbons: {num_quaternary_carbons}")
    print(f"Quaternary carbon indices: {quaternary_carbons}")
