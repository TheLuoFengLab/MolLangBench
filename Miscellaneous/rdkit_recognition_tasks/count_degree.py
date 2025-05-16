from rdkit import Chem
from rdkit.Chem import Draw
import random

def evaluate_random_atom_degree(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    
    atoms = mol.GetAtoms()
    if not atoms:
        raise ValueError("No atoms found in the molecule")
    
    selected_atom = random.choice(atoms)
    selected_idx = selected_atom.GetIdx()
    bonded_atom_indices = [bonded_atom.GetIdx() for bonded_atom in selected_atom.GetNeighbors()]
    degree = selected_atom.GetDegree()
    
    return selected_atom, degree, bonded_atom_indices

if __name__ == "__main__":
    SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"

    selected_atom, degree, bonded_atom_indices = evaluate_random_atom_degree(SMILES)
    print(f"Selected atom: {selected_atom.GetSymbol()}: {selected_atom.GetIdx()}")
    print(f"Degree of the selected atom: {degree}")
    print(f"Bonded atom indices: {bonded_atom_indices}")

