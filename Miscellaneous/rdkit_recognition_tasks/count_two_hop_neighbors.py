from rdkit import Chem
from rdkit.Chem import Draw
import random

def get_n_hop_neighbors(smiles, n):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    atoms = mol.GetAtoms()
    if not atoms:
        raise ValueError("No atoms found in the molecule")

    selected_atom = random.choice(atoms)
    selected_idx = selected_atom.GetIdx()
    
    current_level = {selected_idx}
    visited = set(current_level)
    
    for _ in range(n):
        next_level = set()
        for idx in current_level:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    next_level.add(neighbor_idx)
        
        if not next_level:
            break
        
        visited.update(next_level)
        current_level = next_level
    
    n_hop_neighbors = list(next_level)
    #n_hop_neighbors = list(visited - {selected_idx})

    return selected_idx, len(n_hop_neighbors), n_hop_neighbors
    

if __name__ == "__main__":
    SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"

    selected_atom_idx, num_two_hop_neighbors, two_hop_neighbors = get_n_hop_neighbors(SMILES, 2)
    print(f"Selected atom index: {selected_atom_idx}")
    print(f"Number of two-hop neighbors: {num_two_hop_neighbors}")
    print(f"Two-hop neighbors: {two_hop_neighbors}")