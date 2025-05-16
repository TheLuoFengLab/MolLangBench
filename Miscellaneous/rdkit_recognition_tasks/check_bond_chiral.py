from rdkit import Chem
from rdkit.Chem import Draw
import random
import re

def check_chirality_center(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, 0, False
    
    # Find all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    
    # Convert R/S to 1/2
    config_map = {'R': 1, 'S': 2}
    
    if chiral_centers:
        # 20% chance to choose an atom adjacent to a chiral center in SMILES
        if random.random() < 0.2:
            # Find atoms with @ in their SMILES representation
            pattern = r'\[[^]]*@[^]]*\]'
            matches = re.finditer(pattern, smiles)
            chiral_atoms_in_smiles = []
            for match in matches:
                # Get the atom index from the SMILES string
                atom_str = match.group()
                # Find the corresponding atom in the molecule
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() in atom_str:
                        chiral_atoms_in_smiles.append(atom.GetIdx())
            
            if chiral_atoms_in_smiles:
                # Get all carbon atoms adjacent to chiral centers
                adjacent_carbons = set()
                for chiral_idx in chiral_atoms_in_smiles:
                    atom = mol.GetAtomWithIdx(chiral_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6:  # Only consider carbon atoms
                            adjacent_carbons.add(neighbor.GetIdx())
                
                if adjacent_carbons:
                    # Choose a random adjacent carbon
                    center_idx = random.choice(list(adjacent_carbons))
                    return center_idx, 4, True  # Always carbon for adjacent atoms
        
        # If we didn't choose an adjacent atom, proceed with normal chiral center selection
        # First check for non-carbon chiral centers
        non_carbon_centers = [(idx, config) for idx, config in chiral_centers 
                            if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        
        if non_carbon_centers:
            # Choose a non-carbon chiral center
            center_idx, config = random.choice(non_carbon_centers)
            return center_idx, config_map.get(config, 0), False
        else:
            # Randomly sample a carbon chiral center
            center_idx, config = random.choice(chiral_centers)
            return center_idx, config_map.get(config, 0), True
    else:
        # Find potential chiral centers (carbons with 4 different substituents)
        potential_centers = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 4:  # Tetrahedral carbon
                    # Check if all neighbors are different
                    neighbor_symbols = [n.GetSymbol() for n in neighbors]
                    if len(set(neighbor_symbols)) == 4:
                        potential_centers.append(atom.GetIdx())
        
        if potential_centers:
            # Randomly sample a potential center
            center_idx = random.choice(potential_centers)
            return center_idx, 0, True  # Always carbon for potential centers
        else:
            return None, 0, False

if __name__ == "__main__":
    SMILES = "[C@]12(C)CCC[C@@]34[C@@H]5[C@H]6[C@@]7([C@H](O)C(=C)[C@H]([C@H](OC(c8ccccc8)=O)[C@H]73)C6)C[C@@](OC(=O)c3ccccc3)(N5C1)[C@H]24"

    center_idx, config, is_carbon = check_chirality_center(SMILES)
    if center_idx is not None:
        print(f"Chiral center atom index: {center_idx}")
        print(f"Configuration (1=R, 2=S, 0=no chirality, 4=no chiral but adjacent to chiral center): {config}")
        print(f"Is carbon: {is_carbon}")
        # Print the atom symbol and whether it's a chiral center
        mol = Chem.MolFromSmiles(SMILES)
        atom = mol.GetAtomWithIdx(center_idx)
        print(f"Atom symbol: {atom.GetSymbol()}")
        # Check if this atom is actually a chiral center
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        is_chiral = any(idx == center_idx for idx, _ in chiral_centers)
        print(f"Is chiral center: {is_chiral}")
    else:
        print("No chiral centers or potential chiral centers found")
