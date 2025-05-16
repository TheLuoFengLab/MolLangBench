from rdkit import Chem
from rdkit.Chem import Draw
import random
import re

def check_chirality_center(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, 0, False
    
    # Find all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
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

def check_bond_stereo(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, 0, False
    
    # Find all stereo bonds
    stereo_bonds = []
    for bond in mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            stereo_bonds.append(bond)
    
    if stereo_bonds:
        # First check for non-carbon-carbon stereo bonds
        non_cc_stereo_bonds = []
        for bond in stereo_bonds:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if not (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6):
                non_cc_stereo_bonds.append(bond)
        
        if non_cc_stereo_bonds:
            # Choose a non-carbon-carbon stereo bond
            bond = random.choice(non_cc_stereo_bonds)
            is_cc_bond = False
        else:
            # Randomly choose a stereo bond
            bond = random.choice(stereo_bonds)
            is_cc_bond = True
        
        # Get the bond configuration
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            config = 1  # E configuration
        elif stereo == Chem.BondStereo.STEREOZ:
            config = 2  # Z configuration
        else:
            config = 0  # Unknown configuration
        
        return bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), config, is_cc_bond
    else:
        # Find potential stereo bonds (double bonds with different substituents)
        potential_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                
                # Check if both atoms have at least one other neighbor
                if len(begin_atom.GetNeighbors()) > 1 and len(end_atom.GetNeighbors()) > 1:
                    # Get the other neighbors
                    begin_neighbors = [n for n in begin_atom.GetNeighbors() if n.GetIdx() != end_atom.GetIdx()]
                    end_neighbors = [n for n in end_atom.GetNeighbors() if n.GetIdx() != begin_atom.GetIdx()]
                    
                    # Check if the neighbors are different
                    if (len(begin_neighbors) > 0 and len(end_neighbors) > 0 and
                        begin_neighbors[0].GetSymbol() != end_neighbors[0].GetSymbol()):
                        potential_bonds.append(bond)
        
        if potential_bonds:
            # Randomly choose a potential stereo bond
            bond = random.choice(potential_bonds)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            is_cc_bond = (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6)
            return bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), 0, is_cc_bond
        else:
            return None, None, 0, False

if __name__ == "__main__":
    # Example with E/Z configuration
    SMILES = "C/C=C/C"  # E configuration
    # SMILES = "C/C=C\C"  # Z configuration
    # SMILES = "C/C=N/O"  # Non-carbon-carbon stereo bond

    begin_idx, end_idx, config, is_cc_bond = check_bond_stereo(SMILES)
    if begin_idx is not None:
        print(f"Stereo bond atom indices: {begin_idx}, {end_idx}")
        print(f"Configuration (1=E, 2=Z, 0=no stereo or potential stereo): {config}")
        print(f"Is carbon-carbon bond: {is_cc_bond}")
        # Print the atom symbols
        mol = Chem.MolFromSmiles(SMILES)
        begin_atom = mol.GetAtomWithIdx(begin_idx)
        end_atom = mol.GetAtomWithIdx(end_idx)
        print(f"Bond: {begin_atom.GetSymbol()}-{end_atom.GetSymbol()}")
        # Check if this is actually a stereo bond
        bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
        is_stereo = bond.GetStereo() != Chem.BondStereo.STEREONONE
        print(f"Is stereo bond: {is_stereo}")
    else:
        print("No stereo bonds or potential stereo bonds found")
