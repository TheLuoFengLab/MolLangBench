from rdkit import Chem
from rdkit.Chem import Draw

def count_ring_junction_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    ring_info = mol.GetRingInfo()

    bridgehead_atoms = []
    for atom in mol.GetAtoms():
        # NumAtomRings returns how many rings an atom is in
        if ring_info.NumAtomRings(atom.GetIdx()) > 1:
            bridgehead_atoms.append(atom.GetIdx())
    
    return len(bridgehead_atoms), bridgehead_atoms

if __name__ == "__main__":
    SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"

    num_ring_junction_atoms, ring_junction_atoms = count_ring_junction_atoms(SMILES)
    print(f"Number of ring junction atoms: {num_ring_junction_atoms}")
    print(f"Ring junction atoms: {ring_junction_atoms}")
