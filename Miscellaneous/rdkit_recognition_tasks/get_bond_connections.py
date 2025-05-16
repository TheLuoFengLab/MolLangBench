from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import random

def get_bond_connections(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    bond_info = mol.GetBonds()

    bond_connections = []
    for bond in bond_info:
        temp = dict()
        temp["bond_atoms"] = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        temp["bond_type"] = bond.GetBondType()
        temp['idx_distance'] = abs(bond.GetBeginAtomIdx() - bond.GetEndAtomIdx())
        bond_connections.append(temp)
    
    df = pd.DataFrame(bond_connections)

    num_atoms = mol.GetNumAtoms()
    # we can check if any of pair of atoms are idx distance of 1, but they are not connected
    bond_atoms = df['bond_atoms'].tolist()
    # convert each tuple to a set
    bond_atoms = [set(i) for i in bond_atoms]
    unconnected_atoms = []
    for i in range(num_atoms):
        if {i, i+1} not in bond_atoms:
            temp = dict()
            temp["bond_atoms"] = (i, i+1)
            temp["bond_type"] = 0
            temp['idx_distance'] = 1
            unconnected_atoms.append(temp)

    df_unconnected = pd.DataFrame(unconnected_atoms)
    
    # randomly choose two atoms that are not connected, and add them to the dataframe
    random_unconnected_atoms = []
    for _ in range(2):
        random_pair = random.sample(range(num_atoms), 2)
        random_pair = set(random_pair)
        if random_pair not in bond_atoms:
            temp = dict()
            temp["bond_atoms"] = random_pair
            temp["bond_type"] = 0
            list_random_pair = list(random_pair)
            temp['idx_distance'] = abs(list_random_pair[0] - list_random_pair[1])
            random_unconnected_atoms.append(temp)

    df_random = pd.DataFrame(random_unconnected_atoms)

    df = pd.concat([df, df_random, df_unconnected], ignore_index=True)
    
    return df

if __name__ == "__main__":
    SMILES = "C[C@H]1C[C@H]2OC(=O)C3=CCC[C@@H]([C@]1(C)CCC1=CC(=O)OC1)[C@]32C"

    df = get_bond_connections(SMILES)
    print(df)
