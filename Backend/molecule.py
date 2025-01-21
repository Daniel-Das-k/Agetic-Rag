from rdkit import Chem
from rdkit.Chem import Draw

# SMILES for Paclitaxel
smiles = "C1CCC(CC1)C(=O)OC2=CC=CC=C2"
# Parse the SMILES string
mol = Chem.MolFromSmiles(smiles)

# Generate the depiction
Draw.MolToFile(mol, "Paclitaxel.png")

