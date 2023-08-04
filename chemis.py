import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski, Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem

from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

st.title('RDKit + Py3DMOL 😀')


smiles=st.text_input('SMILES please','CC')
mol = Chem.MolFromSmiles(smiles)

mol_wt = round(Descriptors.MolWt(mol), 4)
st.write(f'Peso molecular: {mol_wt}')
