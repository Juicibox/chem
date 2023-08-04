import streamlit as st

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

st.title('RDKit + Py3DMOL 😀')

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock


compound_smiles=st.text_input('SMILES please','CC')
blks=makeblocks(compound_smiles)
st.writer(blks)
