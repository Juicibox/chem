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
img = Draw.MolToImage(mol)
name = "mole"

AllChem.ComputeGasteigerCharges(mol)
contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
contribs1 = rdMolDescriptors._CalcCrippenContribs(mol)
fig1 = SimilarityMaps.GetSimilarityMapFromWeights(mol, [x for x, y in contribs1], colorMap='jet',contourLines=10)


col1, col2 = st.columns(2)

with col1:

  st.image(img, use_column_width=True, caption=f'Imagen en 2D de la estructura de la molécula {name}')
 
  st.markdown(f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Imagen de las cargas parciales de Gasteiger para la molécula {name}</h2>", unsafe_allow_html=True)

  st.markdown(
      f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Imagen de las contribuciones de Crippen para la solubilidad (logP) para la molécula {name}</h2>",
      unsafe_allow_html=True)

with col2:
  mol_wt = round(Descriptors.MolWt(mol), 4)
  st.write(f'Peso molecular: {mol_wt}')



