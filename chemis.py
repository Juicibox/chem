import streamlit as st
import py3Dmol
import pubchempy as pcp
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



col1, col2 = st.columns(2)

with col1:

  st.image(img, use_column_width=True, caption=f'Imagen en 2D de la estructura de la molécula {name}')
 

with col2:
  mol_wt = round(Descriptors.MolWt(mol), 4)
  st.write(f'Peso molecular: {mol_wt}')
  
  num_atoms = mol.GetNumAtoms()
  st.write(f'Número de átomos: {num_atoms}')

  num_bonds = mol.GetNumBonds()
  st.write(f'Número de enlaces: {num_bonds}')

  num_acceptors = Lipinski.NumHAcceptors(mol)
  st.write(f'Número de aceptores de enlaces de hidrógeno: {num_acceptors}')

  num_donors = Lipinski.NumHDonors(mol)
  st.write(f'Número de donantes de enlaces de hidrógeno: {num_donors}')

  num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
  st.write(f'Número de enlaces giratorios: {num_rotatable_bonds}')

  num_aliphatic_rings = Lipinski.NumAliphaticRings(mol)
  st.write(f'Número de anillos alifáticos: {num_aliphatic_rings}')

  num_aromatic_rings = Lipinski.NumAromaticRings(mol)
  st.write(f'Número de anillos aromáticos: {num_aromatic_rings}')

  fraction_csp3 = Lipinski.FractionCSP3(mol)
  st.write(f'Proporción de átomos de carbono híbrido SP3: {fraction_csp3}')

  inde_flex = round(Crippen.MolMR(mol),4)
  st.write(f'Índice de refracción molar: {inde_flex}')

  mol_log_p = round(Descriptors.MolLogP(mol),4)
  st.write(f'Coeficiente de partición grasa-agua: {mol_log_p}')

  tpsa = round(Descriptors.TPSA(mol), 4)
  st.write(f'Área de superficie del polo topológico: {tpsa}')


