import streamlit as st
import py3Dmol
import matplotlib.pyplot as plt
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski, Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem

from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from chempy import balance_stoichiometry
from stmol import showmol

import pickle
import pandas as pd

with open('gb_model.pkl', 'rb') as f:
    model = pickle.load(f)

def AromaticAtoms(m):
  aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
  aa_count = []
  for i in aromatic_atoms:
    if i == True:
      aa_count.append(i)
  sum_aa_count = sum(aa_count)
  return sum_aa_count

def iupac_to_cid(iupac_name):
    result = pcp.get_compounds(iupac_name, 'name')
    if result:
        return result[0].cid
    else:
        return None

def iupac_to_smiles(iupac_name):
    result = pcp.get_compounds(iupac_name, 'name')
    if result:
        return result[0].isomeric_smiles
    else:
        return None

def visualize_molecule():
    st.title("Propiedades de Solubilidad y Vizualización de Moléculas ⚗️")

    # Casilla de entrada para el nombre en IUPAC
    iupac_name = st.text_input("Ingrese el nombre IUPAC en inglés", value='1,3,7-Trimethylpurine-2,6-dione')

    if iupac_name:
        # Convertir el nombre en IUPAC a SMILES utilizando la función iupac_to_smiles
        smiles = iupac_to_smiles(iupac_name)
        cid = iupac_to_cid(iupac_name)

        if smiles:
            # Visualizar la molécula correspondiente utilizando la función visualize_molecule
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)

            AllChem.ComputeGasteigerCharges(mol)
            contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
            fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
            

            col1, col2 = st.columns(2)

            with col1:

                st.image(img, use_column_width=True, caption=f'Imagen en 2D de la estructura de la molécula {iupac_name}.')
                st.pyplot(fig)
                st.markdown(f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Imagen de las cargas parciales de Gasteiger para la molécula de {iupac_name}.</h2>", unsafe_allow_html=True)

                mole3d = py3Dmol.view(query='cid:' + str(cid))
                mole3d.setStyle({'stick': {'color': 'spectrum'}})

                if st.button("Densidad Electrónica"):
                    mole3d.addSurface('MS', {'opacity': 0.7, 'colorscheme': {'gradient': 'rwb'}})
                showmol(mole3d, height=500, width=800)
                #st.markdown(f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Modelo tridimensional de la molécula de {iupac_name} con su distribución de densidad electrónica </h2>", unsafe_allow_html=True)

            with col2:
                mol_wt = round(Descriptors.MolWt(mol), 4)
                st.write(f"Peso molecular: <span style='font-weight: bold;'> {mol_wt} </span>", unsafe_allow_html=True)

                num_atoms = mol.GetNumAtoms()
                st.write(f"Número de átomos: <span style='font-weight: bold;'>{num_atoms}</span>", unsafe_allow_html=True)

                num_bonds = mol.GetNumBonds()
                st.write(f"Número de enlaces: <span style='font-weight: bold;'>{num_bonds}</span>", unsafe_allow_html=True)

                num_acceptors = Lipinski.NumHAcceptors(mol)
                st.write(f"Número de aceptores de enlaces de hidrógeno: <span style='font-weight: bold;'>{num_acceptors}</span>", unsafe_allow_html=True)

                num_donors = Lipinski.NumHDonors(mol)
                st.write(f"Número de donantes de enlaces de hidrógeno: <span style='font-weight: bold;'>{num_donors}</span>", unsafe_allow_html=True)

                num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
                st.write(f"Número de enlaces giratorios:<span style='font-weight: bold;'> {num_rotatable_bonds}</span>", unsafe_allow_html=True)

                num_aliphatic_rings = Lipinski.NumAliphaticRings(mol)
                st.write(f"Número de anillos alifáticos:<span style='font-weight: bold;'> {num_aliphatic_rings}</span>", unsafe_allow_html=True)

                num_aromatic_rings = Lipinski.NumAromaticRings(mol)
                st.write(f"Número de anillos aromáticos:<span style='font-weight: bold;'> {num_aromatic_rings}</span>", unsafe_allow_html=True)

                fraction_csp3 = Lipinski.FractionCSP3(mol)
                st.write(f"Proporción de átomos de carbono híbrido SP3: <span style='font-weight: bold;'>{fraction_csp3}</span>", unsafe_allow_html=True)

                inde_flex = round(Crippen.MolMR(mol),4)
                st.write(f"Índice de refracción molar:<span style='font-weight: bold;'> {inde_flex}</span>", unsafe_allow_html=True)

                mol_log_p = round(Descriptors.MolLogP(mol),4)
                st.write(f"Coeficiente de partición grasa-agua:<span style='font-weight: bold;'> {mol_log_p}</span>", unsafe_allow_html=True)

                tpsa = round(Descriptors.TPSA(mol), 4)
                st.write(f"Área de superficie del polo topológico:<span style='font-weight: bold;'> {tpsa}</span>", unsafe_allow_html=True)

                bo = Descriptors.NumRotatableBonds(mol)
                atom_heav = Descriptors.HeavyAtomCount(mol)
                atom_aroma = AromaticAtoms(mol)
                aro_portion = atom_aroma / atom_heav

                dat_logs = (mol_log_p, mol_wt, bo, aro_portion)
                dat_logs = pd.DataFrame(dat_logs).T

                st.write("Estimación de la solubilidad mendiante un modelo de machine learning: ")
                logs = model.predict(dat_logs)
                logs = logs[0] 
                st.write(f"Log de la solubilidad acuosa de Delaney (LogS):<span style='font-weight: bold;'> {round(logs, 4)} mol/L </span>", unsafe_allow_html=True)


        else:
            st.error("No se encontró una molécula para el nombre en IUPAC proporcionado.")
    st.markdown(f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Modelo tridimensional de la molécula de {iupac_name} con su distribución de densidad electrónica. </h2>", unsafe_allow_html=True)
def validate_non_empty_input(inputs):
    return all(inputs)

def balanceo_químico():
    st.title("Balanceador de Ecuaciones Químicas 🧪")
    num_reactivos = st.number_input("Ingrese el número de reactivos:", min_value=1, step=1, value=2)

    reactivos = []
    cols = st.columns(num_reactivos)
    for i in range(num_reactivos):
        with cols[i]:
            if i == 0:
                reactivos.append(st.text_input(f"Reactivo {i + 1}", value="H2"))
            elif i == 1:
                reactivos.append(st.text_input(f"Reactivo {i + 1}", value="O2"))

    num_productos = st.number_input("Ingrese el número de productos:", min_value=1, step=1, value=1)
    productos = []
    cols = st.columns(num_productos)
    for n in range(num_productos):
        with cols[n]:
            productos.append(st.text_input(f"Producto {n + 1}", value="H2O"))

    if st.button("Balancear"):
        if validate_non_empty_input(reactivos) and validate_non_empty_input(productos):
            reac, prod = balance_stoichiometry(reactivos, productos)
            rx = dict(reac)
            px = dict(prod)
            reac_str = ' + '.join([f'{v} {k}' for k, v in rx.items()])
            prod_str = ' + '.join([f'{v} {k}' for k, v in px.items()])
            st.markdown(f"La reacción química balanceada queda:<br>"
                        f"<span style='font-size: 20px;'>{reac_str}  :arrow_right: {prod_str}</span>",
                        unsafe_allow_html=True)
        else:
            st.write("Ingrese al menos un reactivo y un producto para balancear la ecuación.")


def main():
    # Menú lateral
    st.sidebar.title("Opciones")
    selected_option = st.sidebar.selectbox("Seleccione una opción", ["Solubilidad", "Balance Químico"])

    if selected_option == "Solubilidad":
        visualize_molecule()
    elif selected_option == "Balance Químico":
        balanceo_químico()

if __name__ == '__main__':
    main()




