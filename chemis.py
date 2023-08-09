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
    st.title("Visualización de moléculas")

    # Casilla de entrada para el nombre en IUPAC
    iupac_name = st.text_input("Ingrese el nombre IUPAC en inglés")

    if iupac_name:
        # Convertir el nombre en IUPAC a SMILES utilizando la función iupac_to_smiles
        smiles = iupac_to_smiles(iupac_name)

        if smiles:
            # Visualizar la molécula correspondiente utilizando la función visualize_molecule
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)

            AllChem.ComputeGasteigerCharges(mol)
            contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
            fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
            

            col1, col2 = st.columns(2)

            with col1:

                st.image(img, use_column_width=True, caption=f'Imagen en 2D de la estructura de la molécula {iupac_name}')
                st.pyplot(fig)
                st.markdown(f"<h2 style='text-align: center; color: #a0a0a0; font-size: 13px;'>Imagen de las cargas parciales de Gasteiger para la molécula {iupac_name}</h2>", unsafe_allow_html=True)
              

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





        else:
            st.error("No se encontró una molécula para el nombre en IUPAC proporcionado.")



def chemical_balancing():
    st.title("Balanceo")
    # Aquí puedes agregar tu lógica de balanceo químico

def main():
    # Menú lateral
    st.sidebar.title("Opciones")
    selected_option = st.sidebar.selectbox("Seleccione una opción", ["Visualización", "Balanceo químico"])

    if selected_option == "Visualización":
        visualize_molecule()
    elif selected_option == "Balanceo":
        chemical_balancing()

if __name__ == '__main__':
    main()




