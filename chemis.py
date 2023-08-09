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
    iupac_name = st.text_input("Ingrese el nombre IUPAC en inglés", value='Glucose')

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

                mole3d = py3Dmol.view(query='cid:' + str(cid))
                mole3d.setStyle({'stick': {'color': 'spectrum'}})

                if st.button("Densidad Electrónica"):
                    mole3d.addSurface('MS', {'opacity': 0.7, 'colorscheme': {'gradient': 'rwb'}})
                showmol(mole3d, height=500, width=800)

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

def validate_non_empty_input(inputs):
    return all(inputs)

def balanceo_químico():
    st.title("Balanceador de ecuaciones químicas 🧪")
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




