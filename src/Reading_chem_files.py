from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import os

'''
This is the file that will allow me to read the mol files.
To process the sdf files, I be sending one mol file at a time to this file from multi_file.py
The process will also be the same for a directory of mol files.
This should return: the image of the structure, name, and the matrices needed to calculate the indices in 
Indeces_Calculation.py
'''

# This function will draw the molecule image and save it to a file.
# This function will be called from the results_frame.py and multi_file.py 
# and would allow the user to see the structure of the molecule and the name of the molecule
# I know the name is inconvenient, but it's because I figured out how to do it on the spot and I don't want to change the name and forget 
# to fix it everywhere else.
# In the multi_file.py, I will create a list: the result ofloading_sdf_file 
# and then I will loop over the list and get the name and the image

def draw_molecule_image(file_path=None, mol=None):
    if mol is None:
        if file_path.endswith(".mol"):
            mol = Chem.MolFromMolFile(file_path)
    
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown Molecule"
    # Save the image to a file
    output_path="data/images"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_path = os.path.join("images", f"{name}.png")

    img = Draw.MolToImage(mol)

    img_width, img_height = img.size
    new_height = img_height + 30
    new_img = Image.new("RGB", (img_width, new_height), "white")
    new_img.paste(img, (0, 0))

    draw = ImageDraw.Draw(new_img)
    font = ImageFont.load_default()
    draw.text((10, img_height + 5), f"{name}", fill="black", font=font)

    # mol is a RDKit molecule object
    return output_path, name, mol

# This function will return the adjacency and distance matrices of the molecule.
def get_molecule_matrices(mol):
    """
    Returns a dictionary containing:
    - adjacency matrix (NumPy array)
    - distance matrix (NumPy array)
    - number of atoms
    - RDKit molecule object
    """

    if mol is None:
        return None

    # Generate adjacency and distance matrices
    adj_matrix = np.array(rdmolops.GetAdjacencyMatrix(mol))
    dist_matrix = np.array(rdmolops.GetDistanceMatrix(mol))

    return {
        "adjacency": adj_matrix,
        "distance": dist_matrix,
    }

def loading_sdf_files(file_path):
    # Load SDF files and return a list of RDKit molecule objects.
    suppl = Chem.SDMolSupplier(file_path)
    mols = [mol for mol in suppl if mol is not None]
    return mols