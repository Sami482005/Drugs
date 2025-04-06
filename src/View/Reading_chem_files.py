from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from PIL import Image, ImageDraw, ImageFont
import numpy as np


def draw_molecule_image(file_path, output_path="molecule_preview.png"):
    mol = None
    if file_path.endswith(".mol"):
        mol = Chem.MolFromMolFile(file_path)
    elif file_path.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(file_path)
        mol = next((m for m in suppl if m is not None), None)
    else:
        print("Unsupported file format. Use .mol or .sdf")
        return None

    if mol is None:
        print("No valid molecule found.")
        return None

    name = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown Molecule"
    img = Draw.MolToImage(mol)

    # Create a new image with space for the name
    img_width, img_height = img.size
    new_height = img_height + 30
    new_img = Image.new("RGB", (img_width, new_height), "white")
    new_img.paste(img, (0, 0))

    draw = ImageDraw.Draw(new_img)
    font = ImageFont.load_default()
    draw.text((10, img_height + 5), f"{name}", fill="black", font=font)

    new_img.save(output_path)
    return output_path

def get_molecule_name(file_path):
    mol = None
    if file_path.endswith(".mol"):
        mol = Chem.MolFromMolFile(file_path)
    elif file_path.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(file_path)
        mol = next((m for m in suppl if m is not None), None)
    
    if mol is None:
        return "Unknown Molecule"

    return mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown Molecule"

def get_molecule_matrices(file_path):
    """
    Returns a dictionary containing:
    - adjacency matrix
    - distance matrix
    - number of atoms
    - RDKit molecule object
    """
    mol = None
    if file_path.endswith(".mol"):
        mol = Chem.MolFromMolFile(file_path)
    elif file_path.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(file_path)
        mol = next((m for m in suppl if m is not None), None)

    if mol is None:
        print("No valid molecule found.")
        return None

    mol = Chem.AddHs(mol)  # Add hydrogens if needed for better topology

    # Get adjacency and distance matrices
    adj_matrix = rdmolops.GetAdjacencyMatrix(mol)
    dist_matrix = rdmolops.GetDistanceMatrix(mol)

    return {
        "adjacency": np.array(adj_matrix),
        "distance": np.array(dist_matrix),
        "num_atoms": mol.GetNumAtoms(),
        "mol": mol
    }