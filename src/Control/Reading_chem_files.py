from rdkit import Chem
from rdkit.Chem import Draw
import sys
from PIL import Image, ImageDraw, ImageFont
import os
import Indeces_Calculation as IC

def draw_molecules(file_path):
    molecules = []
    
    # Determine file type and read molecules
    if file_path.endswith(".mol"):
        mol = Chem.MolFromMolFile(file_path)
        if mol:
            molecules.append(mol)
    elif file_path.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(file_path)
        molecules = [m for m in suppl if m is not None]
    else:
        print("Unsupported file format. Please provide a .mol or .sdf file.")
        return
    
    if not molecules:
        print("No valid molecules found in file.")
        return
    
    images = []
    if os.path.isdir(file_path):
        # If a directory is provided, process all .mol files in the directory
        for filename in os.listdir(file_path):
            if filename.endswith(".mol"):
                mol = Chem.MolFromMolFile(os.path.join(file_path, filename))
                if mol:
                    molecules.append(mol)

    for mol in molecules:
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown Molecule"
        img = Draw.MolToImage(mol)
        
        # Create an image with space for the name and indices
        img_width, img_height = img.size
        new_height = img_height + 50  # Additional space for text
        new_img = Image.new("RGB", (img_width, new_height), "white")
        new_img.paste(img, (0, 0))
        
        draw = ImageDraw.Draw(new_img)
        font = ImageFont.load_default()
        
        
        edge, weiner, petitjean = IC.getEdgeDensity(file_path), IC.getWeinerIndex(file_path), IC.getPetitjeanIndex(file_path)
        
        draw.text((10, img_height + 5), f"Name: {name}", fill="black", font=font)
        draw.text((10, img_height + 25), f"Indices: {edge}, {weiner}, {petitjean}", fill="black", font=font)
        
        images.append(new_img)
    
    # Display all images
    for img in images:
        img.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python draw_molecule.py <path_to_mol_or_sdf>")
    else:
        draw_molecules(sys.argv[1])