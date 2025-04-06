'''
This frame will be the same as the mol file frame
but instead it will allow me to move back and forth between the different structures
found in the sdf file on in the directory of mol files.

'''
import tkinter as tk
from tkinter import Label, Button
from PIL import ImageTk
from rdkit import Chem
import Indeces_Calculation as IC


class SDFDirFrame(tk.Frame):
    def __init__(self, root, sdf_file):
        self.root = root
        self.root.title("Molecule Viewer")

        self.molecules = [mol for mol in Chem.SDMolSupplier(sdf_file) if mol is not None]
        self.index = 0

        self.image_label = Label(root)
        self.image_label.pack()

        self.info_label = Label(root, text="", font=("Arial", 12), justify="left")
        self.info_label.pack(pady=10)

        self.next_button = Button(root, text="Next", command=self.show_next_molecule)
        self.next_button.pack(pady=5)

        self.show_molecule(self.index)

    def show_molecule(self, index):
        # Load the molecule
        # Loop over Reading_chem_files.loading_sdf_file and get the name and the image
        # The before finishing the loop, get the indices (it will call itself the function for the matrices)
        # Print the results into the page.
        # Once the user clicks next, it will go to the next molecule in the list
        # Show computed metrics
        return
        
    def show_next_molecule(self):
        self.index = (self.index + 1) % len(self.molecules)
        self.show_molecule(self.index)

# --- Entry point ---
if __name__ == "__main__":
    root = tk.Tk()
    viewer = SDFDirFrame(root, "your_sdf_file.sdf")  # Replace with your SDF file path
    root.mainloop()
