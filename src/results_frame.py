'''
This will be the displayer of the results for 1 molecule (provided by one mol file) 
A similar page will be created if an sdf or a directory of mol files is uploaded to just 
include an additional next page button that will display the results of the next molecule in the list.
'''
import tkinter as tk
from tkinter import Label, Button
from PIL import Image, ImageTk
import os
import Indeces_Calculation as IC
import Reading_chem_files as RF


class ResultsFrame(tk.Frame):
    def __init__(self, parent, controller, file_path, selected_indices):
        super().__init__(parent)

        canvas = tk.Canvas(self, borderwidth=0)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        super().__init__(parent)
        self.controller = controller
        self.file_path = file_path  
        self.selected_indices = selected_indices  

        # Use the controller to set the title and geometry
        self.controller.parent.title("Mol File Results")
        self.controller.parent.geometry("700x800")

        # Image display
        self.image_label = Label(self)
        self.image_label.pack()

        # Molecule name
        self.molecule_name_label = Label(self, text="", font=("Times New Roman", 12))
        self.molecule_name_label.pack(pady=5)

        # Indices labels
        self.petitjean_label = Label(self, text="PetitJean: ")
        self.petitjean_label.pack()

        self.weiner_label = Label(self, text="Weiner: ")
        self.weiner_label.pack()

        self.edge_density_label = Label(self, text="Edge Density: ")
        self.edge_density_label.pack()

        # Load results after initialization
        self.load_results()

    # Loads and displays the results for the selected file.
    def load_results(self):
        
        # Get molecule image and name
        image_path, molecule_name, mol = RF.draw_molecule_image(self.file_path, None)

        # Compute indices based on selection
        indices = self.calculate_indices(mol)

        if image_path and os.path.exists(image_path):
            img = Image.open(image_path)
            img = img.resize((500, 500), Image.Resampling.LANCZOS)
            img_tk = ImageTk.PhotoImage(img)
            self.image_label.config(image=img_tk)
            self.image_label.image = img_tk

        # Update molecule name and indices
        self.molecule_name_label.config(text=f"Molecule: {molecule_name}")
        self.petitjean_label.config(text=f"PetitJean: {indices.get('Petitjean', 'N/A')}")
        self.weiner_label.config(text=f"Weiner: {indices.get('Weiner', 'N/A')}")
        self.edge_density_label.config(text=f"Edge Density: {indices.get('EdgeDensity', 'N/A')}")

    def calculate_indices(self, mol):
        # Calculates the selected indices and returns a dictionary with results.
        results = {}
        if "Edge Density" in self.selected_indices:
            results["EdgeDensity"] = IC.getEdgeDensity(mol)
        if "Weiner Index" in self.selected_indices:
            results["Weiner"] = IC.getWeinerIndex(mol)
        if "Petitjean Index" in self.selected_indices:
            results["Petitjean"] = IC.getPetitjeanIndex(mol)
        return results
