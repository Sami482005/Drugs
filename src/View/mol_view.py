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
        self.controller = controller
        self.file_path = file_path  # Store the file path
        self.selected_indices = selected_indices  # Store selected indices

        self.controller.title("Mol File Results")
        self.controller.geometry("700x700")

        Label(self, text="Processing MOL File...", font=("Times New Roman", 18, "bold")).pack(pady=10)

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

        # Back button
        Button(self, text="Back", command=self.go_back).pack(pady=20, side="left", anchor="sw")

        # Load results after initialization
        self.load_results()

    # Loads and displays the results for the selected file.
    def load_results(self):
        
        # Get molecule image and name
        image_path = RF.get_image(self.file_path)
        molecule_name = RF.get_molecule_name(self.file_path)

        # Compute indices based on selection
        indices = self.calculate_indices(self.file_path)

        # Display image if available
        if image_path and os.path.exists(image_path):
            img = Image.open(image_path)
            img = img.resize((200, 200), Image.Resampling.LANCZOS)
            img_tk = ImageTk.PhotoImage(img)
            self.image_label.config(image=img_tk)
            self.image_label.image = img_tk

        # Update molecule name and indices
        self.molecule_name_label.config(text=f"Molecule: {molecule_name}")
        self.petitjean_label.config(text=f"PetitJean: {indices.get('Petitjean', 'N/A')}")
        self.weiner_label.config(text=f"Weiner: {indices.get('Weiner', 'N/A')}")
        self.edge_density_label.config(text=f"Edge Density: {indices.get('EdgeDensity', 'N/A')}")

    def calculate_indices(self, file_path):
        # Calculates the selected indices and returns a dictionary with results.
        results = {}
        if "Edge Density" in self.selected_indices:
            results["EdgeDensity"] = IC.getEdgeDensity(file_path)
        if "Weiner Index" in self.selected_indices:
            results["Weiner"] = IC.getWeinerIndex(file_path)
        if "Petitjean Index" in self.selected_indices:
            results["Petitjean"] = IC.getPetitjeanIndex(file_path)
        return results

    def go_back(self):
        # Handles navigation back to the previous frame
        self.controller.show_previous_frame()
