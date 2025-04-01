'''
This will be the displayer of the results for 1 molecule (provided by one mol file) 
A similar page will be created if an sdf or a directory of mol files is uploaded to just 
include an additional next page button that will display the results of the next molecule in the list.
'''
import tkinter as tk
from tkinter import Label, Button
from PIL import Image, ImageTk
import control.Indeces_Calculation as IC
import control.Reading_chem_files as RF

class ResultsFrame(tk.Frame):
    def __init__(self, parent, controller, get_image, get_data):
        super().__init__(parent)
        self.controller = controller
        
        self.controller.title("Results")
        self.controller.geometry("700x700")
        Label(self, text="Results", font=("Times New Roman", 16, "bold")).pack(pady=10)
        
        # Image display
        self.image_label = Label(self)
        self.image_label.pack()
        
        # Image name
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


    # To get the data to fill out the labels    
    def load_results(self):
        image_path = RF.get_image_path()
        molecule_name = RF.get_molecule_name()
        indices = IC.get_indices()
        
        if image_path:
            img = Image.open(image_path)
            img = img.resize((200, 200), Image.Resampling.LANCZOS)
            img_tk = ImageTk.PhotoImage(img)
            self.image_label.config(image=img_tk)
            self.image_label.image = img_tk
        
        self.molecule_name_label.config(text=molecule_name)
        self.petitjean_label.config(text=f"PetitJean: {indices[0]}")
        self.weiner_label.config(text=f"Weiner: {indices[1]}")
        self.edge_density_label.config(text=f"Edge Density: {indices[2]}")
    
    def go_back(self):
        self.controller.show_previous_frame()



    # Function to execute the selected functions on the uploaded file or directory
    # This function will be called when the user clicks the "Go" button
    def execute_selected_functions(self, file_path):
        if self.edge_var.get():
            IC.getEdgeDensity(file_path)
        if self.weiner_var.get():
            IC.getWeinerIndex(file_path)
        if self.petitjean_var.get():
            IC.getPetitjeanIndex(file_path)
    
    # Function to execute the selected functions on all mol files in the uploaded directory
    # This function will be called when the user clicks the "Go" button
    def execute_selected_functions_for_dir(self, dir_path):
        for filename in os.listdir(dir_path):
            if filename.endswith(".mol"):
                file_path = os.path.join(dir_path, filename)
                if self.edge_var.get():
                    IC.getEdgeDensity(file_path)
                if self.weiner_var.get():
                    IC.getWeinerIndex(file_path)
                if self.petitjean_var.get():
                    IC.getPetitjeanIndex(file_path)