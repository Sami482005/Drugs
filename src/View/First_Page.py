'''
 This code will show the first page of the app
 The user should be able to select which of the three indeces they want to calculate: Edge Density,
 Weiner index, Petitjean Index.
 The user should be able to upload an sdf or mol file or a folder that contains several mol files
'''

import tkinter as tk
from tkinter import filedialog, messagebox
import os

import mol_view
import multi_file

class FileUploaderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Main Page")
        self.root.geometry("500x400")
        
        self.selected_path = None  # Stores selected file or directory
        self.selected_indices = [] # Stores the selected indices to calculate
        self.label = tk.Label(root, text="Welcome to the Indices Calculator", fg="blue", font=("Times New Roman", 24))
        self.label.pack(pady=10)

        self.edge_var = tk.BooleanVar()
        self.weiner_var = tk.BooleanVar()
        self.petitjean_var = tk.BooleanVar()

        self.checkbox_frame = tk.Frame(root)
        self.checkbox_frame.pack(pady=5)

        self.edge_cb = tk.Checkbutton(self.checkbox_frame, text="Edge Density", variable=self.edge_var, font=("Times New Roman", 20))
        self.edge_cb.pack(side=tk.LEFT, padx=5)

        self.weiner_cb = tk.Checkbutton(self.checkbox_frame, text="Weiner Index", variable=self.weiner_var, font=("Times New Roman", 20))
        self.weiner_cb.pack(side=tk.LEFT, padx=5)

        self.petitjean_cb = tk.Checkbutton(self.checkbox_frame, text="Petitjean", variable=self.petitjean_var, font=("Times New Roman", 20))
        self.petitjean_cb.pack(side=tk.LEFT, padx=5)

        self.upload_button = tk.Button(root, text="Upload sdf or mol file", font=("Times New Roman", 18), command=self.upload_file)
        self.upload_button.pack(pady=5)

        self.upload_dir_button = tk.Button(root, text="Upload Directory of mol files", font=("Times New Roman", 18), command=self.upload_directory)
        self.upload_dir_button.pack(pady=5)

        self.selected_label = tk.Label(root, text="No file or directory selected", fg="red", font=("Times New Roman", 20))
        self.selected_label.pack(pady=5)

        self.go_button = tk.Button(root, text="Go", command=self.open_correct_frame, font=("Times New Roman", 20))
        self.go_button.pack(side=tk.RIGHT, padx=10, pady=10)

    def upload_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Chemical Files", "*.sdf *.mol")])
        if file_path:
            self.selected_path = file_path
            self.selected_label.config(text=f"Selected file: {os.path.basename(file_path)}", fg="green", font=("Times New Roman", 18))


    def upload_directory(self):
        dir_path = filedialog.askdirectory()
        if dir_path:
            self.selected_path = dir_path
            self.selected_label.config(text=f"Selected directory: {os.path.basename(dir_path)}", fg="green", font=("Times New Roman", 18))

    def open_correct_frame(self):
        # Error Handling
        if not (self.edge_var.get() or self.weiner_var.get() or self.petitjean_var.get()):
            messagebox.showerror("Error", "Please select at least one index to calculate.")
            return

        if not self.selected_path:
            messagebox.showerror("Error", "Please upload a file or directory first.")
            return
        
        # Decide which frame to open based on the uploaded file type.
        if self.selected_path.endswith(".mol"):
            if self.edge_var.get():
                self.selected_indices.append("Edge Density")
            if self.weiner_var.get():
                self.selected_indices.append("Weiner Index")
            if self.petitjean_var.get():
                self.selected_indices.append("Petitjean Index")
            results_frame = mol_view.ResultsFrame(self.root, self.root, self.selected_path, self.selected_indices)
            self.show_frame(results_frame)
        else:
            if self.edge_var.get():
                self.selected_indices.append("Edge Density")
            if self.weiner_var.get():
                self.selected_indices.append("Weiner Index")
            if self.petitjean_var.get():
                self.selected_indices.append("Petitjean Index")
            results_frame = multi_file.SDFDirFrame(self.root, self.root, self.selected_path, self.selected_indices)  # Open SDF/Directory processing frame
            self.show_frame(results_frame)

if __name__ == "__main__":
    root = tk.Tk()
    app = FileUploaderApp(root)
    root.mainloop()