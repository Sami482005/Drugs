'''
 This code will show the first page of the app
 The user should be able to select which of the three indeces they want to calculate: Edge Density,
 Weiner index, Petitjean Index.
 The user should be able to upload an sdf or mol file or a folder that contains several mol files
'''

import tkinter as tk
from tkinter import filedialog, messagebox
import os
import results_frame
import multi_file

class FileUploaderApp(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.title("Indices Calculator")
        self.geometry("500x500")

        tk.Label(self, text="Welcome to the Indices Calculator", fg="blue", font=("Times New Roman", 24)).pack(pady=10)

        # Allow the user to select which indices they want to calculate
        self.edge_var = tk.BooleanVar()
        self.weiner_var = tk.BooleanVar()
        self.petitjean_var = tk.BooleanVar()

        checkbox_frame = tk.Frame(self)
        checkbox_frame.pack(pady=5)

        tk.Checkbutton(checkbox_frame, text="Edge Density", variable=self.edge_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Weiner Index", variable=self.weiner_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Petitjean", variable=self.petitjean_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)

        tk.Button(self, text="Upload sdf or mol file", font=("Times New Roman", 18), command=self.upload_file).pack(pady=5)
        tk.Button(self, text="Upload Directory of mol files", font=("Times New Roman", 18), command=self.upload_directory).pack(pady=5)

        self.selected_label = tk.Label(self, text="No file or directory selected", fg="red", font=("Times New Roman", 20))
        self.selected_label.pack(pady=5)

        tk.Button(self, text="Go", font=("Times New Roman", 20), command=self.open_correct_frame).pack(side=tk.RIGHT, padx=10, pady=10)

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
        self.selected_indices = []
        if self.edge_var.get(): self.selected_indices.append("Edge Density")
        if self.weiner_var.get(): self.selected_indices.append("Weiner Index")
        if self.petitjean_var.get(): self.selected_indices.append("Petitjean Index")
        
        if os.path.endswith(self.selected_path, (".mol")):
            frame = results_frame.ResultsFrame(self.container, self, self.selected_path, self.selected_indices)
            self.show_frame(frame)
            return
        frame = multi_file.SDFDirFrame(self.container, self, self.selected_path, self.selected_indices)
        self.show_frame(frame)


if __name__ == "__main__":
    app = FileUploaderApp()
    app.mainloop()