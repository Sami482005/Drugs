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
    def __init__(self, parent):
        super().__init__(parent)  # Pass the parent to the superclass
        self.parent = parent
        self.parent.title("Indices Calculator")
        self.parent.geometry("500x500")

        self.frames = []  # Stack to keep track of frames

        # Add a welcome label
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

        # By default, always tell the user that no file or directory is selected
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

        # Load the indices selected by the user
        self.selected_indices = []
        if self.edge_var.get(): self.selected_indices.append("Edge Density")
        if self.weiner_var.get(): self.selected_indices.append("Weiner Index")
        if self.petitjean_var.get(): self.selected_indices.append("Petitjean Index")
        
        # There are two frames we prepared. one for only 1 molecule (mol file) and one for multiple molecules
        # (sdf file or directory of mol files)
        if self.selected_path.endswith(".mol"):
            frame = results_frame.ResultsFrame(self, self, self.selected_path, self.selected_indices)
            self.show_frame(frame)
            return
        frame = multi_file.SDFDirFrame(self, self, self.selected_path, self.selected_indices)
        self.show_frame(frame)

    def show_previous_frame(self):
        # Show the previous frame
        if len(self.frames) > 1:
            # Remove the current frame from the stack
            current_frame = self.frames.pop()
            current_frame.pack_forget()  # Hide the current frame

            # Show the previous frame
            previous_frame = self.frames[-1]
            previous_frame.pack(fill="both", expand=True)
        else:
            messagebox.showinfo("Info", "No previous frame to show.")

    def show_frame(self, frame):
        # Push the current frame onto the stack
        if self.frames:
            self.frames[-1].pack_forget()  # Hide the current frame
        self.frames.append(frame)  # Add the new frame to the stack
        frame.pack(fill="both", expand=True)  # Show the new frame
        
if __name__ == "__main__":
    root = tk.Tk()  # Create the main Tkinter window
    app = FileUploaderApp(root)  # Pass the root as parent
    app.pack(fill="both", expand=True)
    root.mainloop()  # Start the Tkinter event loop