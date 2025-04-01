'''
 This code will show the first page of the app
 The user should be able to select which of the three indeces they want to calculate: Edge Density,
 Weiner index, Petitjean Index.
 The user should be able to upload an sdf or mol file or a folder that contains several mol files
'''

import tkinter as tk
from tkinter import filedialog, messagebox
import os

class FileUploaderApp:
    # Initializing the elements of the main page
    def __init__(self, root):
        self.root = root
        self.root.title("Main Page")
        self.root.geometry("500x400")  # Increased frame size
        
        self.label = tk.Label(root, text="Welcome to the Indeces Calculator", fg="blue", font=("Times New Roman", 14))
        self.label.pack(pady=10)
        
        self.edge_var = tk.BooleanVar()
        self.weiner_var = tk.BooleanVar()
        self.petitjean_var = tk.BooleanVar()
        
        self.checkbox_frame = tk.Frame(root)
        self.checkbox_frame.pack(pady=5)
        
        # I need three checkboxes for the three indeces
        self.edge_cb = tk.Checkbutton(self.checkbox_frame, text="Edge Density", variable=self.edge_var)
        self.edge_cb.pack(side=tk.LEFT, padx=5)
        
        self.weiner_cb = tk.Checkbutton(self.checkbox_frame, text="Weiner Index", variable=self.weiner_var)
        self.weiner_cb.pack(side=tk.LEFT, padx=5)
        
        self.petitjean_cb = tk.Checkbutton(self.checkbox_frame, text="Petitjean", variable=self.petitjean_var)
        self.petitjean_cb.pack(side=tk.LEFT, padx=5)
        
        # I need a button to upload a file or a directory
        self.upload_button = tk.Button(root, text="Upload sdf or mol file", command=self.upload_file)
        self.upload_button.pack(pady=5)
        
        self.upload_dir_button = tk.Button(root, text="Upload Directory of mol files", command=self.upload_directory)
        self.upload_dir_button.pack(pady=5)
        
        # Always remind the user that they haven't uploaded a file or directory until they do
        self.selected_label = tk.Label(root, text="No file or directory selected", fg="red")
        self.selected_label.pack(pady=5)
        
        # To move to next page
        self.go_button = tk.Button(root, text="Go", command=self.next_frame)
        self.go_button.pack(side=tk.RIGHT, padx=10, pady=10)
    
    # Function to allow uploading sdf or mol files
    def upload_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Chemical Files", "*.sdf *.mol")])
        if file_path:
            self.selected_label.config(text=f"Selected file: {os.path.basename(file_path)}", fg="green")
            self.execute_selected_functions(file_path) # Calculate the indeces chosen on the selected file
    
    # Function to allow uploading a directory of mol files
    def upload_directory(self):
        dir_path = filedialog.askdirectory()
        if dir_path:
            self.selected_label.config(text=f"Selected directory: {os.path.basename(dir_path)}", fg="green")
            self.execute_selected_functions_for_dir(dir_path) # Calculate the indeces chosen on the mol files in the selected directory
    
    

    # Function to move to the next page
    # Checks that the user has chosen at least one index and uploaded a file or directory
    # If not, it will show an error message
    def next_frame(self):
        if not (self.edge_var.get() or self.weiner_var.get() or self.petitjean_var.get()):
            messagebox.showerror("Error", "Please select at least one index to calculate.")
            return
        
        if not hasattr(self, 'selected_label') or "No file or directory selected" in self.selected_label.cget("text"):
            messagebox.showerror("Error", "Please upload a file or directory first.")
            return
        
        # Proceed to the next frame or functionality
        

if __name__ == "__main__":
    root = tk.Tk()
    app = FileUploaderApp(root)
    root.mainloop()
