import tkinter as tk
from tkinter import filedialog, Checkbutton, IntVar, Label, Button
import os

class IndexSelectorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Index Selector & File Uploader")
        self.root.geometry("1080x720")

        # Checkboxes for indices
        self.index_vars = {
            "Wiener Index": IntVar(),
            "Petitjean Index": IntVar(),
            "TPSA": IntVar()
        }

        Label(root, text="Select Indices:", font=("Arial", 20)).pack(anchor="w", padx=10, pady=5)
        for index, var in self.index_vars.items():
            Checkbutton(root, text=index, variable=var).pack(anchor="w", padx=20)

        # File selection button
        Button(root, text="Select Files", command=self.upload_files).pack(pady=10)
        Button(root, text="Select Directory", command=self.upload_directory).pack(pady=10)

        # Display selected files
        self.file_label = Label(root, text="No files selected", fg="gray")
        self.file_label.pack()

    def upload_files(self):
        files = filedialog.askopenfilenames(title="Select .mol or .sdf files", filetypes=[("MOL/SDF Files", "*.mol;*.sdf")])
        self.display_selected_files(files)

    def upload_directory(self):
        directory = filedialog.askdirectory(title="Select Directory with .mol Files")
        if directory:
            files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".mol")]
            self.display_selected_files(files)

    def display_selected_files(self, files):
        if files:
            file_names = "\n".join(os.path.basename(f) for f in files)
            self.file_label.config(text=f"Selected Files:\n{file_names}", fg="black")
        else:
            self.file_label.config(text="No valid files selected", fg="gray")

if __name__ == "__main__":
    root = tk.Tk()
    app = IndexSelectorApp(root)
    root.mainloop()
