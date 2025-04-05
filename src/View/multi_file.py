'''
This frame will be the same as the mol file frame
but instead it will allow me to move back and forth between the different structures
found in the sdf file on in the directory of mol files.

'''
import tkinter as tk

class SDFDirFrame(tk.Toplevel):
    #Look back at mol_view to properly initiate this one like that one
    def __init__(self, parent):
        super().__init__(parent)
        self.title("SDF/Directory Processor")
        self.geometry("400x300")
        label = tk.Label(self, text="Processing SDF File or Directory...", font=("Times New Roman", 18))
        label.pack(pady=20)
