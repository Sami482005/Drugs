import os  # Paths
import networkx as nx  # Graphs
import numpy as np  # Matrices + Arrays
import matplotlib.pyplot as plt  # Drawing
import tkinter as tk
from tkinter import filedialog, messagebox, Label
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import py3Dmol
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image, ImageTk

def load_molecules(file_path):
    mols = []
    if os.path.isfile(file_path):
        if file_path.endswith(".mol"):
            mol = read_mol_from_file(file_path)
            if mol:
                mols.append((os.path.basename(file_path), mol))
        elif file_path.endswith(".sdf"):
            for i, mol in enumerate(read_mols_from_sdf(file_path)):
                if mol:
                    mols.append((f"Molecule_{i+1}", mol))
    elif os.path.isdir(file_path):
        for file in os.listdir(file_path):
            if file.endswith(".mol"):
                full_path = os.path.join(file_path, file)
                mol = read_mol_from_file(full_path)
                if mol:
                    mols.append((file, mol))
    else:
        raise ValueError("Unsupported input. Provide a .mol, .sdf file, or a directory with .mol files.")
    return mols


def parse_into_sdf(directory):
    output_file = os.path.join(directory, "combined.sdf")
    mol_files = [f for f in os.listdir(directory) if f.endswith(".mol")]
    with open(output_file, 'w') as outfile:
        for mol_file in mol_files:
            with open(os.path.join(directory, mol_file), 'r') as infile:
                content = infile.read().strip()
                outfile.write(content + "\n$$$$\n")
    return output_file
def mol_to_nx(mol):
    G = nx.Graph()

    # Add atoms as nodes with properties
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_num=atom.GetAtomicNum(), symbol=atom.GetSymbol())

    # Add bonds as edges with properties
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=str(bond.GetBondType())
        )

    return G

def read_mols_from_sdf(file_path):
    suppl = SDMolSupplier(file_path)
    return [mol for mol in suppl if mol is not None]

def read_mol_from_file(file_path):
    return Chem.MolFromMolFile(file_path)


def get_molecule_matrices(mol):
    # Get the NetworkX graph from the RDKit molecule
    G = mol_to_nx(mol)
    num_atoms = G.number_of_nodes()
    adj_matrix = nx.to_numpy_array(G, dtype=int)

    # Calculate shortest path distances using NetworkX
    length_dict = dict(nx.all_pairs_shortest_path_length(G))
    dist_matrix = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        for j in range(num_atoms):
            dist_matrix[i][j] = length_dict[i].get(j, -1)  # Use -1 if no path exists

    return {
        "adjacency": adj_matrix,
        "distance": dist_matrix,
        "num_atoms": num_atoms,
        "graph": G
    }

def getEdgeDensity(mol):
    data = get_molecule_matrices(mol)
    adj = data["adjacency"]
    num_atoms = data["num_atoms"]

    num_edges = np.sum(adj) // 2  # Each bond is counted twice
    max_edges = num_atoms * (num_atoms - 1) / 2  # Max possible edges

    return round(num_edges / max_edges, 4)

def getWeinerIndex(mol):
    data = get_molecule_matrices(mol)
    dist = data["distance"]
    
    N = len(dist)
    wiener_index = 0
    for i in range(N):
        for j in range(i + 1, N):
            wiener_index += dist[i][j]
    
    return wiener_index

def getPetitjeanIndex(mol):
    data = get_molecule_matrices(mol)
    dist = data["distance"]
    eccentricities = LongestShortest_paths(dist)
    D = max(eccentricities)
    R = min(eccentricities)
    return round((D - R) / R, 4)

def LongestShortest_paths(matrix):
    N = len(matrix)
    eccentricities = [0] * N
    for i in range(N):
        Max = 0
        for j in range(N):
            if i != j:
                Max = max(Max, matrix[i][j])
        eccentricities[i] = Max
    return eccentricities

class FileUploaderApp(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.parent.title("Indices Calculator")
        self.parent.geometry("500x500")
        self.frames = []
        tk.Label(self, text="Welcome to the Indices Calculator", fg="blue", font=("Times New Roman", 24)).pack(pady=10)
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
        if not (self.edge_var.get() or self.weiner_var.get() or self.petitjean_var.get()):
            messagebox.showerror("Error", "Please select at least one index to calculate.")
            return

        if not self.selected_path:
            messagebox.showerror("Error", "Please upload a file or directory first.")
            return

        self.selected_indices = []
        if self.edge_var.get(): self.selected_indices.append("Edge Density")
        if self.weiner_var.get(): self.selected_indices.append("Weiner Index")
        if self.petitjean_var.get(): self.selected_indices.append("Petitjean Index")

        if self.selected_path.endswith(".mol") or self.selected_path.endswith(".sdf") or os.path.isdir(self.selected_path):
            frame = ResultsFrame(self, self, self.selected_path, self.selected_indices)
            self.show_frame(frame)

    def show_frame(self, frame):
        if self.frames:
            old_frame = self.frames.pop()
            old_frame.destroy()
        self.frames.append(frame)
        frame.pack(fill="both", expand=True)

class ResultsFrame(tk.Frame):
    def __init__(self, parent, controller, file_path, selected_indices):
        for widget in parent.winfo_children():
            widget.destroy()

        super().__init__(parent)
        self.controller = controller
        self.file_path = file_path
        self.selected_indices = selected_indices
        self.controller.parent.title("Mol File Results")
        self.controller.parent.geometry("700x800")

        canvas = tk.Canvas(self)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.pack(fill="both", expand=True)

        molecules = load_molecules(file_path)
        for name, mol in molecules:
            G = mol_to_nx(mol)
            self.display_molecule_results(scrollable_frame, name, G, mol)


    def display_molecule_results(self, frame, name, G, mol):
        # Convert RDKit molecule to image
        img = Draw.MolToImage(mol, size=(600, 600))
        img_tk = ImageTk.PhotoImage(img)

        # Display image in Tkinter
        img_label = tk.Label(frame, image=img_tk)
        img_label.image = img_tk  # Keep a reference to avoid garbage collection
        img_label.pack(pady=10)

        Label(frame, text=f"Molecule: {name}", font=("Times New Roman", 14, "bold")).pack(pady=(5, 10))

        indices = self.calculate_indices(mol)

        if "Petitjean Index" in self.selected_indices:
            Label(frame, text=f"PetitJean: {indices.get('Petitjean', 'N/A')}").pack()
        if "Weiner Index" in self.selected_indices:
            Label(frame, text=f"Weiner: {indices.get('Weiner', 'N/A')}").pack()
        if "Edge Density" in self.selected_indices:
            Label(frame, text=f"Edge Density: {indices.get('EdgeDensity', 'N/A')}").pack()

        tk.Label(frame, text="-" * 60).pack(pady=15)


    
    def calculate_indices(self, mol):
        results = {}
        try:
            if "Edge Density" in self.selected_indices:
                results["EdgeDensity"] = getEdgeDensity(mol)
            if "Weiner Index" in self.selected_indices:
                results["Weiner"] = getWeinerIndex(mol)
            if "Petitjean Index" in self.selected_indices:
                results["Petitjean"] = getPetitjeanIndex(mol)
        except Exception as e:
            print(f"Error calculating indices for molecule: {e}")
        return results

      

if __name__ == "__main__":
    root = tk.Tk()  # Create the main Tkinter window
    app = FileUploaderApp(root)  # Pass the root as parent
    app.pack(fill="both", expand=True)
    root.mainloop()  # Start the Tkinter event loop
